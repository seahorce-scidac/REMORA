#include <REMORA.H>

#include <AMReX_BCRec.H>
#include <AMReX_Box.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_Geometry.H>
#include <AMReX_Interpolater.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_Print.H>

using namespace amrex;

/*
  Coupling reference context (implementation-side):

  1) Legacy state-passing contract:
     Warner et al. (2010), COAWST, Fig. 5 / Block B.
     REMORA receives atmospheric states and computes surface fluxes internally
     (bulk-physics/COARE-style path).

  2) Future direct flux-passing roadmap:
     COAWST's ATM2OCN_FLUXES pathway (documented in COAWST manuals/workshops
     and exercised in Zambon et al., 2014) motivates direct flux exchange
     (tau_x, tau_y, heat/moisture) instead of state-only exchange.
     COAWST code anchors:
       - Master/mct_roms_wrf.h
       - ROMS/Nonlinear/atm2ocn_flux.F
       - ROMS/Nonlinear/bulk_flux.F
*/

namespace {
constexpr int SSTIndex = 0;

Geometry
make_unit_slab_geometry (Box const& domain)
{
    static constexpr Real lo[AMREX_SPACEDIM] = {0.0_rt, 0.0_rt, 0.0_rt};
    static constexpr Real hi[AMREX_SPACEDIM] = {1.0_rt, 1.0_rt, 1.0_rt};
    static constexpr int periodicity[AMREX_SPACEDIM] = {0, 0, 0};
    RealBox rb(lo, hi);
    return Geometry(domain, &rb, CoordSys::cartesian, periodicity);
}

Vector<BCRec>
make_internal_bcs (int ncomp)
{
    Vector<BCRec> bcs(ncomp);
    for (auto& bc : bcs) {
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            bc.setLo(dir, BCType::int_dir);
            bc.setHi(dir, BCType::int_dir);
        }
    }
    return bcs;
}

void
CopyDriverSlabToRemoraLayout (const MultiFab& src,
                              MultiFab& dst,
                              const Geometry& geom,
                              const char* context)
{
    const Box src_domain = enclosedCells(src.boxArray().minimalBox());
    const Box dst_domain = enclosedCells(dst.boxArray().minimalBox());

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        src.boxArray().ixType() == dst.boxArray().ixType(),
        context);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        src_domain.smallEnd(0) == dst_domain.smallEnd(0) &&
        src_domain.smallEnd(1) == dst_domain.smallEnd(1) &&
        src_domain.length(0) == dst_domain.length(0) &&
        src_domain.length(1) == dst_domain.length(1) &&
        src_domain.length(2) == dst_domain.length(2),
        context);

    dst.setVal(0.0_rt);
    dst.ParallelCopy(src, 0, 0, dst.nComp(), IntVect(0), IntVect(0), geom.periodicity());
    dst.FillBoundary(geom.periodicity());
}
}

amrex::Real
REMORA::EvolveOneStep (amrex::Real /*time*/, amrex::Real /*dt_request*/)
{
    Real cur_time = t_new[0];
    const int step = istep[0];

    if (cur_time >= stop_time) {
        return Real(0.0);
    }

    ComputeDt();

    int lev = 0;
    int iteration = 1;
    if (max_level == 0) {
        timeStep(lev, cur_time, iteration);
    } else {
        timeStepML(cur_time, iteration);
    }

    cur_time += dt[0];

    WriteAtIntermediateTime(step, cur_time);

    post_timestep(step, cur_time, dt[0]);

    return dt[0];
}

void
REMORA::ConfigureDriverAtmosToOceanCoupling (bool use_coupling_driver,
                                             bool use_two_way_coupling,
                                             DriverAtmosForcingMode active_mode)
{
    running_with_coupling_driver = use_coupling_driver;
    driver_uses_two_way_coupling = use_two_way_coupling;
    driver_atmos_forcing_mode = active_mode;
}

void
REMORA::SetDriverAtmosToOceanForcingMode (DriverAtmosForcingMode mode)
{
    driver_atmos_forcing_mode = mode;
}

void
REMORA::GetAtmosToOceanRhoLayout (amrex::BoxArray& ba,
                                  amrex::DistributionMapping& dm) const
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        !vec_srflx.empty() && vec_srflx[0] != nullptr,
        "REMORA::GetAtmosToOceanRhoLayout requires post-InitData rho-point forcing storage.");
    ba = vec_srflx[0]->boxArray();
    dm = vec_srflx[0]->DistributionMap();
}

void
REMORA::GetAtmosToOceanUFaceLayout (amrex::BoxArray& ba,
                                    amrex::DistributionMapping& dm) const
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        !vec_sustr.empty() && vec_sustr[0] != nullptr,
        "REMORA::GetAtmosToOceanUFaceLayout requires post-InitData u-face forcing storage.");
    ba = vec_sustr[0]->boxArray();
    dm = vec_sustr[0]->DistributionMap();
}

void
REMORA::GetAtmosToOceanVFaceLayout (amrex::BoxArray& ba,
                                    amrex::DistributionMapping& dm) const
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        !vec_svstr.empty() && vec_svstr[0] != nullptr,
        "REMORA::GetAtmosToOceanVFaceLayout requires post-InitData v-face forcing storage.");
    ba = vec_svstr[0]->boxArray();
    dm = vec_svstr[0]->DistributionMap();
}

/*
 * \brief Extracts SST from the 3D conservative state for the atmospheric driver.
 *
 * Reads Temp_comp at the top water-column cell (k_sfc), converts from
 * Celsius to Kelvin, and copies the result into state[SSTIndex].
 *
 * @param[in,out] state  OCN2ATM slab buffer sized by the driver (one MultiFab
 *                       per ocean-to-atmosphere export layer). state[SSTIndex]
 *                       (index 0) is overwritten with sea-surface temperature
 *                       (SST) sampled from the Temp_comp tracer at the uppermost
 *                       sigma level (k = Nz) and converted from degrees Celsius
 *                       to Kelvin for the atmospheric driver. An empty vector or
 *                       null state[0] is treated as a no-op.
 * @param[in    ] time   Current ocean model time (unused; retained for driver
 *                       interface conformance).
 */
void
REMORA::PackSurfaceState (Vector<MultiFab*>& state, Real /*time*/)
{
    if (state.empty() || state[SSTIndex] == nullptr) { return; }
    const int lev = 0;

    // REMORA stores temperature in Celsius. Surface is at k=N (top of water column).
    const int k_sfc = cons_new[lev]->boxArray().minimalBox().bigEnd(2);

    // Build a temp MultiFab on REMORA's ba2d (k=0) derived from cons_new's BoxArray.
    // Same DistributionMap ensures each box is local; we fill at k=0 from cons at k=k_sfc.
    BoxList bl2d = cons_new[lev]->boxArray().boxList();
    for (auto& b : bl2d) { b.setRange(2, 0); }
    BoxArray ba2d(std::move(bl2d));
    MultiFab tmp(ba2d, cons_new[lev]->DistributionMap(), 1, 0);

    for (MFIter mfi(*cons_new[lev]); mfi.isValid(); ++mfi) {
        auto const& c = cons_new[lev]->const_array(mfi);
        auto         t = tmp.array(mfi);
        Box bx = makeSlab(mfi.validbox(), 2, k_sfc);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int) {
            // Write to k=0 in tmp (ba2d range); convert Celsius → Kelvin.
            t(i, j, 0) = c(i, j, k_sfc, Temp_comp) + 273.15_rt;
        });
    }

    MultiFab& dst = *state[SSTIndex];
    const Box src_domain = tmp.boxArray().minimalBox();
    const Box dst_domain = dst.boxArray().minimalBox();

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        src_domain.smallEnd(2) == src_domain.bigEnd(2) &&
        dst_domain.smallEnd(2) == dst_domain.bigEnd(2),
        "REMORA::PackSurfaceState expects 2D slab source and destination MultiFabs.");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        src_domain.smallEnd(0) == dst_domain.smallEnd(0) &&
        src_domain.smallEnd(1) == dst_domain.smallEnd(1),
        "REMORA::PackSurfaceState requires aligned atmosphere/ocean slab origins.");

    const IntVect src_len = src_domain.length();
    const IntVect dst_len = dst_domain.length();
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        src_len[2] == 1 && dst_len[2] == 1,
        "REMORA::PackSurfaceState expects unit-thickness source and destination slabs.");

    const bool same_xy = (src_len[0] == dst_len[0] && src_len[1] == dst_len[1]);
    const bool dst_finer =
        (dst_len[0] >= src_len[0] && dst_len[1] >= src_len[1] &&
         dst_len[0] % src_len[0] == 0 && dst_len[1] % src_len[1] == 0);
    const bool dst_coarser =
        (src_len[0] >= dst_len[0] && src_len[1] >= dst_len[1] &&
         src_len[0] % dst_len[0] == 0 && src_len[1] % dst_len[1] == 0);

    if (same_xy) {
        dst.ParallelCopy(tmp, 0, 0, 1);
        return;
    }

    if (dst_finer) {
        const IntVect ratio(dst_len[0] / src_len[0], dst_len[1] / src_len[1], 1);
        const auto src_geom = make_unit_slab_geometry(src_domain);
        const auto dst_geom = make_unit_slab_geometry(dst_domain);
        const auto bcs = make_internal_bcs(1);
        amrex::InterpFromCoarseLevel(dst, IntVect(0), IntVect(0),
                                     tmp, 0, 0, 1,
                                     src_geom, dst_geom,
                                     ratio, &pc_interp, bcs, 0);
        return;
    }

    if (dst_coarser) {
        const IntVect ratio(src_len[0] / dst_len[0], src_len[1] / dst_len[1], 1);
        amrex::average_down(tmp, dst, 0, 1, ratio);
        return;
    }

    amrex::Abort("REMORA::PackSurfaceState requires matching horizontal extents and integer-ratio atmosphere/ocean slab resolutions.");
}

/*
 * \brief Receives atmospheric states from the driver and applies unit conversions.
 *
 * Fills REMORA's internal forcing MultiFabs from states and records which
 * lanes were successfully updated in driver_atmos_state_from_driver.
 * Unit conversions applied: Pair Pa to mb; Tair K to Celsius.
 *
 * @param[in] states  ATM2OCN forcing slab buffer from the driver (Warner et al.
 *                    2010, Block B state-passing contract), indexed by AtmosState.
 *                    Expected units per lane:
 *                    Uwind/Vwind: 10-m winds [m/s];
 *                    Pair: mean sea-level pressure [Pa, converted to mb];
 *                    Qair: near-surface specific humidity [kg/kg];
 *                    Tair: 2-m air temperature [K, converted to degC];
 *                    Cloud: cloud fraction [0-1];
 *                    Rain: precipitation rate [kg/m^2/s];
 *                    SWrad/LWrad: downwelling shortwave/longwave radiation [W/m^2].
 *                    Missing lanes (null pointer or index out of range) are skipped;
 *                    driver_atmos_state_from_driver tracks populated lanes for the
 *                    bulk-flux parameterization fallback logic.
 * @param[in] time    Current ocean model time (unused; retained for driver
 *                    interface conformance).
 */
void
REMORA::ApplyAtmosphericStates (const Vector<MultiFab*>& states, Real /*time*/)
{
    running_with_coupling_driver = true;
    driver_atmos_forcing_mode = DriverAtmosForcingMode::State;
    driver_atmos_state_from_driver.fill(false);
    if (finest_level < 0) { return; }

    // Wind (m/s) — no unit conversion
    if (vec_uwind[0] != nullptr) {
        if (states.size() > AtmosState::Uwind && states[AtmosState::Uwind] != nullptr) {
            vec_uwind[0]->ParallelCopy(*states[AtmosState::Uwind], 0, 0, 1);
            vec_uwind[0]->FillBoundary(geom[0].periodicity());
            driver_atmos_state_from_driver[AtmosState::Uwind] = true;
        }
    }
    if (vec_vwind[0] != nullptr) {
        if (states.size() > AtmosState::Vwind && states[AtmosState::Vwind] != nullptr) {
            vec_vwind[0]->ParallelCopy(*states[AtmosState::Vwind], 0, 0, 1);
            vec_vwind[0]->FillBoundary(geom[0].periodicity());
            driver_atmos_state_from_driver[AtmosState::Vwind] = true;
        }
    }

    // Atmospheric pressure: Pa → mb (REMORA bulk flux expects mb)
    if (vec_Pair[0] != nullptr) {
        if (states.size() > AtmosState::Pair && states[AtmosState::Pair] != nullptr) {
            vec_Pair[0]->ParallelCopy(*states[AtmosState::Pair], 0, 0, 1);
            vec_Pair[0]->mult(0.01_rt, 0, 1);
            vec_Pair[0]->FillBoundary(geom[0].periodicity());
            driver_atmos_state_from_driver[AtmosState::Pair] = true;
        }
    }

    // Specific humidity (kg/kg) — no conversion
    if (vec_qair[0] != nullptr) {
        if (states.size() > AtmosState::Qair && states[AtmosState::Qair] != nullptr) {
            vec_qair[0]->ParallelCopy(*states[AtmosState::Qair], 0, 0, 1);
            vec_qair[0]->FillBoundary(geom[0].periodicity());
            driver_atmos_state_from_driver[AtmosState::Qair] = true;
        }
    }

    // Air temperature: K → °C (REMORA stores/uses Celsius internally)
    if (vec_Tair[0] != nullptr) {
        if (states.size() > AtmosState::Tair && states[AtmosState::Tair] != nullptr) {
            vec_Tair[0]->ParallelCopy(*states[AtmosState::Tair], 0, 0, 1);
            vec_Tair[0]->plus(-273.15_rt, 0, 1);
            vec_Tair[0]->FillBoundary(geom[0].periodicity());
            driver_atmos_state_from_driver[AtmosState::Tair] = true;
        }
    }

    // Cloud fraction [0-1], rain, SW/LW radiation — no unit conversion
    if (vec_cloud[0] != nullptr) {
        if (states.size() > AtmosState::Cloud && states[AtmosState::Cloud] != nullptr) {
            vec_cloud[0]->ParallelCopy(*states[AtmosState::Cloud], 0, 0, 1);
            vec_cloud[0]->FillBoundary(geom[0].periodicity());
            driver_atmos_state_from_driver[AtmosState::Cloud] = true;
        }
    }
    if (vec_rain[0] != nullptr) {
        if (states.size() > AtmosState::Rain && states[AtmosState::Rain] != nullptr) {
            vec_rain[0]->ParallelCopy(*states[AtmosState::Rain], 0, 0, 1);
            vec_rain[0]->FillBoundary(geom[0].periodicity());
            driver_atmos_state_from_driver[AtmosState::Rain] = true;
        }
    }
    if (vec_srflx[0] != nullptr) {
        if (states.size() > AtmosState::SWrad && states[AtmosState::SWrad] != nullptr) {
            vec_srflx[0]->ParallelCopy(*states[AtmosState::SWrad], 0, 0, 1);
            vec_srflx[0]->FillBoundary(geom[0].periodicity());
            driver_atmos_state_from_driver[AtmosState::SWrad] = true;
        }
    }
    if (vec_longwave_down[0] != nullptr) {
        if (states.size() > AtmosState::LWrad && states[AtmosState::LWrad] != nullptr) {
            vec_longwave_down[0]->ParallelCopy(*states[AtmosState::LWrad], 0, 0, 1);
            vec_longwave_down[0]->FillBoundary(geom[0].periodicity());
            driver_atmos_state_from_driver[AtmosState::LWrad] = true;
        }
    }

}

void
REMORA::ApplyAtmosphericFluxes (const Vector<MultiFab*>& states, Real /*time*/)
{
    running_with_coupling_driver = true;
    driver_atmos_forcing_mode = DriverAtmosForcingMode::Flux;
    driver_atmos_state_from_driver.fill(false);
    if (finest_level < 0) { return; }

    if (states.size() <= AtmosFluxes::Evap ||
        states[AtmosFluxes::TauX] == nullptr ||
        states[AtmosFluxes::TauY] == nullptr ||
        states[AtmosFluxes::SHflux] == nullptr ||
        states[AtmosFluxes::LHflux] == nullptr ||
        states[AtmosFluxes::SWrad] == nullptr ||
        states[AtmosFluxes::LWrad] == nullptr ||
        states[AtmosFluxes::Rain] == nullptr ||
        states[AtmosFluxes::Evap] == nullptr ||
        vec_sustr[0] == nullptr || vec_svstr[0] == nullptr ||
        vec_stflux[0] == nullptr || vec_mskr[0] == nullptr ||
        vec_msku[0] == nullptr || vec_mskv[0] == nullptr ||
        vec_srflx[0] == nullptr || vec_lrflx[0] == nullptr ||
        vec_lhflx[0] == nullptr || vec_shflx[0] == nullptr ||
        vec_rain[0] == nullptr || vec_evap[0] == nullptr) {
        return;
    }

    const Real Hscale2 = 1.0_rt / (solverChoice.rho0 * Cp);
    const Real rho0 = solverChoice.rho0;

    MultiFab tau_x_tmp(vec_sustr[0]->boxArray(), vec_sustr[0]->DistributionMap(), 1,
                       vec_sustr[0]->nGrowVect());
    MultiFab tau_y_tmp(vec_svstr[0]->boxArray(), vec_svstr[0]->DistributionMap(), 1,
                       vec_svstr[0]->nGrowVect());
    MultiFab shflux_tmp(vec_shflx[0]->boxArray(), vec_shflx[0]->DistributionMap(), 1,
                        vec_shflx[0]->nGrowVect());
    MultiFab lhflux_tmp(vec_lhflx[0]->boxArray(), vec_lhflx[0]->DistributionMap(), 1,
                        vec_lhflx[0]->nGrowVect());
    MultiFab lwflux_tmp(vec_lrflx[0]->boxArray(), vec_lrflx[0]->DistributionMap(), 1,
                        vec_lrflx[0]->nGrowVect());

    CopyDriverSlabToRemoraLayout(*states[AtmosFluxes::TauX], tau_x_tmp, geom[0],
                                 "REMORA::ApplyAtmosphericFluxes expected TauX slab to match REMORA u-face layout.");
    CopyDriverSlabToRemoraLayout(*states[AtmosFluxes::TauY], tau_y_tmp, geom[0],
                                 "REMORA::ApplyAtmosphericFluxes expected TauY slab to match REMORA v-face layout.");
    CopyDriverSlabToRemoraLayout(*states[AtmosFluxes::SHflux], shflux_tmp, geom[0],
                                 "REMORA::ApplyAtmosphericFluxes expected SHflux slab to match REMORA rho layout.");
    CopyDriverSlabToRemoraLayout(*states[AtmosFluxes::LHflux], lhflux_tmp, geom[0],
                                 "REMORA::ApplyAtmosphericFluxes expected LHflux slab to match REMORA rho layout.");
    CopyDriverSlabToRemoraLayout(*states[AtmosFluxes::LWrad], lwflux_tmp, geom[0],
                                 "REMORA::ApplyAtmosphericFluxes expected LWrad slab to match REMORA rho layout.");
    CopyDriverSlabToRemoraLayout(*states[AtmosFluxes::SWrad], *vec_srflx[0], geom[0],
                                 "REMORA::ApplyAtmosphericFluxes expected SWrad slab to match REMORA rho layout.");
    CopyDriverSlabToRemoraLayout(*states[AtmosFluxes::Rain], *vec_rain[0], geom[0],
                                 "REMORA::ApplyAtmosphericFluxes expected rain slab to match REMORA rho layout.");
    CopyDriverSlabToRemoraLayout(*states[AtmosFluxes::Evap], *vec_evap[0], geom[0],
                                 "REMORA::ApplyAtmosphericFluxes expected evap slab to match REMORA rho layout.");
    vec_lrflx[0]->setVal(0.0);
    vec_lhflx[0]->setVal(0.0);
    vec_shflx[0]->setVal(0.0);
    vec_stflux[0]->setVal(0.0);

    for (MFIter mfi(*vec_sustr[0], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        Array4<Real> const& sustr = vec_sustr[0]->array(mfi);
        Array4<const Real> const& msku = vec_msku[0]->const_array(mfi);
        Array4<const Real> const& tau_x = tau_x_tmp.const_array(mfi);
        Box ubx = mfi.grownnodaltilebox(0, IntVect(NGROW,NGROW,0));
        Box ubxD = ubx;
        ubxD.makeSlab(2,0);
        ParallelFor(ubxD, [=] AMREX_GPU_DEVICE (int i, int j, int ) {
            sustr(i,j,0) = tau_x(i,j,0) / rho0 * msku(i,j,0);
        });
    }

    for (MFIter mfi(*vec_svstr[0], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        Array4<Real> const& svstr = vec_svstr[0]->array(mfi);
        Array4<const Real> const& mskv = vec_mskv[0]->const_array(mfi);
        Array4<const Real> const& tau_y = tau_y_tmp.const_array(mfi);
        Box vbx = mfi.grownnodaltilebox(1, IntVect(NGROW,NGROW,0));
        Box vbxD = vbx;
        vbxD.makeSlab(2,0);

        ParallelFor(vbxD, [=] AMREX_GPU_DEVICE (int i, int j, int ) {
            svstr(i,j,0) = tau_y(i,j,0) / rho0 * mskv(i,j,0);
        });
    }

    for (MFIter mfi(*vec_stflux[0], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        Array4<Real> const& stflux = vec_stflux[0]->array(mfi);
        Array4<Real> const& lrflx = vec_lrflx[0]->array(mfi);
        Array4<Real> const& lhflx = vec_lhflx[0]->array(mfi);
        Array4<Real> const& shflx = vec_shflx[0]->array(mfi);
        Array4<const Real> const& mskr = vec_mskr[0]->const_array(mfi);
        Array4<const Real> const& srflx = vec_srflx[0]->const_array(mfi);
        Array4<const Real> const& rain = vec_rain[0]->const_array(mfi);
        Array4<const Real> const& evap = vec_evap[0]->const_array(mfi);
        Array4<const Real> const& shflux = shflux_tmp.const_array(mfi);
        Array4<const Real> const& lhflux = lhflux_tmp.const_array(mfi);
        Array4<const Real> const& lwflux = lwflux_tmp.const_array(mfi);

        Box gbx2 = mfi.growntilebox(IntVect(NGROW,NGROW,0));
        Box gbx2D = gbx2;
        gbx2D.makeSlab(2,0);

        ParallelFor(gbx2D, [=] AMREX_GPU_DEVICE (int i, int j, int ) {
            // ERF exports flux lanes in native surface-flux units; convert once at ingest.
            lrflx(i,j,0) = lwflux(i,j,0) * Hscale2;
            lhflx(i,j,0) = -lhflux(i,j,0) * Hscale2;
            shflx(i,j,0) = -shflux(i,j,0) * Hscale2;
            stflux(i,j,0,Temp_comp) =
                (srflx(i,j,0) * Hscale2 + lrflx(i,j,0) + lhflx(i,j,0) + shflx(i,j,0))
                * mskr(i,j,0);
            stflux(i,j,0,Salt_comp) =
                mskr(i,j,0) * (evap(i,j,0) - rain(i,j,0)) / rhow;
        });
    }

    vec_sustr[0]->FillBoundary(geom[0].periodicity());
    vec_svstr[0]->FillBoundary(geom[0].periodicity());
    vec_srflx[0]->FillBoundary(geom[0].periodicity());
    vec_lrflx[0]->FillBoundary(geom[0].periodicity());
    vec_lhflx[0]->FillBoundary(geom[0].periodicity());
    vec_shflx[0]->FillBoundary(geom[0].periodicity());
    vec_stflux[0]->FillBoundary(geom[0].periodicity());
    vec_rain[0]->FillBoundary(geom[0].periodicity());
    vec_evap[0]->FillBoundary(geom[0].periodicity());
    vec_stflux[0]->FillBoundary(geom[0].periodicity());

    if (amrex::ParallelDescriptor::IOProcessor()) {
        amrex::Print() << "REMORA ApplyAtmosphericFluxes validation:\n"
                       << "  sustr: min=" << vec_sustr[0]->min(0) << " max=" << vec_sustr[0]->max(0) << "\n"
                       << "  svstr: min=" << vec_svstr[0]->min(0) << " max=" << vec_svstr[0]->max(0) << "\n"
                       << "  stflux(Temp): min=" << vec_stflux[0]->min(Temp_comp) << " max=" << vec_stflux[0]->max(Temp_comp) << "\n"
                       << "  stflux(Salt): min=" << vec_stflux[0]->min(Salt_comp) << " max=" << vec_stflux[0]->max(Salt_comp) << "\n"
                       << "  srflx: min=" << vec_srflx[0]->min(0) << " max=" << vec_srflx[0]->max(0) << "\n"
                       << "  lrflx: min=" << vec_lrflx[0]->min(0) << " max=" << vec_lrflx[0]->max(0) << "\n"
                       << "  lhflx: min=" << vec_lhflx[0]->min(0) << " max=" << vec_lhflx[0]->max(0) << "\n"
                       << "  shflx: min=" << vec_shflx[0]->min(0) << " max=" << vec_shflx[0]->max(0) << "\n";
    }
}
