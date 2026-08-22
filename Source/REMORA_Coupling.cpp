#include <REMORA.H>

#include <AMReX_BCRec.H>
#include <AMReX_Box.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_Geometry.H>
#include <AMReX_Interpolater.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_Print.H>
#include <AMReX_Reduce.H>

#include <limits>

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

// -------------------------------------------------------------------------
// NEW: Conservative Sparse Matrix Remap Engine (Reverse: OCN -> ATM)
// -------------------------------------------------------------------------

// Source index region each destination box's stencils actually reference.
//
// The entries in index_mf are *source* indices, and for non-conformal grids they
// routinely fall far outside the destination box's own index range, so staging the
// source on `dst.boxArray()` plus a fixed ghost halo leaves those reads outside
// the valid region - zero at best, out of bounds at worst. Mirror of the same
// helper on the ERF side of the coupling.
//
// The result becomes a BoxArray, which must be globally consistent, hence the
// reduction.
amrex::BoxArray
StagedSourceBoxArray (const amrex::MultiFab& src,
                      const amrex::MultiFab& dst,
                      const amrex::iMultiFab& index_mf,
                      int max_stencil_size)
{
    using namespace amrex;

    const Box src_domain = src.boxArray().minimalBox();
    const int nboxes = static_cast<int>(dst.boxArray().size());
    constexpr int int_big = std::numeric_limits<int>::max();

    Vector<int> lo(2 * nboxes,  int_big);
    Vector<int> hi(2 * nboxes, -int_big);

    for (MFIter mfi(dst); mfi.isValid(); ++mfi) {
        const int b = mfi.index();
        const Box bx = mfi.validbox();
        auto const& idx = index_mf.const_array(mfi);

        ReduceOps<ReduceOpMin, ReduceOpMin, ReduceOpMax, ReduceOpMax> reduce_op;
        ReduceData<int, int, int, int> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

        reduce_op.eval(bx, reduce_data,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
        {
            int i_min = int_big, j_min = int_big, i_max = -int_big, j_max = -int_big;
            for (int m = 0; m < max_stencil_size; ++m) {
                const int si = idx(i, j, k, m * 3);
                const int sj = idx(i, j, k, m * 3 + 1);
                if (si < 0 || sj < 0) { continue; }   // -1 marks an unused slot
                i_min = amrex::min(i_min, si); i_max = amrex::max(i_max, si);
                j_min = amrex::min(j_min, sj); j_max = amrex::max(j_max, sj);
            }
            return {i_min, j_min, i_max, j_max};
        });

        auto const& hv = reduce_data.value(reduce_op);
        lo[2*b  ] = amrex::get<0>(hv); lo[2*b+1] = amrex::get<1>(hv);
        hi[2*b  ] = amrex::get<2>(hv); hi[2*b+1] = amrex::get<3>(hv);
    }

    ParallelDescriptor::ReduceIntMin(lo.dataPtr(), static_cast<int>(lo.size()));
    ParallelDescriptor::ReduceIntMax(hi.dataPtr(), static_cast<int>(hi.size()));

    // The staged boxes must carry the source's index type from the start:
    // intersecting a cell-centered box with a staggered one trips
    // Box::operator&='s sameType assertion.
    const IndexType src_ixtype = src.boxArray().ixType();

    BoxList bl(src_ixtype);
    for (int b = 0; b < nboxes; ++b) {
        Box need;
        if (lo[2*b] > hi[2*b] || lo[2*b+1] > hi[2*b+1]) {
            // No stencil anywhere in this destination box. A BoxArray cannot hold
            // an empty box, so stage a single cell; nothing reads it.
            need = Box(src_domain.smallEnd(), src_domain.smallEnd(), src_ixtype);
        } else {
            need = Box(IntVect(lo[2*b], lo[2*b+1], src_domain.smallEnd(2)),
                       IntVect(hi[2*b], hi[2*b+1], src_domain.bigEnd(2)),
                       src_ixtype);
            need &= src_domain;
        }
        bl.push_back(need);
    }
    return BoxArray(std::move(bl));
}

// fallback_val fills destination cells that no source cell overlaps.
//
// Zero is the wrong default for a physical field: the receiving model cannot tell
// "no data here" from "the ocean says zero", and for an intensive quantity zero is
// not merely inaccurate, it drives the receiver's flux formulas outside their valid
// domain. This mirrors the rationale on the atmosphere->ocean twin at
// ERF/Source/Coupling/ERF_to_REMORA.cpp:122.
//
// The ocean->atmosphere SST lane deliberately leaves this at zero and relies on the
// per-cell coverage flag instead: where REMORA does not cover an ERF cell, ERF keeps
// its own wrflowinp SST rather than consuming a value invented here. That is a
// withholding contract, not a value, so no climatological constant belongs at that
// call site. The parameter exists because the mask blend below needs it and because
// other lanes have a meaningful value to pass.
void
ApplyConservativeRemap (const amrex::MultiFab& src,
                        amrex::MultiFab& dst,
                        const amrex::MultiFab& weight_mf,
                        const amrex::iMultiFab& index_mf,
                        int max_stencil_size,
                        const amrex::MultiFab* dst_mask = nullptr,
                        const amrex::iMultiFab* dst_land_mask = nullptr,
                        amrex::Real fallback_val = amrex::Real(0.0))
{
    using namespace amrex;

    // 1. Data Routing: Route REMORA SST data onto the ERF atmospheric layout,
    // staging exactly the source region the stencils reference.
    MultiFab src_on_dst(StagedSourceBoxArray(src, dst, index_mf, max_stencil_size),
                        dst.DistributionMap(), src.nComp(), 0);
    src_on_dst.setVal(0.0);
    src_on_dst.ParallelCopy(src);

    dst.setVal(0.0);

    // 2. Stencil Application: Execute the local sparse dot product
    for (MFIter mfi(dst, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        Box bx = mfi.tilebox();

        auto const& w_arr   = weight_mf.const_array(mfi);
        auto const& idx_arr = index_mf.const_array(mfi);
        auto const& src_arr = src_on_dst.const_array(mfi);
        auto        dst_arr = dst.array(mfi);
        const bool has_mask = (dst_mask != nullptr);
        auto const& mask_arr = has_mask ? dst_mask->const_array(mfi) : Array4<const Real>{};
        // ERF land mask: 0 = water, anything non-zero is land. Destination cells
        // over land are zeroed out (wet/dry masking only, no vector rotation).
        const bool has_land_mask = (dst_land_mask != nullptr);
        auto const& land_arr = has_land_mask ? dst_land_mask->const_array(mfi) : Array4<const int>{};

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            // No stencil entry at all means no source cell overlapped this
            // destination cell. Hand back the fallback rather than an accumulated
            // zero, and do not apply the masks to it: a masked-out cell still gets
            // evaluated by the receiving model's flux formulas before the mask is
            // applied, so it too must hold an admissible value. Mirrors
            // ERF/Source/Coupling/ERF_to_REMORA.cpp:174.
            if (idx_arr(i, j, k, 0) < 0) {
                dst_arr(i, j, k) = fallback_val;
                return;
            }

            Real sum = 0.0;
            for (int m = 0; m < max_stencil_size; ++m) {
                Real w = w_arr(i, j, k, m);
                if (w > 0.0) {
                    int src_i = idx_arr(i, j, k, m * 3);
                    int src_j = idx_arr(i, j, k, m * 3 + 1);
                    int src_k = idx_arr(i, j, k, m * 3 + 2);

                    sum += w * src_arr(src_i, src_j, src_k);
                }
            }
            // Blend toward the fallback rather than multiplying by the mask.
            // Multiplying drives partially masked cells toward exactly zero, which
            // is right for a flux lane (fallback_val = 0, so this reduces to
            // sum *= mask and is bit-identical) but destructive for an intensive
            // state lane such as SST, where it scales a temperature toward 0 K.
            // Mirrors ERF/Source/Coupling/ERF_to_REMORA.cpp:202.
            if (has_mask) {
                const Real mask = mask_arr(i, j, k);
                sum = mask * sum + (Real(1.0) - mask) * fallback_val;
            }
            // Test against zero, not against 1: ERF stamps lmask = 2 for cells
            // under ImmersedForcing buildings (ERF/Source/ERF_MakeNewArrays.cpp:890),
            // so an == 1 test reads every building as water and hands it a remapped
            // SST. The driver's atmosphere->ocean side already uses the tolerant
            // form (ERFRemoraMultiBlockContainer.cpp:1767), so == 1 also made the
            // two directions disagree about the same cell within one timestep.
            if (has_land_mask && land_arr(i, j, k) != 0) { sum = Real(0.0); }
            dst_arr(i, j, k) = sum;
        });
    }
}
}

amrex::Real
REMORA::EvolveOneStep (amrex::Real /*time*/, amrex::Real /*dt_request*/)
{
    Real cur_time = t_new[0];
    const int step = istep[0];

    if (cur_time >= stop_time) {
        return zero;
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
REMORA::SetLongwaveFromDriver ()
{
    // "constant" is any non-computed type: it makes setup_step pass
    // vec_longwave_down to bulk_fluxes instead of leaving lw_ptr null. The two
    // flags below are what ReadParameters would have derived for that type.
    solverChoice.bulk_flux_type[BulkFlux::LWrad] = BulkForcingType::constant;
    solverChoice.longwave_down   = true;
    solverChoice.longwave_is_net = false;
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

void
REMORA::GetAtmosToOceanPsiCoordinates (const amrex::MultiFab*& x_psi,
                                       const amrex::MultiFab*& y_psi) const
{
    if (vec_xp.empty() || vec_yp.empty() ||
        vec_xp[0] == nullptr || vec_yp[0] == nullptr) {
        x_psi = nullptr;
        y_psi = nullptr;
        return;
    }
    x_psi = vec_xp[0].get();
    y_psi = vec_yp[0].get();
}

void
REMORA::GetAtmosToOceanPsiLonLat (const amrex::MultiFab*& lon_psi,
                                  const amrex::MultiFab*& lat_psi) const
{
    if (vec_lonp.empty() || vec_latp.empty() ||
        vec_lonp[0] == nullptr || vec_latp[0] == nullptr) {
        lon_psi = nullptr;
        lat_psi = nullptr;
        return;
    }
    lon_psi = vec_lonp[0].get();
    lat_psi = vec_latp[0].get();
}

void
REMORA::GetLandSeaMasks (const amrex::MultiFab*& mskr,
                         const amrex::MultiFab*& msku,
                         const amrex::MultiFab*& mskv) const
{
    mskr = (!vec_mskr.empty() && vec_mskr[0]) ? vec_mskr[0].get() : nullptr;
    msku = (!vec_msku.empty() && vec_msku[0]) ? vec_msku[0].get() : nullptr;
    mskv = (!vec_mskv.empty() && vec_mskv[0]) ? vec_mskv[0].get() : nullptr;
}

/*
 * \brief Extracts SST from the 3D conservative state for the atmospheric driver.
 *
 * Reads Temp_comp at the top water-column cell (k_sfc), converts from
 * Celsius to Kelvin, and conservatively remaps the result into state[SSTIndex].
 */
void
REMORA::PackSurfaceState (Vector<MultiFab*>& state,
                          Real /*time*/,
                          const amrex::MultiFab* weight_o2a_mf,
                          const amrex::iMultiFab* index_o2a_mf,
                          int max_stencil_size,
                          const amrex::iMultiFab* dst_land_mask)
{
    if (state.empty() || state[SSTIndex] == nullptr) { return; }
    const int lev = 0;

    // REMORA stores temperature in Celsius. Surface is at k=N (top of water column).
    const int k_sfc = cons_new[lev]->boxArray().minimalBox().bigEnd(2);

    // Build a temp MultiFab on REMORA's ba2d (k=0) derived from cons_new's BoxArray.
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
            t(i, j, 0) = c(i, j, k_sfc, Temp_comp) + Real(273.15);
        });
    }

    MultiFab& dst = *state[SSTIndex];

    if (weight_o2a_mf != nullptr && index_o2a_mf != nullptr) {
        // Execute sparse conservative remap from REMORA SST to ERF layout
        ApplyConservativeRemap(tmp, dst, *weight_o2a_mf, *index_o2a_mf, max_stencil_size,
                               nullptr, dst_land_mask);
    } else {
        // Fallback for un-stenciled or synthetic runs
        dst.setVal(zero);
        dst.ParallelCopy(tmp, 0, 0, 1);
    }
}

/*
 * \brief Receives atmospheric states from the driver and applies unit conversions.
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
            vec_Pair[0]->mult(Real(0.01), 0, 1);
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
            vec_Tair[0]->plus(Real(-273.15), 0, 1);
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

    const Real Hscale2 = one / (solverChoice.rho0 * Cp);
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

    tau_x_tmp.setVal(zero);
    tau_x_tmp.ParallelCopy(*states[AtmosFluxes::TauX], 0, 0, 1);
    tau_x_tmp.FillBoundary(geom[0].periodicity());

    tau_y_tmp.setVal(zero);
    tau_y_tmp.ParallelCopy(*states[AtmosFluxes::TauY], 0, 0, 1);
    tau_y_tmp.FillBoundary(geom[0].periodicity());

    shflux_tmp.setVal(zero);
    shflux_tmp.ParallelCopy(*states[AtmosFluxes::SHflux], 0, 0, 1);
    shflux_tmp.FillBoundary(geom[0].periodicity());

    lhflux_tmp.setVal(zero);
    lhflux_tmp.ParallelCopy(*states[AtmosFluxes::LHflux], 0, 0, 1);
    lhflux_tmp.FillBoundary(geom[0].periodicity());

    lwflux_tmp.setVal(zero);
    lwflux_tmp.ParallelCopy(*states[AtmosFluxes::LWrad], 0, 0, 1);
    lwflux_tmp.FillBoundary(geom[0].periodicity());

    vec_srflx[0]->setVal(zero);
    vec_srflx[0]->ParallelCopy(*states[AtmosFluxes::SWrad], 0, 0, 1);
    vec_srflx[0]->FillBoundary(geom[0].periodicity());

    vec_rain[0]->setVal(zero);
    vec_rain[0]->ParallelCopy(*states[AtmosFluxes::Rain], 0, 0, 1);
    vec_rain[0]->FillBoundary(geom[0].periodicity());

    vec_evap[0]->setVal(zero);
    vec_evap[0]->ParallelCopy(*states[AtmosFluxes::Evap], 0, 0, 1);
    vec_evap[0]->FillBoundary(geom[0].periodicity());

    vec_lrflx[0]->setVal(zero);
    vec_lhflx[0]->setVal(zero);
    vec_shflx[0]->setVal(zero);
    vec_stflux[0]->setVal(zero);

    for (MFIter mfi(*vec_sustr[0], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        Array4<Real> const& sustr = vec_sustr[0]->array(mfi);
        Array4<const Real> const& msku = vec_msku[0]->const_array(mfi);
        Array4<const Real> const& tau_x = tau_x_tmp.const_array(mfi);
        Box ubx = mfi.grownnodaltilebox(0, IntVect(NGROW,NGROW,0));
        Box ubxD = ubx;
        ubxD.makeSlab(2,0);
        ParallelFor(ubxD, [=] AMREX_GPU_DEVICE (int i, int j, int ) {
            sustr(i,j,0) = -tau_x(i,j,0) / rho0 * msku(i,j,0);
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
            svstr(i,j,0) = -tau_y(i,j,0) / rho0 * mskv(i,j,0);
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

    const Real sustr_min = vec_sustr[0]->min(0);
    const Real sustr_max = vec_sustr[0]->max(0);
    const Real svstr_min = vec_svstr[0]->min(0);
    const Real svstr_max = vec_svstr[0]->max(0);
    const Real stflux_temp_min = vec_stflux[0]->min(Temp_comp);
    const Real stflux_temp_max = vec_stflux[0]->max(Temp_comp);
    const Real stflux_salt_min = vec_stflux[0]->min(Salt_comp);
    const Real stflux_salt_max = vec_stflux[0]->max(Salt_comp);
    const Real srflx_min = vec_srflx[0]->min(0);
    const Real srflx_max = vec_srflx[0]->max(0);
    const Real lrflx_min = vec_lrflx[0]->min(0);
    const Real lrflx_max = vec_lrflx[0]->max(0);
    const Real lhflx_min = vec_lhflx[0]->min(0);
    const Real lhflx_max = vec_lhflx[0]->max(0);
    const Real shflx_min = vec_shflx[0]->min(0);
    const Real shflx_max = vec_shflx[0]->max(0);

    amrex::Print() << "REMORA ApplyAtmosphericFluxes validation:\n"
                   << "  sustr: min=" << sustr_min << " max=" << sustr_max << "\n"
                   << "  svstr: min=" << svstr_min << " max=" << svstr_max << "\n"
                   << "  stflux(Temp): min=" << stflux_temp_min << " max=" << stflux_temp_max << "\n"
                   << "  stflux(Salt): min=" << stflux_salt_min << " max=" << stflux_salt_max << "\n"
                   << "  srflx: min=" << srflx_min << " max=" << srflx_max << "\n"
                   << "  lrflx: min=" << lrflx_min << " max=" << lrflx_max << "\n"
                   << "  lhflx: min=" << lhflx_min << " max=" << lhflx_max << "\n"
                   << "  shflx: min=" << shflx_min << " max=" << shflx_max << "\n";
}
