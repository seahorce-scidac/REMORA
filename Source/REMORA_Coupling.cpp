#include <REMORA.H>

#include <AMReX_Box.H>
#include <AMReX_MFIter.H>
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

    if ( (plot_int > 0      && (step+1 - last_plot_file_step) == plot_int         ) ||
         (plot_int_time > 0 && (cur_time >= (last_plot_file_time + plot_int_time))) )
    {
        last_plot_file_step = step+1;
        last_plot_file_time = cur_time;
        WritePlotFile(step+1);
        history_count++;
    }

    if ((check_int > 0 && (step+1 - last_check_file_step) == check_int)
            || (check_int_time > 0 && cur_time >= (last_check_file_time + check_int_time))) {
        last_check_file_step = step+1;
        last_check_file_time = cur_time;
        WriteCheckpointFile();
    }

    post_timestep(step, cur_time, dt[0]);

    return dt[0];
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
    state[SSTIndex]->ParallelCopy(tmp, 0, 0, 1);
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
