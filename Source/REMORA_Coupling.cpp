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

void
REMORA::ApplyAtmosphericStates (const Vector<MultiFab*>& states, Real time)
{
    if (states.empty() || states[0] == nullptr) { return; }
    if (finest_level < 0) { return; }

    // ParallelCopy from driver slab (k=0) into REMORA forcing arrays (also k=0).
    // Uses cross-BoxArray safe ParallelCopy rather than the asserting copy helper.
    auto safe_copy = [&] (int src_idx, MultiFab* dst) {
        if (dst == nullptr) { return; }
        if (src_idx >= static_cast<int>(states.size())) { return; }
        if (states[src_idx] == nullptr) { return; }
        dst->ParallelCopy(*states[src_idx], 0, 0, 1);
    };

    // Wind (m/s) — no unit conversion
    safe_copy(0, vec_uwind[0].get());
    safe_copy(1, vec_vwind[0].get());

    // Atmospheric pressure: Pa → mb (REMORA bulk flux expects mb)
    safe_copy(2, vec_Pair[0].get());
    if (vec_Pair[0]) { vec_Pair[0]->mult(0.01_rt, 0, 1); }

    // Specific humidity (kg/kg) — no conversion
    safe_copy(3, vec_qair[0].get());

    // Air temperature: K → °C (REMORA stores/uses Celsius internally)
    safe_copy(4, vec_Tair[0].get());
    if (vec_Tair[0]) { vec_Tair[0]->plus(-273.15_rt, 0, 1); }

    // Cloud fraction [0-1], rain, SW/LW radiation — no unit conversion
    safe_copy(5, vec_cloud[0].get());
    safe_copy(6, vec_rain[0].get());
    safe_copy(7, vec_srflx[0].get());
    safe_copy(8, vec_longwave_down[0].get());

    if (ParallelDescriptor::IOProcessor()) {
        Print() << "REMORA::ApplyAtmosphericStates time=" << time
                << " sample_Uwind=" << states[0]->min(0) << "\n";
    }
}
