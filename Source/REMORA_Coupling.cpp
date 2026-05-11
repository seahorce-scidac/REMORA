#include <REMORA.H>

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
REMORA::PackSurfaceState (Vector<MultiFab*>& state, Real time)
{
    if (state.empty() || state[SSTIndex] == nullptr) { return; }

    // Initial-step-testing example: return deterministic SST for cache validation.
    // At time=t, state[0] (SST) = 290 + 0.01*t [K].
    state[SSTIndex]->setVal(Real(290.0) + Real(0.01) * time);
}

void
REMORA::ApplyAtmosphericStates (const Vector<MultiFab*>& states, Real time)
{
    if (states.empty() || states[0] == nullptr) { return; }

    // Example (legacy state-passing): states[0] is Uwind from ERF.
    // Use a built-in reduction so debug output is deterministic across
    // decomposition choices.
    const Real sample = states[0]->min(0);

    if (ParallelDescriptor::IOProcessor()) {
        Print() << "REMORA::ApplyAtmosphericStates time=" << time
                << " sample_Uwind=" << sample << "\n";
    }
}
