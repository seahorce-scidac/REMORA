#include <REMORA.H>

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

void
fill_multifab (MultiFab& mf, Real value)
{
    mf.setVal(value);
}

}

void
REMORA::PackSurfaceState (Vector<MultiFab*>& state, Real time)
{
    if (state.empty() || state[SSTIndex] == nullptr) { return; }

    // Initial-step-testing example: return deterministic SST for cache validation.
    // At time=t, state[0] (SST) = 290 + 0.01*t [K].
    fill_multifab(*state[SSTIndex], Real(290.0) + Real(0.01) * time);
}

void
REMORA::ApplyAtmosphericStates (const Vector<MultiFab*>& states, Real time)
{
    if (states.empty() || states[0] == nullptr) { return; }

    // Example (legacy state-passing): states[0] is Uwind from ERF.
    // We read one sample value and print it so reviewers can verify
    // atmosphere->ocean data movement in initial-step-testing runs.
    Real sample = Real(0.0);
    for (MFIter mfi(*states[0], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto& bx = mfi.validbox();
        const auto arr = states[0]->const_array(mfi);
        const auto lo = lbound(bx);
        sample = arr(lo.x, lo.y, lo.z, 0);
        break;
    }

    if (ParallelDescriptor::IOProcessor()) {
        Print() << "REMORA::ApplyAtmosphericStates time=" << time
                << " sample_Uwind=" << sample << "\n";
    }
}
