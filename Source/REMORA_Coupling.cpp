#include <REMORA.H>

#include <AMReX_MFIter.H>
#include <AMReX_Print.H>

using namespace amrex;

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

    // Stage-0 coupling uses synthetic SST so the driver can validate slab ownership.
    fill_multifab(*state[SSTIndex], Real(290.0) + Real(0.01) * time);
}

void
REMORA::ApplyAtmosphericStates (const Vector<MultiFab*>& states, Real time)
{
    if (states.empty() || states[0] == nullptr) { return; }

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
