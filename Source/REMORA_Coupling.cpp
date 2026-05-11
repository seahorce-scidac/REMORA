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

int
bottom_cell_k (const MultiFab& mf)
{
    return mf.boxArray().minimalBox().smallEnd(2);
}

int
top_cell_k (const MultiFab& mf)
{
    return mf.boxArray().minimalBox().bigEnd(2);
}

void
copy_plane_to_plane_xy (MultiFab& dst,
                        int dst_k,
                        const MultiFab& src,
                        int src_k)
{
    AMREX_ALWAYS_ASSERT(dst.nComp() >= 1);
    AMREX_ALWAYS_ASSERT(src.nComp() >= 1);
    AMREX_ALWAYS_ASSERT(dst.boxArray() == src.boxArray());
    AMREX_ALWAYS_ASSERT(dst.DistributionMap() == src.DistributionMap());

    for (MFIter mfi(dst, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        Box bx = makeSlab(mfi.validbox(), 2, dst_k);

        auto const& dst_arr = dst.array(mfi);
        auto const& src_arr = src.const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            dst_arr(i,j,dst_k) = src_arr(i,j,src_k);
        });
    }
}

void
copy_if_present (const Vector<MultiFab*>& states,
                 int src_idx,
                 MultiFab* dst)
{
    if (dst == nullptr) { return; }
    if (src_idx >= static_cast<int>(states.size())) { return; }
    if (states[src_idx] == nullptr) { return; }

    // Initial matched-grid assumption: identical horizontal decomposition.
    // We use an interface-face convention and derive source/destination cell
    // planes from it.
    // For ERF cell-centered states at the lower boundary, interface face is
    // at src.smallEnd(2), so the source cell adjacent to the interface is
    // src_k = src.smallEnd(2) (typically k=0).
    const int dst_k = bottom_cell_k(*dst);
    const int src_face_k = bottom_cell_k(*states[src_idx]);
    const int src_k = src_face_k;
    copy_plane_to_plane_xy(*dst, dst_k, *states[src_idx], src_k);
}

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

    // Legacy state-passing mapping (COAWST Block-B style):
    // 0:Uwind 1:Vwind 2:Pair 3:RH/Hair(qair path) 4:Tair
    // 5:cloud 6:rain 7:SWrad 8:LWrad
    if (finest_level >= 0) {
        copy_if_present(states, 0, vec_uwind[0].get());
        copy_if_present(states, 1, vec_vwind[0].get());
        copy_if_present(states, 2, vec_Pair[0].get());
        copy_if_present(states, 3, vec_qair[0].get());
        copy_if_present(states, 4, vec_Tair[0].get());
        copy_if_present(states, 5, vec_cloud[0].get());
        copy_if_present(states, 6, vec_rain[0].get());
        copy_if_present(states, 7, vec_srflx[0].get());
        copy_if_present(states, 8, vec_longwave_down[0].get());
    }

    // Example (legacy state-passing): states[0] is Uwind from ERF.
    // Use a built-in reduction so debug output is deterministic across
    // decomposition choices.
    const int src_top_k = top_cell_k(*states[0]);
    const Real sample = states[0]->min(0);

    if (ParallelDescriptor::IOProcessor()) {
        Print() << "REMORA::ApplyAtmosphericStates time=" << time
                << " sample_Uwind=" << sample
                << " (src_k=" << src_top_k << ")\n";
    }
}
