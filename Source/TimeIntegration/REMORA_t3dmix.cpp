#include <REMORA.H>

using namespace amrex;

void
REMORA::t3dmix  (const Box& bx,
                const Array4<Real      >& state,
                const Array4<Real      >& state_rhs,
                const Array4<Real const>& diff2,
                const Array4<Real const>& Hz,
                const Array4<Real const>& pm,
                const Array4<Real const>& pn,
                const Real dt_lev, const int ncomp)
{
    //-----------------------------------------------------------------------
    //  Add in harmonic diffusivity s terms.
    //-----------------------------------------------------------------------

    Box xbx(bx); xbx.surroundingNodes(0);
    Box ybx(bx); ybx.surroundingNodes(1);

    FArrayBox fab_FX(xbx,ncomp,The_Async_Arena());
    FArrayBox fab_FE(ybx,ncomp,The_Async_Arena());

    auto FX=fab_FX.array();
    auto FE=fab_FE.array();

    ParallelFor(xbx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
    {
        const Real pmon_u = (pm(i-1,j,0)+pm(i,j,0))/(pn(i-1,j,0)+pn(i,j,0));

        const Real cff = 0.25_rt * (diff2(i,j,n) + diff2(i-1,j,n)) * pmon_u;
        FX(i,j,k,n) = cff * (Hz(i,j,k) + Hz(i-1,j,k)) * (state_rhs(i,j,k,n)-state_rhs(i-1,j,k,n));
    });

    ParallelFor(ybx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
    {
        const Real pnom_v = (pn(i,j-1,0)+pn(i,j,0))/(pm(i,j-1,0)+pm(i,j,0));

        const Real cff = 0.25_rt*(diff2(i,j,n)+diff2(i,j-1,n)) * pnom_v;
        FE(i,j,k,n) = cff * (Hz(i,j,k) + Hz(i,j-1,k)) * (state_rhs(i,j,k,n) - state_rhs(i,j-1,k,n));
    });

    /*
     Time-step harmonic, S-surfaces diffusion term.
    */
    ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
    {
        const Real cff = dt_lev*pm(i,j,0)*pn(i,j,0);

        state(i,j,k,n) += cff * ( (FX(i+1,j  ,k,n)-FX(i,j,k,n))
                                 +(FE(i  ,j+1,k,n)-FE(i,j,k,n)) );
    });
}
