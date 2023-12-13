#include <ROMSX.H>

using namespace amrex;

void
ROMSX::t3dmix  (const Box& bx,
                Array4<Real> state,
                Array4<Real> diff2, Array4<Real> Hz,
                Array4<Real> pm, Array4<Real> pn,
                Array4<Real> pmon_u, Array4<Real> pnom_v,
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
            const Real cff = 0.25 * (diff2(i,j,n) + diff2(i-1,j,n)) * pmon_u(i,j,0);
            FX(i,j,k,n) = cff * (Hz(i,j,k) + Hz(i+1,j,k)) * (state(i,j,k,n)-state(i-1,j,k,n));
        });

    ParallelFor(ybx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
        {
            const Real cff = 0.25*(diff2(i,j,n)+diff2(i,j-1,n)) * pnom_v(i,j,0);
            FE(i,j,k,n) = cff * (Hz(i,j,k) + Hz(i,j-1,k)) * (state(i,j,k,n) - state(i,j-1,k,n));
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
