#include <ROMSX.H>
#include <Utils.H>

using namespace amrex;

void
ROMSX::t3dmix  (const Box& bx,
                 Array4<Real> t,
                 Array4<Real> visc2, Array4<Real> Hz,
                 Array4<Real> pm, Array4<Real> pn,
                 Array4<Real> pmon_u, Array4<Real> pnom_v,
                 int nrhs, int nnew,
                 const amrex::Real dt_lev)
{
    //-----------------------------------------------------------------------
    //  Add in harmonic viscosity s terms.
    //-----------------------------------------------------------------------
    // TODO: Fix comments to match tracer, also fix box dims.
    // And also make sure visc2 set
    // And figure out where this goes in the main flow of the code

    Box gbx2 = bx;
    gbx2.grow(IntVect(2,2,0));
    int ncomp = 1;

    FArrayBox fab_FX(gbx2,1,amrex::The_Async_Arena());
    FArrayBox fab_FE(gbx2,1,amrex::The_Async_Arena());

    auto FX=fab_FX.array();
    auto FE=fab_FE.array();

    amrex::ParallelFor(gbx2,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        FX(i,j,k)=0.0;
        FE(i,j,k)=0.0;
    });

    //K_LOOP : DO k=1,N(ng)
    //DO j=Jstr,Jend
    //DO i=Istr,Iend+1
    amrex::ParallelFor(bx, ncomp,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                const amrex::Real cff = 0.25 * (diff2(i,j,0) + diff2(i-1,j,0)) * pmon_u(i,j);
                FX(i,j,k) = cff * (Hz(i,j,k)+Hz(i+1,j,k))*(t(i,j,k,nrhs)-t(i-1,j,k,nrhs));
            });

    //K_LOOP : DO k=1,N(ng)
    //DO j=Jstr,Jend+1
    //DO i=Istr,Iend
    amrex::ParallelFor(bx, ncomp,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                const amrex::Real cff=0.25*(diff2(i,j,0)+diff2(i,j-1,0)) * pnom_v(i,j);
                FE(i,j) = cff * (Hz(i,j,k) + Hz(i,j-1,k)) * (t(i,j,k,nrhs) - t(i,j-1,k,nrhs));
            });
/*
 Time-step harmonic, S-surfaces diffusion term. Notice that momentum
 at this stage is HzU and HzV and has m2/s units. Add contribution for
 barotropic forcing terms.
 */
    //K_LOOP : DO k=1,N(ng)
    //DO j=Jstr,Jend
    //DO i=IstrU,Iend
    amrex::ParallelFor(bx, ncomp,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                const amrex::Real cff=dt_lev*pm(i,j)*pn(i,j);
                const amrex::Real cff1=cff*(FX(i+1,j  )-FX(i,j));
                const amrex::Real cff2=cff*(FE(i  ,j+1)-FE(i,j));
                const amrex::Real cff3=cff1+cff2;
                const amrex::Real t(i,j,k,nnew)=t(i,j,k,nnew)+cff3;
            });

}
