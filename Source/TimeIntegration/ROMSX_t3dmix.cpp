#include <ROMSX.H>
#include <Utils.H>

using namespace amrex;

void
ROMSX::t3dmix  (const Box& bx,
                 Array4<Real> t,
                 Array4<Real> diff2, Array4<Real> Hz,
                 Array4<Real> pm, Array4<Real> pn,
                 Array4<Real> pmon_u, Array4<Real> pnom_v,
                 int nrhs, int nnew,
                 const amrex::Real dt_lev)
{
    //-----------------------------------------------------------------------
    //  Add in harmonic diffusivity s terms.
    //-----------------------------------------------------------------------

    Box gbx2 = bx;
    gbx2.grow(IntVect(NGROW,NGROW,0));

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
    amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                const amrex::Real cff = 0.25 * (diff2(i,j,0) + diff2(i-1,j,0)) * pmon_u(i,j,0);
                FX(i,j,k) = cff * (Hz(i,j,k)+Hz(i+1,j,k))*(t(i,j,k,nrhs)-t(i-1,j,k,nrhs));
            });

    //K_LOOP : DO k=1,N(ng)
    //DO j=Jstr,Jend+1
    //DO i=Istr,Iend
    amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                const amrex::Real cff=0.25*(diff2(i,j,0)+diff2(i,j-1,0)) * pnom_v(i,j,0);
                FE(i,j,k) = cff * (Hz(i,j,k) + Hz(i,j-1,k)) * (t(i,j,k,nrhs) - t(i,j-1,k,nrhs));
            });
    /*
     Time-step harmonic, S-surfaces diffusion term.
    */
    //K_LOOP : DO k=1,N(ng)
    //DO j=Jstr,Jend
    //DO i=IstrU,Iend
    amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                const amrex::Real cff=dt_lev*pm(i,j,0)*pn(i,j,0);
                const amrex::Real cff1=cff*(FX(i+1,j  ,k)-FX(i,j,k));
                const amrex::Real cff2=cff*(FE(i  ,j+1,k)-FE(i,j,k));
                const amrex::Real cff3=cff1+cff2;
                t(i,j,k,nnew)=t(i,j,k,nnew)+cff3;
            });

}
