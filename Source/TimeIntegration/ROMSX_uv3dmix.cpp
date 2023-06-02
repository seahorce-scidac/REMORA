#include <ROMSX.H>
#include <Utils.H>

using namespace amrex;

void
ROMSX::uv3dmix  (const Box& bx,
                 Array4<Real> u  , Array4<Real> v,
                 Array4<Real> rufrc, Array4<Real> rvfrc,
                 Array4<Real> visc2_p,
                 Array4<Real> visc2_r,
                 Array4<Real> Hz,
                 Array4<Real> om_r, Array4<Real> on_r,
                 Array4<Real> om_p, Array4<Real> on_p,
                 Array4<Real> pm, Array4<Real> pn,
                 int nrhs, int nnew,
                 const amrex::Real dt_lev)
{
    //-----------------------------------------------------------------------
    //  Add in harmonic viscosity s terms.
    //-----------------------------------------------------------------------

    Box gbx2 = bx;
    gbx2.grow(IntVect(NGROW,NGROW,0));
    int ncomp = 1;

    FArrayBox fab_UFx(gbx2,1,amrex::The_Async_Arena());
    FArrayBox fab_UFe(gbx2,1,amrex::The_Async_Arena());
    FArrayBox fab_VFx(gbx2,1,amrex::The_Async_Arena());
    FArrayBox fab_VFe(gbx2,1,amrex::The_Async_Arena());

    auto UFx=fab_UFx.array();
    auto UFe=fab_UFe.array();
    auto VFx=fab_VFx.array();
    auto VFe=fab_VFe.array();

    amrex::ParallelFor(gbx2,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        UFx(i,j,k)=0.0;
        UFe(i,j,k)=0.0;
        VFx(i,j,k)=0.0;
        VFe(i,j,k)=0.0;
    });

    //K_LOOP : DO k=1,N(ng)
    //DO j=JstrV-1,Jend
    //DO i=IstrU-1,Iend
    // Should these be on ubox/vbox?
    amrex::ParallelFor(bx, ncomp,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                // cff depends on k, but UFx and VFe will only be affected by the last cell?
                const amrex::Real cff = 0.5*Hz(i,j,k) * (pm(i,j,0) / pn(i,j,0) *
                        ((pn(i,  j,0) + pn(i+1,j,0)) * u(i+1,j,k,nrhs)-
                         (pn(i-1,j,0) + pn(i,  j,0)) * u(i  ,j,k,nrhs))-
                        pn(i,j,0) / pm(i,j,0) *
                        ((pm(i,j  ,0)+pm(i,j+1,0))*v(i,j+1,k,nrhs)-
                         (pm(i,j-1,0)+pm(i,j  ,0))*v(i,j  ,k,nrhs)));
                // ifndef VISC_3DCOEF branch
                UFx(i,j,k) = on_r(i,j,0)*on_r(i,j,0)*visc2_r(i,j,0)*cff;
                VFe(i,j,k) = om_r(i,j,0)*om_r(i,j,0)*visc2_r(i,j,0)*cff;
            });

    //K_LOOP : DO k=1,N(ng)
    //DO j=Jstr,Jend+1
    //DO i=Istr,Iend+1
    amrex::ParallelFor(bx, ncomp,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                const amrex::Real cff = 0.125 * (Hz(i-1,j  ,k)+Hz(i,j  ,k)+
                                      Hz(i-1,j-1,k)+Hz(i,j-1,k))*
                            (pm(i,j,0)/pn(i,j,0)*
                             ((pn(i  ,j-1,0)+pn(i  ,j,0))*v(i  ,j,k,nrhs)-
                              (pn(i-1,j-1,0)+pn(i-1,j,0))*v(i-1,j,k,nrhs)));
                // ifndef VISC_3DCOEF branch
                UFe(i,j,k) = om_p(i,j,0)*om_p(i,j,0)*visc2_p(i,j,0)*cff;
                VFx(i,j,k) = on_p(i,j,0)*on_p(i,j,0)*visc2_p(i,j,0)*cff;
            });
/*
 Time-step harmonic, S-surfaces viscosity term. Notice that momentum
 at this stage is HzU and HzV and has m2/s units. Add contribution for
 barotropic forcing terms.
 */
    //K_LOOP : DO k=1,N(ng)
    //DO j=Jstr,Jend
    //DO i=IstrU,Iend
           amrex::PrintToFile("rufrc").SetPrecision(18)<<FArrayBox(rufrc)<<std::endl;
       amrex::PrintToFile("rvfrc").SetPrecision(18)<<FArrayBox(rvfrc)<<std::endl;
    amrex::ParallelFor(bx, ncomp,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                const amrex::Real cff=dt_lev*0.25*(pm(i-1,j,0)+pm(i,j,0))*(pn(i-1,j,0)+pn(i,j,0));
                const amrex::Real cff1=0.5*(pn(i-1,j,0)+pn(i,j,0))*(UFx(i,j  ,k)-UFx(i-1,j,k));
                const amrex::Real cff2=0.5*(pm(i-1,j,0)+pm(i,j,0))*(UFe(i,j+1,k)-UFe(i  ,j,k));
                const amrex::Real cff3=cff*(cff1+cff2);
                rufrc(i,j,0)=rufrc(i,j,0)+cff1+cff2;
                u(i,j,k,nnew)=u(i,j,k,nnew)+cff3;
/*#ifdef DIAGNOSTICS_UV
                DiaRUfrc(i,j,3,M2hvis)=DiaRUfrc(i,j,3,M2hvis)+cff1+cff2
                DiaRUfrc(i,j,3,M2xvis)=DiaRUfrc(i,j,3,M2xvis)+cff1
                DiaRUfrc(i,j,3,M2yvis)=DiaRUfrc(i,j,3,M2yvis)+cff2
                DiaU3wrk(i,j,k,M3hvis)=cff3
                DiaU3wrk(i,j,k,M3xvis)=cff*cff1
                DiaU3wrk(i,j,k,M3yvis)=cff*cff2
#endif*/
            });

    //K_LOOP : DO k=1,N(ng)
    //DO j=JstrV,Jend
    //DO i=Istr,Iend
    amrex::ParallelFor(bx, ncomp,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                const amrex::Real cff=dt_lev*0.25*(pm(i,j,0)+pm(i,j-1,0))*(pn(i,j,0)+pn(i,j-1,0));
                const amrex::Real cff1=0.5*(pn(i,j-1,0)+pn(i,j,0))*(VFx(i+1,j,k)-VFx(i,j  ,k));
                const amrex::Real cff2=0.5*(pm(i,j-1,0)+pm(i,j,0))*(VFe(i  ,j,k)-VFe(i,j-1,k));
                const amrex::Real cff3=cff*(cff1-cff2);
                rvfrc(i,j,0)=rvfrc(i,j,0)+cff1-cff2;
                v(i,j,k,nnew)=v(i,j,k,nnew)+cff3;
/*#ifdef DIAGNOSTICS_UV
                DiaRVfrc(i,j,3,M2hvis)=DiaRVfrc(i,j,3,M2hvis)+cff1-cff2
                DiaRVfrc(i,j,3,M2xvis)=DiaRVfrc(i,j,3,M2xvis)+cff1
                DiaRVfrc(i,j,3,M2yvis)=DiaRVfrc(i,j,3,M2yvis)-cff2
                DiaV3wrk(i,j,k,M3hvis)=cff3
                DiaV3wrk(i,j,k,M3xvis)= cff*cff1
                DiaV3wrk(i,j,k,M3yvis)=-cff*cff2
#endif*/
            });


}
