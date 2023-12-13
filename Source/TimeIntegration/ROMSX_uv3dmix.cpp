#include <ROMSX.H>

using namespace amrex;

void
ROMSX::uv3dmix  (const Box& xbx, const Box& ybx,
                 Array4<Real> u  , Array4<Real> v,
                 Array4<Real> uold  , Array4<Real> vold,
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

    FArrayBox fab_UFx(growLo(xbx,0,1),1,amrex::The_Async_Arena()); fab_UFx.template setVal<RunOn::Device>(0.);
    FArrayBox fab_UFe(growHi(xbx,1,1),1,amrex::The_Async_Arena()); fab_UFe.template setVal<RunOn::Device>(0.);
    FArrayBox fab_VFe(growLo(ybx,1,1),1,amrex::The_Async_Arena()); fab_VFe.template setVal<RunOn::Device>(0.);
    FArrayBox fab_VFx(growHi(ybx,0,1),1,amrex::The_Async_Arena()); fab_VFx.template setVal<RunOn::Device>(0.);

    auto UFx=fab_UFx.array();
    auto UFe=fab_UFe.array();
    auto VFx=fab_VFx.array();
    auto VFe=fab_VFe.array();

    // Think of the cell-centered box as                              [ 0:nx-1, 0:ny-1] (0,0,0) cc
    //
    // xbx is the x-face-centered box on which we update u  [ 0:nx  , 0:ny-1] (1,0,0) x-faces
    //   to do so requires                        UFx on    [-1:nx  , 0:ny-1] (0,0,0) cc
    //      which requires                        uold on   [-1:nx+1, 0:ny-1] (1,0,0) x-faces
    //        and requires                        vold on   [-1:nx  , 0:ny  ] (0,1,0) y-faces
    //   to do so requires                        UFe on    [ 0:nx  , 0:ny  ] (1,1,0) xy-nodes
    //      which requires                        uold on   [ 0:nx  ,-1:ny  ] (1,0,0) x-faces
    //       and  requires                        vold on   [ -1:nx , 0:ny  ] (0,1,0) y-faces
    ParallelFor(growLo(xbx,0,1),
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        const amrex::Real cff = 0.5*Hz(i,j,k) * (pm(i,j,0) / pn(i,j,0) *
                                    ((pn(i,  j,0) + pn(i+1,j,0)) * uold(i+1,j,k,nrhs)-
                                    (pn(i-1,j,0) + pn(i,  j,0)) * uold(i  ,j,k,nrhs))-
                                pn(i,j,0) / pm(i,j,0) *
                                    ((pm(i,j  ,0)+pm(i,j+1,0))*vold(i,j+1,k,nrhs)-
                                    (pm(i,j-1,0)+pm(i,j  ,0))*vold(i,j  ,k,nrhs)));
        UFx(i,j,k) = on_r(i,j,0)*on_r(i,j,0)*visc2_r(i,j,0)*cff;
    });
    ParallelFor(growHi(xbx,1,1),
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        const amrex::Real cff = 0.125 * (Hz(i-1,j  ,k)+Hz(i,j ,k)+
                              Hz(i-1,j-1,k)+Hz(i,j-1,k))*
                    (pm(i,j,0)/pn(i,j,0)*
                     ((pn(i  ,j-1,0)+pn(i  ,j,0))*vold(i  ,j,k,nrhs)-
                      (pn(i-1,j-1,0)+pn(i-1,j,0))*vold(i-1,j,k,nrhs))+
                     pn(i,j,0)/pm(i,j,0)*
                     ((pm(i-1,j  ,0)+pm(i,j  ,0))*uold(i,j  ,k,nrhs)-
                      (pm(i-1,j-1,0)+pm(i,j-1,0))*uold(i,j-1,k,nrhs)));
        UFe(i,j,k) = om_p(i,j,0)*om_p(i,j,0)*visc2_p(i,j,0)*cff;
    });

    ParallelFor(xbx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        const amrex::Real cff=dt_lev*0.25*(pm(i-1,j,0)+pm(i,j,0))*(pn(i-1,j,0)+pn(i,j,0));
        const amrex::Real cff1=0.5*(pn(i-1,j,0)+pn(i,j,0))*(UFx(i,j  ,k)-UFx(i-1,j,k));
        const amrex::Real cff2=0.5*(pm(i-1,j,0)+pm(i,j,0))*(UFe(i,j+1,k)-UFe(i  ,j,k));
        const amrex::Real cff3=cff*(cff1+cff2);
        amrex::Gpu::Atomic::Add(&(rufrc(i,j,0)), cff1+cff2);
        u(i,j,k,nnew)=u(i,j,k,nnew)+cff3;
    });


    // Think of the cell-centered box as                              [ 0:nx-1, 0:ny-1] (0,0,0) cc
    //
    // ybx is the y-face-centered box on which we update v  [ 0:nx-1, 0:ny  ] (1,0,0) x-faces
    //   to do so requires                        VFe on    [ 0:nx-1,-1:ny  ] (0,0,0) cc
    //      which requires                        uold on   [ 0:nx  ,-1:ny  ] (1,0,0) x-faces
    //        and requires                        vold on   [ 0:nx-1,-1:ny+1] (0,1,0) y-faces
    //   to do so requires                        VFx on    [ 0:nx  , 0:ny  ] (1,1,0) xy-nodes
    //      which requires                        uold on   [ 0:nx  ,-1:ny  ] (1,0,0) x-faces
    //       and  requires                        vold on   [-1:nx  , 0:ny  ] (0,1,0) y-faces
    ParallelFor(growLo(ybx,1,1),
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        // cff depends on k, but UFx and VFe will only be affected by the last cell?
        const amrex::Real cff = 0.5*Hz(i,j,k) * (pm(i,j,0) / pn(i,j,0) *
                ((pn(i,  j,0) + pn(i+1,j,0)) * uold(i+1,j,k,nrhs)-
                 (pn(i-1,j,0) + pn(i,  j,0)) * uold(i  ,j,k,nrhs))-
                pn(i,j,0) / pm(i,j,0) *
                ((pm(i,j  ,0)+pm(i,j+1,0))*vold(i,j+1,k,nrhs)-
                 (pm(i,j-1,0)+pm(i,j  ,0))*vold(i,j  ,k,nrhs)));
        VFe(i,j,k) = om_r(i,j,0)*om_r(i,j,0)*visc2_r(i,j,0)*cff;
    });
    ParallelFor(growHi(ybx,0,1),
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        const amrex::Real cff = 0.125 * (Hz(i-1,j  ,k)+Hz(i,j ,k)+
                              Hz(i-1,j-1,k)+Hz(i,j-1,k))*
                    (pm(i,j,0)/pn(i,j,0)*
                     ((pn(i  ,j-1,0)+pn(i  ,j,0))*vold(i  ,j,k,nrhs)-
                      (pn(i-1,j-1,0)+pn(i-1,j,0))*vold(i-1,j,k,nrhs))+
                     pn(i,j,0)/pm(i,j,0)*
                     ((pm(i-1,j  ,0)+pm(i,j  ,0))*uold(i,j  ,k,nrhs)-
                      (pm(i-1,j-1,0)+pm(i,j-1,0))*uold(i,j-1,k,nrhs)));
        VFx(i,j,k) = on_p(i,j,0)*on_p(i,j,0)*visc2_p(i,j,0)*cff;
    });

    ParallelFor(ybx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        const amrex::Real cff=dt_lev*0.25*(pm(i,j,0)+pm(i,j-1,0))*(pn(i,j,0)+pn(i,j-1,0));
        const amrex::Real cff1=0.5*(pn(i,j-1,0)+pn(i,j,0))*(VFx(i+1,j,k)-VFx(i,j  ,k));
        const amrex::Real cff2=0.5*(pm(i,j-1,0)+pm(i,j,0))*(VFe(i  ,j,k)-VFe(i,j-1,k));
        const amrex::Real cff3=cff*(cff1-cff2);
        amrex::Gpu::Atomic::Add(&(rvfrc(i,j,0)), cff1-cff2);
        v(i,j,k,nnew)=v(i,j,k,nnew)+cff3;
    });

}
