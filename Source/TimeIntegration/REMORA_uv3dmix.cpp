#include <REMORA.H>

using namespace amrex;

void
REMORA::uv3dmix  (const Box& xbx, const Box& ybx,
                 const Array4<Real      >& u,
                 const Array4<Real      >& v,
                 const Array4<Real const>& uold,
                 const Array4<Real const>& vold,
                 const Array4<Real      >& rufrc,
                 const Array4<Real      >& rvfrc,
                 const Array4<Real const>& visc2_p,
                 const Array4<Real const>& visc2_r,
                 const Array4<Real const>& Hz,
                 const Array4<Real const>& pm,
                 const Array4<Real const>& pn,
                 const Array4<Real const>& mskp,
                 int nrhs, int nnew,
                 const Real dt_lev)
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

    auto N = xbx.hiVect()[2] - ybx.loVect()[2];

    // Think of the cell-centered box as                              [ 0:nx-1, 0:ny-1] (0,0,0) cc
    //
    // xbx is the x-face-centered box on which we update u  [ 0:nx  , 0:ny-1] (1,0,0) x-faces
    //   to do so requires                        UFx on    [-1:nx  , 0:ny-1] (0,0,0) cc
    //      which requires                        uold on   [-1:nx+1, 0:ny-1] (1,0,0) x-faces
    //        and requires                        vold on   [-1:nx  , 0:ny  ] (0,1,0) y-faces
    //   to do so requires                        UFe on    [ 0:nx  , 0:ny  ] (1,1,0) xy-nodes
    //      which requires                        uold on   [ 0:nx  ,-1:ny  ] (1,0,0) x-faces
    //       and  requires                        vold on   [ -1:nx , 0:ny  ] (0,1,0) y-faces

    ParallelFor(growLo(xbx,0,1), [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        const Real cff = 0.5_rt*Hz(i,j,k) * ( pm(i,j,0) / pn(i,j,0) *
                                    ( (pn(i  ,j,0) + pn(i+1,j,0)) * uold(i+1,j,k,nrhs)-
                                      (pn(i-1,j,0) + pn(i  ,j,0)) * uold(i  ,j,k,nrhs) )-
                                                pn(i,j,0) / pm(i,j,0) *
                                    ( (pm(i,j  ,0) + pm(i,j+1,0)) * vold(i,j+1,k,nrhs)-
                                      (pm(i,j-1,0) + pm(i,j  ,0)) * vold(i,j  ,k,nrhs) ) );

        Real on_r = 1.0_rt / pn(i,j,0);
        UFx(i,j,k) = on_r * on_r * visc2_r(i,j,0) * cff;
    });

    ParallelFor(growHi(xbx,1,1), [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        const Real pmon_p = (pm(i-1,j-1,0)+pm(i-1,j,0)+pm(i,j-1,0)+pm(i,j,0)) /
                            (pn(i-1,j-1,0)+pn(i-1,j,0)+pn(i,j-1,0)+pn(i,j,0));
        const Real pnom_p = (pn(i-1,j-1,0)+pn(i-1,j,0)+pn(i,j-1,0)+pn(i,j,0)) /
                            (pm(i-1,j-1,0)+pm(i-1,j,0)+pm(i,j-1,0)+pm(i,j,0));
        const Real cff = mskp(i,j,0) * 0.125_rt *
                    (Hz(i-1,j  ,k) + Hz(i,j ,k)+ Hz(i-1,j-1,k) + Hz(i,j-1,k))*
                    (pmon_p*
                     ((pn(i  ,j-1,0)+pn(i  ,j,0))*vold(i  ,j,k,nrhs)-
                      (pn(i-1,j-1,0)+pn(i-1,j,0))*vold(i-1,j,k,nrhs))+
                     pnom_p*
                     ((pm(i-1,j  ,0)+pm(i,j  ,0))*uold(i,j  ,k,nrhs)-
                      (pm(i-1,j-1,0)+pm(i,j-1,0))*uold(i,j-1,k,nrhs)));

        const Real om_p =  4.0_rt / (pm(i-1,j-1,0)+pm(i-1,j,0)+pm(i,j-1,0)+pm(i,j,0));
        UFe(i,j,k) = om_p*om_p*visc2_p(i,j,0)*cff;
    });

    ParallelFor(makeSlab(xbx,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int)
    {
        for (int k=0; k<=N; k++) {
            const Real cff=dt_lev*0.25_rt*(pm(i-1,j,0)+pm(i,j,0))*(pn(i-1,j,0)+pn(i,j,0));
            const Real cff1=0.5_rt*(pn(i-1,j,0)+pn(i,j,0))*(UFx(i,j  ,k)-UFx(i-1,j,k));
            const Real cff2=0.5_rt*(pm(i-1,j,0)+pm(i,j,0))*(UFe(i,j+1,k)-UFe(i  ,j,k));
            const Real cff3=cff*(cff1+cff2);
            u(i,j,k,nnew)=u(i,j,k,nnew)+cff3;
            rufrc(i,j,0) += cff1+cff2;
        }
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

    ParallelFor(growLo(ybx,1,1), [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        // cff depends on k, but UFx and VFe will only be affected by the last cell?
        const Real cff = 0.5_rt*Hz(i,j,k) * (pm(i,j,0) / pn(i,j,0) *
                ((pn(i,  j,0) + pn(i+1,j,0)) * uold(i+1,j,k,nrhs)-
                 (pn(i-1,j,0) + pn(i,  j,0)) * uold(i  ,j,k,nrhs))-
                pn(i,j,0) / pm(i,j,0) *
                ((pm(i,j  ,0)+pm(i,j+1,0))*vold(i,j+1,k,nrhs)-
                 (pm(i,j-1,0)+pm(i,j  ,0))*vold(i,j  ,k,nrhs)));

        Real om_r = 1.0_rt / pm(i,j,0);
        VFe(i,j,k) = om_r * om_r * visc2_r(i,j,0) * cff;
    });

    ParallelFor(growHi(ybx,0,1), [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        const Real pmon_p = (pm(i-1,j-1,0)+pm(i-1,j,0)+pm(i,j-1,0)+pm(i,j,0)) /
                            (pn(i-1,j-1,0)+pn(i-1,j,0)+pn(i,j-1,0)+pn(i,j,0));
        const Real pnom_p = (pn(i-1,j-1,0)+pn(i-1,j,0)+pn(i,j-1,0)+pn(i,j,0)) /
                            (pm(i-1,j-1,0)+pm(i-1,j,0)+pm(i,j-1,0)+pm(i,j,0));
        const Real cff = mskp(i,j,0) * 0.125_rt * (Hz(i-1,j  ,k)+Hz(i,j ,k)+
                              Hz(i-1,j-1,k)+Hz(i,j-1,k))*
                    (pmon_p*
                     ((pn(i  ,j-1,0)+pn(i  ,j,0))*vold(i  ,j,k,nrhs)-
                      (pn(i-1,j-1,0)+pn(i-1,j,0))*vold(i-1,j,k,nrhs))+
                     pnom_p*
                     ((pm(i-1,j  ,0)+pm(i,j  ,0))*uold(i,j  ,k,nrhs)-
                      (pm(i-1,j-1,0)+pm(i,j-1,0))*uold(i,j-1,k,nrhs)));

        const Real on_p =  4.0_rt / (pn(i-1,j-1,0)+pn(i-1,j,0)+pn(i,j-1,0)+pn(i,j,0));
        VFx(i,j,k) = on_p*on_p*visc2_p(i,j,0)*cff;
    });

    ParallelFor(makeSlab(ybx,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int )
    {
        for (int k=0; k<=N; k++) {
            const Real cff=dt_lev*0.25_rt*(pm(i,j,0)+pm(i,j-1,0))*(pn(i,j,0)+pn(i,j-1,0));
            const Real cff1=0.5_rt*(pn(i,j-1,0)+pn(i,j,0))*(VFx(i+1,j,k)-VFx(i,j  ,k));
            const Real cff2=0.5_rt*(pm(i,j-1,0)+pm(i,j,0))*(VFe(i  ,j,k)-VFe(i,j-1,k));
            const Real cff3=cff*(cff1-cff2);

            v(i,j,k,nnew)=v(i,j,k,nnew)+cff3;
            rvfrc(i,j,0) += cff1-cff2;
        }
    });

}
