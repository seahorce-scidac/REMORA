#include <REMORA.H>

using namespace amrex;

/**
 * rhs_uv_2d
 *
 * @param[in   ] xbx Box for operations on x-velocity
 * @param[in   ] ybx Box for operations on y-velocity
 * @param[in   ] ubar
 * @param[in   ] vbar
 * @param[  out] rhs_ubar
 * @param[  out] rhs_vbar
 * @param[in   ] DUon
 * @param[in   ] DVom
 * @param[in   ] krhs
 */

void
REMORA::rhs_uv_2d (const Box& xbx, const Box& ybx,
                  const Array4<Real const>& ubar,
                  const Array4<Real const>& vbar,
                  const Array4<Real      >& rhs_ubar  ,
                  const Array4<Real      >& rhs_vbar,
                  const Array4<Real const>& DUon,
                  const Array4<Real const>& DVom,
                  const int krhs)
{
    //
    // Scratch space
    //
    FArrayBox fab_UFx(growLo(xbx,0,1),1,amrex::The_Async_Arena()); fab_UFx.template setVal<RunOn::Device>(0.);
    FArrayBox fab_UFe(growHi(xbx,1,1),1,amrex::The_Async_Arena()); fab_UFe.template setVal<RunOn::Device>(0.);
    FArrayBox fab_VFe(growLo(ybx,1,1),1,amrex::The_Async_Arena()); fab_VFe.template setVal<RunOn::Device>(0.);
    FArrayBox fab_VFx(growHi(ybx,0,1),1,amrex::The_Async_Arena()); fab_VFx.template setVal<RunOn::Device>(0.);

    auto UFx=fab_UFx.array();
    auto UFe=fab_UFe.array();
    auto VFx=fab_VFx.array();
    auto VFe=fab_VFe.array();

    // *************************************************************
    // UPDATING U
    // *************************************************************

    // Think of the cell-centered box as                              [ 0:nx-1, 0:ny-1] (0,0,0) cc
    //
    // xbx is the x-face-centered box on which we update ubar (with rhs_ubar)  [ 0:nx  , 0:ny-1] (1,0,0) x-faces
    //   to do so requires                                           UFx on    [-1:nx  , 0:ny-1] (0,0,0) cc
    //      which requires                                           ubar on   [-2:nx+2, 0:ny-1] (1,0,0) x-faces
    //       and  requires                                           DUon on   [-2:nx+2, 0:ny-1] (1,0,0) x-faces
    //   to do so requires                                           UFe on    [ 0:nx  , 0:ny  ] (1,1,0) xy-nodes
    //      which requires                                           ubar on   [ 0:nx  ,-2:ny+1] (1,0,0) x-faces
    //       and  requires                                           DVom on   [-2:nx+1, 0:ny-1] (0,1,0) y-faces

    //
    // Define UFx, the x-fluxes at cell centers for updating u
    // (Note that grow arguments are (bx, dir, ng)
    //
    ParallelFor(growLo(xbx,0,1),
    [=] AMREX_GPU_DEVICE (int i, int j, int )
    {
        Real uxx_i   = ubar(i-1,j,0,krhs)-2.0*ubar(i  ,j,0,krhs)+ubar(i+1,j,0,krhs);
        Real uxx_ip1 = ubar(i  ,j,0,krhs)-2.0*ubar(i+1,j,0,krhs)+ubar(i+2,j,0,krhs);
        Real uxx_avg = uxx_i + uxx_ip1;

        Real Huxx_i   = DUon(i-1,j,0)-2.0*DUon(i  ,j,0)+DUon(i+1,j,0);
        Real Huxx_ip1 = DUon(i  ,j,0)-2.0*DUon(i+1,j,0)+DUon(i+2,j,0);

        Real cff=1.0/6.0;
        Real ubar_avg = ubar(i  ,j,0,krhs)+ubar(i+1,j,0,krhs);

        UFx(i,j,0)=0.25*(ubar_avg-cff*uxx_avg) * (DUon(i,j,0)+ DUon(i+1,j,0)-cff*(Huxx_i+ Huxx_ip1));
    });

    //
    // Define UFe, the y-fluxes at nodes for updating u
    //
    ParallelFor(growHi(xbx,1,1),
    [=] AMREX_GPU_DEVICE (int i, int j, int )
    {
        //should not include grow cells
        Real uee_j   = ubar(i,j-1,0,krhs)-2.0*ubar(i,j  ,0,krhs)+ubar(i,j+1,0,krhs);
        Real uee_jm1 = ubar(i,j-2,0,krhs)-2.0*ubar(i,j-1,0,krhs)+ubar(i,j  ,0,krhs);
        Real uee_avg = uee_j + uee_jm1;

        Real Hvxx_i   = DVom(i-1,j,0)-2.0*DVom(i  ,j,0)+DVom(i+1,j,0);
        Real Hvxx_im1 = DVom(i-2,j,0)-2.0*DVom(i-1,j,0)+DVom(i  ,j,0);
        Real cff=1.0/6.0;
        Real cff1=ubar(i,j  ,0,krhs)+ubar(i,j-1,0,krhs);
        Real cff2=DVom(i,j,0)+DVom(i-1,j,0);

        UFe(i,j,0)=0.25*(cff1-uee_avg*cff)*
          (cff2-cff*(Hvxx_i+Hvxx_im1));
    });

    //
    //  Add in horizontal advection.
    //
    ParallelFor(makeSlab(xbx,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int)
    {
         Real cff1=UFx(i,j  ,0)-UFx(i-1,j,0);
         Real cff2=UFe(i,j+1,0)-UFe(i  ,j,0);

         rhs_ubar(i,j,0) -= (cff1 + cff2);
    });

    // *************************************************************
    // UPDATING V
    // *************************************************************

    // Think of the cell-centered box as                              [ 0:nx-1, 0:ny-1] (0,0,0) cc
    //
    // ybx is the y-face-centered box on which we update vbar (with rhs_vbar)  [ 0:nx-1, 0:ny  ] (1,0,0)  y-faces
    //   to do so requires                                  VFe on    [ 0:nx-1,-1:ny  ] (1,1,0)  xy-nodes
    //      which requires                                    vbar on [ 0:nx-1,-2:ny+2] (1,0,0)  y-faces
    //       and  requires                                    DVom on [ 0:nx-1,-2:ny+2] (0,1,0)  x-faces
    //   to do so requires                                  VFx on    [ 0:nx  , 0:ny  ] (0,0,0)  cc
    //      which requires                                    vbar on [-2:nx+1, 0:ny  ] (1,0,0)  y-faces
    //       and  requires                                    DUon on [ 0:nx-1,-2:ny+1] (0,0,0)  y-faces

    // Grow ybx by one in high x-direction
    ParallelFor(growHi(ybx,0,1),
    [=] AMREX_GPU_DEVICE (int i, int j, int )
    {
        Real cff=1.0/6.0;
        Real vxx_i   = vbar(i-1,j,0,krhs)-2.0*vbar(i  ,j,0,krhs)+vbar(i+1,j,0,krhs);
        Real vxx_im1 = vbar(i-2,j,0,krhs)-2.0*vbar(i-1,j,0,krhs)+vbar(i  ,j,0,krhs);
        Real vxx_avg = vxx_i + vxx_im1;
        //auto vxx_im1 = (i == gbx1.smallEnd(0)) ? vxx(i-1,j,k) :
        //    (vbar(i-2,j,0,krhs)-2.0*vbar(i-1,j,0,krhs)+vbar(i,j,0,krhs));
        //neglecting terms about periodicity since testing only periodic for now
        Real Huee_j   = DUon(i,j-1,0)-2.0*DUon(i,j  ,0)+DUon(i,j+1,0);
        Real Huee_jm1 = DUon(i,j-2,0)-2.0*DUon(i,j-1,0)+DUon(i,j  ,0);
        Real cff1=vbar(i  ,j,0,krhs)+vbar(i-1,j,0,krhs);
        Real cff2=DUon(i,j,0)+DUon(i,j-1,0);

        //auto Huee_jm1 = (j == gbx1.smallEnd(1)) ? Huee(i,j-1,k) :
        //    (DUon(i,j-2,k)-2.0*DUon(i,j-1,k)+DUon(i,j,k));

        VFx(i,j,0)=0.25*(cff1-vxx_avg*cff)* (cff2-cff*(Huee_j+ Huee_jm1));
    });

    ParallelFor(growLo(ybx,1,1),
    [=] AMREX_GPU_DEVICE (int i, int j, int)
    {
        Real vee_j    = vbar(i,j-1,0,krhs)-2.0*vbar(i,j  ,0,krhs)+vbar(i,j+1,0,krhs);
        Real vee_jp1  = vbar(i,j  ,0,krhs)-2.0*vbar(i,j+1,0,krhs)+vbar(i,j+2,0,krhs);
        Real vee_avg  = vee_j + vee_jp1;

        Real Hvee_j   = DVom(i,j-1,0)-2.0*DVom(i,j  ,0)+DVom(i,j+1,0);
        Real Hvee_jp1 = DVom(i,j  ,0)-2.0*DVom(i,j+1,0)+DVom(i,j+2,0);

        Real cff=1.0/6.0;
        Real cff1=vbar(i,j  ,0,krhs)+vbar(i,j+1,0,krhs);

        VFe(i,j,0) = 0.25 * (cff1-vee_avg*cff) * (DVom(i,j  ,0)+ DVom(i,j+1,0) -
                                           cff  * (Hvee_j+ Hvee_jp1));
    });

    ParallelFor(makeSlab(ybx,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int)
    {
        Real cff1=VFx(i+1,j,0)-VFx(i  ,j,0);
        Real cff2=VFe(i  ,j,0)-VFe(i,j-1,0);

        rhs_vbar(i,j,0) -= (cff1 + cff2);
    });

}
