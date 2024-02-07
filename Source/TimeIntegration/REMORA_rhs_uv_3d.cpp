#include <REMORA.H>

using namespace amrex;

/**
 * rhs_uv_2d
 *
 * @param[in   ] xbx Box for operations on x-velocity
 * @param[in   ] ybx Box for operations on y-velocity
 * @param[in   ] uold
 * @param[in   ] vold
 * @param[  out] ru
 * @param[  out] rv
 * @param[  out] rufrc
 * @param[  out] rvfrc
 * @param[in   ] sustr
 * @param[in   ] svstr
 * @param[in   ] bustr
 * @param[in   ] bvstr
 * @param[in   ] Huon
 * @param[in   ] Hvom
 * @param[in   ] pm
 * @param[in   ] pn
 * @param[in   ] W
 * @param[inout] FC
 * @param[in   ] nrhs
 * @param[in   ] N
 */

void
REMORA::rhs_uv_3d (const Box& xbx, const Box& ybx,
                  const Array4<Real const>& uold  ,
                  const Array4<Real const>& vold,
                  const Array4<Real      >& ru,
                  const Array4<Real      >& rv,
                  const Array4<Real      >& rufrc,
                  const Array4<Real      >& rvfrc,
                  const Array4<Real const>& sustr,
                  const Array4<Real const>& svstr,
                  const Array4<Real const>& bustr,
                  const Array4<Real const>& bvstr,
                  const Array4<Real const>& Huon,
                  const Array4<Real const>& Hvom,
                  const Array4<Real const>& pm,
                  const Array4<Real const>& pn,
                  const Array4<Real const>& W   ,
                  const Array4<Real      >& FC,
                  int nrhs, int N)
{
    const Box& domain = geom[0].Domain();
    const auto dlo = amrex::lbound(domain);
    const auto dhi = amrex::ubound(domain);

    int ncomp = 1;
    Vector<BCRec> bcrs_x(ncomp);
    Vector<BCRec> bcrs_y(ncomp);
    amrex::setBC(xbx,domain,BCVars::xvel_bc,0,1,domain_bcs_type,bcrs_x);
    amrex::setBC(ybx,domain,BCVars::yvel_bc,0,1,domain_bcs_type,bcrs_y);
    auto bcr_x = bcrs_x[0];
    auto bcr_y = bcrs_y[0];

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

    //check this////////////
    const Real Gadv = -0.25_rt;

    // *************************************************************
    // UPDATING U
    // *************************************************************

    // Think of the cell-centered box as                              [ 0:nx-1, 0:ny-1] (0,0,0) cc
    //
    // xbx is the x-face-centered box on which we update u (with ru)  [ 0:nx  , 0:ny-1] (1,0,0) x-faces
    //   to do so requires                                  UFx on    [-1:nx  , 0:ny-1] (0,0,0) cc
    //      which requires                                  uold on   [-2:nx+2, 0:ny-1] (1,0,0) x-faces
    //       and  requires                                  Huon on   [-2:nx+2, 0:ny-1] (0,0,0) x-faces
    //   to do so requires                                  UFe on    [ 0:nx  , 0:ny  ] (1,1,0) xy-nodes
    //      which requires                                  uold on   [ 0:nx  ,-2:ny+1] (1,0,0) x-faces
    //       and  requires                                  Hvom on   [-2:nx+1, 0:ny-1] (0,1,0) y-faces

    //
    // Define UFx, the x-fluxes at cell centers for updating u
    // (Note that grow arguments are (bx, dir, ng)
    //
    ParallelFor(growLo(xbx,0,1), [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        Real cff1 = uold(i,j,k,nrhs)+uold(i+1,j,k,nrhs);

        Real uxx_i   = uold(i-1,j,k,nrhs)-2.0_rt*uold(i  ,j,k,nrhs)+uold(i+1,j,k,nrhs);
        Real uxx_ip1 = uold(i  ,j,k,nrhs)-2.0_rt*uold(i+1,j,k,nrhs)+uold(i+2,j,k,nrhs);
        // Upwinding

        Real Huxx_i   = Huon(i-1,j,k)-2.0_rt*Huon(i  ,j,k)+Huon(i+1,j,k);
        Real Huxx_ip1 = Huon(i  ,j,k)-2.0_rt*Huon(i+1,j,k)+Huon(i+2,j,k);

        if (i == dlo.x && bcr_x.lo(0) == REMORABCType::ext_dir) {
            uxx_i = uxx_ip1;
            Huxx_i = Huxx_ip1;
        }
        else if (i == dhi.x && bcr_x.hi(0) == REMORABCType::ext_dir) {
            uxx_ip1 = uxx_i;
            Huxx_ip1 = Huxx_i;
        }

        Real cff = (cff1 > 0.0_rt) ? uxx_i : uxx_ip1;

        Real Huon_avg = (Huon(i,j,k) + Huon(i+1,j,k));

        UFx(i,j,k) = 0.25_rt*(cff1+Gadv*cff) * ( Huon_avg + 0.5_rt*Gadv*(Huxx_i + Huxx_ip1) );
    });

    //
    // Define UFe, the y-fluxes at nodes for updating u
    //
    ParallelFor(growHi(xbx,1,1), [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        Real cff1 = uold(i,j,k,nrhs) + uold(i  ,j-1,k,nrhs);
        Real cff2 = Hvom(i,j,k)      + Hvom(i-1,j  ,k);

        Real uee_jm1 = uold(i,j-2,k,nrhs) - 2.0_rt*uold(i,j-1,k,nrhs) + uold(i  ,j,k,nrhs);
        Real uee_j   = uold(i,j-1,k,nrhs) - 2.0_rt*uold(i,j  ,k,nrhs) + uold(i,j+1,k,nrhs);

        if (j == dlo.y and bcr_y.lo(1) == REMORABCType::ext_dir) {
            uee_jm1 = uee_j;
        } else if (j == dhi.y+1 and bcr_y.hi(1) == REMORABCType::ext_dir) {
            uee_j = uee_jm1;
        }

        // Upwinding
        Real cff = (cff2 > 0.0_rt) ?  uee_jm1 : uee_j;

        Real Hvxx_i   = Hvom(i-1,j,k)-2.0_rt*Hvom(i  ,j,k)+Hvom(i+1,j,k);
        Real Hvxx_im1 = Hvom(i-2,j,k)-2.0_rt*Hvom(i-1,j,k)+Hvom(i  ,j,k);

        UFe(i,j,k) = 0.25_rt * (cff1+Gadv*cff)* (cff2+Gadv*0.5_rt*(Hvxx_i + Hvxx_im1));
    });

    //
    // Define the RHS for u by differencing fluxes
    //
    ParallelFor(xbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
          ru(i,j,k,nrhs) -= ( (UFx(i,j,k)-UFx(i-1,j,k)) + (UFe(i,j+1,k)-UFe(i  ,j,k)) );
    });

    ParallelFor(xbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
          //-----------------------------------------------------------------------
          //  Add in vertical advection.
          //-----------------------------------------------------------------------
          Real cff1=9.0_rt/16.0_rt;
          Real cff2=1.0_rt/16.0_rt;

          if (k>=1 && k<=N-2)
          {
                  FC(i,j,k)=( cff1*(uold(i  ,j,k  ,nrhs)+ uold(i,j,k+1,nrhs))
                             -cff2*(uold(i  ,j,k-1,nrhs)+ uold(i,j,k+2,nrhs)) )*
                                ( cff1*(   W(i  ,j,k)+ W(i-1,j,k))
                                 -cff2*(   W(i+1,j,k)+ W(i-2,j,k)) );
          }
          else // this needs to be split up so that the following can be concurrent
          {
              FC(i,j,N)=0.0_rt;

              FC(i,j,N-1)=( cff1*(uold(i  ,j,N-1,nrhs)+ uold(i,j,N  ,nrhs))
                           -cff2*(uold(i  ,j,N-2,nrhs)+ uold(i,j,N  ,nrhs)) )*
                              ( cff1*(   W(i  ,j,N-1)+ W(i-1,j,N-1))
                               -cff2*(   W(i+1,j,N-1)+ W(i-2,j,N-1)) );

              FC(i,j,0)=( cff1*(uold(i  ,j,0,nrhs)+ uold(i,j,1,nrhs))
                         -cff2*(uold(i  ,j,0,nrhs)+ uold(i,j,2,nrhs)) )*
                            ( cff1*(   W(i  ,j,0)+ W(i-1,j,0))
                             -cff2*(   W(i+1,j,0)+ W(i-2,j,0)) );

              //              FC(i,0,-1)=0.0_rt;
          }
    });

    ParallelFor(xbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        Real cff = (k >= 1) ? FC(i,j,k)-FC(i,j,k-1) : FC(i,j,k);

        ru(i,j,k,nrhs) -= cff;
    });

    Gpu::synchronize();

    AMREX_ASSERT(xbx.smallEnd(2) == 0 && xbx.bigEnd(2) == N);
    ParallelFor(makeSlab(xbx,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int)
    {
       for (int k = 0; k <= N; ++k)
       {
          rufrc(i,j,0) += ru(i,j,k,nrhs);

          Real om_u = 2.0_rt / (pm(i-1,j,0)+pm(i,j,0));
          Real on_u = 2.0_rt / (pn(i-1,j,0)+pn(i,j,0));
          Real cff  = om_u * on_u;

          Real cff1 = (k == N) ?  sustr(i,j,0)*cff : 0.0_rt;
          Real cff2 = (k == 0) ? -bustr(i,j,0)*cff : 0.0_rt;

          rufrc(i,j,0) += cff1+cff2;
       }
    });

    // *************************************************************
    // UPDATING V
    // *************************************************************

    // Think of the cell-centered box as                              [ 0:nx-1, 0:ny-1] (0,0,0) cc
    //
    // ybx is the y-face-centered box on which we update v (with rv)  [ 0:nx-1, 0:ny  ] (1,0,0)  y-faces
    //   to do so requires                                  VFe on    [ 0:nx-1,-1:ny  ] (1,1,0)  xy-nodes
    //      which requires                                    vold on [ 0:nx-1,-2:ny+2] (1,0,0)  y-faces
    //       and  requires                                    Hvom on [ 0:nx-1,-2:ny+2] (0,1,0)  x-faces
    //   to do so requires                                  VFx on    [ 0:nx  , 0:ny  ] (0,0,0)  cc
    //      which requires                                    vold on [-2:nx+1, -:ny  ] (1,0,0)  y-faces
    //       and  requires                                    Hvom on [ 0:nx-1,-2:ny+1] (0,0,0)  y-faces

    // Grow ybx by one in high x-direction
    ParallelFor(growHi(ybx,0,1), [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        Real cff1 = vold(i,j,k,nrhs) + vold(i-1,j  ,k,nrhs);
        Real cff2 = Huon(i,j,k)      + Huon(i  ,j-1,k);

        Real vxx_im1 = vold(i-2,j,k,nrhs)-2.0_rt*vold(i-1,j,k,nrhs)+vold(i  ,j,k,nrhs);
        Real vxx_i   = vold(i-1,j,k,nrhs)-2.0_rt*vold(i  ,j,k,nrhs)+vold(i+1,j,k,nrhs);

        if (i == dlo.x and bcr_x.lo(0) == REMORABCType::ext_dir) {
            vxx_i = vxx_im1;
        } else if (i == dhi.x+1 and bcr_x.hi(0) == REMORABCType::ext_dir) {
            vxx_im1 = vxx_i;
        }

        // Upwinding
        Real cff = (cff2 > 0.0_rt) ? vxx_im1 : vxx_i;


        Real Huee_j   = Huon(i,j-1,k)-2.0_rt*Huon(i,j  ,k)+Huon(i,j+1,k);
        Real Huee_jm1 = Huon(i,j-2,k)-2.0_rt*Huon(i,j-1,k)+Huon(i,j  ,k);

        VFx(i,j,k) = 0.25_rt*(cff1+Gadv*cff)* (cff2+Gadv*0.5_rt*(Huee_j + Huee_jm1));
    });

    // Grow ybx by one in low y-direction
    ParallelFor(growLo(ybx,1,1), [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        Real cff1=vold(i,j,k,nrhs)+vold(i,j+1,k,nrhs);

        // Upwinding
        Real vee_j    = vold(i,j-1,k,nrhs)-2.0_rt*vold(i,j  ,k,nrhs)+ vold(i,j+1,k,nrhs);
        Real vee_jp1  = vold(i,j  ,k,nrhs)-2.0_rt*vold(i,j+1,k,nrhs)+ vold(i,j+2,k,nrhs);

        Real Hvee_j   = Hvom(i,j-1,k)-2.0_rt*Hvom(i,j  ,k)+Hvom(i,j+1,k);
        Real Hvee_jp1 = Hvom(i,j  ,k)-2.0_rt*Hvom(i,j+1,k)+Hvom(i,j+2,k);

        if (j == dlo.y and bcr_y.lo(1) == REMORABCType::ext_dir) {
            vee_j = vee_jp1;
            Hvee_j = Hvee_jp1;
        }
        else if (j == dhi.y and bcr_y.hi(1) == REMORABCType::ext_dir) {
            vee_jp1 = vee_j;
            Hvee_jp1 = Hvee_j;
        }
        Real cff = (cff1 > 0.0_rt) ? vee_j : vee_jp1;

        VFe(i,j,k) = 0.25_rt * (cff1+Gadv*cff) * ( Hvom(i,j  ,k)+ Hvom(i,j+1,k) + 0.5_rt * Gadv * (Hvee_j + Hvee_jp1) );
    });

    ParallelFor(ybx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
          rv(i,j,k,nrhs) -= ( (VFx(i+1,j,k)-VFx(i,j,k)) + (VFe(i,j,k)-VFe(i,j-1,k)) );
    });

    ParallelFor(ybx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
          Real cff1=9.0_rt/16.0_rt;
          Real cff2=1.0_rt/16.0_rt;
          if (k>=1 && k<=N-2)
          {
              FC(i,j,k)=( cff1*(vold(i,j,k  ,nrhs)+ vold(i,j,k+1,nrhs))
                         -cff2*(vold(i,j,k-1,nrhs)+ vold(i,j,k+2,nrhs)) )*
                            ( cff1*(W(i,j  ,k)+ W(i,j-1,k))
                             -cff2*(W(i,j+1,k)+ W(i,j-2,k)) );
          }
          else // this needs to be split up so that the following can be concurrent
          {
              FC(i,j,N)=0.0_rt;
              FC(i,j,N-1)=( cff1*(vold(i,j,N-1,nrhs)+ vold(i,j,N  ,nrhs))
                           -cff2*(vold(i,j,N-2,nrhs)+ vold(i,j,N  ,nrhs)) )*
                              ( cff1*(W(i,j  ,N-1)+ W(i,j-1,N-1))
                               -cff2*(W(i,j+1,N-1)+ W(i,j-2,N-1)) );
              FC(i,j,0)=( cff1*(vold(i,j,0,nrhs)+ vold(i,j,1,nrhs))
                             -cff2*(vold(i,j,0,nrhs)+ vold(i,j,2,nrhs)) )*
                            ( cff1*(W(i,j  ,0)+ W(i,j-1,0))
                             -cff2*(W(i,j+1,0)+ W(i,j-2,0)) );
              //              FC(i,0,-1)=0.0_rt;
          }
    }); Gpu::synchronize();

    ParallelFor(ybx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        Real cff = (k >= 1) ? FC(i,j,k)-FC(i,j,k-1) : FC(i,j,k);

        rv(i,j,k,nrhs) -= cff;
    });

    Gpu::synchronize();

    AMREX_ASSERT(ybx.smallEnd(2) == 0 && ybx.bigEnd(2) == N);
    ParallelFor(makeSlab(ybx,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int)
    {
       for (int k = 0; k <= N; ++k)
       {
          rvfrc(i,j,0) += rv(i,j,k,nrhs);

          Real om_v = 2.0_rt / (pm(i,j-1,0)+pm(i,j,0));
          Real on_v = 2.0_rt / (pn(i,j-1,0)+pn(i,j,0));
          Real cff = om_v * on_v;

          Real cff1 = (k == N) ?  svstr(i,j,0)*cff : 0.0_rt;
          Real cff2 = (k == 0) ? -bvstr(i,j,0)*cff : 0.0_rt;

          rvfrc(i,j,0) += cff1+cff2;
       }
    });
}
