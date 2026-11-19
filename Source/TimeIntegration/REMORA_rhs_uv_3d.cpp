#include <REMORA.H>

using namespace amrex;

/**
 * @param[in   ] lev            level to operate on
 * @param[in   ] xbx     Box for operations on x-velocity
 * @param[in   ] ybx     Box for operations on y-velocity
 * @param[in   ] uold    u-velocity at last time step
 * @param[in   ] vold    v-velocity at last time step
 * @param[inout] ru      u-velocity RHS
 * @param[inout] rv      v-velocity RHS
 * @param[inout] rufrc   forcing for u-velocity RHS
 * @param[inout] rvfrc   forcing for v-velocity RHS
 * @param[in   ] sustr   u-direction surface momentum flux
 * @param[in   ] svstr   v-direction surface momentum flux
 * @param[in   ] bustr   u-direction bottom stress
 * @param[in   ] bvstr   v-direction bottom stress
 * @param[in   ] Huon    u-volume flux
 * @param[in   ] Hvom    v-volume flux
 * @param[in   ] omn_u   om_u*on_u at u points
 * @param[in   ] omn_v   om_v*on_v at v points
 * @param[in   ] W       vertical velocity
 * @param[in   ] nrhs    index of component for RHS
 * @param[in   ] N       number of vertical levels
 */
void
REMORA::rhs_uv_3d (int lev,
                  const Box& xbx, const Box& ybx,
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
                  const Array4<Real const>& omn_u,
                  const Array4<Real const>& omn_v,
                  const Array4<Real const>& W   ,
                  int nrhs, int N)
{
    BL_PROFILE("REMORA::rhs_uv_3d()");
    const MultiFab* bottom_transport = GetBottomExchangeTransport(lev);
    const Box& domain = geom[lev].Domain();
    const auto dlo = amrex::lbound(domain);
    const auto dhi = amrex::ubound(domain);

    GeometryData const& geomdata = geom[0].data();
    bool is_periodic_in_x = geomdata.isPeriodic(0);
    bool is_periodic_in_y = geomdata.isPeriodic(1);

    //
    // Scratch space, each defined on exactly the box the flux kernels below write
    //
    FArrayBox fab_UFx(growLo(xbx,0,1),1,amrex::The_Async_Arena());
    FArrayBox fab_UFe(growHi(xbx,1,1),1,amrex::The_Async_Arena());
    FArrayBox fab_VFe(growLo(ybx,1,1),1,amrex::The_Async_Arena());
    FArrayBox fab_VFx(growHi(ybx,0,1),1,amrex::The_Async_Arena());

    auto UFx=fab_UFx.array();
    auto UFe=fab_UFe.array();
    auto VFx=fab_VFx.array();
    auto VFe=fab_VFe.array();

    auto uv_hadv_scheme = solverChoice.uv_Hadv_scheme;
    // The flux scratch is not zeroed: each scheme below writes every flux it reads
    AMREX_ALWAYS_ASSERT(uv_hadv_scheme == AdvectionScheme::upstream3 ||
                        uv_hadv_scheme == AdvectionScheme::centered2);

    //check this////////////
    const Real Gadv = Real(-0.25);

    if (uv_hadv_scheme == AdvectionScheme::upstream3) {
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

            Real uxx_i   = uold(i-1,j,k,nrhs)-two*uold(i  ,j,k,nrhs)+uold(i+1,j,k,nrhs);
            Real uxx_ip1 = uold(i  ,j,k,nrhs)-two*uold(i+1,j,k,nrhs)+uold(i+2,j,k,nrhs);
            // Upwinding

            Real Huxx_i   = Huon(i-1,j,k)-two*Huon(i  ,j,k)+Huon(i+1,j,k);
            Real Huxx_ip1 = Huon(i  ,j,k)-two*Huon(i+1,j,k)+Huon(i+2,j,k);

            if (i == dlo.x && !is_periodic_in_x) {
                uxx_i = uxx_ip1;
                Huxx_i = Huxx_ip1;
            }
            else if (i == dhi.x && !is_periodic_in_x) {
                uxx_ip1 = uxx_i;
                Huxx_ip1 = Huxx_i;
            }

            Real cff = (cff1 > zero) ? uxx_i : uxx_ip1;

            Real Huon_avg = (Huon(i,j,k) + Huon(i+1,j,k));

            UFx(i,j,k) = Real(0.25)*(cff1+Gadv*cff) * ( Huon_avg + Real(0.5)*Gadv*(Huxx_i + Huxx_ip1) );
        });

        //
        // Define UFe, the y-fluxes at nodes for updating u
        //
        ParallelFor(growHi(xbx,1,1), [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real cff1 = uold(i,j,k,nrhs) + uold(i  ,j-1,k,nrhs);
            Real cff2 = Hvom(i,j,k)      + Hvom(i-1,j  ,k);

            Real uee_jm1 = uold(i,j-2,k,nrhs) - two*uold(i,j-1,k,nrhs) + uold(i  ,j,k,nrhs);
            Real uee_j   = uold(i,j-1,k,nrhs) - two*uold(i,j  ,k,nrhs) + uold(i,j+1,k,nrhs);

            if (j == dlo.y and !is_periodic_in_y) {
                uee_jm1 = uee_j;
            } else if (j == dhi.y+1 and !is_periodic_in_y) {
                uee_j = uee_jm1;
            }

            // Upwinding
            Real cff = (cff2 > zero) ?  uee_jm1 : uee_j;

            Real Hvxx_i   = Hvom(i-1,j,k)-Real(2.0)*Hvom(i  ,j,k)+Hvom(i+1,j,k);
            Real Hvxx_im1 = Hvom(i-2,j,k)-Real(2.0)*Hvom(i-1,j,k)+Hvom(i  ,j,k);

            UFe(i,j,k) = Real(0.25) * (cff1+Gadv*cff)* (cff2+Gadv*Real(0.5)*(Hvxx_i + Hvxx_im1));
        });
    } else if (uv_hadv_scheme == AdvectionScheme::centered2) {
        ParallelFor(growLo(xbx,0,1), [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            UFx(i,j,k) = Real(0.25) * (uold(i,j,k,nrhs) + uold(i+1,j,k,nrhs)) * (Huon(i,j,k)+Huon(i+1,j,k));
        });
        ParallelFor(growHi(xbx,1,1), [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            UFe(i,j,k) = Real(0.25) * (uold(i,j-1,k,nrhs) + uold(i,j,k,nrhs)) * (Hvom(i-1,j,k)+Hvom(i,j,k));
        });
    }

    if (uv_hadv_scheme == AdvectionScheme::upstream3) {
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

            Real vxx_im1 = vold(i-2,j,k,nrhs)-Real(2.0)*vold(i-1,j,k,nrhs)+vold(i  ,j,k,nrhs);
            Real vxx_i   = vold(i-1,j,k,nrhs)-Real(2.0)*vold(i  ,j,k,nrhs)+vold(i+1,j,k,nrhs);

            if (i == dlo.x and !is_periodic_in_x) {
                vxx_im1 = vxx_i;
            } else if (i == dhi.x+1 and !is_periodic_in_x) {
                vxx_i = vxx_im1;
            }

            // Upwinding
            Real cff = (cff2 > zero) ? vxx_im1 : vxx_i;


            Real Huee_j   = Huon(i,j-1,k)-Real(2.0)*Huon(i,j  ,k)+Huon(i,j+1,k);
            Real Huee_jm1 = Huon(i,j-2,k)-Real(2.0)*Huon(i,j-1,k)+Huon(i,j  ,k);

            VFx(i,j,k) = Real(0.25)*(cff1+Gadv*cff)* (cff2+Gadv*Real(0.5)*(Huee_j + Huee_jm1));
        });

        // Grow ybx by one in low y-direction
        ParallelFor(growLo(ybx,1,1), [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real cff1=vold(i,j,k,nrhs)+vold(i,j+1,k,nrhs);

            // Upwinding
            Real vee_j    = vold(i,j-1,k,nrhs)-Real(2.0)*vold(i,j  ,k,nrhs)+ vold(i,j+1,k,nrhs);
            Real vee_jp1  = vold(i,j  ,k,nrhs)-Real(2.0)*vold(i,j+1,k,nrhs)+ vold(i,j+2,k,nrhs);

            Real Hvee_j   = Hvom(i,j-1,k)-Real(2.0)*Hvom(i,j  ,k)+Hvom(i,j+1,k);
            Real Hvee_jp1 = Hvom(i,j  ,k)-Real(2.0)*Hvom(i,j+1,k)+Hvom(i,j+2,k);

            if (j == dlo.y and !is_periodic_in_y) {
                vee_j = vee_jp1;
                Hvee_j = Hvee_jp1;
            }
            else if (j == dhi.y and !is_periodic_in_y) {
                vee_jp1 = vee_j;
                Hvee_jp1 = Hvee_j;
            }
            Real cff = (cff1 > zero) ? vee_j : vee_jp1;

            VFe(i,j,k) = Real(0.25) * (cff1+Gadv*cff) * ( Hvom(i,j  ,k)+ Hvom(i,j+1,k) + Real(0.5) * Gadv * (Hvee_j + Hvee_jp1) );
        });
    } else if (uv_hadv_scheme == AdvectionScheme::centered2) {
        ParallelFor(growHi(ybx,0,1), [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            VFx(i,j,k) = Real(0.25) * (vold(i-1,j,k,nrhs) + vold(i,j,k,nrhs)) * (Huon(i,j-1,k)+Huon(i,j,k));
        });
        ParallelFor(growLo(ybx,1,1), [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            VFe(i,j,k) = Real(0.25) * (vold(i,j,k,nrhs) + vold(i,j+1,k,nrhs)) * (Hvom(i,j,k)+Hvom(i,j+1,k));
        });
    }

    //
    // Update ru and rv level by level from the bottom up: subtract the horizontal advection
    // (the divergence of the fluxes above) and the vertical advection, and add the result
    // and the surface and bottom stresses to the depth integrals rufrc and rvfrc. The
    // vertical flux FC at the top face of a level is computed there and kept in a 2D buffer
    // as the flux at the bottom face of the next level. Only outflow carries horizontal
    // momentum through the bottom; the surface flux is zero.
    //
    AMREX_ASSERT(xbx.smallEnd(2) == 0 && xbx.bigEnd(2) == N);
    AMREX_ASSERT(ybx.smallEnd(2) == 0 && ybx.bigEnd(2) == N);
    const Box xbx2d = makeSlab(xbx,2,0);
    const Box ybx2d = makeSlab(ybx,2,0);
    FArrayBox fab_FCu(xbx2d,1,amrex::The_Async_Arena());
    FArrayBox fab_FCv(ybx2d,1,amrex::The_Async_Arena());
    auto FCu = fab_FCu.array();
    auto FCv = fab_FCv.array();

    const bool bottom_exchange = bottom_transport != nullptr;
    const int upstream3 = (uv_hadv_scheme == AdvectionScheme::upstream3) ? 1 : 0;
    AnyCTO(TypeList<CompileTimeOptions<0,1>>{}, {upstream3},
    [&] (auto const& fu, auto const& fv) {
        ParallelFor(xbx2d, N+1, fu, ybx2d, N+1, fv);
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int, int k, auto up3)
    {
        constexpr int UP3 = decltype(up3)::value;
        const Real cff1=Real(9.0)/Real(16.0);
        const Real cff2=one/Real(16.0);

        // Vertical flux at w-level kf
        auto FC = [&] (int kf) -> Real
        {
            if (kf == 0) {
                return bottom_exchange ? Real(0.5) *
                    (amrex::min(W(i,j,0),Real(0.0)) +
                     amrex::min(W(i-1,j,0),Real(0.0))) * uold(i,j,0,nrhs) : zero;
            } else if (kf == N+1) {
                return zero;
            } else if (UP3 == 1) {
                if (kf == 1) {
                    return ( cff1*(uold(i  ,j,0,nrhs)+ uold(i,j,1,nrhs))
                            -cff2*(uold(i  ,j,0,nrhs)+ uold(i,j,2,nrhs)) )*
                               ( cff1*(   W(i  ,j,1)+ W(i-1,j,1))
                                -cff2*(   W(i+1,j,1)+ W(i-2,j,1)) );
                } else if (kf == N) {
                    return ( cff1*(uold(i  ,j,N-1,nrhs)+ uold(i,j,N  ,nrhs))
                            -cff2*(uold(i  ,j,N-2,nrhs)+ uold(i,j,N  ,nrhs)) )*
                               ( cff1*(   W(i  ,j,N)+ W(i-1,j,N))
                                -cff2*(   W(i+1,j,N)+ W(i-2,j,N)) );
                } else {
                    return ( cff1*(uold(i  ,j,kf-1,nrhs)+ uold(i,j,kf  ,nrhs))
                            -cff2*(uold(i  ,j,kf-2,nrhs)+ uold(i,j,kf+1,nrhs)) )*
                               ( cff1*(   W(i  ,j,kf)+ W(i-1,j,kf))
                                -cff2*(   W(i+1,j,kf)+ W(i-2,j,kf)) );
                }
            } else {
                return Real(0.25) * (uold(i,j,kf-1,nrhs)+uold(i,j,kf,nrhs)) * (W(i,j,kf) + W(i-1,j,kf));
            }
        };

        const Real FC_lo = (k == 0) ? FC(0) : FCu(i,j,0);
        const Real FC_hi = FC(k+1);
        FCu(i,j,0) = FC_hi;

        Real r = ru(i,j,k,nrhs) - ( (UFx(i,j,k)-UFx(i-1,j,k)) + (UFe(i,j+1,k)-UFe(i  ,j,k)) );
        r -= FC_hi-FC_lo;
        ru(i,j,k,nrhs) = r;

        rufrc(i,j,0) += r;

        Real cff  = omn_u(i,j,0);

        Real cff3 = (k == N) ?  sustr(i,j,0)*cff : zero;
        Real cff4 = (k == 0) ? -bustr(i,j,0)*cff : zero;

        rufrc(i,j,0) += cff3+cff4;
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int, int k, auto up3)
    {
        constexpr int UP3 = decltype(up3)::value;
        const Real cff1=Real(9.0)/Real(16.0);
        const Real cff2=one/Real(16.0);

        // Vertical flux at w-level kf
        auto FC = [&] (int kf) -> Real
        {
            if (kf == 0) {
                return bottom_exchange ? Real(0.5) *
                    (amrex::min(W(i,j,0),Real(0.0)) +
                     amrex::min(W(i,j-1,0),Real(0.0))) * vold(i,j,0,nrhs) : zero;
            } else if (kf == N+1) {
                return zero;
            } else if (UP3 == 1) {
                if (kf == 1) {
                    return ( cff1*(vold(i,j,0,nrhs)+ vold(i,j,1,nrhs))
                            -cff2*(vold(i,j,0,nrhs)+ vold(i,j,2,nrhs)) )*
                           ( cff1*(W(i,j  ,1)+ W(i,j-1,1))
                            -cff2*(W(i,j+1,1)+ W(i,j-2,1)) );
                } else if (kf == N) {
                    return ( cff1*(vold(i,j,N-1,nrhs)+ vold(i,j,N  ,nrhs))
                            -cff2*(vold(i,j,N-2,nrhs)+ vold(i,j,N  ,nrhs)) )*
                           ( cff1*(W(i,j  ,N)+ W(i,j-1,N))
                            -cff2*(W(i,j+1,N)+ W(i,j-2,N)) );
                } else {
                    return ( cff1*(vold(i,j,kf-1,nrhs)+ vold(i,j,kf  ,nrhs))
                            -cff2*(vold(i,j,kf-2,nrhs)+ vold(i,j,kf+1,nrhs)) )*
                           ( cff1*(W(i,j  ,kf)+ W(i,j-1,kf))
                            -cff2*(W(i,j+1,kf)+ W(i,j-2,kf)) );
                }
            } else {
                return Real(0.25) * (vold(i,j,kf-1,nrhs)+vold(i,j,kf,nrhs)) * (W(i,j,kf) + W(i,j-1,kf));
            }
        };

        const Real FC_lo = (k == 0) ? FC(0) : FCv(i,j,0);
        const Real FC_hi = FC(k+1);
        FCv(i,j,0) = FC_hi;

        Real r = rv(i,j,k,nrhs) - ( (VFx(i+1,j,k)-VFx(i,j,k)) + (VFe(i,j,k)-VFe(i,j-1,k)) );
        r -= FC_hi-FC_lo;
        rv(i,j,k,nrhs) = r;

        rvfrc(i,j,0) += r;

        Real cff = omn_v(i,j,0);

        Real cff3 = (k == N) ?  svstr(i,j,0)*cff : zero;
        Real cff4 = (k == 0) ? -bvstr(i,j,0)*cff : zero;

        rvfrc(i,j,0) += cff3+cff4;
    });
}
