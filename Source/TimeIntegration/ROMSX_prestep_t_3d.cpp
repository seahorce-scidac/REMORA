#include <ROMSX.H>

using namespace amrex;

//
// prestep_uv_3d
//

void
ROMSX::prestep_t_3d (const Box& tbx, const Box& gbx,
                      Array4<Real> tempold,
                      Array4<Real> temp,
                      Array4<Real> tempcache,
                      Array4<Real> ru,
                      Array4<Real> Hz, Array4<Real> Akt,
                      Array4<Real> Huon, Array4<Real> Hvom,
                      Array4<Real> pm,   Array4<Real> pn,
                      Array4<Real> W  ,  Array4<Real> DC,
                      Array4<Real> FC ,  Array4<Real> tempstore,
                      Array4<Real> z_r,  Array4<Real> z_w, Array4<Real> h,
                      int iic, int ntfirst, int nnew, int nstp, int nrhs, int N,
                      Real lambda, Real dt_lev)
{
    //copy the tilebox
    Box gbx1 = tbx;
    Box gbx2 = tbx;
    Box tbxp1 = tbx;
    Box tbxp2 = tbx;

    tbxp1.grow(IntVect(NGROW-1,NGROW-1,0));
    tbxp2.grow(IntVect(NGROW,NGROW,0));
    FArrayBox fab_FX(tbxp2,1,amrex::The_Async_Arena()); //3D
    FArrayBox fab_FE(tbxp2,1,amrex::The_Async_Arena()); //3D
    FArrayBox fab_curv(tbxp2,1,amrex::The_Async_Arena()); //fab_curv.setVal(0.0);
    FArrayBox fab_grad(tbxp2,1,amrex::The_Async_Arena()); //fab_curv.setVal(0.0);

    auto FX=fab_FX.array();
    auto FE=fab_FE.array();
    auto curv=fab_curv.array();
    auto grad=fab_grad.array();
    ParallelFor(tbxp2,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        grad(i,j,k)=0.0;

        curv(i,j,k)=0.0;

        FX(i,j,k)=0.0;
        FE(i,j,k)=0.0;
    });


    Box ubx = surroundingNodes(tbx,0);
    Box vbx = surroundingNodes(tbx,1);

    Box utbxp1 = surroundingNodes(tbxp1,0);
    Box vtbxp1 = surroundingNodes(tbxp1,1);

    Box gbx3uneven_init(IntVect(AMREX_D_DECL(tbx.smallEnd(0)-3,tbx.smallEnd(1)-3,tbx.smallEnd(2))),
                   IntVect(AMREX_D_DECL(tbx.bigEnd(0)+2,tbx.bigEnd(1)+2,tbx.bigEnd(2))));
    BoxArray ba_gbx3uneven = intersect(BoxArray(gbx3uneven_init), gbx);
    AMREX_ASSERT((ba_gbx3uneven.size() == 1));
    Box gbx3uneven = ba_gbx3uneven[0];

    gbx2.grow(IntVect(NGROW,NGROW,0));
    BoxArray ba_gbx2 = intersect(BoxArray(gbx2), gbx);
    AMREX_ASSERT((ba_gbx2.size() == 1));
    gbx2 = ba_gbx2[0];

    gbx1.grow(IntVect(NGROW-1,NGROW-1,0));
    BoxArray ba_gbx1 = intersect(BoxArray(gbx1), gbx);
    AMREX_ASSERT((ba_gbx1.size() == 1));
    gbx1 = ba_gbx1[0];

    //------------------------------------------------------------------------
    //  Vertically integrate horizontal mass flux divergence.
    //------------------------------------------------------------------------
    //
    //Should really use gbx3uneven
    Box gbx3unevenD = gbx3uneven;
    gbx3unevenD.makeSlab(2,0);
    Box gbx1D = gbx1;
    gbx1D.makeSlab(2,0);

    ParallelFor(gbx1D,
    [=] AMREX_GPU_DEVICE (int i, int j, int )
    {
        //  Starting with zero vertical velocity at the bottom, integrate
        //  from the bottom (k=0) to the free-surface (k=N).  The w(:,:,N(ng))
        //  contains the vertical velocity at the free-surface, d(zeta)/d(t).
        //  Notice that barotropic mass flux divergence is not used directly.
        //
        //W(i,j,-1)=0.0;
        int k=0;
        W(i,j,k) = - (Huon(i+1,j,k)-Huon(i,j,k)) - (Hvom(i,j+1,k)-Hvom(i,j,k));
        for(k=1;k<=N;k++) {
            W(i,j,k) = W(i,j,k-1) - (Huon(i+1,j,k)-Huon(i,j,k)) - (Hvom(i,j+1,k)-Hvom(i,j,k));
        }
    });
    ParallelFor(gbx1,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        //  Starting with zero vertical velocity at the bottom, integrate
        //  from the bottom (k=0) to the free-surface (k=N).  The w(:,:,N(ng))
        //  contains the vertical velocity at the free-surface, d(zeta)/d(t).
        //  Notice that barotropic mass flux divergence is not used directly.
        //
        Real wrk_i=W(i,j,N)/(z_w(i,j,N)+h(i,j,0,0));

        if(k!=N) {
            W(i,j,k) = W(i,j,k)- wrk_i*(z_w(i,j,k)+h(i,j,0,0));
        }
    });

    ParallelFor(gbx1,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        if (k == N) {
            W(i,j,N) = 0.0;
        }
    });

    //From ini_fields and .in file
    //fab_Akt.setVal(1e-6);
    FArrayBox fab_stflux(tbxp2,1,amrex::The_Async_Arena());
    auto stflux= fab_stflux.array();
    FArrayBox fab_btflux(tbxp2,1,amrex::The_Async_Arena());
    auto btflux= fab_btflux.array();

    //From ini_fields and .in file
    //fab_stflux.setVal(0.0);
    //also set btflux=0 (as in ana_btflux.H)

    ParallelFor(tbxp2,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        stflux(i,j,k)=0.0;
        btflux(i,j,k)=0.0;
    });

    //Use FC and DC as intermediate arrays for FX and FE
    //First pass do centered 2d terms

    if (solverChoice.flat_bathymetry) {
    ParallelFor(tbxp1,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        FX(i,j,k)=Box(tempold).contains(i-1,j,k) ? Huon(i,j,k)*
                    0.5*(tempold(i-1,j,k)+
                         tempold(i  ,j,k)) : 1e34;
    });
    ParallelFor(tbxp1,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        FE(i,j,k)=Box(tempold).contains(i,j-1,k) ? Hvom(i,j,k)*
                    0.5*(tempold(i,j-1,k)+
                         tempold(i,j,k)) : 1e34;
    });
    }
    else {
    ParallelFor(utbxp1,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        //should be t index 3
        FX(i,j,k)=tempold(i,j,k,nrhs)-tempold(i-1,j,k,nrhs);
    });
    ParallelFor(vtbxp1,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        //should be t index 3
        FE(i,j,k)=tempold(i,j,k,nrhs)-tempold(i,j-1,k,nrhs);
    });

    Real cffa=1.0/6.0;
    Real cffb=1.0/3.0;
    if (solverChoice.Hadv_scheme == AdvectionScheme::upstream3) {
        ParallelFor(tbxp1,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            //Upstream3
            curv(i,j,k)=-FX(i,j,k)+FX(i+1,j,k);
        });
        //HACK to avoid using the wrong index of t (using upstream3)
        ParallelFor(tbxp1,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real max_Huon=max(Huon(i,j,k),0.0); //FArrayBox(Huon).max<RunOn::Device>();
            Real min_Huon=min(Huon(i,j,k),0.0); //FArrayBox(Huon).min<RunOn::Device>();
            FX(i,j,k)=Huon(i,j,k)*0.5*(tempold(i,j,k)+tempold(i-1,j,k))-
                cffa*(curv(i,j,k)*min_Huon+ curv(i-1,j,k)*max_Huon);
        });
    }
    else if (solverChoice.Hadv_scheme == AdvectionScheme::centered4) {
        ParallelFor(tbxp1,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            //Centered4
            grad(i,j,k)=0.5*(FX(i,j,k)+FX(i+1,j,k));
        });
        ParallelFor(ubx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            FX(i,j,k)=Huon(i,j,k)*0.5*(tempold(i,j,k)+tempold(i-1,j,k)-
                                       cffb*(grad(i,j,k)-grad(i-1,j,k)));
        });
    }
    else {
        Error("Not a valid horizontal advection scheme");
    }
    if (solverChoice.Hadv_scheme == AdvectionScheme::upstream3) {
        ParallelFor(tbxp1,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            curv(i,j,k)=-FE(i,j,k)+FE(i,j+1,k);
        });
        //HACK to avoid using the wrong index of t (using upstream3)
        ParallelFor(tbxp1,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real max_Hvom=max(Hvom(i,j,k),0.0); //FArrayBox(Huon).max<RunOn::Device>();
            Real min_Hvom=min(Hvom(i,j,k),0.0); //FArrayBox(Huon).min<RunOn::Device>();
            FE(i,j,k)=Hvom(i,j,k)*0.5*(tempold(i,j,k)+tempold(i,j-1,k))-
                cffa*(curv(i,j,k)*min_Hvom+ curv(i,j-1,k)*max_Hvom);
        });
    }
    else if (solverChoice.Hadv_scheme == AdvectionScheme::centered4) {
        ParallelFor(tbxp1,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            grad(i,j,k)=0.5*(FE(i,j,k)+FE(i,j+1,k));
        });
        ParallelFor(vbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            FE(i,j,k)=Hvom(i,j,k)*0.5*(tempold(i,j,k)+tempold(i,j-1,k)-
                                       cffb*(grad(i,j,k)- grad(i,j-1,k)));
        });
    }
    else {
        Error("Not a valid horizontal advection scheme");
    }
    }

    //Intermediate tracer at 3
    //
    //  Time-step horizontal advection (m Tunits).
    //

    Real cff1 = 0.0, cff2 = 0.0, cff;

    Real GammaT = 1.0/6.0;

    if (iic==ntfirst)
    {
        cff=0.5*dt_lev;
        cff1=1.0;
        cff2=0.0;
    } else {
        cff=(1-GammaT)*dt_lev;
        cff1=0.5+GammaT;
        cff2=0.5-GammaT;
    }

    ParallelFor(tbx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        tempstore(i,j,k)=Hz(i,j,k)*(cff1*tempold(i,j,k)+
                                    cff2*tempcache(i,j,k))-
                                    cff*pm(i,j,0)*pn(i,j,0)*
                                    (FX(i+1,j,k)-FX(i,j,k)+
                                     FE(i,j+1,k)-FE(i,j,k));
         /*
         tempstore(i,j,k,3)=Hz(i,j,k)*(cff1*tempold(i,j,k,nstp)+
                                       cff2*temp(i,j,k,nnew))-
                            cff*pm(i,j,0)*pn(i,j,0)*
                            (FC(i+1,j)-FC(i,j)+
                             DC(i,j+1)-DC(i,j));*/
    });

    //
    // Time-step vertical advection of tracers (Tunits). Impose artificial
    // continuity equation.
    //
    ParallelFor(tbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
              //-----------------------------------------------------------------------
              //  Add in vertical advection.
              //-----------------------------------------------------------------------

              Real c1=0.5;
              Real c2=7.0/12.0;
              Real c3=1.0/12.0;

              if (k>=1 && k<=N-2)
              {
                      FC(i,j,k)=( c2*(tempold(i  ,j,k  ,nrhs)+ tempold(i,j,k+1,nrhs))
                                 -c3*(tempold(i  ,j,k-1,nrhs)+ tempold(i,j,k+2,nrhs)) )*
                                    ( W(i,j,k));
              }
              else // this needs to be split up so that the following can be concurrent
              {
                  FC(i,j,N)=0.0;

                  FC(i,j,N-1) = ( c2*tempold(i,j,N-1,nrhs)+ c1*tempold(i,j,N,nrhs)-c3*tempold(i,j,N-2,nrhs) )
                              * W(i,j,N-1);

                  FC(i,j,  0) = ( c2*tempold(i,j,  1,nrhs)+ c1*tempold(i,j,0,nrhs)-c3*tempold(i,j,2,nrhs) )
                              * W(i,j,0);
              }

    });

    ParallelFor(tbxp1,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        if(k-1>=0) {
        DC(i,j,k)=1.0/(Hz(i,j,k)-
                        cff*pm(i,j,0)*pn(i,j,0)*
                        (Huon(i+1,j,k)-Huon(i,j,k)+
                         Hvom(i,j+1,k)-Hvom(i,j,k)+
                         (W(i,j,k)-W(i,j,k-1))));
        } else {
        DC(i,j,k)=1.0/(Hz(i,j,k)-
                        cff*pm(i,j,0)*pn(i,j,0)*
                        (Huon(i+1,j,k)-Huon(i,j,k)+
                         Hvom(i,j+1,k)-Hvom(i,j,k)+
                         (W(i,j,k))));
        }

         /*
         tempstore(i,j,k,3)=Hz(i,j,k)*(cff1*tempold(i,j,k,nstp)+
                                       cff2*temp(i,j,k,nnew))-
                            cff*pm(i,j,0)*pn(i,j,0)*
                            (FC(i+1,j)-FC(i,j)+
                             DC(i,j+1)-DC(i,j));*/
    });

    ParallelFor(tbx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        Real c1 = cff*pm(i,j,0)*pn(i,j,0);

        Real c4 = (k>0) ? FC(i,j,k)-FC(i,j,k-1) : FC(i,j,k);

        tempstore(i,j,k) = DC(i,j,k)*(tempstore(i,j,k)-c1*c4);
    });

    //-----------------------------------------------------------------------
    //  Start computation of tracers at n+1 time-step, t(i,j,k,nnew,itrc).
    //-----------------------------------------------------------------------
    //
    //  Compute vertical diffusive fluxes "FC" of the tracer fields at
    //  ru is passed in here, but not actually used since ioff=0 and joff=0
    update_vel_3d(tbx, gbx, 0, 0, temp, tempold, ru, Hz, Akt, DC, FC,
                  stflux, btflux, z_r, pm, pn, iic, iic, nnew, nstp, nrhs, N, lambda, dt_lev);
}
