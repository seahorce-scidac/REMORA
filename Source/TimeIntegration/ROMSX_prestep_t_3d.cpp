#include <ROMSX.H>
#include <Utils.H>

using namespace amrex;

//
// prestep_uv_3d
//

void
ROMSX::prestep_t_3d (const Box& bx,
                      Array4<Real> uold  , Array4<Real> vold,
                      Array4<Real> u , Array4<Real> v,
                      Array4<Real> tempold  , Array4<Real> salstold,
                      Array4<Real> temp , Array4<Real> salt,
                      Array4<Real> ru, Array4<Real> rv,
                      Array4<Real> Hz, Array4<Real> Akv,
                      Array4<Real> on_u, Array4<Real> om_v,
                      Array4<Real> Huon, Array4<Real> Hvom,
                      Array4<Real> pm, Array4<Real> pn,
                      Array4<Real> W   , Array4<Real> DC,
                      Array4<Real> FC  , Array4<Real> tempstore, Array4<Real> saltstore,
                      Array4<Real> FX, Array4<Real> FE,
                      Array4<Real> z_r,
                      int iic, int ntfirst, int nnew, int nstp, int nrhs, int N,
                      Real lambda, Real dt_lev)
{
    //copy the tilebox
    Box gbx1 = bx;
    Box gbx11 = bx;
    Box gbx2 = bx;

    Box ubx = surroundingNodes(bx,0);
    Box vbx = surroundingNodes(bx,1);

    Box gbx3uneven(IntVect(AMREX_D_DECL(bx.smallEnd(0)-3,bx.smallEnd(1)-3,bx.smallEnd(2))),
                   IntVect(AMREX_D_DECL(bx.bigEnd(0)+2,bx.bigEnd(1)+2,bx.bigEnd(2))));
    Box gbx2uneven(IntVect(AMREX_D_DECL(bx.smallEnd(0)-2,bx.smallEnd(1)-2,bx.smallEnd(2))),
                   IntVect(AMREX_D_DECL(bx.bigEnd(0)+1,bx.bigEnd(1)+1,bx.bigEnd(2))));
    //
    //------------------------------------------------------------------------
    //  Vertically integrate horizontal mass flux divergence.
    //------------------------------------------------------------------------
    //
    //Should really use gbx3uneven
    amrex::ParallelFor(Box(W),
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        //  Starting with zero vertical velocity at the bottom, integrate
        //  from the bottom (k=0) to the free-surface (k=N).  The w(:,:,N(ng))
        //  contains the vertical velocity at the free-surface, d(zeta)/d(t).
        //  Notice that barotropic mass flux divergence is not used directly.
        //
        W(i,j,k)=0.0;
    });
    amrex::ParallelFor(ubx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        //  Starting with zero vertical velocity at the bottom, integrate
        //  from the bottom (k=0) to the free-surface (k=N).  The w(:,:,N(ng))
        //  contains the vertical velocity at the free-surface, d(zeta)/d(t).
        //  Notice that barotropic mass flux divergence is not used directly.
        //
        if(k==0) {
            W(i,j,k)=0.0;
        } else {
            W(i,j,k) = - (Huon(i+1,j,k)-Huon(i,j,k));
        }
    });
    amrex::ParallelFor(vbx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        //  Starting with zero vertical velocity at the bottom, integrate
        //  from the bottom (k=0) to the free-surface (k=N).  The w(:,:,N(ng))
        //  contains the vertical velocity at the free-surface, d(zeta)/d(t).
        //  Notice that barotropic mass flux divergence is not used directly.
        //
        if(k==0) {
            W(i,j,k)=W(i,j,k);
        } else {
            W(i,j,k) = W(i,j,k)- (Hvom(i,j+1,k)-Hvom(i,j,k));
        }
    });
    amrex::ParallelFor(Box(W),
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        //  Starting with zero vertical velocity at the bottom, integrate
        //  from the bottom (k=0) to the free-surface (k=N).  The w(:,:,N(ng))
        //  contains the vertical velocity at the free-surface, d(zeta)/d(t).
        //  Notice that barotropic mass flux divergence is not used directly.
        //
        if(k==0) {
            W(i,j,k)=W(i,j,k);
        } else {
            W(i,j,k) = W(i,j,k) + W(i,j,k-1);
        }
    });
    //make only gbx be grown to match multifabs
    gbx2.grow(IntVect(2,2,0));
    gbx1.grow(IntVect(1,1,0));
    gbx11.grow(IntVect(1,1,1));
    FArrayBox fab_Akt(gbx2,1,amrex::The_Async_Arena());
    auto Akt= fab_Akt.array();

    //From ini_fields and .in file
    //fab_Akt.setVal(1e-6);
    FArrayBox fab_stflux(gbx2,1,amrex::The_Async_Arena());
    auto stflux= fab_stflux.array();
    FArrayBox fab_btflux(gbx2,1,amrex::The_Async_Arena());
    auto btflux= fab_btflux.array();

    //From ini_fields and .in file
    //fab_stflux.setVal(0.0);
    //also set btflux=0 (as in ana_btflux.H)

    amrex::ParallelFor(gbx2,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        Akt(i,j,k)=1e-6;
        stflux(i,j,k)=0.0;
        btflux(i,j,k)=0.0;
    });
    amrex::AllPrint() << "Box(Huon) " << Box(Huon) << std::endl;
    amrex::AllPrint() << "Box(Hvom) " << Box(Hvom) << std::endl;

    //Use FC and DC as intermediate arrays for FX and FE
    //First pass do centered 2d terms
    Print()<<(Box(Huon))<<std::endl;
    Print()<<Box(ubx)<<std::endl;
    Print()<<Box(FX)<<std::endl;
    Print()<<(Box(tempold))<<std::endl;

    amrex::ParallelFor(Box(FX),
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        FX(i,j,k)=Huon(i,j,k)*
                    0.5*(tempold(i-1,j,k)+
                         tempold(i  ,j,k));
    });
    amrex::ParallelFor(Box(FE),
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        FE(i,j,k)=Hvom(i,j,k)*
                    0.5*(tempold(i,j-1,k)+
                         tempold(i,j,k));
    });

    //Intermediate tracer at 3
    //
    //  Time-step horizontal advection (m Tunits).
    //

    Real cff1 = 0.0, cff2 = 0.0, cff;

    int indx=0; //nrhs-3
    Real GammaT = 1.0/6.0;

    if (iic==ntfirst)
    {
        cff=0.5*dt_lev;
        cff1=1.0;
        cff2=0.0;
    } else {
        cff=(1-GammaT)*dt_lev;
        cff1=0.5+GammaT;
        cff2=0.5+GammaT;
    }
    amrex::ParallelFor(gbx1,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        tempstore(i,j,k)=Hz(i,j,k)*(cff1*tempold(i,j,k)+
                                    cff2*temp(i,j,k))-
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
    amrex::ParallelFor(Box(FC),
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
              //-----------------------------------------------------------------------
              //  Add in vertical advection.
              //-----------------------------------------------------------------------

              Real cff1=0.5;
              Real cff2=7.0/12.0;
              Real cff3=1.0/12.0;

              if (k>=1 && k<=N-2)
              {
                      FC(i,j,k)=( cff2*(tempold(i  ,j,k  ,nrhs)+ tempold(i,j,k+1,nrhs))
                                 -cff3*(tempold(i  ,j,k-1,nrhs)+ tempold(i,j,k+2,nrhs)) )*
                                    ( W(i,j,k));
              }
              else // this needs to be split up so that the following can be concurrent
              {
                  FC(i,j,N)=0.0;

                  FC(i,j,N-1)=( cff2*tempold(i  ,j,N-1,nrhs)+ cff1*tempold(i,j,N  ,nrhs)
                               -cff3*tempold(i  ,j,N-2,nrhs) )*
                                  ( W(i  ,j,N-1));

                  FC(i,j,0)=( cff2*tempold(i  ,j,1,nrhs)+ cff1*tempold(i,j,0,nrhs)
                             -cff3*tempold(i  ,j,2,nrhs) )*
                                ( W(i  ,j,0));

                  //              FC(i,0,-1)=0.0;
              }

    });
/*
    Real cff;

    int indx=0; //nrhs-3
    Real GammaT = 1.0/6.0;

    if (iic==ntfirst)
    {
        cff=0.5*dt_lev;
    } else {
        cff=(1-GammaT)*dt_lev;
        }*/
/*
    Print()<<"boxes gbx1 dc hz pm pn huon hvom w"<<std::endl;
    Print()<<gbx1<<std::endl;
    Print()<<Box(DC)<<std::endl;
    Print()<<Box(Hz)<<std::endl;
    Print()<<Box(pm)<<std::endl;
    Print()<<Box(pn)<<std::endl;
    Print()<<Box(Huon)<<std::endl;
    Print()<<Box(Hvom)<<std::endl;
    Print()<<Box(W)<<std::endl;
*/
    amrex::ParallelFor(gbx1,
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
    //    //Print()<<cff<<std::endl;
    //    exit(1);
    //Print()<<FArrayBox(Hz)<<std::endl;
    //Print()<<FArrayBox(pm)<<std::endl;
    //Print()<<FArrayBox(pn)<<std::endl;
    //Print()<<FArrayBox(Huon)<<std::endl;
    //Print()<<FArrayBox(Hvom)<<std::endl;
    //Print()<<FArrayBox(W)<<std::endl;
    //Print()<<FArrayBox(DC)<<std::endl;
    //Print()<<FArrayBox(uold)<<std::endl;
    //Print()<<FArrayBox(u)<<std::endl;
    //Print()<<FArrayBox(tempold)<<std::endl;
    //Print()<<FArrayBox(tempstore)<<std::endl;
    //Print()<<FArrayBox(temp)<<std::endl;
    //    exit(1);
    amrex::ParallelFor(gbx1,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        Real cff1=cff*pm(i,j,0)*pn(i,j,0);
        Real cff4;
        if(k-1>=0) {
            cff4=FC(i,j,k)-FC(i,j,k-1);
        } else {
            cff4=FC(i,j,k);
        }
        tempstore(i,j,k)=DC(i,j,k)*(tempstore(i,j,k)-cff1*cff4);
//      temp(i,j,k)=tempold(i,j,k);
    });

    //-----------------------------------------------------------------------
    //  Start computation of tracers at n+1 time-step, t(i,j,k,nnew,itrc).
    //-----------------------------------------------------------------------
    //
    //  Compute vertical diffusive fluxes "FC" of the tracer fields at
    update_vel_3d(gbx1, 0, 0, tempstore, temp, ru, Hz, Akt, DC, FC,
                  stflux, btflux, z_r, pm, pn, iic, iic, nnew, nstp, nrhs, N, lambda, dt_lev);
    //Print()<<FArrayBox(tempold)<<std::endl;
    //Print()<<FArrayBox(tempstore)<<std::endl;
    //Print()<<FArrayBox(temp)<<std::endl;
}
