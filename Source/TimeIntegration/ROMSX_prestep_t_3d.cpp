#include <ROMSX.H>
#include <Utils.H>

using namespace amrex;

//
// prestep_uv_3d
//

void
ROMSX::prestep_t_3d (const Box& bx,
                      Array4<Real> uold  , Array4<Real> vold,
                      Array4<Real> u_arr , Array4<Real> v_arr,
                      Array4<Real> tempold  , Array4<Real> salstold,
                      Array4<Real> temp_arr , Array4<Real> salt_arr,
                      Array4<Real> ru_arr, Array4<Real> rv_arr,
                      Array4<Real> Hz_arr, Array4<Real> Akv_arr,
                      Array4<Real> on_u, Array4<Real> om_v,
                      Array4<Real> Huon, Array4<Real> Hvom,
                      Array4<Real> pm_arr, Array4<Real> pn_arr,
                      Array4<Real> W   , Array4<Real> DC_arr,
                      Array4<Real> FC_arr  , Array4<Real> tempstore, Array4<Real> saltstore,
                      Array4<Real> FX_arr, Array4<Real> FE_arr,
                      Array4<Real> z_r_arr,
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
    //make only gbx be grown to match multifabs
    gbx2.grow(IntVect(2,2,0));
    gbx1.grow(IntVect(1,1,0));
    gbx11.grow(IntVect(1,1,1));
    FArrayBox fab_Akt(gbx2,1,amrex::The_Async_Arena());
    auto Akt_arr= fab_Akt.array();

    //From ini_fields and .in file
    //fab_Akt.setVal(1e-6);
    FArrayBox fab_stflux(gbx2,1,amrex::The_Async_Arena());
    auto stflux_arr= fab_stflux.array();

    //From ini_fields and .in file
    //fab_stflux.setVal(0.0);

    amrex::ParallelFor(gbx2,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        Akt_arr(i,j,k)=1e-6;
        stflux_arr(i,j,k)=0.0;
    });
    amrex::AllPrint() << "Box(Huon) " << Box(Huon) << std::endl;
    amrex::AllPrint() << "Box(Hvom) " << Box(Hvom) << std::endl;

    //Use FC and DC as intermediate arrays for FX and FE
    //First pass do centered 2d terms
    amrex::ParallelFor(ubx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        FX_arr(i,j,k)=Huon(i,j,k)*
                    0.5*(tempold(i-1,j,k)+
                         tempold(i  ,j,k));
    });
    amrex::ParallelFor(vbx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        FE_arr(i,j,k)=Hvom(i,j,k)*
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
        tempstore(i,j,k)=Hz_arr(i,j,k)*(cff1*tempold(i,j,k)+
                                       cff2*temp_arr(i,j,k))-
                        cff*pm_arr(i,j,0)*pn_arr(i,j,0)*
                        (FX_arr(i+1,j,k)-FX_arr(i,j,k)+
                         FE_arr(i,j+1,k)-FE_arr(i,j,k));
         /*
         tempstore(i,j,k,3)=Hz(i,j,k)*(cff1*tempold(i,j,k,nstp)+
                                       cff2*temp(i,j,k,nnew))-
                            cff*pm_arr(i,j,0)*pn_arr(i,j,0)*
                            (FC(i+1,j)-FC(i,j)+
                             DC(i,j+1)-DC(i,j));*/
    });
    //
    // Time-step vertical advection of tracers (Tunits). Impose artificial
    // continuity equation.
    //
    amrex::ParallelFor(Box(FC_arr),
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
                      FC_arr(i,j,k)=( cff2*(tempold(i  ,j,k  ,nrhs)+ tempold(i,j,k+1,nrhs))
                                     -cff3*(tempold(i  ,j,k-1,nrhs)+ tempold(i,j,k+2,nrhs)) )*
                                    ( W(i,j,k));
              }
              else // this needs to be split up so that the following can be concurrent
              {
                  FC_arr(i,j,N)=0.0;

                  FC_arr(i,j,N-1)=( cff2*tempold(i  ,j,N-1,nrhs)+ cff1*tempold(i,j,N  ,nrhs)
                                    -cff3*tempold(i  ,j,N-2,nrhs) )*
                                  ( W(i  ,j,N-1));

                  FC_arr(i,j,0)=( cff2*tempold(i  ,j,1,nrhs)+ cff1*tempold(i,j,0,nrhs)
                                 -cff3*tempold(i  ,j,2,nrhs) )*
                                ( W(i  ,j,0));

                  //              FC_arr(i,0,-1)=0.0;
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
    amrex::ParallelFor(gbx1,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        if(k-1>=0) {
        DC_arr(i,j,k)=1.0/(Hz_arr(i,j,k)-
                        cff*pm_arr(i,j,0)*pn_arr(i,j,0)*
                        (Huon(i+1,j,k)-Huon(i,j,k)+
                         Hvom(i,j+1,k)-Hvom(i,j,k)+
                         (W(i,j,k)-W(i,j,k-1))));
        } else {
        DC_arr(i,j,k)=1.0/(Hz_arr(i,j,k)-
                        cff*pm_arr(i,j,0)*pn_arr(i,j,0)*
                        (Huon(i+1,j,k)-Huon(i,j,k)+
                         Hvom(i,j+1,k)-Hvom(i,j,k)+
                         (W(i,j,k))));
        }

         /*
         tempstore(i,j,k,3)=Hz(i,j,k)*(cff1*tempold(i,j,k,nstp)+
                                       cff2*temp(i,j,k,nnew))-
                            cff*pm_arr(i,j,0)*pn_arr(i,j,0)*
                            (FC(i+1,j)-FC(i,j)+
                             DC(i,j+1)-DC(i,j));*/
    });
    //    Print()<<cff<<std::endl;
    //    exit(1);
    Print()<<FArrayBox(Hz_arr)<<std::endl;
    Print()<<FArrayBox(pm_arr)<<std::endl;
    Print()<<FArrayBox(pn_arr)<<std::endl;
    Print()<<FArrayBox(Huon)<<std::endl;
    Print()<<FArrayBox(Hvom)<<std::endl;
    Print()<<FArrayBox(W)<<std::endl;
    Print()<<FArrayBox(DC_arr)<<std::endl;
    Print()<<FArrayBox(uold)<<std::endl;
    Print()<<FArrayBox(u_arr)<<std::endl;
    Print()<<FArrayBox(tempold)<<std::endl;
    Print()<<FArrayBox(tempstore)<<std::endl;
    Print()<<FArrayBox(temp_arr)<<std::endl;
    //    exit(1);
    amrex::ParallelFor(gbx1,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        Real cff1=cff*pm_arr(i,j,0)*pn_arr(i,j,0);
        Real cff4;
        if(k-1>=0) {
            cff4=FC_arr(i,j,k)-FC_arr(i,j,k-1);
        } else {
            cff4=FC_arr(i,j,k);
        }
        tempstore(i,j,k)=DC_arr(i,j,k)*(tempstore(i,j,k)-cff1*cff4);
	temp_arr(i,j,k)=tempold(i,j,k);
    });
    Print()<<FArrayBox(tempold)<<std::endl;
    Print()<<FArrayBox(tempstore)<<std::endl;
    Print()<<FArrayBox(temp_arr)<<std::endl;

    //-----------------------------------------------------------------------
    //  Start computation of tracers at n+1 time-step, t(i,j,k,nnew,itrc).
    //-----------------------------------------------------------------------
    //
    //  Compute vertical diffusive fluxes "FC" of the tracer fields at
    update_vel_3d(gbx1, 0, 0, tempstore, temp_arr, ru_arr, Hz_arr, Akt_arr, DC_arr, FC_arr,
                  stflux_arr, z_r_arr, pm_arr, pn_arr, iic, iic, nnew, nstp, nrhs, N, lambda, dt_lev);
    Print()<<FArrayBox(tempold)<<std::endl;
    Print()<<FArrayBox(tempstore)<<std::endl;
    Print()<<FArrayBox(temp_arr)<<std::endl;
}
