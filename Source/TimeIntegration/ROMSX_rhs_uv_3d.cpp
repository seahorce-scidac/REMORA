#include <ROMSX.H>

using namespace amrex;

//
// rhs_3d
//

void
ROMSX::rhs_3d (const Box& bx, const Box& gbx,
               Array4<Real> uold  , Array4<Real> vold,
               Array4<Real> ru, Array4<Real> rv,
               Array4<Real> rufrc, Array4<Real> rvfrc,
               Array4<Real> sustr, Array4<Real> svstr,
               Array4<Real> bustr, Array4<Real> bvstr,
               Array4<Real> Huon, Array4<Real> Hvom,
               Array4<Real> on_u, Array4<Real> om_v,
               Array4<Real> om_u, Array4<Real> on_v,
               Array4<Real> W   , Array4<Real> FC,
               int nrhs, int N)
{
    if (verbose > 2) {
        amrex::PrintToFile("ru_begin_rhs3d").SetPrecision(18)<<FArrayBox(ru)<<std::endl;
        amrex::PrintToFile("rv_begin_rhs3d").SetPrecision(18)<<FArrayBox(rv)<<std::endl;
        amrex::PrintToFile("rufrc_begin_rhs3d").SetPrecision(18)<<FArrayBox(rufrc)<<std::endl;
        amrex::PrintToFile("rvfrc_begin_rhs3d").SetPrecision(18)<<FArrayBox(rvfrc)<<std::endl;
    }
    //copy the tilebox
    Box tbxp1 = bx;
    Box tbxp2 = bx;
    //make only gbx be grown to match multifabs
    tbxp1.grow(IntVect(NGROW-1,NGROW-1,0));
    tbxp2.grow(IntVect(NGROW,NGROW,0));

    BoxArray ba_gbx1 = intersect(BoxArray(tbxp1), gbx);
    AMREX_ASSERT((ba_gbx1.size() == 1));
    Box gbx1 = ba_gbx1[0];

    BoxArray ba_gbx2 = intersect(BoxArray(tbxp2), gbx);
    AMREX_ASSERT((ba_gbx2.size() == 1));
    Box gbx2 = ba_gbx2[0];

    Box ubx = surroundingNodes(bx,0);
    Box vbx = surroundingNodes(bx,1);


    Box bxD = bx;
    bxD.makeSlab(2,0);
    Box gbx1D = gbx1;
    gbx1D.makeSlab(2,0);

    //
    // Scratch space
    //
    FArrayBox fab_Huee(tbxp2,1,amrex::The_Async_Arena()); //fab_Huee.setVal(0.0);
    FArrayBox fab_uee(tbxp2,1,amrex::The_Async_Arena()); //fab_uee.setVal(0.0);

    FArrayBox fab_Hvee(tbxp2,1,amrex::The_Async_Arena()); //fab_Hvee.setVal(0.0);
    FArrayBox fab_vee(tbxp2,1,amrex::The_Async_Arena()); //fab_vee.setVal(0.0);

    FArrayBox fab_Hvxx(tbxp2,1,amrex::The_Async_Arena()); //fab_Hvxx.setVal(0.0);
    FArrayBox fab_uxx(tbxp2,1,amrex::The_Async_Arena()); //fab_uxx.setVal(0.0);

    FArrayBox fab_Huxx(tbxp2,1,amrex::The_Async_Arena()); //fab_Huxx.setVal(0.0);
    FArrayBox fab_vxx(tbxp2,1,amrex::The_Async_Arena()); //fab_vxx.setVal(0.0);

    FArrayBox fab_UFx(tbxp2,1,amrex::The_Async_Arena()); //fab_UFx.setVal(0.0);
    FArrayBox fab_UFe(tbxp2,1,amrex::The_Async_Arena()); //fab_UFe.setVal(0.0);
    FArrayBox fab_VFx(tbxp2,1,amrex::The_Async_Arena()); //fab_VFx.setVal(0.0);
    FArrayBox fab_VFe(tbxp2,1,amrex::The_Async_Arena()); //fab_VFe.setVal(0.0);

    auto Huxx=fab_Huxx.array();
    auto Hvxx=fab_Hvxx.array();
    auto Huee=fab_Huee.array();
    auto Hvee=fab_Hvee.array();
    auto uxx=fab_uxx.array();
    auto uee=fab_uee.array();
    auto vxx=fab_vxx.array();
    auto vee=fab_vee.array();
    auto UFx=fab_UFx.array();
    auto UFe=fab_UFe.array();
    auto VFx=fab_VFx.array();
    auto VFe=fab_VFe.array();

    //check this////////////
    const Real Gadv = -0.25;

    amrex::ParallelFor(tbxp2,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        Huee(i,j,k)=0.0;
        uee(i,j,k)=0.0;
        Hvee(i,j,k)=0.0;
        vee(i,j,k)=0.0;

        Huxx(i,j,k)=0.0;
        uxx(i,j,k)=0.0;
        Hvxx(i,j,k)=0.0;
        vxx(i,j,k)=0.0;

        UFx(i,j,k)=0.0;
        UFe(i,j,k)=0.0;
        VFx(i,j,k)=0.0;
        VFe(i,j,k)=0.0;
    });

    amrex::ParallelFor(tbxp1,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        //should not include grow cells
        uxx(i,j,k)=uold(i-1,j,k,nrhs)-2.0*uold(i,j,k,nrhs)+uold(i+1,j,k,nrhs);

        //neglecting terms about periodicity since testing only periodic for now
        Huxx(i,j,k)=Huon(i-1,j,k)-2.0*Huon(i,j,k)+Huon(i+1,j,k);
    });

    amrex::ParallelFor(tbxp1,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        Real cff;
        Real cff1=uold(i  ,j,k,nrhs)+uold(i+1,j,k,nrhs);

        if (cff1 > 0.0) {
          cff=uxx(i,j,k);
        } else {
          cff=uxx(i+1,j,k);
        }

        UFx(i,j,k)=0.25*(cff1+Gadv*cff) * (Huon(i,j,k)+ Huon(i+1,j,k)+
                         Gadv*0.5*(Huxx(i,j,k)+ Huxx(i+1,j,k)));

        //should not include grow cells
        uee(i,j,k)=uold(i,j-1,k,nrhs)-2.0*uold(i,j,k,nrhs)+uold(i,j+1,k,nrhs);
    });

    amrex::ParallelFor(tbxp1,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        /////////////MIGHT NEED NEW LOOP HERE
        //neglecting terms about periodicity since testing only periodic for now
        Hvxx(i,j,k)=Hvom(i-1,j,k)-2.0*Hvom(i,j,k)+Hvom(i+1,j,k);
    });

    amrex::ParallelFor(tbxp1,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
            Real cff;
            Real cff1=uold(i,j  ,k,nrhs)+uold(i,j-1,k,nrhs);
            Real cff2=Hvom(i,j,k)+Hvom(i-1,j,k);
            if (cff2>0.0) {
              cff=uee(i,j-1,k);
            } else {
              cff=uee(i,j,k);
            }

            UFe(i,j,k)=0.25*(cff1+Gadv*cff)*
              (cff2+Gadv*0.5*(Hvxx(i  ,j,k)+Hvxx(i-1,j,k)));

            vxx(i,j,k)=vold(i-1,j,k,nrhs)-2.0*vold(i,j,k,nrhs)+
              vold(i+1,j,k,nrhs);
            //neglecting terms about periodicity since testing only periodic for now
            Huee(i,j,k)=Huon(i,j-1,k)-2.0*Huon(i,j,k)+Huon(i,j+1,k);
            cff1=vold(i  ,j,k,nrhs)+vold(i-1,j,k,nrhs);
            cff2=Huon(i,j,k)+Huon(i,j-1,k);
            if (cff2>0.0) {
              cff = (i == tbxp1.smallEnd(0)) ? vxx(i-1,j,k) :
                  (vold(i-2,j,k,nrhs)-2.0*vold(i-1,j,k,nrhs)+vold(i,j,k,nrhs));
            } else {
              cff=vxx(i,j,k);
            }
            auto Huee_jm1 = (j == tbxp1.smallEnd(1)) ? Huee(i,j-1,k) :
                (Huon(i,j-2,k)-2.0*Huon(i,j-1,k)+Huon(i,j,k));
            VFx(i,j,k)=0.25*(cff1+Gadv*cff)* (cff2+Gadv*0.5*(Huee(i,j,k)+ Huee_jm1));
            vee(i,j,k)=vold(i,j-1,k,nrhs)-2.0*vold(i,j,k,nrhs)+
              vold(i,j+1,k,nrhs);
            Hvee(i,j,k)=Hvom(i,j-1,k)-2.0*Hvom(i,j,k)+Hvom(i,j+1,k);
        });

        amrex::ParallelFor(tbxp1,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            //neglecting terms about periodicity since testing only periodic for now
            Real cff;
            Real cff1=vold(i,j  ,k,nrhs)+vold(i,j+1,k,nrhs);
            if (cff1>0.0) {
              cff=vee(i,j,k);
            } else {
              cff=vee(i,j+1,k);
            }

            VFe(i,j,k) = 0.25 * (cff1+Gadv*cff) * (Hvom(i,j  ,k)+ Hvom(i,j+1,k) +
                         0.5  *            Gadv * (Hvee(i,j  ,k)+ Hvee(i,j+1,k)));
        });
        //This uses W being an extra grow cell sized
        amrex::ParallelFor(gbx1,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
              //
              //  Add in horizontal advection.
              //

              Real cff1=UFx(i,j  ,k)-UFx(i-1,j,k);
              Real cff2=UFe(i,j+1,k)-UFe(i  ,j,k);
              Real cff=cff1+cff2;

              ru(i,j,k,nrhs) -= cff;

              cff1=VFx(i+1,j,k)-VFx(i  ,j,k);
              cff2=VFe(i  ,j,k)-VFe(i,j-1,k);
              cff=cff1+cff2;
              rv(i,j,k,nrhs) -= cff;
        }); Gpu::synchronize();
        //This uses W being an extra grow cell sized
        amrex::ParallelFor(gbx1,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
              //-----------------------------------------------------------------------
              //  Add in vertical advection.
              //-----------------------------------------------------------------------
              Real cff1=9.0/16.0;
              Real cff2=1.0/16.0;
              //              if(i>=0)
              if (k>=1 && k<=N-2)
              {
                      FC(i,j,k)=( cff1*(uold(i  ,j,k  ,nrhs)+ uold(i,j,k+1,nrhs))
                                 -cff2*(uold(i  ,j,k-1,nrhs)+ uold(i,j,k+2,nrhs)) )*
                                    ( cff1*(   W(i  ,j,k)+ W(i-1,j,k))
                                     -cff2*(   W(i+1,j,k)+ W(i-2,j,k)) );
              }
              else // this needs to be split up so that the following can be concurrent
              {
                  FC(i,j,N)=0.0;

                  FC(i,j,N-1)=( cff1*(uold(i  ,j,N-1,nrhs)+ uold(i,j,N  ,nrhs))
                               -cff2*(uold(i  ,j,N-2,nrhs)+ uold(i,j,N  ,nrhs)) )*
                                  ( cff1*(   W(i  ,j,N-1)+ W(i-1,j,N-1))
                                   -cff2*(   W(i+1,j,N-1)+ W(i-2,j,N-1)) );

                  FC(i,j,0)=( cff1*(uold(i  ,j,0,nrhs)+ uold(i,j,1,nrhs))
                             -cff2*(uold(i  ,j,0,nrhs)+ uold(i,j,2,nrhs)) )*
                                ( cff1*(   W(i  ,j,0)+ W(i-1,j,0))
                                 -cff2*(   W(i+1,j,0)+ W(i-2,j,0)) );

                  //              FC(i,0,-1)=0.0;
              }
        }); Gpu::synchronize();
        //This uses W being an extra grow cell sized
        amrex::ParallelFor(gbx1,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
              Real cff;
              if(k-1>=0) {
                  cff=FC(i,j,k)-FC(i,j,k-1);
              } else {
                  cff=FC(i,j,k);
              }

              ru(i,j,k,nrhs) -= cff;
        }); Gpu::synchronize();
        //This uses W being an extra grow cell sized
        amrex::ParallelFor(gbx1,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
              Real cff1=9.0/16.0;
              Real cff2=1.0/16.0;
              if (k>=1 && k<=N-2)
              {
                  FC(i,j,k)=( cff1*(vold(i,j,k  ,nrhs)+ vold(i,j,k+1,nrhs))
                             -cff2*(vold(i,j,k-1,nrhs)+ vold(i,j,k+2,nrhs)) )*
                                ( cff1*(W(i,j  ,k)+ W(i,j-1,k))
                                 -cff2*(W(i,j+1,k)+ W(i,j-2,k)) );
              }
              else // this needs to be split up so that the following can be concurrent
              {
                  FC(i,j,N)=0.0;
                  FC(i,j,N-1)=( cff1*(vold(i,j,N-1,nrhs)+ vold(i,j,N  ,nrhs))
                               -cff2*(vold(i,j,N-2,nrhs)+ vold(i,j,N  ,nrhs)) )*
                                  ( cff1*(W(i,j  ,N-1)+ W(i,j-1,N-1))
                                   -cff2*(W(i,j+1,N-1)+ W(i,j-2,N-1)) );
                  FC(i,j,0)=( cff1*(vold(i,j,0,nrhs)+ vold(i,j,1,nrhs))
                                 -cff2*(vold(i,j,0,nrhs)+ vold(i,j,2,nrhs)) )*
                                ( cff1*(W(i,j  ,0)+ W(i,j-1,0))
                                 -cff2*(W(i,j+1,0)+ W(i,j-2,0)) );
                  //              FC(i,0,-1)=0.0;
              }
        }); Gpu::synchronize();
        //This uses W being an extra grow cell sized
        amrex::ParallelFor(gbx1,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
              Real cff1=9.0/16.0;
              Real cff2=1.0/16.0;
            Real cff;
              if(k-1>=0) {
                  cff=FC(i,j,k)-FC(i,j,k-1);
              } else {
                  cff=FC(i,j,k);
              }
              rv(i,j,k,nrhs) -= cff;
        }); Gpu::synchronize();
        //This uses W being an extra grow cell sized
        AMREX_ASSERT(gbx1.smallEnd(2) == 0 && gbx2.bigEnd(2) == N);
        amrex::ParallelFor(gbx1D,
        [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
           for (int k = 0; k <= N; ++k) {
              Real cff1=9.0/16.0;
              Real cff2=1.0/16.0;
              Real cff;
              //Recursive summation:
              rufrc(i,j,0) += ru(i,j,k,nrhs);
              rvfrc(i,j,0) += rv(i,j,k,nrhs);
// This toggles whether to upate forcing terms on slabbed box or not. Slabbing it changes plotfile to machine precision
#if 1
              //These forcing terms should possibly be updated on a slabbed box
              cff=om_u(i,j,0)*on_u(i,j,0);
              if(k==N) // this is consistent with update_vel_3d
                  cff1=sustr(i,j,0)*cff;
              else
                  cff1=0.0;
              if(k==0) //should this be k==-1?
                  cff2=-bustr(i,j,0)*cff;
              else
                  cff2=0.0;
              if (verbose > 2) {
                printf("%d %d %d  %15.15g %15.15g %15.15g  rufrc rhs3d\n", i,j,k, rufrc(i,j,0),cff1,cff2);
              }
              rufrc(i,j,0) += cff1+cff2;

              //These forcing terms should possibly be updated on a slabbed box
              cff=om_v(i,j,0)*on_v(i,j,0);
              if(k==N) // this is consistent with update_vel_3d
                  cff1=svstr(i,j,0)*cff;
              else
                  cff1=0.0;
              if(k==0) //should this be k==-1?
                  cff2=-bvstr(i,j,0)*cff;
              else
                  cff2=0.0;
              rvfrc(i,j,0) += cff1+cff2;
           }
#else
           }
        });
        amrex::ParallelFor(gbx1D,
        [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
              Real cff=om_u(i,j,0)*on_u(i,j,0);
              Real cff1=sustr(i,j,0)*cff;
              Real cff2=-bustr(i,j,0)*cff;

              rufrc(i,j,0) += cff1+cff2;

              cff=om_v(i,j,0)*on_v(i,j,0);
              cff1=svstr(i,j,0)*cff;
              cff2=-bvstr(i,j,0)*cff;

              rvfrc(i,j,0)+=cff1+cff2;
#endif
        });
        if (verbose > 2) {
        amrex::PrintToFile("rufrc_rhs3d").SetPrecision(18)<<FArrayBox(rufrc)<<std::endl;
        amrex::PrintToFile("rvfrc_rhs3d").SetPrecision(18)<<FArrayBox(rvfrc)<<std::endl;
        amrex::PrintToFile("sustr_rhs3d").SetPrecision(18)<<FArrayBox(sustr)<<std::endl;
        amrex::PrintToFile("svstr_rhs3d").SetPrecision(18)<<FArrayBox(svstr)<<std::endl;
        amrex::PrintToFile("bustr_rhs3d").SetPrecision(18)<<FArrayBox(bustr)<<std::endl;
        amrex::PrintToFile("bvstr_rhs3d").SetPrecision(18)<<FArrayBox(bvstr)<<std::endl;
        amrex::PrintToFile("ru_rhs3d").SetPrecision(18)<<FArrayBox(ru)<<std::endl;
        amrex::PrintToFile("rv_rhs3d").SetPrecision(18)<<FArrayBox(rv)<<std::endl;
        }
}
