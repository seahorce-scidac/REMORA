#include <ROMSX.H>
#include <Utils.H>

using namespace amrex;

//
// rhs_3d
//
void
ROMSX::rhs_3d (const Box& bx,
               Array4<Real> uold  , Array4<Real> vold,
               Array4<Real> ru_arr, Array4<Real> rv_arr,
               Array4<Real> Huon, Array4<Real> Hvom,
               Array4<Real> W   , Array4<Real> FC_arr,
               int nrhs, int N)
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

    //
    // Scratch space
    //
    FArrayBox fab_Huee(gbx2,1,amrex::The_Async_Arena()); fab_Huee.setVal(0.0);
    FArrayBox fab_uee(gbx2,1,amrex::The_Async_Arena()); fab_uee.setVal(0.0);

    FArrayBox fab_Hvee(gbx2,1,amrex::The_Async_Arena()); fab_Hvee.setVal(0.0);
    FArrayBox fab_vee(gbx2,1,amrex::The_Async_Arena()); fab_vee.setVal(0.0);

    FArrayBox fab_Hvxx(gbx2,1,amrex::The_Async_Arena()); fab_Hvxx.setVal(0.0);
    FArrayBox fab_uxx(gbx2,1,amrex::The_Async_Arena()); fab_uxx.setVal(0.0);

    FArrayBox fab_Huxx(gbx2,1,amrex::The_Async_Arena()); fab_Huxx.setVal(0.0);
    FArrayBox fab_vxx(gbx2,1,amrex::The_Async_Arena()); fab_vxx.setVal(0.0);

    FArrayBox fab_UFx(gbx2,1,amrex::The_Async_Arena()); fab_UFx.setVal(0.0);
    FArrayBox fab_UFe(gbx2,1,amrex::The_Async_Arena()); fab_UFe.setVal(0.0);
    FArrayBox fab_VFx(gbx2,1,amrex::The_Async_Arena()); fab_VFx.setVal(0.0);
    FArrayBox fab_VFe(gbx2,1,amrex::The_Async_Arena()); fab_VFe.setVal(0.0);

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

    amrex::ParallelFor(gbx1,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        //should not include grow cells
        uxx(i,j,k)=uold(i-1,j,k,nrhs)-2.0*uold(i,j,k,nrhs)+uold(i+1,j,k,nrhs);

        //neglecting terms about periodicity since testing only periodic for now
        Huxx(i,j,k)=Huon(i-1,j,k)-2.0*Huon(i,j,k)+Huon(i+1,j,k);
    });

    amrex::ParallelFor(gbx1,
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

    amrex::ParallelFor(gbx1,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        /////////////MIGHT NEED NEW LOOP HERE
        //neglecting terms about periodicity since testing only periodic for now
        Hvxx(i,j,k)=Hvom(i-1,j,k)-2.0*Hvom(i,j,k)+Hvom(i+1,j,k);
    });

    amrex::ParallelFor(gbx1,
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
              cff=vxx(i-1,j,k);
            } else {
              cff=vxx(i,j,k);
            }
            VFx(i,j,k)=0.25*(cff1+Gadv*cff)* (cff2+Gadv*0.5*(Huee(i,j,k)+ Huee(i,j-1,k)));
            vee(i,j,k)=vold(i,j-1,k,nrhs)-2.0*vold(i,j,k,nrhs)+
              vold(i,j+1,k,nrhs);
            Hvee(i,j,k)=Hvom(i,j-1,k)-2.0*Hvom(i,j,k)+Hvom(i,j+1,k);
        });

        amrex::ParallelFor(gbx1,
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

        amrex::ParallelFor(gbx1,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
              //
              //  Add in horizontal advection.
              //

              Real cff1=UFx(i,j  ,k)-UFx(i-1,j,k);
              Real cff2=UFe(i,j+1,k)-UFe(i  ,j,k);
              Real cff=cff1+cff2;

              ru_arr(i,j,k,nrhs) -= cff;

              cff1=VFx(i+1,j,k)-VFx(i  ,j,k);
              cff2=VFe(i  ,j,k)-VFe(i,j-1,k);
              cff=cff1+cff2;
              rv_arr(i,j,k,nrhs) -= cff;

              //-----------------------------------------------------------------------
              //  Add in vertical advection.
              //-----------------------------------------------------------------------
              cff1=9.0/16.0;
              cff2=1.0/16.0;
              //              if(i>=0)
              if (k>=1 && k<=N-2)
              {
                      FC_arr(i,j,k)=( cff1*(uold(i  ,j,k  ,nrhs)+ uold(i,j,k+1,nrhs))
                                     -cff2*(uold(i  ,j,k-1,nrhs)+ uold(i,j,k+2,nrhs)) )*
                                    ( cff1*(   W(i  ,j,k)+ W(i-1,j,k))
                                     -cff2*(   W(i+1,j,k)+ W(i-2,j,k)) );
              }
              else // this needs to be split up so that the following can be concurrent
              {
                  FC_arr(i,j,N)=0.0;

                  FC_arr(i,j,N-1)=( cff1*(uold(i  ,j,N-1,nrhs)+ uold(i,j,N  ,nrhs))
                                   -cff2*(uold(i  ,j,N-2,nrhs)+ uold(i,j,N  ,nrhs)) )*
                                  ( cff1*(   W(i  ,j,N-1)+ W(i-1,j,N-1))
                                   -cff2*(   W(i+1,j,N-1)+ W(i-2,j,N-1)) );

                  FC_arr(i,j,0)=( cff1*(uold(i  ,j,1,nrhs)+ uold(i,j,2,nrhs))
                                 -cff2*(uold(i  ,j,1,nrhs)+ uold(i,j,3,nrhs)) )*
                                ( cff1*(   W(i  ,j,1)+ W(i-1,j,1))
                                 -cff2*(   W(i+1,j,1)+ W(i-2,j,1)) );

                  //              FC_arr(i,0,-1)=0.0;
              }

              if(k-1>=0) {
                  cff=FC_arr(i,j,k)-FC_arr(i,j,k-1);
              } else {
                  cff=FC_arr(i,j,k);
              }

              ru_arr(i,j,k,nrhs) -= cff;

              if (k>=1 && k<=N-2)
              {
                  FC_arr(i,j,k)=( cff1*(vold(i,j,k  ,nrhs)+ vold(i,j,k+1,nrhs))
                                 -cff2*(vold(i,j,k-1,nrhs)+ vold(i,j,k+2,nrhs)) )*
                                ( cff1*(W(i,j  ,k)+ W(i,j-1,k))
                                 -cff2*(W(i,j+1,k)+ W(i,j-2,k)) );
              }
              else // this needs to be split up so that the following can be concurent
              {
                  FC_arr(i,j,N)=0.0;
                  FC_arr(i,j,N-1)=( cff1*(vold(i,j,N-1,nrhs)+ vold(i,j,N  ,nrhs))
                                   -cff2*(vold(i,j,N-2,nrhs)+ vold(i,j,N  ,nrhs)) )*
                                  ( cff1*(W(i,j  ,N-1)+ W(i,j-1,N-1))
                                   -cff2*(W(i,j+1,N-1)+ W(i,j-2,N-1)) );

                  FC_arr(i,j,0)=( cff1*(vold(i,j,1,nrhs)+ vold(i,j,2,nrhs))
                                 -cff2*(vold(i,j,1,nrhs)+ vold(i,j,3,nrhs)) )*
                                ( cff1*(W(i,j  ,1)+ W(i,j-1,1))
                                 -cff2*(W(i,j+1,1)+ W(i,j-2,1)) );

                  //              FC_arr(i,0,-1)=0.0;
              }

              if(k-1>=0) {
                  cff=FC_arr(i,j,k)-FC_arr(i,j,k-1);
              } else {
                  cff=FC_arr(i,j,k);
              }
              rv_arr(i,j,k,nrhs) -= cff;

    });
}
