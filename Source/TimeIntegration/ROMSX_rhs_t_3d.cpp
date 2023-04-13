#include <ROMSX.H>
#include <Utils.H>

using namespace amrex;

//
// rhs_3d
//

void
ROMSX::rhs_t_3d (const Box& bx,
                 Array4<Real> told  , Array4<Real> t,
                 Array4<Real> Huon, Array4<Real> Hvom,
                 Array4<Real> pn, Array4<Real> pm,
                 Array4<Real> W   , Array4<Real> FC_arr,
                 int nrhs, int nnew, int N, Real dt_lev)
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
    FArrayBox fab_grad(gbx2,1,amrex::The_Async_Arena()); //fab_grad.setVal(0.0);
    FArrayBox fab_uee(gbx2,1,amrex::The_Async_Arena()); //fab_uee.setVal(0.0);

    FArrayBox fab_uxx(gbx2,1,amrex::The_Async_Arena()); //fab_uxx.setVal(0.0);

    FArrayBox fab_curv(gbx2,1,amrex::The_Async_Arena()); //fab_curv.setVal(0.0);

    FArrayBox fab_FX(gbx2,1,amrex::The_Async_Arena()); //fab_FX.setVal(0.0);
    FArrayBox fab_FE(gbx2,1,amrex::The_Async_Arena()); //fab_FE.setVal(0.0);

    auto curv=fab_curv.array();
    auto grad=fab_grad.array();
    auto uxx=fab_uxx.array();
    auto uee=fab_uee.array();
    auto FX=fab_FX.array();
    auto FE=fab_FE.array();

    //check this////////////
    const Real Gadv = -0.25;

    amrex::ParallelFor(gbx2,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        grad(i,j,k)=0.0;
        uee(i,j,k)=0.0;

        curv(i,j,k)=0.0;
        uxx(i,j,k)=0.0;

        FX(i,j,k)=0.0;
        FE(i,j,k)=0.0;
    });

    amrex::ParallelFor(gbx1,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        //should be t index 3
        FX(i,j,k)=told(i,j,k,nrhs)-told(i-1,j,k,nrhs);
    });
    amrex::ParallelFor(gbx1,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        //Upstream3
        curv(i,j,k)=-FX(i,j,k)+FX(i+1,j,k);
        //Centered4
        grad(i,j,k)=0.5*(FX(i,j,k)+FX(i+1,j,k));
    });

    Real cffa=1.0/6.0;
    Real cffb=1.0/3.0;
    amrex::ParallelFor(gbx1,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
#if 0
        FX(i,j,k)=Huon(i,j,k)*0.5*(told(i,j,k)+told(i-1,j,k))+
                  cffa*(curv(i,j,k)*min(Huon)+ curv(i-1,j,k)*max(Huon));
#else
        FX(i,j,k)=Huon(i,j,k)*0.5*(told(i,j,k)+told(i-1,j,k))+
                  cffb*(grad(i,j,k)+ grad(i-1,j,k));
#endif
    });

    amrex::ParallelFor(gbx1,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        //should be t index 3
        FE(i,j,k)=told(i,j,k,nrhs)-told(i,j-1,k,nrhs);
    });
    amrex::ParallelFor(gbx1,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        //Upstream3
        curv(i,j,k)=-FE(i,j,k)+FE(i,j+1,k);
        //Centered4
        grad(i,j,k)=0.5*(FX(i,j,k)+FX(i,j+1,k));
    });

    cffa=1.0/6.0;
    cffb=1.0/3.0;
    amrex::ParallelFor(gbx1,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
#if 0
        FE(i,j,k)=Hvom(i,j,k)*0.5*(told(i,j,k)+told(i,j-1,k))+
                  cffa*(curv(i,j,k)*min(Hvom)+ curv(i,j-1,k)*max(Hvom));
#else
        FE(i,j,k)=Hvom(i,j,k)*0.5*(told(i,j,k)+told(i,j-1,k))+
                  cffb*(grad(i,j,k)+ grad(i,j-1,k));
#endif
    });

        amrex::ParallelFor(gbx1,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
              //
              //  Add in horizontal advection.
              //
              Real cff = dt_lev*pm(i,j,0)*pn(i,j,0);
              Real cff1=cff*(FX(i+1,j,k)-FX(i,j,k));
              Real cff2=cff*(FE(i,j+1,k)-FE(i,j,k));
              Real cff3=cff1+cff2;

              t(i,j,k,nnew) -= cff3;
#if 0
              //-----------------------------------------------------------------------
              //  Add in vertical advection.
              //-----------------------------------------------------------------------
              cff1=9.0/16.0;
              cff2=1.0/16.0;
              //              if(i>=0)
              if (k>=1 && k<=N-2)
              {
                      FC_arr(i,j,k)=( cff1*(told(i  ,j,k  ,nrhs)+ told(i,j,k+1,nrhs))
                                     -cff2*(told(i  ,j,k-1,nrhs)+ told(i,j,k+2,nrhs)) )*
                                    ( cff1*(   W(i  ,j,k)+ W(i-1,j,k))
                                     -cff2*(   W(i+1,j,k)+ W(i-2,j,k)) );
              }
              else // this needs to be split up so that the following can be concurrent
              {
                  FC_arr(i,j,N)=0.0;

                  FC_arr(i,j,N-1)=( cff1*(told(i  ,j,N-1,nrhs)+ told(i,j,N  ,nrhs))
                                   -cff2*(told(i  ,j,N-2,nrhs)+ told(i,j,N  ,nrhs)) )*
                                  ( cff1*(   W(i  ,j,N-1)+ W(i-1,j,N-1))
                                   -cff2*(   W(i+1,j,N-1)+ W(i-2,j,N-1)) );

                  FC_arr(i,j,0)=( cff1*(told(i  ,j,1,nrhs)+ told(i,j,2,nrhs))
                                 -cff2*(told(i  ,j,1,nrhs)+ told(i,j,3,nrhs)) )*
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
#endif
    });
}
