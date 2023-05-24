#include <ROMSX.H>
#include <Utils.H>

using namespace amrex;

//
// rhs_3d
//

void
ROMSX::rhs_t_3d (const Box& bx,
                 Array4<Real> told  , Array4<Real> t, Array4<Real> tempstore,
                 Array4<Real> Huon, Array4<Real> Hvom,
                 Array4<Real> oHz,
                 Array4<Real> pn, Array4<Real> pm,
                 Array4<Real> W   , Array4<Real> FC,
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
    gbx2.grow(IntVect(NGROW,NGROW,0));
    gbx1.grow(IntVect(NGROW-1,NGROW-1,0));
    gbx11.grow(IntVect(NGROW-1,NGROW-1,NGROW-1));

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

    amrex::ParallelFor(ubx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        //should be t index 3
        FX(i,j,k)=tempstore(i,j,k,nrhs)-tempstore(i-1,j,k,nrhs);
    });
    amrex::ParallelFor(ubx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        //Upstream3
        curv(i,j,k)=-FX(i,j,k)+FX(i+1,j,k);
        //Centered4
        grad(i,j,k)=0.5*(FX(i,j,k)+FX(i+1,j,k));
    });

    Real cffa=1.0/6.0;
    Real cffb=1.0/3.0;
    //HACK to avoid using the wrong index of t (using upstream3)
    Real max_Huon=FArrayBox(Huon).max();
    Real min_Huon=FArrayBox(Huon).min();
    amrex::ParallelFor(ubx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
#if 1
        FX(i,j,k)=Huon(i,j,k)*0.5*(tempstore(i,j,k)+tempstore(i-1,j,k))+
                  cffa*(curv(i,j,k)*min_Huon+ curv(i-1,j,k)*max_Huon);
#else
        FX(i,j,k)=Huon(i,j,k)*0.5*(tempstore(i,j,k)+tempstore(i-1,j,k))+
                  cffb*(grad(i,j,k)+ grad(i-1,j,k));
#endif
    });

    amrex::ParallelFor(vbx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        //should be t index 3
        FE(i,j,k)=tempstore(i,j,k,nrhs)-tempstore(i,j-1,k,nrhs);
    });
    amrex::ParallelFor(vbx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        //Upstream3
        curv(i,j,k)=-FE(i,j,k)+FE(i,j+1,k);
        //Centered4
        grad(i,j,k)=0.5*(FE(i,j,k)+FE(i,j+1,k));
    });

    cffa=1.0/6.0;
    cffb=1.0/3.0;
    //HACK to avoid using the wrong index of t (using upstream3)
    Real max_Hvom=FArrayBox(Hvom).max();
    Real min_Hvom=FArrayBox(Hvom).min();
    amrex::ParallelFor(vbx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
#if 1
        FE(i,j,k)=Hvom(i,j,k)*0.5*(tempstore(i,j,k)+tempstore(i,j-1,k))+
                  cffa*(curv(i,j,k)*min_Hvom+ curv(i,j-1,k)*max_Hvom);
#else
        FE(i,j,k)=Hvom(i,j,k)*0.5*(tempstore(i,j,k)+tempstore(i,j-1,k))+
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
    });

        //-----------------------------------------------------------------------
        //  Time-step vertical advection term.
        //-----------------------------------------------------------------------
        //Check which type of differences:
        //
        //  Fourth-order, central differences vertical advective flux
        //  (Tunits m3/s).
        //
    amrex::ParallelFor(bx,
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
                      FC(i,j,k)=( cff2*(tempstore(i  ,j,k  )+ tempstore(i,j,k+1))
                                 -cff3*(tempstore(i  ,j,k-1)+ tempstore(i,j,k+2)) )*
                                    ( W(i,j,k));
              }
              else // this needs to be split up so that the following can be concurrent
              {
                  FC(i,j,N)=0.0;

                  FC(i,j,N-1)=( cff2*tempstore(i  ,j,N-1)+ cff1*tempstore(i,j,N  )
                               -cff3*tempstore(i  ,j,N-2) )*
                                  ( W(i  ,j,N-1));

                  FC(i,j,0)=( cff2*tempstore(i  ,j,1)+ cff1*tempstore(i,j,0)
                             -cff3*tempstore(i  ,j,2) )*
                                ( W(i  ,j,0));

                  //              FC(i,0,-1)=0.0;
              }

    });
    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        Real cff1=dt_lev*pm(i,j,0)*pn(i,j,0);
        Real cff4;
        if(k-1>=0) {
            cff4=FC(i,j,k)-FC(i,j,k-1);
        } else {
            cff4=FC(i,j,k);
        }
        t(i,j,k)=oHz(i,j,k)*(t(i,j,k)-cff1*cff4);
        //    if(i==2&&j==2&&k==2) {
        //    Print()<<i<<j<<k<<t(i,j,k)<<"\t"<<oHz(i,j,k)<<"\t"<<cff1<<"\t"<<cff4<<std::endl;
        //    Abort("any nans?");
        //}
    });
}
