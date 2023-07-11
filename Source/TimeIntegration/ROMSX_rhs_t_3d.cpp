#include <ROMSX.H>
#include <Utils.H>

using namespace amrex;

//
// rhs_3d
//

void
ROMSX::rhs_t_3d (const Box& bx, const Box& gbx,
                 Array4<Real> told  , Array4<Real> t, Array4<Real> tempstore,
                 Array4<Real> Huon, Array4<Real> Hvom,
                 Array4<Real> Hz, Array4<Real> oHz,
                 Array4<Real> pn, Array4<Real> pm,
                 Array4<Real> W   , Array4<Real> FC,
                 int nrhs, int nnew, int N, Real dt_lev)
{
    //copy the tilebox
    Box tbxp1 = bx;
    Box tbxp2 = bx;

    Box ubx = surroundingNodes(bx,0);
    Box vbx = surroundingNodes(bx,1);

    //make only gbx be grown to match multifabs
    tbxp2.grow(IntVect(NGROW,NGROW,0));
    tbxp1.grow(IntVect(NGROW-1,NGROW-1,0));

    BoxArray ba_gbx1 = intersect(BoxArray(tbxp1),gbx);
    AMREX_ASSERT((ba_gbx1.size() == 1));
    Box gbx1 = ba_gbx1[0];

    BoxArray ba_gbx2 = intersect(BoxArray(tbxp2),gbx);
    AMREX_ASSERT((ba_gbx2.size() == 1));
    Box gbx2 = ba_gbx2[0];

    //
    // Scratch space
    //
    FArrayBox fab_grad(tbxp2,1,amrex::The_Async_Arena()); //fab_grad.setVal(0.0);
    FArrayBox fab_uee(tbxp2,1,amrex::The_Async_Arena()); //fab_uee.setVal(0.0);

    FArrayBox fab_uxx(tbxp2,1,amrex::The_Async_Arena()); //fab_uxx.setVal(0.0);

    FArrayBox fab_curv(tbxp2,1,amrex::The_Async_Arena()); //fab_curv.setVal(0.0);

    FArrayBox fab_FX(tbxp2,1,amrex::The_Async_Arena()); //fab_FX.setVal(0.0);
    FArrayBox fab_FE(tbxp2,1,amrex::The_Async_Arena()); //fab_FE.setVal(0.0);

    auto curv=fab_curv.array();
    auto grad=fab_grad.array();
    auto uxx=fab_uxx.array();
    auto uee=fab_uee.array();
    auto FX=fab_FX.array();
    auto FE=fab_FE.array();

    //check this////////////
    const Real Gadv = -0.25;

    amrex::ParallelFor(tbxp2,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        grad(i,j,k)=0.0;
        uee(i,j,k)=0.0;

        curv(i,j,k)=0.0;
        uxx(i,j,k)=0.0;

        FX(i,j,k)=0.0;
        FE(i,j,k)=0.0;
    });
    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        oHz(i,j,k) = 1.0/ Hz(i,j,k);
    });

    if (verbose > 0) {
        Print() << "tempstore bx " << Box(tempstore) << std::endl;
    }
    amrex::ParallelFor(tbxp1,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        //should be t index 3
        FX(i,j,k)=tempstore(i,j,k,nrhs)-tempstore(i-1,j,k,nrhs);
        if ((verbose >= 2) && printinloop) {
            printf("FX ts2 %d %d %d %15.15g %15.15g\n",i,j,k,tempstore(i,j,k,nrhs),tempstore(i-1,j,k,nrhs));
        }
    });
    //PrintToFile("FX_init").SetPrecision(18) << FArrayBox(FX) << std::endl;
    Real cffa=1.0/6.0;
    Real cffb=1.0/3.0;
    if (solverChoice.Hadv_scheme == AdvectionScheme::upstream3) {
        amrex::ParallelFor(tbxp1,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            //Upstream3
            curv(i,j,k)=-FX(i,j,k)+FX(i+1,j,k);
        });
        //HACK to avoid using the wrong index of t (using upstream3)
        Real max_Huon=FArrayBox(Huon).max();
        Real min_Huon=FArrayBox(Huon).min();
        amrex::ParallelFor(tbxp1,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            FX(i,j,k)=Huon(i,j,k)*0.5*(tempstore(i,j,k)+tempstore(i-1,j,k))+
                      cffa*(curv(i,j,k)*min_Huon+ curv(i-1,j,k)*max_Huon);
        });
    }
    else if (solverChoice.Hadv_scheme == AdvectionScheme::centered4) {
        amrex::ParallelFor(tbxp1,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            //Centered4
            grad(i,j,k)=0.5*(FX(i,j,k)+FX(i+1,j,k));
            if ((verbose >= 2) && printinloop) {
                printf("grad FX %d %d %d %25.25g %25.25g %25.25g\n",i,j,k,grad(i,j,k),FX(i,j,k),FX(i+1,j,k));
            }
        });
        amrex::ParallelFor(tbxp1,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            FX(i,j,k)=Huon(i,j,k)*0.5*(tempstore(i,j,k)+tempstore(i-1,j,k))+
                      cffb*(grad(i,j,k)+ grad(i-1,j,k));
            if ((verbose >= 2) && printinloop) {
                printf("FX Huon ts2 grad2 %d %d %d %25.25g %25.25g %25.25g %25.25g %25.25g %25.25g\n",i,j,k,FX(i,j,k),Huon(i,j,k),tempstore(i,j,k),tempstore(i-1,j,k),grad(i,j,k),grad(i-1,j,k));
            }
        });
    }
    else {
        Error("Not a valid horizontal advection scheme");
    }
    if (verbose >= 2)
        PrintToFile("FX_set1").SetPrecision(18) << FArrayBox(FX) << std::endl;

    amrex::ParallelFor(tbxp1,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        //should be t index 3
        FE(i,j,k)=tempstore(i,j,k,nrhs)-tempstore(i,j-1,k,nrhs);
    });

    cffa=1.0/6.0;
    cffb=1.0/3.0;
    if (solverChoice.Hadv_scheme == AdvectionScheme::upstream3) {
        amrex::ParallelFor(tbxp1,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            curv(i,j,k)=-FE(i,j,k)+FE(i,j+1,k);
        });
        //HACK to avoid using the wrong index of t (using upstream3)
        Real max_Hvom=FArrayBox(Hvom).max();
        Real min_Hvom=FArrayBox(Hvom).min();
        amrex::ParallelFor(tbxp1,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            FE(i,j,k)=Hvom(i,j,k)*0.5*(tempstore(i,j,k)+tempstore(i,j-1,k))+
                      cffa*(curv(i,j,k)*min_Hvom+ curv(i,j-1,k)*max_Hvom);
        });
    }
    else if (solverChoice.Hadv_scheme == AdvectionScheme::centered4) {
        amrex::ParallelFor(tbxp1,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            grad(i,j,k)=0.5*(FE(i,j,k)+FE(i,j+1,k));
        });
        amrex::ParallelFor(tbxp1,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            FE(i,j,k)=Hvom(i,j,k)*0.5*(tempstore(i,j,k)+tempstore(i,j-1,k))+
                      cffb*(grad(i,j,k)+ grad(i,j-1,k));
        });
    }
    else {
        Error("Not a valid horizontal advection scheme");
    }

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
              if ((verbose >= 2) && printinloop)
                printf("update %d %d %d %15.15g %15.15g %15.15g %15.15g %15.15g\n",i,j,k,t(i,j,k,nnew),cff,cff1,cff2,cff3);

              t(i,j,k,nnew) -= cff3;
        });
        if (verbose >= 2) {
            PrintToFile("FX").SetPrecision(18) << FArrayBox(FX) << std::endl;
            PrintToFile("FE").SetPrecision(18) << FArrayBox(FE) << std::endl;
            PrintToFile("t_int").SetPrecision(18) << FArrayBox(t) << std::endl;
        }

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
        if ((verbose >= 2) && printinloop)
            printf("update2 %d %d %d %15.15g %15.15g %15.15g %15.15g\n",i,j,k,oHz(i,j,k), t(i,j,k), cff1, cff4);
        t(i,j,k)=oHz(i,j,k)*(t(i,j,k)-cff1*cff4);
        //    if(i==2&&j==2&&k==2) {
        //    Print()<<i<<j<<k<<t(i,j,k)<<"\t"<<oHz(i,j,k)<<"\t"<<cff1<<"\t"<<cff4<<std::endl;
        //    Abort("any nans?");
        //}
    });
}
