#include <ROMSX.H>

using namespace amrex;

/**
 * rhs_t_3d
 *
 * @param[in   ] bx
 * @param[in   ] xbx
 * @param[in   ] ybx
 * @param[inout] t
 * @param[in   ] tempstore
 * @param[in   ] Huon
 * @param[in   ] Hvom
 * @param[in   ] Hz
 * @param[in   ] oHz
 * @param[in   ] pn
 * @param[in   ] pm
 * @param[in   ] W
 * @param[inout] FC
 * @param[in   ] nrhs
 * @param[in   ] nnew
 * @param[in   ] N
 * @param[in   ] dt_lev
 */

void
ROMSX::rhs_t_3d (const Box& bx, const Box& gbx,
                 Array4<Real> t, Array4<Real> tempstore,
                 Array4<Real> Huon, Array4<Real> Hvom,
                 Array4<Real> Hz, Array4<Real> oHz,
                 Array4<Real> pn, Array4<Real> pm,
                 Array4<Real> W   , Array4<Real> FC,
                 int nrhs, int nnew, int N, Real dt_lev)
{
    //copy the tilebox
    Box tbxp1 = bx;
    Box tbxp2 = bx;

    //make only gbx be grown to match multifabs
    tbxp2.grow(IntVect(NGROW,NGROW,0));
    tbxp1.grow(IntVect(NGROW-1,NGROW-1,0));

    // Because grad, curv, FX, FE, are all local, do surroundinNodes
    Box utbxp1 = surroundingNodes(tbxp1, 0);
    Box vtbxp1 = surroundingNodes(tbxp1, 1);
    Box ubx = surroundingNodes(bx, 0);
    Box vbx = surroundingNodes(bx, 1);

    BoxArray ba_gbx1 = intersect(BoxArray(tbxp1),gbx);
    AMREX_ASSERT((ba_gbx1.size() == 1));

    //
    // Scratch space
    //
    FArrayBox fab_grad(tbxp2,1,amrex::The_Async_Arena()); //fab_grad.setVal(0.0);
    FArrayBox fab_curv(tbxp2,1,amrex::The_Async_Arena()); //fab_curv.setVal(0.0);

    FArrayBox fab_FX(tbxp2,1,amrex::The_Async_Arena()); //fab_FX.setVal(0.0);
    FArrayBox fab_FE(tbxp2,1,amrex::The_Async_Arena()); //fab_FE.setVal(0.0);

    auto curv=fab_curv.array();
    auto grad=fab_grad.array();

    auto FX=fab_FX.array();
    auto FE=fab_FE.array();

    fab_grad.template setVal<RunOn::Device>(0.);
    fab_curv.template setVal<RunOn::Device>(0.);

    fab_FX.template setVal<RunOn::Device>(0.);
    fab_FE.template setVal<RunOn::Device>(0.);

    ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        oHz(i,j,k) = 1.0/ Hz(i,j,k);
    });

    ParallelFor(utbxp1,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        //should be t index 3
        FX(i,j,k)=tempstore(i,j,k,nrhs)-tempstore(i-1,j,k,nrhs);
    });

    Real cffa=1.0/6.0;
    Real cffb=1.0/3.0;

    if(solverChoice.flat_bathymetry) {

        if (solverChoice.Hadv_scheme == AdvectionScheme::upstream3) {

            ParallelFor(tbxp1,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                //Upstream3
                curv(i,j,k)=-FX(i,j,k)+FX(i+1,j,k);
            });
            //HACK to avoid using the wrong index of t (using upstream3)
            Real max_Huon=FArrayBox(Huon).max<RunOn::Device>();
            Real min_Huon=FArrayBox(Huon).min<RunOn::Device>();
            ParallelFor(ubx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                FX(i,j,k)=Huon(i,j,k)*0.5*(tempstore(i,j,k)+tempstore(i-1,j,k))+
                    cffa*(curv(i,j,k)*min_Huon+ curv(i-1,j,k)*max_Huon);
            });

        } else if (solverChoice.Hadv_scheme == AdvectionScheme::centered4) {

            ParallelFor(tbxp1,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                //Centered4
                grad(i,j,k)=0.5*(FX(i,j,k)+FX(i+1,j,k));
            });
            ParallelFor(ubx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                FX(i,j,k)=Huon(i,j,k)*0.5*(tempstore(i,j,k)+tempstore(i-1,j,k))+
                    cffb*(grad(i,j,k)+ grad(i-1,j,k));
            });

        } else {
            Error("Not a valid horizontal advection scheme");
        }

    } else {

        if (solverChoice.Hadv_scheme == AdvectionScheme::upstream3) {
            ParallelFor(tbxp1,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
            //Upstream3
                curv(i,j,k)=-FX(i,j,k)+FX(i+1,j,k);
            });
            //HACK to avoid using the wrong index of t (using upstream3)
            ParallelFor(ubx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                Real max_Huon=max(Huon(i,j,k),0.0); //FArrayBox(Huon).max<RunOn::Device>();
                Real min_Huon=min(Huon(i,j,k),0.0); //FArrayBox(Huon).min<RunOn::Device>();
                FX(i,j,k)=Huon(i,j,k)*0.5*(tempstore(i,j,k)+tempstore(i-1,j,k))-
                    cffa*(curv(i,j,k)*min_Huon+ curv(i-1,j,k)*max_Huon);
            });

        } else if (solverChoice.Hadv_scheme == AdvectionScheme::centered4) {

            ParallelFor(tbxp1,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                //Centered4
                grad(i,j,k)=0.5*(FX(i,j,k)+FX(i+1,j,k));
            });
            ParallelFor(ubx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                FX(i,j,k)=Huon(i,j,k)*0.5*(tempstore(i,j,k)+tempstore(i-1,j,k)-
                                           cffb*(grad(i,j,k)- grad(i-1,j,k)));
            });

        } else {
            Error("Not a valid horizontal advection scheme");
        }
    } // flat bathymetry?

    ParallelFor(vtbxp1,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        //should be t index 3
        FE(i,j,k)=tempstore(i,j,k,nrhs)-tempstore(i,j-1,k,nrhs);
    });

    cffa=1.0/6.0;
    cffb=1.0/3.0;
    if (solverChoice.flat_bathymetry) {
        if (solverChoice.Hadv_scheme == AdvectionScheme::upstream3) {
            ParallelFor(tbxp1,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                curv(i,j,k)=-FE(i,j,k)+FE(i,j+1,k);
            });
            //HACK to avoid using the wrong index of t (using upstream3)
            Real max_Hvom=FArrayBox(Hvom).max<RunOn::Device>();
            Real min_Hvom=FArrayBox(Hvom).min<RunOn::Device>();
            ParallelFor(vbx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                FE(i,j,k)=Hvom(i,j,k)*0.5*(tempstore(i,j,k)+tempstore(i,j-1,k))+
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
                FE(i,j,k)=Hvom(i,j,k)*0.5*(tempstore(i,j,k)+tempstore(i,j-1,k))+
                    cffb*(grad(i,j,k)+ grad(i,j-1,k));
            });
        }
        else {
            Error("Not a valid horizontal advection scheme");
        }
    }
    else {
        if (solverChoice.Hadv_scheme == AdvectionScheme::upstream3) {
            ParallelFor(tbxp1,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                curv(i,j,k)=-FE(i,j,k)+FE(i,j+1,k);
            });
            //HACK to avoid using the wrong index of t (using upstream3)
            ParallelFor(vbx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                Real max_Hvom=max(Hvom(i,j,k),0.0); //FArrayBox(Huon).max<RunOn::Device>();
                Real min_Hvom=min(Hvom(i,j,k),0.0); //FArrayBox(Huon).min<RunOn::Device>();
                FE(i,j,k)=Hvom(i,j,k)*0.5*(tempstore(i,j,k)+tempstore(i,j-1,k))-
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
                FE(i,j,k)=Hvom(i,j,k)*0.5*(tempstore(i,j,k)+tempstore(i,j-1,k)-
                                           cffb*(grad(i,j,k)- grad(i,j-1,k)));
            });
        }
        else {
            Error("Not a valid horizontal advection scheme");
        }
    }

        ParallelFor(bx,
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
    ParallelFor(bx,
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
                                 -cff3*(tempstore(i  ,j,k-1)+ tempstore(i,j,k+2)) ) * ( W(i,j,k));

              } else {

                  FC(i,j,N)=0.0;

                  FC(i,j,N-1)=( cff2*tempstore(i  ,j,N-1)+ cff1*tempstore(i,j,N  )
                               -cff3*tempstore(i  ,j,N-2) ) * ( W(i  ,j,N-1));

                  FC(i,j,0)=( cff2*tempstore(i  ,j,1)+ cff1*tempstore(i,j,0)
                             -cff3*tempstore(i  ,j,2) ) * ( W(i  ,j,0));
              }

    });
    ParallelFor(bx,
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
    });
}
