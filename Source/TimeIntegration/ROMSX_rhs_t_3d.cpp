#include <ROMSX.H>

using namespace amrex;

/**
 * rhs_t_3d
 *
 * @param[in   ] gbx
 * @param[inout] t
 * @param[in   ] sstore
 * @param[in   ] Huon
 * @param[in   ] Hvom
 * @param[in   ] Hz
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
                 const Array4<Real      >& t,
                 const Array4<Real const>& sstore,
                 const Array4<Real const>& Huon,
                 const Array4<Real const>& Hvom,
                 const Array4<Real const>& Hz,
                 const Array4<Real const>& pn,
                 const Array4<Real const>& pm,
                 const Array4<Real const>& W ,
                 const Array4<Real      >& FC,
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
    FArrayBox fab_grad(tbxp2,1,amrex::The_Async_Arena());
    FArrayBox fab_curv(tbxp2,1,amrex::The_Async_Arena());

    FArrayBox fab_FX(tbxp2,1,amrex::The_Async_Arena());
    FArrayBox fab_FE(tbxp2,1,amrex::The_Async_Arena());

    auto curv=fab_curv.array();
    auto grad=fab_grad.array();

    auto FX=fab_FX.array();
    auto FE=fab_FE.array();

    fab_grad.template setVal<RunOn::Device>(0.);
    fab_curv.template setVal<RunOn::Device>(0.);

    fab_FX.template setVal<RunOn::Device>(0.);
    fab_FE.template setVal<RunOn::Device>(0.);

    ParallelFor(utbxp1, [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        //should be t index 3
        FX(i,j,k)=sstore(i,j,k,nrhs)-sstore(i-1,j,k,nrhs);
    });

    Real cffa=1.0/6.0;
    Real cffb=1.0/3.0;

    if(solverChoice.flat_bathymetry) {

        if (solverChoice.Hadv_scheme == AdvectionScheme::upstream3) {

            ParallelFor(tbxp1, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                //Upstream3
                curv(i,j,k)=-FX(i,j,k)+FX(i+1,j,k);
            });

            ParallelFor(ubx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                Real max_Huon = std::max(Huon(i,j,k),0.0);
                Real min_Huon = std::min(Huon(i,j,k),0.0);
                FX(i,j,k)=Huon(i,j,k)*0.5*(sstore(i,j,k)+sstore(i-1,j,k))+
                    cffa*(curv(i,j,k)*min_Huon+ curv(i-1,j,k)*max_Huon);
            });

        } else if (solverChoice.Hadv_scheme == AdvectionScheme::centered4) {

            ParallelFor(tbxp1, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                //Centered4
                grad(i,j,k)=0.5*(FX(i,j,k)+FX(i+1,j,k));
            });

            ParallelFor(ubx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                FX(i,j,k)=Huon(i,j,k)*0.5*(sstore(i,j,k)+sstore(i-1,j,k))+
                    cffb*(grad(i,j,k)+ grad(i-1,j,k));
            });

        } else {
            Error("Not a valid horizontal advection scheme");
        }

    } else {

        if (solverChoice.Hadv_scheme == AdvectionScheme::upstream3) {

            ParallelFor(tbxp1, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                //Upstream3
                curv(i,j,k)=-FX(i,j,k)+FX(i+1,j,k);
            });

            //HACK to avoid using the wrong index of t (using upstream3)
            ParallelFor(ubx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                Real max_Huon = std::max(Huon(i,j,k),0.0);
                Real min_Huon = std::min(Huon(i,j,k),0.0);
                FX(i,j,k)=Huon(i,j,k)*0.5*(sstore(i,j,k)+sstore(i-1,j,k))-
                    cffa*(curv(i,j,k)*min_Huon+ curv(i-1,j,k)*max_Huon);
            });

        } else if (solverChoice.Hadv_scheme == AdvectionScheme::centered4) {

            ParallelFor(tbxp1, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                //Centered4
                grad(i,j,k)=0.5*(FX(i,j,k)+FX(i+1,j,k));
            });

            ParallelFor(ubx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                FX(i,j,k)=Huon(i,j,k)*0.5*(sstore(i,j,k)+sstore(i-1,j,k)-
                                           cffb*(grad(i,j,k)- grad(i-1,j,k)));
            });

        } else {
            Error("Not a valid horizontal advection scheme");
        }
    } // flat bathymetry?

    ParallelFor(vtbxp1, [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        //should be t index 3
        FE(i,j,k)=sstore(i,j,k,nrhs)-sstore(i,j-1,k,nrhs);
    });

    cffa=1.0/6.0;
    cffb=1.0/3.0;
    if (solverChoice.flat_bathymetry) {

        if (solverChoice.Hadv_scheme == AdvectionScheme::upstream3) {

            ParallelFor(tbxp1, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                curv(i,j,k)=-FE(i,j,k)+FE(i,j+1,k);
            });

            ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                Real max_Hvom = std::max(Hvom(i,j,k),0.0);
                Real min_Hvom = std::min(Hvom(i,j,k),0.0);

                FE(i,j,k)=Hvom(i,j,k)*0.5*(sstore(i,j,k)+sstore(i,j-1,k))+
                    cffa*(curv(i,j,k)*min_Hvom+ curv(i,j-1,k)*max_Hvom);
            });

        } else if (solverChoice.Hadv_scheme == AdvectionScheme::centered4) {

            ParallelFor(tbxp1, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                grad(i,j,k)=0.5*(FE(i,j,k)+FE(i,j+1,k));
            });

            ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                FE(i,j,k)=Hvom(i,j,k)*0.5*(sstore(i,j,k)+sstore(i,j-1,k))+
                    cffb*(grad(i,j,k)+ grad(i,j-1,k));
            });

        } else {
            Error("Not a valid horizontal advection scheme");
        }

    } else {

        if (solverChoice.Hadv_scheme == AdvectionScheme::upstream3) {
            ParallelFor(tbxp1, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                curv(i,j,k)=-FE(i,j,k)+FE(i,j+1,k);
            });


            ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                Real max_Hvom = std::max(Hvom(i,j,k),0.0);
                Real min_Hvom = std::min(Hvom(i,j,k),0.0);

                FE(i,j,k)=Hvom(i,j,k)*0.5*(sstore(i,j,k)+sstore(i,j-1,k))-
                    cffa*(curv(i,j,k)*min_Hvom+ curv(i,j-1,k)*max_Hvom);
            });

        } else if (solverChoice.Hadv_scheme == AdvectionScheme::centered4) {

            ParallelFor(tbxp1, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                grad(i,j,k)=0.5*(FE(i,j,k)+FE(i,j+1,k));
            });

            ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                FE(i,j,k)=Hvom(i,j,k)*0.5*(sstore(i,j,k)+sstore(i,j-1,k)-
                                           cffb*(grad(i,j,k)- grad(i,j-1,k)));
            });

        } else {
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

    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
              //-----------------------------------------------------------------------
              //  Add in vertical advection.
              //-----------------------------------------------------------------------

              Real cff1=0.5;
              Real cff2=7.0/12.0;
              Real cff3=1.0/12.0;

              if (k>=1 && k<=N-2)
              {
                      FC(i,j,k)=( cff2*(sstore(i  ,j,k  )+ sstore(i,j,k+1))
                                 -cff3*(sstore(i  ,j,k-1)+ sstore(i,j,k+2)) ) * ( W(i,j,k));

              } else {

                  FC(i,j,N)=0.0;

                  FC(i,j,N-1)=( cff2*sstore(i  ,j,N-1)+ cff1*sstore(i,j,N  )
                               -cff3*sstore(i  ,j,N-2) ) * ( W(i  ,j,N-1));

                  FC(i,j,0)=( cff2*sstore(i  ,j,1)+ cff1*sstore(i,j,0)
                             -cff3*sstore(i  ,j,2) ) * ( W(i  ,j,0));
              }

    });

    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        Real cff1=dt_lev*pm(i,j,0)*pn(i,j,0);
        Real cff4;
        if(k-1>=0) {
            cff4=FC(i,j,k)-FC(i,j,k-1);
        } else {
            cff4=FC(i,j,k);
        }

        t(i,j,k) = (t(i,j,k)-cff1*cff4) / Hz(i,j,k);
    });
}
