#include <ROMSX.H>
#include <Utils.H>

using namespace amrex;

//
// Start 3d step
//
void
ROMSX::advance_3d (int lev,
                   MultiFab& mf_u , MultiFab& mf_v ,
                   std::unique_ptr<MultiFab>& mf_ru,
                   std::unique_ptr<MultiFab>& mf_rv,
                   std::unique_ptr<MultiFab>& mf_DU_avg1,
                   std::unique_ptr<MultiFab>& mf_DU_avg2,
                   std::unique_ptr<MultiFab>& mf_DV_avg1,
                   std::unique_ptr<MultiFab>& mf_DV_avg2,
                   std::unique_ptr<MultiFab>& mf_ubar,
                   std::unique_ptr<MultiFab>& mf_vbar,
                   MultiFab& mf_AK, MultiFab& mf_DC,
                   MultiFab& mf_Hzk,
                   std::unique_ptr<MultiFab>& mf_Akv,
                   std::unique_ptr<MultiFab>& mf_Hz,
                   Real dt_lev)
{
    // Need to include uv3dmix

    auto geomdata  = Geom(lev).data();
    const auto dxi = Geom(lev).InvCellSizeArray();

    const int Mm = Geom(lev).Domain().size()[1];

    const int ncomp = 1;
    const int nrhs  = ncomp-1;
    const int nnew  = ncomp-1;

    int iic = istep[lev];
    int ntfirst = 0;

    auto N = Geom(lev).Domain().size()[2]-1; // Number of vertical "levs" aka, NZ

    for ( MFIter mfi(mf_u, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<Real> const& u = mf_u.array(mfi);
        Array4<Real> const& v = mf_v.array(mfi);

        Array4<Real> const& ru_arr = mf_ru->array(mfi);
        Array4<Real> const& rv_arr = mf_rv->array(mfi);

        Array4<Real> const& AK = mf_AK.array(mfi);
        Array4<Real> const& DC = mf_DC.array(mfi);

        Array4<Real> const& Hzk_arr = mf_Hzk.array(mfi);
        Array4<Real> const& Akv_arr = mf_Akv->array(mfi);
        Array4<Real> const& Hz_arr  = mf_Hz->array(mfi);

        Box bx = mfi.tilebox();
        //copy the tilebox
        Box gbx1 = bx;
        Box gbx11 = bx;
        Box gbx2 = bx;
        //make only gbx be grown to match multifabs
        gbx2.grow(IntVect(2,2,0));
        gbx1.grow(IntVect(1,1,0));
        gbx11.grow(IntVect(1,1,1));
        Box gbx=gbx2;

        Box ubx = surroundingNodes(bx,0);
        Box vbx = surroundingNodes(bx,1);
        amrex::Print() << " BX " <<  bx << std::endl;
        amrex::Print() << "UBX " << ubx << std::endl;
        amrex::Print() << "VBX " << vbx << std::endl;

        FArrayBox fab_FC(gbx2,1,amrex::The_Async_Arena);
        FArrayBox fab_BC(gbx2,1,amrex::The_Async_Arena);
        FArrayBox fab_CF(gbx2,1,amrex::The_Async_Arena);
        FArrayBox fab_oHz(gbx11,1,amrex::The_Async_Arena);
        FArrayBox fab_pn(gbx2,1,amrex::The_Async_Arena);
        FArrayBox fab_pm(gbx2,1,amrex::The_Async_Arena);
        FArrayBox fab_on_u(gbx2,1,amrex::The_Async_Arena);
        FArrayBox fab_om_v(gbx2,1,amrex::The_Async_Arena);
        FArrayBox fab_fomn(gbx2,1,amrex::The_Async_Arena);
#if 0
        FArrayBox fab_Huon(gbx2,1,amrex::The_Async_Arena);
        FArrayBox fab_Hvom(gbx2,1,amrex::The_Async_Arena);
        //rhs3d work arrays
        FArrayBox fab_Huxx(gbx2,1,amrex::The_Async_Arena);
        FArrayBox fab_Huee(gbx2,1,amrex::The_Async_Arena);
        FArrayBox fab_Hvxx(gbx2,1,amrex::The_Async_Arena);
        FArrayBox fab_Hvee(gbx2,1,amrex::The_Async_Arena);
        FArrayBox fab_uxx(gbx2,1,amrex::The_Async_Arena);
        FArrayBox fab_uee(gbx2,1,amrex::The_Async_Arena);
        FArrayBox fab_vxx(gbx2,1,amrex::The_Async_Arena);
        FArrayBox fab_vee(gbx2,1,amrex::The_Async_Arena);
        FArrayBox fab_UFx(gbx2,1,amrex::The_Async_Arena);
        FArrayBox fab_UFe(gbx2,1,amrex::The_Async_Arena);
        FArrayBox fab_VFx(gbx2,1,amrex::The_Async_Arena);
        FArrayBox fab_VFe(gbx2,1,amrex::The_Async_Arena);
#endif

        auto FC=fab_FC.array();
        auto BC=fab_BC.array();
        auto CF=fab_CF.array();
        auto oHz_arr=fab_oHz.array();
        auto pn=fab_pn.array();
        auto pm=fab_pm.array();
        auto on_u=fab_on_u.array();
        auto om_v=fab_om_v.array();
        auto fomn=fab_fomn.array();
#if 0
        auto Huon=fab_Huon.array();
        auto Hvom=fab_Hvom.array();
        auto Huxx=fab_Huxx.array();
        auto Huee=fab_Huee.array();
        auto Hvxx=fab_Hvxx.array();
        auto Hvee=fab_Hvee.array();
        auto uxx=fab_uxx.array();
        auto uee=fab_uee.array();
        auto vxx=fab_vxx.array();
        auto vee=fab_vee.array();
        auto UFx=fab_UFx.array();
        auto UFe=fab_UFe.array();
        auto VFx=fab_VFx.array();
        auto VFe=fab_VFe.array();
#endif
        //From ana_grid.h and metrics.F

        //
        // Update to u
        //
        amrex::ParallelFor(gbx2,
        [=] AMREX_GPU_DEVICE (int i, int j, int  )
        {

            const auto prob_lo         = geomdata.ProbLo();
            const auto dx              = geomdata.CellSize();

            pm(i,j,0)=dxi[0];
            pn(i,j,0)=dxi[1];
            //defined UPWELLING
            Real f0=-8.26e-5;
            Real beta=0.0;
            Real Esize=1000*(Mm);
            Real y = prob_lo[1] + (j + 0.5) * dx[1];
            Real f=fomn(i,j,0)=f0+beta*(y-.5*Esize);
            fomn(i,j,0)=f*(1.0/(pm(i,j,0)*pn(i,j,0)));
        });

        amrex::ParallelFor(gbx2,
        [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            om_v(i,j,0)=1.0/dxi[0];
            on_u(i,j,0)=1.0/dxi[1];
        });

        // Begin step3d_uv.F
        amrex::ParallelFor(gbx1,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            AK(i,j,k)=0.5*(Akv_arr(i-1,j,k)+Akv_arr(i  ,j,k));
        });

        amrex::ParallelFor(gbx1,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Hzk_arr(i,j,k)=0.5*(Hz_arr(i-1,j,k)+Hz_arr(i  ,j,k));
        });

        amrex::ParallelFor(gbx1,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            oHz_arr(i,j,k) = 1.0/Hzk_arr(i,j,k);
        });

        amrex::ParallelFor(gbx1,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                Real cff;
                if(iic==ntfirst) {
                  cff=0.25*dt_lev;
                } else if(iic==ntfirst+1) {
                  cff=0.25*dt_lev*3.0/2.0;
                } else {
                  cff=0.25*dt_lev*23.0/12.0;
                }

                DC(i,j,k)=cff*(pm(i,j,0)+pm(i-1,j,0))*(pn(i,j,0)+pn(i-1,j,0));

                u(i,j,k) += DC(i,j,k)*ru_arr(i,j,k,nrhs);
                v(i,j,k) += DC(i,j,k)*rv_arr(i,j,k,nrhs);

                //ifdef SPLINES_VVISC is true
                u(i,j,k) *= oHz_arr(i,j,k);

                //if(j>0&&j<Mm-1)
                v(i,j,k) *= oHz_arr(i,j,k);
            });
        // End previous

       // vertical viscosity term for tracers is identical except AK=Akt
       // Begin vertical viscosity term
       //should be gbx1, but need to fix some bounds inside this loop:
       amrex::ParallelFor(gbx1,
       [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                //
               //  Use conservative, parabolic spline reconstruction of vertical
                //  viscosity derivatives.  Then, time step vertical viscosity term
                //  implicitly by solving a tridiagonal system.
                //
                //  implicitly by solving a tridiagonal system.
                //
               Real cff;
               Real cff1=1.0/6.0;

               if(k>=1)
                   FC(i,j,k)=cff1*Hzk_arr(i,j,k  )-dt_lev*AK(i,j,k-1)*oHz_arr(i,j,k  );
               else
                   FC(i,j,k)=cff1*Hzk_arr(i,j,k  );
               if(k<=N-1)
                {
                    CF(i,j,k)=cff1*Hzk_arr(i,j,k+1)-dt_lev*AK(i,j,k+1)*oHz_arr(i,j,k+1);
                }

               //
               //  LU decomposition and forward substitution.
               //
               cff1=1.0/3.0;
               if(k==0)
               {
                   BC(i,j,k)=cff1*(Hzk_arr(i,j,k)+Hzk_arr(i,j,k+1))+
                       dt_lev*AK(i,j,k)*(oHz_arr(i,j,k)+oHz_arr(i,j,k+1));
                   cff=1.0/(BC(i,j,k)-FC(i,j,k)*0.0);
                   CF(i,j,k)=cff*CF(i,j,k);
                   DC(i,j,k)=cff*(u(i,j,k+1,nnew)-u(i,j,k,nnew)-
                                  FC(i,j,k)*0.0);
               }
               if(k+1<=N&&k>=1)
               {
                   BC(i,j,k)=cff1*(Hzk_arr(i,j,k)+Hzk_arr(i,j,k+1))+
                       dt_lev*AK(i,j,k)*(oHz_arr(i,j,k)+oHz_arr(i,j,k+1));
                   cff=1.0/(BC(i,j,k)-FC(i,j,k)*CF(i,j,k-1));
                   CF(i,j,k)=cff*CF(i,j,k);
                   DC(i,j,k)=cff*(u(i,j,k+1,nnew)-u(i,j,k,nnew)-
                                  FC(i,j,k)*DC(i,j,k-1));
               }

            });
       amrex::ParallelFor(gbx1,
       [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
               //
               //  Backward substitution.
               //
               DC(i,j,N)=0.0;

               if(N-k+1<=N&&N-k>=0) //-N,1,-1 => kidx =N-k+1
               {
                   //              amrex::Print()<<"index k: "<<k<<"corresponds to : "<<N-k<<"prev: "<<DC(i,j,N-k+1)<<std::endl;
                   if(N-k+1<0||N-k+2<0)
                       amrex::Abort("-1 here");
                   DC(i,j,N-k)=DC(i,j,N-k)-CF(i,j,N-k)*DC(i,j,N-k+1);
                   //amrex::Print()<<"index k: "<<k<<"corresponds to : "<<N-k<<"cur:  "<<DC(i,j,N-k)<<std::endl;
                   //              DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1);
               }
            });

       amrex::ParallelFor(gbx1,
       [=] AMREX_GPU_DEVICE (int i, int j, int k)
       {
           DC(i,j,k) *= AK(i,j,k);
       });

       amrex::ParallelFor(gbx1,
       [=] AMREX_GPU_DEVICE (int i, int j, int k)
       {
           Real cff;
           if(k-1>=0) {
               cff = dt_lev*oHz_arr(i,j,k)*(DC(i,j,k)-DC(i,j,k-1));
           } else {
               cff = dt_lev*oHz_arr(i,j,k)*(DC(i,j,k));
           }

           u(i,j,k) += cff;
        });

        amrex::ParallelFor(gbx1,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            AK(i,j,k)=0.5*(Akv_arr(i,j-1,k)+Akv_arr(i,j,k));
        });

        amrex::ParallelFor(gbx1,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Hzk_arr(i,j,k)=0.5*(Hz_arr(i,j-1,k)+Hz_arr(i,j,k));
        });

        amrex::ParallelFor(gbx1,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            oHz_arr(i,j,k) = 1.0/Hzk_arr(i,j,k);
        });

       //
       // Begin vertical velocity term
       //
       amrex::ParallelFor(gbx1, ncomp,
       [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
       {
                //
               //  Use conservative, parabolic spline reconstruction of vertical
                //  viscosity derivatives.  Then, time step vertical viscosity term
                //  implicitly by solving a tridiagonal system.
                //
                //  implicitly by solving a tridiagonal system.
                //
               Real cff;
               Real cff1=1.0/6.0;

               if(k>=1)
                   FC(i,j,k)=cff1*Hzk_arr(i,j,k  )-dt_lev*AK(i,j,k-1)*oHz_arr(i,j,k  );
               else
                   FC(i,j,k)=cff1*Hzk_arr(i,j,k  );
               if(k<=N-1)
                {
                    CF(i,j,k)=cff1*Hzk_arr(i,j,k+1)-dt_lev*AK(i,j,k+1)*oHz_arr(i,j,k+1);
                }

               //
               //  LU decomposition and forward substitution.
               //
               cff1=1.0/3.0;
               if(k==0)
               {
                   BC(i,j,k)=cff1*(Hzk_arr(i,j,k)+Hzk_arr(i,j,k+1))+
                       dt_lev*AK(i,j,k)*(oHz_arr(i,j,k)+oHz_arr(i,j,k+1));
                   cff=1.0/(BC(i,j,k)-FC(i,j,k)*0.0);
                   CF(i,j,k)=cff*CF(i,j,k);
                   DC(i,j,k)=cff*(v(i,j,k+1,nnew)-v(i,j,k,nnew)-
                                  FC(i,j,k)*0.0);
               }
               if(k+1<=N&&k>=1)
               {
                   BC(i,j,k)=cff1*(Hzk_arr(i,j,k)+Hzk_arr(i,j,k+1))+
                       dt_lev*AK(i,j,k)*(oHz_arr(i,j,k)+oHz_arr(i,j,k+1));
                   cff=1.0/(BC(i,j,k)-FC(i,j,k)*CF(i,j,k-1));
                   CF(i,j,k)=cff*CF(i,j,k);
                   DC(i,j,k)=cff*(v(i,j,k+1,nnew)-v(i,j,k,nnew)-
                                  FC(i,j,k)*DC(i,j,k-1));
               }

       });
       amrex::ParallelFor(vbx,
       [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
               //
               //  Backward substitution.
               //
               DC(i,j,N)=0.0;

               if(N-k+1<=N&&N-k>=0) //-N,1,-1 => kidx =N-k+1
               {
                   if(N-k+1<0||N-k+2<0) amrex::Abort("-1 here");
                   DC(i,j,N-k) -= CF(i,j,N-k)*DC(i,j,N-k+1);
               }
              //   DC(i,k) -= CF(i,k)*DC(i,k+1);
            });

       amrex::ParallelFor(vbx,
       [=] AMREX_GPU_DEVICE (int i, int j, int k)
       {
           DC(i,j,k) *= AK(i,j,k);
       });

       amrex::ParallelFor(vbx,
       [=] AMREX_GPU_DEVICE (int i, int j, int k)
       {
           Real cff;
           if (k-1>=0) {
               cff=dt_lev*oHz_arr(i,j,k)*(DC(i,j,k)-DC(i,j,k-1));
           } else {
               cff=dt_lev*oHz_arr(i,j,k)*(DC(i,j,k));
           }
           //if(j>0&&j<Mm-1)
           v(i,j,k) += cff;
       });
    } // MFiter
}
