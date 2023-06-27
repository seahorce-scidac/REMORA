#include <ROMSX.H>
#include <Utils.H>

using namespace amrex;

//
// Start 3d step
//

void
ROMSX::advance_3d (int lev,
                   MultiFab& mf_u , MultiFab& mf_v ,
                   MultiFab& mf_tempold , MultiFab& mf_saltold ,
                   MultiFab& mf_temp , MultiFab& mf_salt ,
                   std::unique_ptr<MultiFab>& mf_tempstore,
                   std::unique_ptr<MultiFab>& mf_saltstore,
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
                   std::unique_ptr<MultiFab>& mf_Huon,
                   std::unique_ptr<MultiFab>& mf_Hvom,
                   const int ncomp, const int N, Real dt_lev)
{

    auto geomdata  = Geom(lev).data();
    const auto dxi = Geom(lev).InvCellSizeArray();

    const int Mm = Geom(lev).Domain().size()[1];

    const int nrhs  = ncomp-1;
    const int nnew  = ncomp-1;

    int iic = istep[lev];
    int ntfirst = 0;

    for ( MFIter mfi(mf_u, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<Real> const& u = mf_u.array(mfi);
        Array4<Real> const& v = mf_v.array(mfi);

        Array4<Real> const& tempold = (mf_tempold).array(mfi);
        Array4<Real> const& saltold = (mf_saltold).array(mfi);

        Array4<Real> const& temp = (mf_temp).array(mfi);
        Array4<Real> const& salt = (mf_salt).array(mfi);

        Array4<Real> const& tempstore = mf_tempstore->array(mfi);
        Array4<Real> const& saltstore = mf_saltstore->array(mfi);

        Array4<Real> const& ru = mf_ru->array(mfi);
        Array4<Real> const& rv = mf_rv->array(mfi);

        Array4<Real> const& AK = mf_AK.array(mfi);
        Array4<Real> const& DC = mf_DC.array(mfi);

        Array4<Real> const& Hzk = mf_Hzk.array(mfi);
        Array4<Real> const& Akv = mf_Akv->array(mfi);
        Array4<Real> const& Hz  = mf_Hz->array(mfi);

        Array4<Real> const& DU_avg1  = mf_DU_avg1->array(mfi);
        Array4<Real> const& DV_avg1  = mf_DV_avg1->array(mfi);

        Array4<Real> const& DU_avg2  = mf_DU_avg2->array(mfi);
        Array4<Real> const& DV_avg2  = mf_DV_avg2->array(mfi);

        Array4<Real> const& ubar = mf_ubar->array(mfi);
        Array4<Real> const& vbar = mf_vbar->array(mfi);

        Array4<Real> const& Huon = mf_Huon->array(mfi);
        Array4<Real> const& Hvom = mf_Hvom->array(mfi);

        Box bx = mfi.tilebox();
        Box gbx = mfi.growntilebox();
        //copy the tilebox
        Box gbx1 = bx;
        Box gbx11 = bx;
        Box gbx2 = bx;
        Box gbx21 = bx;
        //make only gbx be grown to match multifabs
        gbx21.grow(IntVect(NGROW,NGROW,NGROW-1));
        gbx2.grow(IntVect(NGROW,NGROW,0));
        gbx1.grow(IntVect(NGROW-1,NGROW-1,0));
        gbx11.grow(IntVect(NGROW-1,NGROW-1,NGROW-1));

        Box ubx = surroundingNodes(bx,0);
        Box vbx = surroundingNodes(bx,1);
        Box ubx2 = surroundingNodes(ubx,0);
        Box vbx2 = surroundingNodes(vbx,1);
        amrex::Print() << " BX " <<  bx << std::endl;
        amrex::Print() << "UBX " << ubx << std::endl;
        amrex::Print() << "VBX " << vbx << std::endl;

        FArrayBox fab_FC(gbx2,1,amrex::The_Arena());
        FArrayBox fab_BC(gbx2,1,amrex::The_Arena());
        FArrayBox fab_CF(gbx21,1,amrex::The_Arena());
        FArrayBox fab_oHz(gbx11,1,amrex::The_Arena());
        FArrayBox fab_pn(gbx2,1,amrex::The_Arena());
        FArrayBox fab_pm(gbx2,1,amrex::The_Arena());
        FArrayBox fab_fomn(gbx2,1,amrex::The_Arena());
        FArrayBox fab_Akt(gbx2,1,amrex::The_Arena());
        FArrayBox fab_W(gbx2,1,amrex::The_Arena());

        FArrayBox fab_on_u(gbx2,1,amrex::The_Arena());
        FArrayBox fab_om_v(gbx2,1,amrex::The_Arena());
        auto on_u=fab_on_u.array();
        auto om_v=fab_om_v.array();

        auto FC = fab_FC.array();
        auto BC = fab_BC.array();
        auto CF = fab_CF.array();
        auto oHz= fab_oHz.array();
        auto pn=fab_pn.array();
        auto pm=fab_pm.array();
        auto fomn=fab_fomn.array();
        auto Akt= fab_Akt.array();
        auto W= fab_W.array();
        //From ini_fields and .in file
        //fab_Akt.setVal(1e-6);
        //From ana_grid.h and metrics.F
        //
       //amrex::ParallelFor(gbx2,
       // [=] AMREX_GPU_DEVICE (int i, int j, int k)
       // {
       //     printf("%d %d %d %25.25g uvel begin advance3d\n", i,j,k,u(i,j,k));
       //     printf("%d %d %d %25.25g temp begin advance3d\n", i,j,k,temp(i,j,k));
       //     printf("%d %d %d %25.25g tempstore begin advance3d\n", i,j,k,tempstore(i,j,k));
       // });

        //
        // Update to u and v
        //
        amrex::ParallelFor(gbx2,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
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
            Akt(i,j,k)=1e-6;
        });
       amrex::ParallelFor(gbx2,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
          om_v(i,j,0)=1.0/dxi[0];
          on_u(i,j,0)=1.0/dxi[1];
        });

        Real cff;
        if (iic==ntfirst) {
          cff=0.25*dt_lev;
        } else if (iic==ntfirst+1) {
          cff=0.25*dt_lev*3.0/2.0;
        } else {
          cff=0.25*dt_lev*23.0/12.0;
        }

        amrex::ParallelFor(gbx2,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
            //printf("%d %d %d %d %25.25g %25.25g %25.25g u(k+1) u(k) ru(k)\n",i,j,k,0,u(i,j,k+1,nnew),u(i,j,k,nnew),ru(i,j,k,nrhs));
                u(i,j,k) += gbx2.contains(i-1,j,0) ? cff * (pm(i,j,0)+pm(i-1,j,0)) * (pn(i,j,0)+pn(i-1,j,0)) * ru(i,j,k,nrhs) : cff * (2.0 * pm(i,j,0)) * (2.0 * pn(i,j,0)) * ru(i,j,k,nrhs) ;
         //if (k+1 <= N)
            //printf("%d %d %d %d %25.25g %25.25g u(k+1) u(k)\n",i,j,k,0,u(i,j,k+1,nnew),u(i,j,k,nnew));
         v(i,j,k) += gbx2.contains(i,j-1,0) ? cff * (pm(i,j,0)+pm(i,j-1,0)) * (pn(i,j,0)+pn(i,j-1,0)) * rv(i,j,k,nrhs) : cff * (2.0 * pm(i,j,0)) * (2.0 * pn(i,j,0)) * rv(i,j,k,nrhs);

                //ifdef SPLINES_VVISC is true
                u(i,j,k) *= gbx2.contains(i-1,j,0) ? 2.0 / (Hz(i-1,j,k) + Hz(i,j,k)) :  1.0 / (Hz(i,j,k));
         //if (k+1 <= N)
            //printf("%d %d %d %d %25.25g %25.25g %25.25g %25.25g u(k+1) u(k) Hz(i-1) Hz(i) \n",i,j,k,0,u(i,j,k+1,nnew),u(i,j,k,nnew), Hz(i-1,j,k), Hz(i,j,k));
                //if(j>0&&j<Mm-1)
                v(i,j,k) *= gbx2.contains(i,j-1,0) ? 2.0 / (Hz(i,j-1,k) + Hz(i,j,k)) : 1.0 / (Hz(i,j,k));
            });
        // End previous
        {
        amrex::Gpu::synchronize();
        amrex::Gpu::LaunchSafeGuard lsg(false);
       //amrex::ParallelFor(gbx2,
       // [=] AMREX_GPU_DEVICE (int i, int j, int k)
       // {
       //     printf("%d %d %d %25.25g uvel after uv in advance3d\n", i,j,k,u(i,j,k));
       // });
       // NOTE: DC is only used as scratch in vert_visc_3d -- no need to pass or return a value
       // NOTE: may not actually need to set these to zero
        amrex::ParallelFor(gbx21,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                DC(i,j,k) = 0.0;
                CF(i,j,k) = 0.0;
            });
#ifdef AMREX_USE_GPU
    Gpu::synchronize();
    amrex::PrintToFile("DC3_cuda")<<FArrayBox(DC)<<std::endl;
#else
    amrex::PrintToFile("DC3_nocuda")<<FArrayBox(DC)<<std::endl;
#endif
       vert_visc_3d(gbx1,bx,1,0,u,Hz,Hzk,oHz,AK,Akv,BC,DC,FC,CF,nnew,N,dt_lev);

        amrex::ParallelFor(gbx21,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                DC(i,j,k) = 0.0;
                CF(i,j,k) = 0.0;
            });
       vert_visc_3d(gbx1,bx,0,1,v,Hz,Hzk,oHz,AK,Akv,BC,DC,FC,CF,nnew,N,dt_lev);
        }
        amrex::ParallelFor(gbx21,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                DC(i,j,k) = 0.0;
                CF(i,j,k) = 0.0;
            });
       vert_mean_3d(bx,1,0,u,Hz,Hzk,DU_avg1,oHz,Akv,BC,DC,FC,CF,pm,nnew,N,dt_lev);

        amrex::ParallelFor(gbx21,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                DC(i,j,k) = 0.0;
                CF(i,j,k) = 0.0;
            });
       vert_mean_3d(bx,0,1,v,Hz,Hzk,DV_avg1,oHz,Akv,BC,DC,FC,CF,pn,nnew,N,dt_lev);

    //amrex::ParallelFor(ubx,
    //[=] AMREX_GPU_DEVICE (int i, int j, int k)
    //{
        //printf("%d %d %d %25.25g Huon before massflux\n", i, j, k, Huon(i,j,k));
    //    });

       update_massflux_3d(Box(Huon),bx,1,0,u,Huon,Hz,on_u,DU_avg1,DU_avg2,DC,FC,CF,nnew);
        amrex::ParallelFor(Box(ubar),
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                ubar(i,j,k,0) = DC(i,j,-1)*DU_avg1(i,j,0);
                ubar(i,j,k,1) = ubar(i,j,k,0);
            });
        update_massflux_3d(Box(Hvom),bx,0,1,v,Hvom,Hz,om_v,DV_avg1,DV_avg2,DC,FC,CF,nnew);
        amrex::ParallelFor(Box(vbar),
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                vbar(i,j,k,0) = DC(i,j,-1)*DV_avg1(i,j,0);
                vbar(i,j,k,1) = vbar(i,j,k,0);
            });
    //amrex::ParallelFor(ubx,
    //[=] AMREX_GPU_DEVICE (int i, int j, int k)
    //{
    //    printf("%d %d %d %25.25g Huon after massflux\n", i, j, k, Huon(i,j,k));
    //    });
    bool test_functionality=true;
    if(test_functionality) {
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
    amrex::ParallelFor(gbx1,
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

    amrex::ParallelFor(gbx1,
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
       //
       //-----------------------------------------------------------------------
       // rhs_3d
       //-----------------------------------------------------------------------
       //

       //amrex::ParallelFor(gbx2,
       // [=] AMREX_GPU_DEVICE (int i, int j, int k)
       // {
       //     printf("%d %d %d %25.25g uvel before rhst3d advance3d\n", i,j,k,u(i,j,k));
       //     printf("%d %d %d %25.25g temp before rhst3d advance3d\n", i,j,k,temp(i,j,k));
       //     printf("%d %d %d %25.25g tempstore before rhst3d advance3d\n", i,j,k,tempstore(i,j,k));
       // });
       rhs_t_3d(bx, tempold, temp, tempstore, Huon, Hvom, oHz, pn, pm, W, FC, nrhs, nnew, N,dt_lev);

       rhs_t_3d(bx, saltold, salt, saltstore, Huon, Hvom, oHz, pn, pm, W, FC, nrhs, nnew, N,dt_lev);
       //Print()<<FArrayBox(salt)<<std::endl;
    }
    }

    mf_temp.FillBoundary(geom[lev].periodicity());
    mf_salt.FillBoundary(geom[lev].periodicity());

    for ( MFIter mfi(mf_u, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<Real> const& u = mf_u.array(mfi);
        Array4<Real> const& v = mf_v.array(mfi);

        Array4<Real> const& tempold = (mf_tempold).array(mfi);
        Array4<Real> const& saltold = (mf_saltold).array(mfi);

        Array4<Real> const& temp = (mf_temp).array(mfi);
        Array4<Real> const& salt = (mf_salt).array(mfi);

        Array4<Real> const& tempstore = mf_tempstore->array(mfi);
        Array4<Real> const& saltstore = mf_saltstore->array(mfi);

        Array4<Real> const& ru = mf_ru->array(mfi);
        Array4<Real> const& rv = mf_rv->array(mfi);

        Array4<Real> const& AK = mf_AK.array(mfi);
        Array4<Real> const& DC = mf_DC.array(mfi);

        Array4<Real> const& Hzk = mf_Hzk.array(mfi);
        Array4<Real> const& Akv = mf_Akv->array(mfi);
        Array4<Real> const& Hz  = mf_Hz->array(mfi);

        Array4<Real> const& DU_avg1  = mf_DU_avg1->array(mfi);
        Array4<Real> const& DV_avg1  = mf_DV_avg1->array(mfi);

        Array4<Real> const& DU_avg2  = mf_DU_avg2->array(mfi);
        Array4<Real> const& DV_avg2  = mf_DV_avg2->array(mfi);

        Array4<Real> const& Huon = mf_Huon->array(mfi);
        Array4<Real> const& Hvom = mf_Hvom->array(mfi);

        Box bx = mfi.tilebox();
        Box gbx = mfi.growntilebox();
        //copy the tilebox
        Box gbx1 = bx;
        Box gbx11 = bx;
        Box gbx2 = bx;
        Box gbx21 = bx;
        //make only gbx be grown to match multifabs
        gbx21.grow(IntVect(NGROW,NGROW,NGROW-1));
        gbx2.grow(IntVect(NGROW,NGROW,0));
        gbx1.grow(IntVect(NGROW-1,NGROW-1,0));
        gbx11.grow(IntVect(NGROW-1,NGROW-1,NGROW-1));

        Box ubx = surroundingNodes(bx,0);
        Box vbx = surroundingNodes(bx,1);
        Box ubx2 = surroundingNodes(ubx,0);
        Box vbx2 = surroundingNodes(vbx,1);
        amrex::Print() << " BX " <<  bx << std::endl;
        amrex::Print() << "UBX " << ubx << std::endl;
        amrex::Print() << "VBX " << vbx << std::endl;

        FArrayBox fab_FC(gbx2,1,amrex::The_Arena());
        FArrayBox fab_BC(gbx2,1,amrex::The_Arena());
        FArrayBox fab_CF(gbx21,1,amrex::The_Arena());
        FArrayBox fab_oHz(gbx11,1,amrex::The_Arena());
        FArrayBox fab_pn(gbx2,1,amrex::The_Arena());
        FArrayBox fab_pm(gbx2,1,amrex::The_Arena());
        FArrayBox fab_fomn(gbx2,1,amrex::The_Arena());
        FArrayBox fab_Akt(gbx2,1,amrex::The_Arena());
        FArrayBox fab_W(gbx2,1,amrex::The_Arena());

        FArrayBox fab_on_u(gbx2,1,amrex::The_Arena());
        FArrayBox fab_om_v(gbx2,1,amrex::The_Arena());
        auto on_u=fab_on_u.array();
        auto om_v=fab_om_v.array();

        auto FC = fab_FC.array();
        auto BC = fab_BC.array();
        auto CF = fab_CF.array();
        auto oHz= fab_oHz.array();
        auto pn=fab_pn.array();
        auto pm=fab_pm.array();
        auto fomn=fab_fomn.array();
        auto Akt= fab_Akt.array();
        auto W= fab_W.array();

        //From ini_fields and .in file
        //fab_Akt.setVal(1e-6);
        //From ana_grid.h and metrics.F

        //
        // Update to u and v
        //
        amrex::ParallelFor(gbx2,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
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
            Akt(i,j,k)=1e-6;
        });
       amrex::ParallelFor(gbx2,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
          om_v(i,j,0)=1.0/dxi[0];
          on_u(i,j,0)=1.0/dxi[1];
        });
       //amrex::ParallelFor(gbx2,
       // [=] AMREX_GPU_DEVICE (int i, int j, int k)
       // {
       //     printf("%d %d %d %25.25g uvel before visc advance3d\n", i,j,k,u(i,j,k));
       //     printf("%d %d %d %25.25g temp before visc advance3d\n", i,j,k,temp(i,j,k));
       //     printf("%d %d %d %25.25g tempstore before visc advance3d\n", i,j,k,tempstore(i,j,k));
       // });
       ////Print()<<FArrayBox(temp)<<std::endl;
       vert_visc_3d(gbx1,bx,0,0,temp,Hz,Hzk,oHz,AK,Akt,BC,DC,FC,CF,nnew,N,dt_lev);
       //Print()<<FArrayBox(temp)<<std::endl;
       //Print()<<FArrayBox(salt)<<std::endl;
       vert_visc_3d(gbx1,bx,0,0,salt,Hz,Hzk,oHz,AK,Akt,BC,DC,FC,CF,nnew,N,dt_lev);
       //Print()<<FArrayBox(salt)<<std::endl;
       //if(iic==ntfirst+2)
       //exit(1);
       //amrex::ParallelFor(gbx2,
       // [=] AMREX_GPU_DEVICE (int i, int j, int k)
       // {
       //     printf("%d %d %d %25.25g uvel end advance3d\n", i,j,k,u(i,j,k));
       //     printf("%d %d %d %25.25g temp end advance3d\n", i,j,k,temp(i,j,k));
       //     printf("%d %d %d %25.25g tempstore end advance3d\n", i,j,k,tempstore(i,j,k));
       // });

    } // MFiter
}
