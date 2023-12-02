#include <ROMSX.H>

using namespace amrex;

//
// Start 3d step
//

void
ROMSX::advance_3d (int lev,
                   MultiFab& mf_u , MultiFab& mf_v ,
                   MultiFab& mf_tempold, MultiFab& mf_saltold,
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
                   std::unique_ptr<MultiFab>& mf_Akt,
                   std::unique_ptr<MultiFab>& mf_Hz,
                   std::unique_ptr<MultiFab>& mf_Huon,
                   std::unique_ptr<MultiFab>& mf_Hvom,
                   std::unique_ptr<MultiFab>& mf_z_w,
                   std::unique_ptr<MultiFab>& mf_h,
                   const int ncomp, const int N, Real dt_lev)
{
    const int Mm = Geom(lev).Domain().size()[1];

    const int nrhs  = ncomp-1;
    const int nnew  = ncomp-1;

    const GpuArray<Real,AMREX_SPACEDIM> dx  = Geom(lev).CellSizeArray();
    const GpuArray<Real,AMREX_SPACEDIM> dxi = Geom(lev).InvCellSizeArray();

    auto const& prob_lo = Geom(lev).ProbLoArray();

    int iic = istep[lev];
    int ntfirst = 0;

    // Because zeta may have changed
    stretch_transform(lev);

    for ( MFIter mfi(mf_temp, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<Real> const& u = mf_u.array(mfi);
        Array4<Real> const& v = mf_v.array(mfi);

        Array4<Real> const& ru = mf_ru->array(mfi);
        Array4<Real> const& rv = mf_rv->array(mfi);

        Array4<Real> const& AK = mf_AK.array(mfi);
        Array4<Real> const& DC = mf_DC.array(mfi);

        Array4<Real> const& Hzk = mf_Hzk.array(mfi);
        Array4<Real> const& Akv = mf_Akv->array(mfi);

        Array4<Real const> const& Hz  = mf_Hz->const_array(mfi);

        Array4<Real const> const& DU_avg1  = mf_DU_avg1->const_array(mfi);
        Array4<Real const> const& DV_avg1  = mf_DV_avg1->const_array(mfi);

        Array4<Real> const& DU_avg2  = mf_DU_avg2->array(mfi);
        Array4<Real> const& DV_avg2  = mf_DV_avg2->array(mfi);

        Array4<Real> const& ubar = mf_ubar->array(mfi);
        Array4<Real> const& vbar = mf_vbar->array(mfi);

        Array4<Real> const& Huon = mf_Huon->array(mfi);
        Array4<Real> const& Hvom = mf_Hvom->array(mfi);

        Box bx = mfi.tilebox();
        Box gbx1 = mfi.growntilebox(IntVect(NGROW-1,NGROW-1,0));
        Box gbx2 = mfi.growntilebox(IntVect(NGROW,NGROW,0));
        Box gbx11 = mfi.growntilebox(IntVect(NGROW-1,NGROW-1,NGROW-1));
        Box gbx21 = mfi.growntilebox(IntVect(NGROW,NGROW,NGROW-1));

        Box gbx2D = gbx2;
        gbx2D.makeSlab(2,0);

        Box tbxp1 = bx;
        Box tbxp11 = bx;
        Box tbxp2 = bx;
        tbxp1.grow(IntVect(NGROW-1,NGROW-1,0));
        tbxp2.grow(IntVect(NGROW,NGROW,0));
        tbxp11.grow(IntVect(NGROW-1,NGROW-1,NGROW-1));

        FArrayBox fab_FC(gbx2,1,amrex::The_Async_Arena());
        FArrayBox fab_BC(gbx2,1,amrex::The_Async_Arena());
        FArrayBox fab_CF(gbx21,1,amrex::The_Async_Arena());
        FArrayBox fab_oHz(gbx11,1,amrex::The_Async_Arena());
        FArrayBox fab_pn(tbxp2,1,amrex::The_Async_Arena());
        FArrayBox fab_pm(tbxp2,1,amrex::The_Async_Arena());
        FArrayBox fab_fomn(tbxp2,1,amrex::The_Async_Arena());
        FArrayBox fab_W(tbxp2,1,amrex::The_Async_Arena());

        FArrayBox fab_on_u(tbxp2,1,amrex::The_Async_Arena());
        FArrayBox fab_om_v(tbxp2,1,amrex::The_Async_Arena());
        auto on_u = fab_on_u.array();
        auto om_v = fab_om_v.array();

        auto FC = fab_FC.array();
        auto BC = fab_BC.array();
        auto CF = fab_CF.array();
        auto oHz= fab_oHz.array();
        auto pn=fab_pn.array();
        auto pm=fab_pm.array();
        auto fomn=fab_fomn.array();

        //
        // Update to u and v
        //
        ParallelFor(amrex::makeSlab(tbxp2,2,0),
        [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
            pm(i,j,0)=dxi[0];
            pn(i,j,0)=dxi[1];

            //defined UPWELLING
            Real f0=-8.26e-5;
            Real beta=0.0;
            Real Esize=1000*(Mm);
            Real y = prob_lo[1] + (j + 0.5) * dx[1];
            Real f=f0+beta*(y-.5*Esize);
            fomn(i,j,0)=f*(1.0/(pm(i,j,0)*pn(i,j,0)));
        });

        fab_on_u.template setVal<RunOn::Device>(dx[1],amrex::makeSlab(tbxp2,2,0));
        fab_om_v.template setVal<RunOn::Device>(dx[0],amrex::makeSlab(tbxp2,2,0));

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
                u(i,j,k) += tbxp2.contains(i-1,j,0) ? cff * (pm(i,j,0)+pm(i-1,j,0)) * (pn(i,j,0)+pn(i-1,j,0)) * ru(i,j,k,nrhs) :
                                                      cff * (2.0 * pm(i,j,0)) * (2.0 * pn(i,j,0)) * ru(i,j,k,nrhs) ;
                v(i,j,k) += tbxp2.contains(i,j-1,0) ? cff * (pm(i,j,0)+pm(i,j-1,0)) * (pn(i,j,0)+pn(i,j-1,0)) * rv(i,j,k,nrhs) :
                                                      cff * (2.0 * pm(i,j,0)) * (2.0 * pn(i,j,0)) * rv(i,j,k,nrhs);
                u(i,j,k) *= tbxp2.contains(i-1,j,0) ? 2.0 / (Hz(i-1,j,k) + Hz(i,j,k)) :  1.0 / (Hz(i,j,k));
                v(i,j,k) *= tbxp2.contains(i,j-1,0) ? 2.0 / (Hz(i,j-1,k) + Hz(i,j,k)) : 1.0 / (Hz(i,j,k));
            });

        {
        amrex::Gpu::synchronize();
        amrex::Gpu::LaunchSafeGuard lsg(true);

        // NOTE: DC is only used as scratch in vert_visc_3d -- no need to pass or return a value
        // NOTE: may not actually need to set these to zero

        mf_DC[mfi].template setVal<RunOn::Device>(0.,gbx21);
        fab_CF.template setVal<RunOn::Device>(0.,gbx21);

#ifdef AMREX_USE_GPU
    Gpu::synchronize();
#endif
        vert_visc_3d(gbx1,1,0,u,Hz,Hzk,oHz,AK,Akv,BC,DC,FC,CF,nnew,N,dt_lev);

        // Reset to zero
        mf_DC[mfi].template setVal<RunOn::Device>(0.,gbx21);
        fab_CF.template setVal<RunOn::Device>(0.,gbx21);

        vert_visc_3d(gbx1,0,1,v,Hz,Hzk,oHz,AK,Akv,BC,DC,FC,CF,nnew,N,dt_lev);
        }

        // Reset to zero
        mf_DC[mfi].template setVal<RunOn::Device>(0.,gbx21);
        fab_CF.template setVal<RunOn::Device>(0.,gbx21);

        vert_mean_3d(gbx1,1,0,u,Hz,DU_avg1,DC,CF,pm,nnew,N);

        // Reset to zero
        mf_DC[mfi].template setVal<RunOn::Device>(0.,gbx21);
        fab_CF.template setVal<RunOn::Device>(0.,gbx21);

        vert_mean_3d(gbx1,0,1,v,Hz,DV_avg1,DC,CF,pn,nnew,N);

        update_massflux_3d(gbx2,1,0,u,Huon,Hz,on_u,DU_avg1,DU_avg2,DC,FC,CF,nnew);

        ParallelFor(gbx2D, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                ubar(i,j,k,0) = DC(i,j,-1)*DU_avg1(i,j,0);
                ubar(i,j,k,1) = ubar(i,j,k,0);
            });

        update_massflux_3d(gbx2,0,1,v,Hvom,Hz,om_v,DV_avg1,DV_avg2,DC,FC,CF,nnew);

        ParallelFor(gbx2D, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                vbar(i,j,k,0) = DC(i,j,-1)*DV_avg1(i,j,0);
                vbar(i,j,k,1) = vbar(i,j,k,0);
            });
    }

    mf_Huon->FillBoundary(geom[lev].periodicity());
    mf_Hvom->FillBoundary(geom[lev].periodicity());

    for ( MFIter mfi(mf_temp, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<Real> const& tempold = (mf_tempold).array(mfi);
        Array4<Real> const& saltold = (mf_saltold).array(mfi);

        Array4<Real> const& temp = (mf_temp).array(mfi);
        Array4<Real> const& salt = (mf_salt).array(mfi);

        Array4<Real> const& tempstore = mf_tempstore->array(mfi);
        Array4<Real> const& saltstore = mf_saltstore->array(mfi);

        Array4<Real> const& Hz  = mf_Hz->array(mfi);

        Array4<Real> const& Huon = mf_Huon->array(mfi);
        Array4<Real> const& Hvom = mf_Hvom->array(mfi);

        Array4<Real> const& z_w= (mf_z_w)->array(mfi);
        Array4<Real> const& h= (mf_h)->array(mfi);

        Box bx = mfi.tilebox();
        Box gbx = mfi.growntilebox();
        Box gbx1 = mfi.growntilebox(IntVect(NGROW-1,NGROW-1,0));
        Box gbx2 = mfi.growntilebox(IntVect(NGROW,NGROW,0));
        Box gbx11 = mfi.growntilebox(IntVect(NGROW-1,NGROW-1,NGROW-1));
        Box gbx21 = mfi.growntilebox(IntVect(NGROW,NGROW,NGROW-1));

        Box tbxp1 = bx;
        Box tbxp11 = bx;
        Box tbxp2 = bx;
        tbxp1.grow(IntVect(NGROW-1,NGROW-1,0));
        tbxp2.grow(IntVect(NGROW,NGROW,0));
        tbxp11.grow(IntVect(NGROW-1,NGROW-1,NGROW-1));

        FArrayBox fab_FC(gbx2,1,amrex::The_Async_Arena());
        FArrayBox fab_BC(gbx2,1,amrex::The_Async_Arena());
        FArrayBox fab_CF(gbx21,1,amrex::The_Async_Arena());
        FArrayBox fab_oHz(gbx11,1,amrex::The_Async_Arena());
        FArrayBox fab_pn(tbxp2,1,amrex::The_Async_Arena());
        FArrayBox fab_pm(tbxp2,1,amrex::The_Async_Arena());
        FArrayBox fab_fomn(tbxp2,1,amrex::The_Async_Arena());
        FArrayBox fab_W(tbxp2,1,amrex::The_Async_Arena());

        FArrayBox fab_on_u(tbxp2,1,amrex::The_Async_Arena());
        FArrayBox fab_om_v(tbxp2,1,amrex::The_Async_Arena());

        auto FC  = fab_FC.array();
        auto oHz = fab_oHz.array();
        auto pn  = fab_pn.array();
        auto pm  = fab_pm.array();
        auto fomn= fab_fomn.array();
        auto W= fab_W.array();

        // Update to u and v
        //
        ParallelFor(amrex::makeSlab(tbxp2,2,0),
        [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
            pm(i,j,0)=dxi[0];
            pn(i,j,0)=dxi[1];

            //defined UPWELLING
            Real f0=-8.26e-5;
            Real beta=0.0;
            Real Esize=1000*(Mm);
            Real y = prob_lo[1] + (j + 0.5) * dx[1];
            Real f=f0+beta*(y-.5*Esize);
            fomn(i,j,0)=f*(1.0/(pm(i,j,0)*pn(i,j,0)));
        });

        fab_on_u.template setVal<RunOn::Device>(dx[1],makeSlab(tbxp2,2,0));
        fab_om_v.template setVal<RunOn::Device>(dx[0],makeSlab(tbxp2,2,0));

        bool test_functionality=true;
        if (test_functionality) {
        //
        //------------------------------------------------------------------------
        //  Vertically integrate horizontal mass flux divergence.
        //------------------------------------------------------------------------
        //
        //Should really use gbx3uneven
        Box gbx1D = gbx1;
        gbx1D.makeSlab(2,0);

        ParallelFor(gbx1D,
        [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
        //  Starting with zero vertical velocity at the bottom, integrate
        //  from the bottom (k=0) to the free-surface (k=N).  The w(:,:,N(ng))
        //  contains the vertical velocity at the free-surface, d(zeta)/d(t).
        //  Notice that barotropic mass flux divergence is not used directly.
        //
        //W(i,j,-1)=0.0;
        int k=0;
        W(i,j,k) = - (Huon(i+1,j,k)-Huon(i,j,k)) - (Hvom(i,j+1,k)-Hvom(i,j,k));
        for(k=1;k<=N;k++) {
            W(i,j,k) = W(i,j,k-1) - (Huon(i+1,j,k)-Huon(i,j,k)) - (Hvom(i,j+1,k)-Hvom(i,j,k));
        }
    });

    ParallelFor(gbx1,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        //  Starting with zero vertical velocity at the bottom, integrate
        //  from the bottom (k=0) to the free-surface (k=N).  The w(:,:,N(ng))
        //  contains the vertical velocity at the free-surface, d(zeta)/d(t).
        //  Notice that barotropic mass flux divergence is not used directly.
        //
        Real wrk_i=W(i,j,N)/(z_w(i,j,N)+h(i,j,0,0));

        if(k!=N) {
            W(i,j,k) = W(i,j,k)- wrk_i*(z_w(i,j,k)+h(i,j,0,0));
        }
    });

    // probably not the most efficient way
    ParallelFor(gbx1,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        if (k == N) {
            W(i,j,N) = 0.0;
        }
    });

       //
       //-----------------------------------------------------------------------
       // rhs_t_3d
       //-----------------------------------------------------------------------
       //
       rhs_t_3d(bx, gbx, tempold, temp, tempstore, Huon, Hvom, Hz, oHz, pn, pm, W, FC, nrhs, nnew, N,dt_lev);
       rhs_t_3d(bx, gbx, saltold, salt, saltstore, Huon, Hvom, Hz, oHz, pn, pm, W, FC, nrhs, nnew, N,dt_lev);
    }

    } // mfi
    mf_temp.FillBoundary(geom[lev].periodicity());
    mf_salt.FillBoundary(geom[lev].periodicity());

    for ( MFIter mfi(mf_temp, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<Real> const& temp = (mf_temp).array(mfi);
        Array4<Real> const& salt = (mf_salt).array(mfi);

        Array4<Real> const& AK = mf_AK.array(mfi);
        Array4<Real> const& DC = mf_DC.array(mfi);

        Array4<Real> const& Hzk = mf_Hzk.array(mfi);
        Array4<Real      > const& Akt = mf_Akt->array(mfi);
        Array4<Real const> const& Hz  = mf_Hz->const_array(mfi);

        Box bx = mfi.tilebox();
        Box gbx1 = mfi.growntilebox(IntVect(NGROW-1,NGROW-1,0));

        // Copy the tilebox
        Box tbxp1 = bx;
        Box tbxp11 = bx;
        Box tbxp2 = bx;
        Box tbxp21 = bx;
        //make only gbx be grown to match multifabs
        tbxp21.grow(IntVect(NGROW,NGROW,NGROW-1));
        tbxp2.grow(IntVect(NGROW,NGROW,0));
        tbxp1.grow(IntVect(NGROW-1,NGROW-1,0));
        tbxp11.grow(IntVect(NGROW-1,NGROW-1,NGROW-1));

        FArrayBox fab_FC(tbxp2,1,amrex::The_Async_Arena());
        FArrayBox fab_BC(tbxp2,1,amrex::The_Async_Arena());
        FArrayBox fab_CF(tbxp21,1,amrex::The_Async_Arena());
        FArrayBox fab_oHz(tbxp11,1,amrex::The_Async_Arena());
        FArrayBox fab_pn(tbxp2,1,amrex::The_Async_Arena());
        FArrayBox fab_pm(tbxp2,1,amrex::The_Async_Arena());
        FArrayBox fab_fomn(tbxp2,1,amrex::The_Async_Arena());
        FArrayBox fab_W(tbxp2,1,amrex::The_Async_Arena());

        FArrayBox fab_on_u(tbxp2,1,amrex::The_Async_Arena());
        FArrayBox fab_om_v(tbxp2,1,amrex::The_Async_Arena());

        auto FC = fab_FC.array();
        auto BC = fab_BC.array();
        auto CF = fab_CF.array();
        auto oHz= fab_oHz.array();
        auto pn=fab_pn.array();
        auto pm=fab_pm.array();
        auto fomn=fab_fomn.array();

        //From ini_fields and .in file
        //fab_Akt.setVal(1e-6);
        //From ana_grid.h and metrics.F

        //
        // Update to u and v
        //
        ParallelFor(makeSlab(tbxp2,2,0),
        [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
            pm(i,j,0)=dxi[0];
            pn(i,j,0)=dxi[1];

            //defined UPWELLING
            Real f0=-8.26e-5;
            Real beta=0.0;
            Real Esize=1000*(Mm);
            Real y = prob_lo[1] + (j + 0.5) * dx[1];
            Real f=f0+beta*(y-.5*Esize);
            fomn(i,j,0)=f*(1.0/(pm(i,j,0)*pn(i,j,0)));
        });

        fab_on_u.template setVal<RunOn::Device>(dx[1],makeSlab(tbxp2,2,0));
        fab_om_v.template setVal<RunOn::Device>(dx[0],makeSlab(tbxp2,2,0));

        vert_visc_3d(gbx1,0,0,temp,Hz,Hzk,oHz,AK,Akt,BC,DC,FC,CF,nnew,N,dt_lev);
        vert_visc_3d(gbx1,0,0,salt,Hz,Hzk,oHz,AK,Akt,BC,DC,FC,CF,nnew,N,dt_lev);
    } // MFiter
}
