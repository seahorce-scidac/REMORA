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
                   std::unique_ptr<MultiFab>& /*mf_ubar*/,
                   std::unique_ptr<MultiFab>& /*mf_vbar*/,
                   MultiFab& mf_AK, MultiFab& mf_DC,
                   MultiFab& mf_Hzk,
                   std::unique_ptr<MultiFab>& mf_Akv,
                   std::unique_ptr<MultiFab>& mf_Hz,
                   std::unique_ptr<MultiFab>& mf_Huon,
                   std::unique_ptr<MultiFab>& mf_Hvom,
                   const int ncomp, const int N, Real dt_lev)
{
    // Need to include uv3dmix

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
        gbx21.grow(IntVect(2,2,1));
        gbx2.grow(IntVect(2,2,0));
        gbx1.grow(IntVect(1,1,0));
        gbx11.grow(IntVect(1,1,1));

        Box ubx = surroundingNodes(bx,0);
        Box vbx = surroundingNodes(bx,1);
	Box ubx2 = surroundingNodes(ubx,0);
        Box vbx2 = surroundingNodes(vbx,1);
        amrex::Print() << " BX " <<  bx << std::endl;
        amrex::Print() << "UBX " << ubx << std::endl;
        amrex::Print() << "VBX " << vbx << std::endl;

        FArrayBox fab_FC(gbx2,1,amrex::The_Async_Arena());
        FArrayBox fab_BC(gbx2,1,amrex::The_Async_Arena());
        FArrayBox fab_CF(gbx21,1,amrex::The_Async_Arena());
        FArrayBox fab_oHz(gbx11,1,amrex::The_Async_Arena());
        FArrayBox fab_pn(gbx2,1,amrex::The_Async_Arena());
        FArrayBox fab_pm(gbx2,1,amrex::The_Async_Arena());
        FArrayBox fab_fomn(gbx2,1,amrex::The_Async_Arena());
        FArrayBox fab_Akt(gbx2,1,amrex::The_Async_Arena());
        FArrayBox fab_W(gbx2,1,amrex::The_Async_Arena());

        FArrayBox fab_on_u(gbx2,1,amrex::The_Async_Arena());
        FArrayBox fab_om_v(gbx2,1,amrex::The_Async_Arena());
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

        Real cff;
        if (iic==ntfirst) {
          cff=0.25*dt_lev;
        } else if (iic==ntfirst+1) {
          cff=0.25*dt_lev*3.0/2.0;
        } else {
          cff=0.25*dt_lev*23.0/12.0;
        }

        amrex::ParallelFor(gbx1,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                u(i,j,k) += cff * (pm(i,j,0)+pm(i-1,j,0)) * (pn(i,j,0)+pn(i-1,j,0)) * ru(i,j,k,nrhs);
                v(i,j,k) += cff * (pm(i,j,0)+pm(i,j-1,0)) * (pn(i,j,0)+pn(i,j-1,0)) * rv(i,j,k,nrhs);

                //ifdef SPLINES_VVISC is true
                u(i,j,k) *= 2.0 / (Hz(i-1,j,k) + Hz(i,j,k));

                //if(j>0&&j<Mm-1)
                v(i,j,k) *= 2.0 / (Hz(i,j-1,k) + Hz(i,j,k));
            });
        // End previous

       // NOTE: DC is only used as scratch in vert_visc_3d -- no need to pass or return a value
       vert_visc_3d(ubx,1,0,u,Hz,Hzk,oHz,AK,Akv,BC,DC,FC,CF,nnew,N,dt_lev);

       vert_visc_3d(vbx,0,1,v,Hz,Hzk,oHz,AK,Akv,BC,DC,FC,CF,nnew,N,dt_lev);

#if 0
       mf_DC[mfi].setVal(0.0,gbx21);
       fab_CF.setVal(0.0,gbx21);
       vert_mean_3d(bx,1,0,u,Hz,Hzk,DU_avg1,oHz,Akv,BC,DC,FC,CF,pm,nnew,N,dt_lev);

       mf_DC[mfi].setVal(0.0,gbx21);
       fab_CF.setVal(0.0,gbx21);
       vert_mean_3d(bx,0,1,v,Hz,Hzk,DV_avg1,oHz,Akv,BC,DC,FC,CF,pn,nnew,N,dt_lev);
#endif
       /*
       /////////////////////////// doesn't match boxes well, want the whole of Huon updated
       update_massflux_3d(Box(Huon),1,0,u,Huon,Hz,on_u,DU_avg1,DU_avg2,DC,FC,CF,nnew);
       update_massflux_3d(Box(Hvom),0,1,v,Hvom,Hz,om_v,DV_avg1,DV_avg2,DC,FC,CF,nnew);
       */
       /////////////////////////// doesn't match boxes well, want the whole of Huon updated
       update_massflux_3d(ubx,1,0,u,Huon,Hz,on_u,DU_avg1,DU_avg2,DC,FC,CF,nnew);
       update_massflux_3d(vbx,0,1,v,Hvom,Hz,om_v,DV_avg1,DV_avg2,DC,FC,CF,nnew);
#if 0
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
    /*    Print()<<FArrayBox(W)<<std::endl;
    Print()<<FArrayBox(Hvom)<<std::endl;
    Print()<<FArrayBox(Huon)<<std::endl;*/
    amrex::ParallelFor(ubx,
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
    /*    Print()<<FArrayBox(W)<<std::endl;
    Print()<<FArrayBox(Hvom)<<std::endl;
    Print()<<FArrayBox(Huon)<<std::endl;*/
    amrex::ParallelFor(vbx,
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
    /*    Print()<<FArrayBox(W)<<std::endl;
    Print()<<FArrayBox(Hvom)<<std::endl;
    Print()<<FArrayBox(Huon)<<std::endl;*/
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
    /*    Print()<<FArrayBox(W)<<std::endl;
    Print()<<FArrayBox(Hvom)<<std::endl;
    Print()<<FArrayBox(Huon)<<std::endl;
    exit(1);*/
       //
       //-----------------------------------------------------------------------
       // rhs_3d
       //-----------------------------------------------------------------------
       //
    /*
       Print()<<FArrayBox(tempold)<<std::endl;
       Print()<<FArrayBox(temp)<<std::endl;
       Print()<<FArrayBox(tempstore)<<std::endl;
    */
       rhs_t_3d(bx, tempold, temp, tempstore, Huon, Hvom, pn, pm, W, FC, nrhs, nnew, N,dt_lev);
       /*
Print()<<FArrayBox(tempold)<<std::endl;
       Print()<<FArrayBox(temp)<<std::endl;
       Print()<<FArrayBox(tempstore)<<std::endl;
       Print()<<FArrayBox(salt)<<std::endl;
       */
       rhs_t_3d(bx, saltold, salt, saltstore, Huon, Hvom, pn, pm, W, FC, nrhs, nnew, N,dt_lev);
       //Print()<<FArrayBox(salt)<<std::endl;
#endif
       //Print()<<FArrayBox(temp)<<std::endl;
       vert_visc_3d(gbx1,0,0,temp,Hz,Hzk,oHz,AK,Akt,BC,DC,FC,CF,nnew,N,dt_lev);
       //Print()<<FArrayBox(temp)<<std::endl;
       //Print()<<FArrayBox(salt)<<std::endl;
       vert_visc_3d(gbx1,0,0,salt,Hz,Hzk,oHz,AK,Akt,BC,DC,FC,CF,nnew,N,dt_lev);
       //Print()<<FArrayBox(salt)<<std::endl;
       //if(iic==ntfirst+2)
       //exit(1);

    } // MFiter
}
