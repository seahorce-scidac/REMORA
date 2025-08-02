#include <REMORA.H>

using namespace amrex;

/**
 * @param[in   ] lev     level to operate on
 * @param[in   ] time    time at start of step
 * @param[in   ] dt_lev  time step at level
 */
void
REMORA::setup_step (int lev, Real time, Real dt_lev)
{
    BL_PROFILE("REMORA::setup_step()");

    MultiFab& S_old = *cons_old[lev];
    MultiFab& S_new = *cons_new[lev];

    MultiFab& U_old = *xvel_old[lev];
    MultiFab& V_old = *yvel_old[lev];
    MultiFab& W_old = *zvel_old[lev];

    MultiFab& U_new = *xvel_new[lev];
    MultiFab& V_new = *yvel_new[lev];
    MultiFab& W_new = *zvel_new[lev];

    int nvars = S_old.nComp();

    // Fill ghost cells/faces at old time
    FillPatchNoBC(lev, time, *cons_old[lev], cons_old, BdyVars::t);
    FillPatchNoBC(lev, time, *xvel_old[lev], xvel_old, BdyVars::u);
    FillPatchNoBC(lev, time, *yvel_old[lev], yvel_old, BdyVars::v);
    FillPatch(lev, time, *zvel_old[lev], zvel_old, BCVars::zvel_bc, BdyVars::null);

    //////////    //pre_step3d corrections to boundaries

    const BoxArray&            ba = S_old.boxArray();
    const DistributionMapping& dm = S_old.DistributionMap();

    const int nrhs  = 0;
    int nstp  = 0;

    // Place-holder for source array -- for now just set to 0
    MultiFab source(ba,dm,nvars,1);
    source.setVal(0.0_rt);

    //-----------------------------------------------------------------------
    //  Time step momentum equation
    //-----------------------------------------------------------------------

    //Only used locally, probably should be rearranged into FArrayBox declaration
    MultiFab mf_AK(ba,dm,1,IntVect(NGROW,NGROW,0)); //2d missing j coordinate
    MultiFab mf_DC(ba,dm,1,IntVect(NGROW,NGROW,NGROW-1)); //2d missing j coordinate
    MultiFab mf_Hzk(ba,dm,1,IntVect(NGROW,NGROW,NGROW-1)); //2d missing j coordinate
    MultiFab mf_logdrg_tmp(ba,dm,1,IntVect(NGROW,NGROW,0));

    MultiFab* mf_z_r = vec_z_r[lev].get();
    MultiFab* mf_z_w = vec_z_w[lev].get();
    MultiFab* mf_h   = vec_hOfTheConfusingName[lev].get();
    MultiFab* mf_pm  = vec_pm[lev].get();
    MultiFab* mf_pn  =   vec_pn[lev].get();
    MultiFab* mf_fcor  = vec_fcor[lev].get();

    MultiFab* mf_gls = vec_gls[lev].get();
    MultiFab* mf_tke = vec_tke[lev].get();

    //Consider passing these into the advance function or renaming relevant things

    MultiFab mf_rho(ba,dm,1,IntVect(NGROW,NGROW,0));
    std::unique_ptr<MultiFab>& mf_rhoS = vec_rhoS[lev];
    std::unique_ptr<MultiFab>& mf_rhoA = vec_rhoA[lev];
    std::unique_ptr<MultiFab>& mf_bvf = vec_bvf[lev];
    std::unique_ptr<MultiFab>& mf_ru = vec_ru[lev];
    std::unique_ptr<MultiFab>& mf_rv = vec_rv[lev];
    std::unique_ptr<MultiFab>& mf_rufrc = vec_rufrc[lev];
    std::unique_ptr<MultiFab>& mf_rvfrc = vec_rvfrc[lev];
    std::unique_ptr<MultiFab>& mf_sustr = vec_sustr[lev];
    std::unique_ptr<MultiFab>& mf_svstr = vec_svstr[lev];
    std::unique_ptr<MultiFab>& mf_rdrag = vec_rdrag[lev];
    std::unique_ptr<MultiFab>& mf_rdrag2 = vec_rdrag2[lev];
    std::unique_ptr<MultiFab>& mf_ZoBot = vec_ZoBot[lev];
    std::unique_ptr<MultiFab>& mf_bustr = vec_bustr[lev];
    std::unique_ptr<MultiFab>& mf_bvstr = vec_bvstr[lev];

    std::unique_ptr<MultiFab>& mf_mskr = vec_mskr[lev];
    std::unique_ptr<MultiFab>& mf_msku = vec_msku[lev];
    std::unique_ptr<MultiFab>& mf_mskv = vec_mskv[lev];
    std::unique_ptr<MultiFab>& mf_mskp = vec_mskp[lev];

    MultiFab mf_rw(ba,dm,1,IntVect(NGROW,NGROW,0));

    std::unique_ptr<MultiFab>& mf_visc2_p = vec_visc2_p[lev];
    std::unique_ptr<MultiFab>& mf_visc2_r = vec_visc2_r[lev];

    // We need to set these because otherwise in the first call to remora_advance we may
    //    read uninitialized data on ghost values in setting the bc's on the velocities
    mf_rho.setVal(0.e34_rt,IntVect(AMREX_D_DECL(NGROW-1,NGROW-1,0)));
    mf_rhoS->setVal(0.e34_rt,IntVect(AMREX_D_DECL(NGROW-1,NGROW-1,0)));
    mf_rhoA->setVal(0.e34_rt,IntVect(AMREX_D_DECL(NGROW-1,NGROW-1,0)));

    mf_DC.setVal(0);

    FillPatchNoBC(lev, time, *cons_new[lev], cons_new, BdyVars::t);
    FillPatchNoBC(lev, time, *xvel_new[lev], xvel_new, BdyVars::u);
    FillPatchNoBC(lev, time, *yvel_new[lev], yvel_new, BdyVars::v);

    mf_rw.setVal(0.0_rt);
    mf_rufrc->setVal(0);
    mf_rvfrc->setVal(0);

    int iic = istep[lev];
    int ntfirst = 0;
    if(iic==ntfirst) {
        MultiFab::Copy(S_new,S_old,0,0,S_new.nComp(),S_new.nGrowVect());
        MultiFab::Copy(U_new,U_old,0,0,U_new.nComp(),U_new.nGrowVect());
        MultiFab::Copy(V_new,V_old,0,0,V_new.nComp(),V_new.nGrowVect());
        MultiFab::Copy(W_new,W_old,0,0,W_new.nComp(),W_new.nGrowVect());
    }

    // If we're not doing bulk fluxes, set surface momentum fluxes directly.
    // Otherwise, calculate them from winds, so those need to be set
    if (!solverChoice.bulk_fluxes) {
        set_smflux(lev);
    } else {
        set_wind(lev);
    }

    auto N = Geom(lev).Domain().size()[2]-1; // Number of vertical "levs" aka, NZ

    const auto dlo = amrex::lbound(Geom(lev).Domain());
    const auto dhi = amrex::ubound(Geom(lev).Domain());

    for ( MFIter mfi(S_new, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<Real const> const& h     = vec_hOfTheConfusingName[lev]->const_array(mfi);
        Array4<Real const> const& Hz    = vec_Hz[lev]->const_array(mfi);
        Array4<Real      > const& Huon  = vec_Huon[lev]->array(mfi);
        Array4<Real      > const& Hvom  = vec_Hvom[lev]->array(mfi);

        Array4<Real const> const& z_w   = mf_z_w->const_array(mfi);
        Array4<Real const> const& z_r   = mf_z_r->const_array(mfi);
        Array4<Real const> const& uold  = U_old.const_array(mfi);
        Array4<Real const> const& vold  = V_old.const_array(mfi);
        Array4<Real      > const& rho   = mf_rho.array(mfi);
        Array4<Real      > const& rhoA  = mf_rhoA->array(mfi);
        Array4<Real      > const& rhoS  = mf_rhoS->array(mfi);
        Array4<Real      > const& bvf   = mf_bvf->array(mfi);
        Array4<Real      > const& alpha = (solverChoice.bulk_fluxes) ? vec_alpha[lev]->array(mfi) : Array4<Real>();
        Array4<Real      > const& beta  = (solverChoice.bulk_fluxes) ? vec_beta[lev]->array(mfi)  : Array4<Real>();

        Array4<Real const> const& pm = mf_pm->const_array(mfi);
        Array4<Real const> const& pn = mf_pn->const_array(mfi);

        Array4<Real const> const& mskr = mf_mskr->const_array(mfi);

        Box  bx = mfi.tilebox();
        Box gbx1 = mfi.growntilebox(IntVect(NGROW-1,NGROW-1,0));
        Box gbx2 = mfi.growntilebox(IntVect(NGROW,NGROW,0));
        Box ugbx2 = mfi.grownnodaltilebox(0,IntVect(NGROW,NGROW,0));
        Box vgbx2 = mfi.grownnodaltilebox(1,IntVect(NGROW,NGROW,0));

        Box bxD = bx;
        bxD.makeSlab(2,0);
        Box gbx1D = gbx1;
        gbx1D.makeSlab(2,0);
        Box gbx2D = gbx2;
        gbx2D.makeSlab(2,0);

        FArrayBox fab_FC(gbx2,1,amrex::The_Async_Arena()); //3D
        FArrayBox fab_FX(gbx2,1,amrex::The_Async_Arena()); //3D
        FArrayBox fab_FE(gbx2,1,amrex::The_Async_Arena()); //3D
        FArrayBox fab_BC(gbx2,1,amrex::The_Async_Arena());
        FArrayBox fab_CF(gbx2,1,amrex::The_Async_Arena());

        //
        //-----------------------------------------------------------------------
        //  Compute horizontal mass fluxes, Hz*u/n and Hz*v/m (set_massflux_3d)
        //-----------------------------------------------------------------------
        //
        ParallelFor(ugbx2, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real on_u = 2.0_rt / (pn(i-1,j,0)+pn(i,j,0));
            Huon(i,j,k)=0.5_rt*(Hz(i,j,k)+Hz(i-1,j,k))*uold(i,j,k)* on_u;
        });

        ParallelFor(vgbx2, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real om_v= 2.0_rt / (pm(i,j-1,0)+pm(i,j,0));
            Hvom(i,j,k)=0.5_rt*(Hz(i,j,k)+Hz(i,j-1,k))*vold(i,j,k)* om_v;
        });

        Array4<Real const> const& state_old = S_old.const_array(mfi);
        rho_eos(gbx2,state_old,rho,rhoA,rhoS,bvf,alpha,beta,Hz,z_w,z_r,h,mskr,N);
    }

    const Real Cdb_min = solverChoice.Cdb_min;
    const Real Cdb_max = solverChoice.Cdb_max;

    if (solverChoice.bulk_fluxes) {
        bulk_fluxes(lev, cons_old[lev],vec_uwind[lev].get(),vec_vwind[lev].get(),
                    vec_evap[lev].get(),
                    vec_sustr[lev].get(),vec_svstr[lev].get(),vec_stflux[lev].get(),
                    vec_lrflx[lev].get(),vec_lhflx[lev].get(),vec_shflx[lev].get(),N);
        vec_evap[lev]->FillBoundary(geom[lev].periodicity());
    }

    if (solverChoice.do_temp_flux) {
        for ( MFIter mfi(S_new, TilingIfNotGPU()); mfi.isValid(); ++mfi )
        {
            Array4<Real      > const& stflx =  vec_stflx[lev]->array(mfi);
            Array4<Real      > const& btflx =  vec_btflx[lev]->array(mfi);
            Array4<Real const> const& stflux = vec_stflux[lev]->array(mfi);
            Array4<Real const> const& btflux = vec_btflux[lev]->array(mfi);
            Box gbx2 = mfi.growntilebox(IntVect(NGROW,NGROW,0));
            Box gbx2D = gbx2;
            gbx2D.makeSlab(2,0);
            ParallelFor(gbx2D, [=] AMREX_GPU_DEVICE (int i, int j, int ) {
                stflx(i,j,0,Temp_comp) = stflux(i,j,0,Temp_comp);
                btflx(i,j,0,Temp_comp) = btflux(i,j,0,Temp_comp);
            });
        }
    }
    if (solverChoice.do_salt_flux) {
        for ( MFIter mfi(S_new, TilingIfNotGPU()); mfi.isValid(); ++mfi )
        {
            Array4<Real      > const& stflx =  vec_stflx[lev]->array(mfi);
            Array4<Real      > const& btflx =  vec_btflx[lev]->array(mfi);
            Array4<Real const> const& stflux = vec_stflux[lev]->const_array(mfi);
            Array4<Real const> const& salt_old = S_old.const_array(mfi,Salt_comp);
            Box gbx2 = mfi.growntilebox(IntVect(NGROW,NGROW,0));
            Box gbx2D = gbx2;
            gbx2D.makeSlab(2,0);
            ParallelFor(gbx2D, [=] AMREX_GPU_DEVICE (int i, int j, int ) {
                stflx(i,j,0,Salt_comp) = stflux(i,j,0,Salt_comp) * salt_old(i,j,N);
                // The fact that this is btflx on the RHS matches what's in ROMS even though
                // it's weird -- if it's non-zero, does that mean that it will run away since
                // it's always getting multiplied by salt?
                btflx(i,j,0,Salt_comp) = btflx(i,j,0,Salt_comp) * salt_old(i,j,0);
            });
        }
    }

    for ( MFIter mfi(S_new, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<Real      > const& bustr = mf_bustr->array(mfi);
        Array4<Real      > const& bvstr = mf_bvstr->array(mfi);
        Array4<Real const> const& uold  = U_old.const_array(mfi);
        Array4<Real const> const& vold  = V_old.const_array(mfi);
        Array4<Real      > const& logdrg_tmp = mf_logdrg_tmp.array(mfi);
        Array4<Real const> const& z_r   = mf_z_r->const_array(mfi);
        Array4<Real const> const& z_w   = mf_z_w->const_array(mfi);

        Box gbx2 = mfi.growntilebox(IntVect(NGROW,NGROW,0));
        Box gbx2D = gbx2;
        gbx2D.makeSlab(2,0);
        Box ubx1 = mfi.grownnodaltilebox(0,IntVect(NGROW-1,NGROW-1,0));
        Box ubx1D = ubx1;
        ubx1D.makeSlab(2,0);
        Box vbx1 = mfi.grownnodaltilebox(1,IntVect(NGROW-1,NGROW-1,0));
        Box vbx1D = vbx1;
        vbx1D.makeSlab(2,0);
        // Set bottom stress as defined in set_vbx.F
        if (solverChoice.bottom_stress_type == BottomStressType::linear) {
            Array4<Real const> const& rdrag = mf_rdrag->const_array(mfi);
            ParallelFor(ubx1D, [=] AMREX_GPU_DEVICE (int i, int j, int )
            {
                bustr(i,j,0) = 0.5_rt * (rdrag(i-1,j,0)+rdrag(i,j,0))*(uold(i,j,0));
            });
            ParallelFor(vbx1D, [=] AMREX_GPU_DEVICE (int i, int j, int )
            {
                bvstr(i,j,0) = 0.5_rt * (rdrag(i,j-1,0)+rdrag(i,j,0))*(vold(i,j,0));
            });
        } else if (solverChoice.bottom_stress_type == BottomStressType::quadratic) {
            Array4<Real const> const& rdrag2 = mf_rdrag2->const_array(mfi);
            ParallelFor(ubx1D, [=] AMREX_GPU_DEVICE (int i, int j, int )
            {
                Real avg_v = 0.25_rt * (vold(i,j,0) + vold(i,j+1,0) + vold(i-1,j,0) + vold(i-1,j+1,0));
                Real vel_mag = std::sqrt(uold(i,j,0)*uold(i,j,0) + avg_v * avg_v);
                bustr(i,j,0) = 0.5_rt * (rdrag2(i-1,j,0) + rdrag2(i,j,0)) * uold(i,j,0) * vel_mag;
            });
            ParallelFor(vbx1D, [=] AMREX_GPU_DEVICE (int i, int j, int )
            {
                Real avg_u = 0.25_rt * (uold(i,j,0) + uold(i+1,j,0) + uold(i,j-1,0) + uold(i+1,j-1,0));
                Real vel_mag = std::sqrt(avg_u * avg_u + vold(i,j,0) * vold(i,j,0));
                bvstr(i,j,0) = 0.5_rt * (rdrag2(i,j-1,0) + rdrag2(i,j,0)) * vold(i,j,0) * vel_mag;
            });
        } else if (solverChoice.bottom_stress_type == BottomStressType::logarithmic) {
            Array4<Real const> const& ZoBot = mf_ZoBot->const_array(mfi);
            ParallelFor(gbx2D, [=] AMREX_GPU_DEVICE (int i, int j, int )
            {
                Real logz = 1.0_rt / (std::log((z_r(i,j,0) - z_w(i,j,0)) / ZoBot(i,j,0)));
                Real cff = vonKar * vonKar * logz * logz;
                logdrg_tmp(i,j,0) = std::min(Cdb_max,std::max(Cdb_min,cff));
            });
            ParallelFor(ubx1D, [=] AMREX_GPU_DEVICE (int i, int j, int )
            {
                Real avg_v = 0.25_rt * (vold(i,j,0) + vold(i,j+1,0) + vold(i-1,j,0) + vold(i-1,j+1,0));
                Real vel_mag = std::sqrt(uold(i,j,0)*uold(i,j,0) + avg_v * avg_v);
                bustr(i,j,0) = 0.5_rt * (logdrg_tmp(i-1,j,0)+logdrg_tmp(i,j,0)) * uold(i,j,0) * vel_mag;
            });
            ParallelFor(vbx1D, [=] AMREX_GPU_DEVICE (int i, int j, int )
            {
                Real avg_u = 0.25_rt * (uold(i,j,0) + uold(i+1,j,0) + uold(i,j-1,0) + uold(i+1,j-1,0));
                Real vel_mag = std::sqrt(avg_u * avg_u + vold(i,j,0) * vold(i,j,0));
                bvstr(i,j,0) = 0.5_rt * (logdrg_tmp(i,j-1,0) + logdrg_tmp(i,j,0)) * vold(i,j,0) * vel_mag;
            });
        }
    }
    FillPatch(lev, time, *vec_bustr[lev].get(), GetVecOfPtrs(vec_bustr), BCVars::u2d_simple_bc, BdyVars::null,0,true,false);
    FillPatch(lev, time, *vec_bvstr[lev].get(), GetVecOfPtrs(vec_bvstr), BCVars::v2d_simple_bc, BdyVars::null,0,true,false);

    if (solverChoice.vert_mixing_type == VertMixingType::analytic) {
        // Update Akv if using analytic mixing
        set_analytic_vmix(lev);
    }

    set_zeta_to_Ztavg(lev);

    MultiFab mf_W(convert(ba,IntVect(0,0,1)),dm,1,IntVect(NGROW+1,NGROW+1,0));
    mf_W.setVal(0.0_rt);


    if (solverChoice.use_prestep) {
        const int nnew  = 0;
        prestep(lev, U_old, V_old, U_new, V_new,
                mf_ru.get(), mf_rv.get(),
                S_old, S_new, mf_W,
                mf_DC, mf_z_r, mf_z_w, mf_h, mf_pm, mf_pn,
                mf_sustr.get(), mf_svstr.get(), mf_bustr.get(), mf_bvstr.get(),
                mf_msku.get(), mf_mskv.get(),
                iic, ntfirst, nnew, nstp, nrhs, N, dt_lev);
    }

    // We use FillBoundary not FillPatch here since mf_W is single-level scratch space
    mf_W.FillBoundary(geom[lev].periodicity());
    (*physbcs[lev])(mf_W,*mf_mskr.get(),0,1,mf_W.nGrowVect(),t_new[lev],BCVars::zvel_bc);

#ifdef REMORA_USE_NETCDF
    // Get u and v climatology if we're going to do nudging
    if (solverChoice.do_m3_clim_nudg) {
        u_clim_data_from_file->update_interpolated_to_time(t_new[lev]);
        v_clim_data_from_file->update_interpolated_to_time(t_new[lev]);
    }
#endif

    for ( MFIter mfi(S_old, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<Real const> const& Hz    = vec_Hz[lev]->const_array(mfi);
        Array4<Real> const& Huon  = vec_Huon[lev]->array(mfi);
        Array4<Real> const& Hvom  = vec_Hvom[lev]->array(mfi);
        Array4<Real> const& z_r   = (mf_z_r)->array(mfi);
        Array4<Real> const& z_w   = (mf_z_w)->array(mfi);
        Array4<Real const> const& uold  = U_old.const_array(mfi);
        Array4<Real const> const& vold  = V_old.const_array(mfi);
        Array4<Real> const& u    = U_new.array(mfi);
        Array4<Real> const& v    = V_new.array(mfi);
        Array4<Real> const& rho = (mf_rho).array(mfi);
        Array4<Real> const& ru = (mf_ru)->array(mfi);
        Array4<Real> const& rv = (mf_rv)->array(mfi);
        Array4<Real> const& rufrc = (mf_rufrc)->array(mfi);
        Array4<Real> const& rvfrc = (mf_rvfrc)->array(mfi);
        Array4<Real> const& W = (mf_W).array(mfi);
        Array4<Real> const& sustr = (mf_sustr)->array(mfi);
        Array4<Real> const& svstr = (mf_svstr)->array(mfi);
        Array4<Real> const& bustr = (mf_bustr)->array(mfi);
        Array4<Real> const& bvstr = (mf_bvstr)->array(mfi);
        Array4<Real> const& visc2_p = (mf_visc2_p)->array(mfi);
        Array4<Real> const& visc2_r = (mf_visc2_r)->array(mfi);

        Array4<Real const> const& pm = mf_pm->const_array(mfi);
        Array4<Real const> const& pn = mf_pn->const_array(mfi);
        Array4<Real const> const& fcor = mf_fcor->const_array(mfi);

        Array4<Real const> const& msku = mf_msku->const_array(mfi);
        Array4<Real const> const& mskv = mf_mskv->const_array(mfi);
        Array4<Real const> const& mskp = mf_mskp->const_array(mfi);

        Box bx = mfi.tilebox();

        Box tbxp1 = bx;
        Box tbxp2 = bx;
        Box xbx = mfi.nodaltilebox(0);
        Box ybx = mfi.nodaltilebox(1);
        Box xbx_adj = mfi.nodaltilebox(0);
        Box ybx_adj = mfi.nodaltilebox(1);

        auto xbx_lo = lbound(xbx_adj);
        auto xbx_hi = ubound(xbx_adj);
        auto ybx_lo = lbound(ybx_adj);
        auto ybx_hi = ubound(ybx_adj);

        if (xbx_lo.x == dlo.x) {
            xbx_adj.growLo(0,-1);
        } else if (xbx_hi.x == dhi.x) {
            xbx_adj.growHi(0,-1);
        }

        if (ybx_lo.y == dlo.y) {
            ybx_adj.growLo(1,-1);
        } else if (ybx_hi.y == dhi.y) {
            ybx_adj.growHi(1,-1);
        }

        Box gbx1 = mfi.growntilebox(IntVect(NGROW-1,NGROW-1,0));
        Box gbx2 = mfi.growntilebox(IntVect(NGROW,NGROW,0));

        Box utbx = mfi.nodaltilebox(0);
        Box vtbx = mfi.nodaltilebox(1);

        tbxp1.grow(IntVect(NGROW-1,NGROW-1,0));
        tbxp2.grow(IntVect(NGROW,NGROW,0));

        Box bxD = bx;
        bxD.makeSlab(2,0);
        Box gbx1D = gbx1;
        gbx1D.makeSlab(2,0);
        Box gbx2D = gbx2;
        gbx2D.makeSlab(2,0);

        Box tbxp1D = tbxp1;
        tbxp1D.makeSlab(2,0);
        Box tbxp2D = tbxp2;
        tbxp2D.makeSlab(2,0);

        FArrayBox fab_FC(surroundingNodes(tbxp2,2),1,amrex::The_Async_Arena()); //3D
        FArrayBox fab_FX(gbx2,1,amrex::The_Async_Arena()); //3D
        FArrayBox fab_FE(gbx2,1,amrex::The_Async_Arena()); //3D
        FArrayBox fab_BC(gbx2,1,amrex::The_Async_Arena());
        FArrayBox fab_CF(gbx2,1,amrex::The_Async_Arena());

        FArrayBox fab_fomn(tbxp2D,1,amrex::The_Async_Arena());

        auto FC=fab_FC.array();

        auto fomn=fab_fomn.array();

        if (solverChoice.use_coriolis) {
            ParallelFor(tbxp2D, [=] AMREX_GPU_DEVICE (int i, int j, int  )
            {
                fomn(i,j,0) = fcor(i,j,0)*(1.0_rt/(pm(i,j,0)*pn(i,j,0)));
            });
        }

        ParallelFor(gbx2, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            FC(i,j,k)=0.0_rt;
        });

        prsgrd(tbxp1,gbx1,utbx,vtbx,ru,rv,pn,pm,rho,FC,Hz,z_r,z_w,msku,mskv,nrhs,N);

        // Apply mixing to temperature and, if use_salt, salt
        int ncomp = solverChoice.use_salt ? 2 : 1;
        Array4<Real> const&     s_arr = S_new.array(mfi);
        Array4<Real> const& s_arr_rhs = S_old.array(mfi);
        Array4<Real> const& diff2_arr = vec_diff2[lev]->array(mfi);

        t3dmix(bx, s_arr, s_arr_rhs, diff2_arr, Hz, pm, pn, msku, mskv, dt_lev, ncomp);

        Array4<Real> const& diff2_arr_scalar = vec_diff2[lev]->array(mfi,Scalar_comp);
        t3dmix(bx, S_new.array(mfi,Scalar_comp), S_old.array(mfi,Scalar_comp), diff2_arr_scalar, Hz, pm, pn, msku, mskv, dt_lev, 1);

        if (solverChoice.use_coriolis) {
            //-----------------------------------------------------------------------
            // coriolis
            //-----------------------------------------------------------------------
            //
            // ru, rv updated
            // In ROMS, coriolis is the first (un-ifdefed) thing to happen in rhs3d_tile, which gets called after t3dmix
            coriolis(xbx, ybx, uold, vold, ru, rv, Hz, fomn, nrhs, nrhs);
        }

#ifdef REMORA_USE_NETCDF
        if (solverChoice.do_m3_clim_nudg) {
            Array4<const Real> const& uclim = u_clim_data_from_file->mf_interpolated->const_array(mfi);
            Array4<const Real> const& vclim = v_clim_data_from_file->mf_interpolated->const_array(mfi);
            Array4<const Real> const& u_nudg_coeff = vec_nudg_coeff[BdyVars::u][lev]->const_array(mfi);
            Array4<const Real> const& v_nudg_coeff = vec_nudg_coeff[BdyVars::v][lev]->const_array(mfi);
            // These boxes are set to match ROMS
            apply_clim_nudg(xbx_adj, 1, 0, ru, uold, uclim, u_nudg_coeff, Hz, pm, pn);
            apply_clim_nudg(ybx_adj, 0, 1, rv, vold, vclim, v_nudg_coeff, Hz, pm, pn);
        }
#endif

        ////rufrc from 3d is set to ru, then the wind stress (and bottom stress) is added, then the mixing is added
        //rufrc=ru+sustr*om_u*on_u

        rhs_uv_3d(lev, xbx, ybx, uold, vold, ru, rv, rufrc, rvfrc,
                  sustr, svstr, bustr, bvstr, Huon, Hvom,
                  pm, pn, W, FC, nrhs, N);

        if(solverChoice.use_uv3dmix) {
            const int nnew = 0;
            uv3dmix(xbx, ybx, u, v, uold, vold, rufrc, rvfrc, visc2_p, visc2_r, Hz, pm, pn, mskp, nrhs, nnew, dt_lev);
        }
    } // MFIter

    int nnew = (iic +1)% 2;
    nstp = iic % 2;
    if (solverChoice.vert_mixing_type == VertMixingType::GLS) {
        gls_prestep(lev, mf_gls, mf_tke, mf_W, mf_msku.get(), mf_mskv.get(),
                nstp, nnew, iic, ntfirst, N, dt_lev);
    }
    nstp = 0;

    // Commenting out for now, but not sure it's necessary
    //FillPatch(lev, time, *cons_old[lev], cons_old, BCVars::cons_bc, BdyVars::t);
    //FillPatch(lev, time, *cons_new[lev], cons_new, BCVars::cons_bc, BdyVars::t);
    FillPatch(lev, time, *vec_sstore[lev], GetVecOfPtrs(vec_sstore), BCVars::cons_bc, BdyVars::t,0,true,true,0,0,dt_lev,*cons_old[lev]);

    // Don't actually want to apply boundary conditions here
    vec_Huon[lev]->FillBoundary(geom[lev].periodicity());
    vec_Hvom[lev]->FillBoundary(geom[lev].periodicity());
}
