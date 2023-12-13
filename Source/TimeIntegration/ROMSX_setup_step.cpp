#include <ROMSX.H>

using namespace amrex;

// set up a time step for a single level
void
ROMSX::setup_step (int lev, Real time, Real dt_lev)
{
    BL_PROFILE("ROMSX::setup_step()");

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
    FillPatch(lev, time, cons_old[lev], cons_old);
    FillPatch(lev, time, xvel_old[lev], xvel_old);
    FillPatch(lev, time, yvel_old[lev], yvel_old);
    FillPatch(lev, time, zvel_old[lev], zvel_old);

    //////////    //pre_step3d corrections to boundaries

    const BoxArray&            ba = S_old.boxArray();
    const DistributionMapping& dm = S_old.DistributionMap();

    const int ncomp = 1;
    const int nrhs  = ncomp-1;
    const int nnew  = ncomp-1;
    const int nstp  = ncomp-1;

    // Place-holder for source array -- for now just set to 0
    MultiFab source(ba,dm,nvars,1);
    source.setVal(0.0);

    //-----------------------------------------------------------------------
    //  Time step momentum equation
    //-----------------------------------------------------------------------

    //Only used locally, probably should be rearranged into FArrayBox declaration
    MultiFab mf_AK(ba,dm,1,IntVect(NGROW,NGROW,0)); //2d missing j coordinate
    MultiFab mf_DC(ba,dm,1,IntVect(NGROW,NGROW,NGROW-1)); //2d missing j coordinate
    MultiFab mf_Hzk(ba,dm,1,IntVect(NGROW,NGROW,NGROW-1)); //2d missing j coordinate

    std::unique_ptr<MultiFab>& mf_z_r = vec_z_r[lev];
    std::unique_ptr<MultiFab>& mf_z_w = vec_z_w[lev];
    std::unique_ptr<MultiFab>& mf_h = vec_hOfTheConfusingName[lev];

    //Consider passing these into the advance function or renaming relevant things

    MultiFab mf_u(U_new, amrex::make_alias, 0, 1);
    MultiFab mf_v(V_new, amrex::make_alias, 0, 1);
    MultiFab mf_uold(U_old, amrex::make_alias, 0, 1);
    MultiFab mf_vold(V_old, amrex::make_alias, 0, 1);
    MultiFab mf_w(ba,dm,1,IntVect(NGROW,NGROW,0));
    MultiFab mf_rho(ba,dm,1,IntVect(NGROW,NGROW,0));
    std::unique_ptr<MultiFab>& mf_rhoS = vec_rhoS[lev];
    std::unique_ptr<MultiFab>& mf_rhoA = vec_rhoA[lev];
    std::unique_ptr<MultiFab>& mf_ru = vec_ru[lev];
    std::unique_ptr<MultiFab>& mf_rv = vec_rv[lev];
    std::unique_ptr<MultiFab>& mf_rufrc = vec_rufrc[lev];
    std::unique_ptr<MultiFab>& mf_rvfrc = vec_rvfrc[lev];
    std::unique_ptr<MultiFab>& mf_sustr = vec_sustr[lev];
    std::unique_ptr<MultiFab>& mf_svstr = vec_svstr[lev];
    std::unique_ptr<MultiFab>& mf_rdrag = vec_rdrag[lev];
    std::unique_ptr<MultiFab>& mf_bustr = vec_bustr[lev];
    std::unique_ptr<MultiFab>& mf_bvstr = vec_bvstr[lev];
    MultiFab mf_temp(S_new, amrex::make_alias, Temp_comp, 1);
#ifdef ROMSX_USE_SALINITY
    MultiFab mf_salt(S_new, amrex::make_alias, Salt_comp, 1);
#else
    MultiFab mf_salt(S_new, amrex::make_alias, Temp_comp, 1);
#endif
    MultiFab mf_tempold(S_old, amrex::make_alias, Temp_comp, 1);
#ifdef ROMSX_USE_SALINITY
    MultiFab mf_saltold(S_old, amrex::make_alias, Salt_comp, 1);
#else
    MultiFab mf_saltold(S_old, amrex::make_alias, Temp_comp, 1);
#endif
    MultiFab mf_rw(ba,dm,1,IntVect(NGROW,NGROW,0));
    MultiFab mf_W(ba,dm,1,IntVect(NGROW+1,NGROW+1,0));

    std::unique_ptr<MultiFab>& mf_visc2_p = vec_visc2_p[lev];
    std::unique_ptr<MultiFab>& mf_visc2_r = vec_visc2_r[lev];
    std::unique_ptr<MultiFab>& mf_diff2_temp = vec_diff2_temp[lev];
    std::unique_ptr<MultiFab>& mf_diff2_salt = vec_diff2_salt[lev];

    // We need to set these because otherwise in the first call to romsx_advance we may
    //    read uninitialized data on ghost values in setting the bc's on the velocities
    mf_rho.setVal(0.e34,IntVect(AMREX_D_DECL(NGROW-1,NGROW-1,0)));
    mf_rhoS->setVal(0.e34,IntVect(AMREX_D_DECL(NGROW-1,NGROW-1,0)));
    mf_rhoA->setVal(0.e34,IntVect(AMREX_D_DECL(NGROW-1,NGROW-1,0)));

    mf_w.setVal(0);
    mf_DC.setVal(0);
    mf_w.setVal(0.e34,IntVect(AMREX_D_DECL(NGROW-1,NGROW-1,0)));

    MultiFab::Copy(mf_u,U_new,0,0,U_new.nComp(),IntVect(AMREX_D_DECL(NGROW,NGROW,0)));
    MultiFab::Copy(mf_v,V_new,0,0,V_new.nComp(),IntVect(AMREX_D_DECL(NGROW,NGROW,0)));
    MultiFab::Copy(mf_uold,U_old,0,0,U_old.nComp(),IntVect(AMREX_D_DECL(NGROW,NGROW,0)));
    MultiFab::Copy(mf_vold,V_old,0,0,V_old.nComp(),IntVect(AMREX_D_DECL(NGROW,NGROW,0)));
    MultiFab::Copy(mf_w,W_new,0,0,W_new.nComp(),IntVect(AMREX_D_DECL(NGROW,NGROW,0)));
    MultiFab::Copy(mf_W,S_old,Omega_comp,0,mf_W.nComp(),IntVect(AMREX_D_DECL(NGROW,NGROW,0)));

    mf_u.FillBoundary(geom[lev].periodicity());
    mf_v.FillBoundary(geom[lev].periodicity());
    mf_uold.FillBoundary(geom[lev].periodicity());
    mf_vold.FillBoundary(geom[lev].periodicity());
    mf_w.FillBoundary(geom[lev].periodicity());
    mf_W.FillBoundary(geom[lev].periodicity());
    mf_tempold.FillBoundary(geom[lev].periodicity());
    mf_temp.FillBoundary(geom[lev].periodicity());
    mf_saltold.FillBoundary(geom[lev].periodicity());
    mf_salt.FillBoundary(geom[lev].periodicity());

    mf_rw.setVal(0.0);
    mf_W.setVal(0.0);
    U_old.FillBoundary(geom[lev].periodicity());
    V_old.FillBoundary(geom[lev].periodicity());
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
    set_smflux(lev,t_old[lev]);

    auto N = Geom(lev).Domain().size()[2]-1; // Number of vertical "levs" aka, NZ

    const auto prob_lo          = Geom(lev).ProbLoArray();
    const auto dxi              = Geom(lev).InvCellSizeArray();
    const auto dx               = Geom(lev).CellSizeArray();
    const int Mm = Geom(lev).Domain().size()[1];

    //MFIter::allowMultipleMFIters(true);
    for ( MFIter mfi(mf_temp, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<Real const> const& h = (vec_hOfTheConfusingName[lev])->const_array(mfi);
        Array4<Real const> const& Hz  = (vec_Hz[lev])->const_array(mfi);
        Array4<Real      > const& Huon  = (vec_Huon[lev])->array(mfi);
        Array4<Real      > const& Hvom  = (vec_Hvom[lev])->array(mfi);
        Array4<Real const> const& z_w = (mf_z_w)->const_array(mfi);
        Array4<Real      > const& uold = (mf_uold).array(mfi);
        Array4<Real      > const& vold = (mf_vold).array(mfi);
        Array4<Real      > const& rho = (mf_rho).array(mfi);
        Array4<Real      > const& rhoA = (mf_rhoA)->array(mfi);
        Array4<Real      > const& rhoS = (mf_rhoS)->array(mfi);
        Array4<Real const> const& tempold = (mf_tempold).const_array(mfi);
        Array4<Real const> const& saltold = (mf_saltold).const_array(mfi);
        Array4<Real      > const& rdrag = (mf_rdrag)->array(mfi);
        Array4<Real      > const& bustr = (mf_bustr)->array(mfi);
        Array4<Real      > const& bvstr = (mf_bvstr)->array(mfi);

        Box  bx = mfi.tilebox();
        Box ubx = Box(uold);
        Box vbx = Box(vold);
        Box gbx1 = mfi.growntilebox(IntVect(NGROW-1,NGROW-1,0));
        Box gbx2 = mfi.growntilebox(IntVect(NGROW,NGROW,0));

        //TODO: adjust for tiling
        //Box gbx3uneven(IntVect(AMREX_D_DECL(bx.smallEnd(0)-3,bx.smallEnd(1)-3,bx.smallEnd(2))),
        //               IntVect(AMREX_D_DECL(bx.bigEnd(0)+2,bx.bigEnd(1)+2,bx.bigEnd(2))));
        //Box gbx2uneven(IntVect(AMREX_D_DECL(bx.smallEnd(0)-2,bx.smallEnd(1)-2,bx.smallEnd(2))),
        //               IntVect(AMREX_D_DECL(bx.bigEnd(0)+1,bx.bigEnd(1)+1,bx.bigEnd(2))));

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

        FArrayBox fab_pn(gbx2D,1,amrex::The_Async_Arena());
        FArrayBox fab_pm(gbx2D,1,amrex::The_Async_Arena());

        FArrayBox fab_om_r(gbx2D,1,amrex::The_Async_Arena());
        FArrayBox fab_on_r(gbx2D,1,amrex::The_Async_Arena());
        FArrayBox fab_om_p(gbx2D,1,amrex::The_Async_Arena());
        FArrayBox fab_on_p(gbx2D,1,amrex::The_Async_Arena());

        FArrayBox fab_pmon_u(gbx2D,1,amrex::The_Async_Arena());
        FArrayBox fab_pnom_u(gbx2D,1,amrex::The_Async_Arena());
        FArrayBox fab_pmon_v(gbx2D,1,amrex::The_Async_Arena());
        FArrayBox fab_pnom_v(gbx2D,1,amrex::The_Async_Arena());

        FArrayBox fab_on_u(makeSlab(ubx,2,0),1,amrex::The_Async_Arena());
        FArrayBox fab_om_v(makeSlab(vbx,2,0),1,amrex::The_Async_Arena());
        FArrayBox fab_om_u(gbx2D,1,amrex::The_Async_Arena());
        FArrayBox fab_on_v(gbx2D,1,amrex::The_Async_Arena());

        auto pn=fab_pn.array();
        auto pm=fab_pm.array();
        auto on_u=fab_on_u.array();
        auto om_v=fab_om_v.array();
        auto om_u=fab_om_u.array();
        auto on_v=fab_on_v.array();
        auto om_r=fab_om_r.array();
        auto on_r=fab_on_r.array();
        auto om_p=fab_om_p.array();
        auto on_p=fab_on_p.array();
        auto pmon_u=fab_pmon_u.array();
        auto pnom_u=fab_pnom_u.array();
        auto pmon_v=fab_pmon_v.array();
        auto pnom_v=fab_pnom_v.array();

        amrex::ParallelFor(ubx, [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
          on_u(i,j,0)=1.0/dxi[1]; // 2/(pm(i,j-1)+pm(i,j))
        });
        amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
          om_v(i,j,0)=1.0/dxi[0]; // 2/(pm(i,j-1)+pm(i,j))
        });

        amrex::ParallelFor(gbx2D, [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
          //Note: are the comment definitons right? Don't seem to match metrics.f90
          om_r(i,j,0)=1.0/dxi[0]; // 1/pm(i,j)
          on_r(i,j,0)=1.0/dxi[1]; // 1/pn(i,j)
          //todo: om_p on_p
          om_p(i,j,0)=1.0/dxi[0]; // 4/(pm(i-1,j-1)+pm(i-1,j)+pm(i,j-1)+pm(i,j))
          on_p(i,j,0)=1.0/dxi[1]; // 4/(pn(i-1,j-1)+pn(i-1,j)+pn(i,j-1)+pn(i,j))
          on_v(i,j,0)=1.0/dxi[1]; // 2/(pn(i-1,j)+pn(i,j))
          om_u(i,j,0)=1.0/dxi[0]; // 2/(pm(i-1,j)+pm(i,j))
          pmon_u(i,j,0)=1.0;        // (pm(i-1,j)+pm(i,j))/(pn(i-1,j)+pn(i,j))
          pnom_u(i,j,0)=1.0;        // (pn(i-1,j)+pn(i,j))/(pm(i-1,j)+pm(i,j))
          pmon_v(i,j,0)=1.0;        // (pm(i,j-1)+pm(i,j))/(pn(i,j-1)+pn(i,j))
          pnom_v(i,j,0)=1.0;        // (pn(i,j-1)+pn(i,j))/(pm(i,j-1)+pm(i,j))
        });

        amrex::ParallelFor(gbx2,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
          Huon(i,j,k,0)=0.0;
          Hvom(i,j,k,0)=0.0;
        });

        // Set bottom stress as defined in set_vbx.F
        amrex::ParallelFor(gbx1D,
        [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
            bustr(i,j,0) = 0.5 * (rdrag(i-1,j,0)+rdrag(i,j,0))*(uold(i,j,0,nrhs));
            bvstr(i,j,0) = 0.5 * (rdrag(i,j-1,0)+rdrag(i,j,0))*(vold(i,j,0,nrhs));
        });

        // Updates Huon and Hvom

        set_massflux_3d(uold,Huon,on_u,vold,Hvom,om_v,Hz,nnew);

        rho_eos(gbx2,tempold,saltold,rho,rhoA,rhoS,Hz,z_w,h,nrhs,N);
    }

    if(solverChoice.use_prestep) {
        prestep(lev, mf_uold, mf_vold,
                mf_u, mf_v,
                mf_ru, mf_rv, mf_tempold, mf_saltold,
                mf_temp, mf_salt, mf_W,
                mf_DC, mf_z_r, mf_z_w, mf_h, mf_sustr, mf_svstr, mf_bustr,
                mf_bvstr, iic, ntfirst, nnew, nstp, nrhs, N, dt_lev);
    }


    mf_W.FillBoundary(geom[lev].periodicity());

    for ( MFIter mfi(mf_temp, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<Real> const& Hz  = (vec_Hz[lev])->array(mfi);
        Array4<Real> const& Huon  = (vec_Huon[lev])->array(mfi);
        Array4<Real> const& Hvom  = (vec_Hvom[lev])->array(mfi);
        Array4<Real> const& z_r = (mf_z_r)->array(mfi);
        Array4<Real> const& z_w = (mf_z_w)->array(mfi);
        Array4<Real> const& uold = (mf_uold).array(mfi);
        Array4<Real> const& vold = (mf_vold).array(mfi);
        Array4<Real> const& u = (mf_u).array(mfi);
        Array4<Real> const& v = (mf_v).array(mfi);
        Array4<Real> const& rho = (mf_rho).array(mfi);
        Array4<Real> const& temp = (mf_temp).array(mfi);
        Array4<Real> const& salt = (mf_salt).array(mfi);
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
        Array4<Real> const& diff2_salt = (mf_diff2_salt)->array(mfi);
        Array4<Real> const& diff2_temp = (mf_diff2_temp)->array(mfi);

        Array4<Real> const& zeta = (vec_zeta[lev])->array(mfi);
        Array4<Real> const& Zt_avg1 = (vec_Zt_avg1[lev])->array(mfi);

        Box bx = mfi.tilebox();

        Box tbxp1 = bx;
        Box tbxp2 = bx;
        Box xbx = mfi.nodaltilebox(0);
        Box ybx = mfi.nodaltilebox(1);
        Box gbx1 = mfi.growntilebox(IntVect(NGROW-1,NGROW-1,0));
        Box gbx2 = mfi.growntilebox(IntVect(NGROW,NGROW,0));

        Box ubx = Box(uold);
        Box vbx = Box(vold);

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
        //gbx1D.grow(IntVect(NGROW-1,NGROW-1,0));
        //gbx2D.grow(IntVect(NGROW,NGROW,0));

        FArrayBox fab_FC(tbxp2,1,amrex::The_Async_Arena()); //3D
        FArrayBox fab_FX(gbx2,1,amrex::The_Async_Arena()); //3D
        FArrayBox fab_FE(gbx2,1,amrex::The_Async_Arena()); //3D
        FArrayBox fab_BC(gbx2,1,amrex::The_Async_Arena());
        FArrayBox fab_CF(gbx2,1,amrex::The_Async_Arena());
        FArrayBox fab_pn(tbxp2D,1,amrex::The_Async_Arena());
        FArrayBox fab_pm(tbxp2D,1,amrex::The_Async_Arena());
        FArrayBox fab_om_u(tbxp2D,1,amrex::The_Async_Arena());
        FArrayBox fab_on_v(tbxp2D,1,amrex::The_Async_Arena());
        FArrayBox fab_om_r(tbxp2D,1,amrex::The_Async_Arena());
        FArrayBox fab_on_r(tbxp2D,1,amrex::The_Async_Arena());
        FArrayBox fab_om_p(tbxp2D,1,amrex::The_Async_Arena());
        FArrayBox fab_on_p(tbxp2D,1,amrex::The_Async_Arena());
        FArrayBox fab_pmon_u(tbxp2D,1,amrex::The_Async_Arena());
        FArrayBox fab_pnom_u(tbxp2D,1,amrex::The_Async_Arena());
        FArrayBox fab_pmon_v(tbxp2D,1,amrex::The_Async_Arena());
        FArrayBox fab_pnom_v(tbxp2D,1,amrex::The_Async_Arena());
        FArrayBox fab_fomn(tbxp2D,1,amrex::The_Async_Arena());

        FArrayBox fab_on_u(makeSlab(ubx,2,0),1,amrex::The_Async_Arena());
        FArrayBox fab_om_v(makeSlab(vbx,2,0),1,amrex::The_Async_Arena());

        auto FC=fab_FC.array();
        auto pn=fab_pn.array();
        auto pm=fab_pm.array();
        auto on_u=fab_on_u.array();
        auto om_v=fab_om_v.array();
        auto om_u=fab_om_u.array();
        auto on_v=fab_on_v.array();
        auto om_r=fab_om_r.array();
        auto on_r=fab_on_r.array();
        auto om_p=fab_om_p.array();
        auto on_p=fab_on_p.array();
        auto pmon_u=fab_pmon_u.array();
        auto pnom_u=fab_pnom_u.array();
        auto pmon_v=fab_pmon_v.array();
        auto pnom_v=fab_pnom_v.array();
        auto fomn=fab_fomn.array();

        Real coriolis_f0 = solverChoice.coriolis_f0;
        Real coriolis_beta = solverChoice.coriolis_beta;

        ParallelFor(tbxp2D,
        [=] AMREX_GPU_DEVICE (int i, int j, int  )
        {
            const auto dx              = geom[lev].CellSize();

            pm(i,j,0) = dxi[0];
            pn(i,j,0) = dxi[1];
            Real Esize=geom[lev].ProbHi()[1] - geom[lev].ProbLo()[1];
            Real y = prob_lo[1] + (j + 0.5) * dx[1];
            Real f=coriolis_f0 + coriolis_beta*(y-.5*Esize);
            fomn(i,j,0)=f*(1.0/(pm(i,j,0)*pn(i,j,0)));
        });
        amrex::ParallelFor(ubx, [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
            on_u(i,j,0)=1.0/dxi[1]; // 2/(pm(i,j-1)+pm(i,j))
        });
        amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
            om_v(i,j,0)=1.0/dxi[0]; // 2/(pm(i,j-1)+pm(i,j))
        });

        amrex::ParallelFor(tbxp2D,
        [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
          //Note: are the comment definitons right? Don't seem to match metrics.f90
          om_r(i,j,0)=1.0/dxi[0]; // 1/pm(i,j)
          on_r(i,j,0)=1.0/dxi[1]; // 1/pn(i,j)
          //todo: om_p on_p
          om_p(i,j,0)=1.0/dxi[0]; // 4/(pm(i-1,j-1)+pm(i-1,j)+pm(i,j-1)+pm(i,j))
          on_p(i,j,0)=1.0/dxi[1]; // 4/(pn(i-1,j-1)+pn(i-1,j)+pn(i,j-1)+pn(i,j))
          on_v(i,j,0)=1.0/dxi[1]; // 2/(pn(i-1,j)+pn(i,j))
          om_u(i,j,0)=1.0/dxi[0]; // 2/(pm(i-1,j)+pm(i,j))
          pmon_u(i,j,0)=1.0;        // (pm(i-1,j)+pm(i,j))/(pn(i-1,j)+pn(i,j))
          pnom_u(i,j,0)=1.0;        // (pn(i-1,j)+pn(i,j))/(pm(i-1,j)+pm(i,j))
          pmon_v(i,j,0)=1.0;        // (pm(i,j-1)+pm(i,j))/(pn(i,j-1)+pn(i,j))
          pnom_v(i,j,0)=1.0;        // (pn(i,j-1)+pn(i,j))/(pm(i,j-1)+pm(i,j))
        });

        amrex::ParallelFor(gbx2,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            FC(i,j,k)=0.0;
        });

        prsgrd(tbxp1,gbx1,utbx,vtbx,ru,rv,on_u,om_v,rho,FC,Hz,z_r,z_w,nrhs,N);

        t3dmix(bx, temp, diff2_temp, Hz, pm, pn, pmon_u, pnom_v, nrhs, nnew, dt_lev);
        t3dmix(bx, salt, diff2_salt, Hz, pm, pn, pmon_u, pnom_v, nrhs, nnew, dt_lev);

        if (solverChoice.use_coriolis) {
            //-----------------------------------------------------------------------
            // coriolis
            //-----------------------------------------------------------------------
            //
            // ru, rv updated
            // In ROMS, coriolis is the first (un-ifdefed) thing to happen in rhs3d_tile, which gets called after t3dmix
            coriolis(xbx, ybx, uold, vold, ru, rv, Hz, fomn, nrhs, nrhs);
        }

        //
        //-----------------------------------------------------------------------
        //

        ////rufrc from 3d is set to ru, then the wind stress (and bottom stress) is added, then the mixing is added
        //rufrc=ru+sustr*om_u*on_u

        rhs_uv_3d(xbx, ybx, uold, vold, ru, rv, rufrc, rvfrc,
                  sustr, svstr, bustr, bvstr, Huon, Hvom,
                  on_u, om_v, om_u, on_v, W, FC, nrhs, N);

        if(solverChoice.use_uv3dmix) {
            uv3dmix(xbx, ybx, u, v, uold, vold, rufrc, rvfrc, visc2_p, visc2_r, Hz, om_r, on_r, om_p, on_p, pm, pn, nrhs, nnew, dt_lev);
        }

        // Set first two components of zeta to time-averaged values before barotropic update
        amrex::ParallelFor(gbx2D,
        [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            zeta(i,j,0,0) = Zt_avg1(i,j,0);
            zeta(i,j,0,1) = Zt_avg1(i,j,0);
        });
    } // MFIter

    // Update Akv with new depth. NOTE: this happens before set_zeta in ROMS
    set_vmix(lev);

    mf_temp.FillBoundary(geom[lev].periodicity());
    mf_salt.FillBoundary(geom[lev].periodicity());

    mf_tempold.FillBoundary(geom[lev].periodicity());
    mf_saltold.FillBoundary(geom[lev].periodicity());

    vec_t3[lev]->FillBoundary(geom[lev].periodicity());
    vec_s3[lev]->FillBoundary(geom[lev].periodicity());

    vec_Huon[lev]->FillBoundary(geom[lev].periodicity());
    vec_Hvom[lev]->FillBoundary(geom[lev].periodicity());
}
