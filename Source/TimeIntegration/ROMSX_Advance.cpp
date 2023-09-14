#include <ROMSX.H>
#include <Utils.H>

using namespace amrex;

// advance a single level for a single time step
 void
ROMSX::Advance (int lev, Real time, Real dt_lev, int /*iteration*/, int /*ncycle*/)
{
    BL_PROFILE("ROMSX::Advance()");

    MultiFab& S_old = vars_old[lev][Vars::cons];
    MultiFab& S_new = vars_new[lev][Vars::cons];

    MultiFab& U_old = vars_old[lev][Vars::xvel];
    MultiFab& V_old = vars_old[lev][Vars::yvel];
    MultiFab& W_old = vars_old[lev][Vars::zvel];

    MultiFab& U_new = vars_new[lev][Vars::xvel];
    MultiFab& V_new = vars_new[lev][Vars::yvel];
    MultiFab& W_new = vars_new[lev][Vars::zvel];

    U_old.FillBoundary(geom[lev].periodicity());
    V_old.FillBoundary(geom[lev].periodicity());
    W_old.FillBoundary(geom[lev].periodicity());
    MultiFab::Copy(S_new,S_old,0,0,S_new.nComp(),S_new.nGrowVect());
    MultiFab::Copy(U_new,U_old,0,0,U_new.nComp(),U_new.nGrowVect());
    MultiFab::Copy(V_new,V_old,0,0,V_new.nComp(),V_new.nGrowVect());
    MultiFab::Copy(W_new,W_old,0,0,W_new.nComp(),W_new.nGrowVect());

    //////////    //pre_step3d corrections to boundaries

    const BoxArray&            ba = S_old.boxArray();
    const DistributionMapping& dm = S_old.DistributionMap();

    int nvars = S_old.nComp();

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
    std::unique_ptr<MultiFab>& mf_Hz = vec_Hz[lev];
    std::unique_ptr<MultiFab>& mf_z_r = vec_z_r[lev];
    std::unique_ptr<MultiFab>& mf_z_w = vec_z_w[lev];
    std::unique_ptr<MultiFab>& mf_h = vec_hOfTheConfusingName[lev];
    //Consider passing these into the advance function or renaming relevant things
    /*
    MultiFab mf_u(ba,dm,1,IntVect(NGROW,NGROW,0));
    MultiFab mf_v(ba,dm,1,IntVect(NGROW,NGROW,0));
    MultiFab mf_uold(ba,dm,1,IntVect(NGROW,NGROW,0));
    MultiFab mf_vold(ba,dm,1,IntVect(NGROW,NGROW,0));*/
    MultiFab mf_u(U_new, amrex::make_alias, 0, 1);
    MultiFab mf_v(V_new, amrex::make_alias, 0, 1);
    MultiFab mf_uold(U_old, amrex::make_alias, 0, 1);
    MultiFab mf_vold(V_old, amrex::make_alias, 0, 1);
    MultiFab mf_w(ba,dm,1,IntVect(NGROW,NGROW,0));
    MultiFab mf_pden(ba,dm,1,IntVect(NGROW,NGROW,0));
    MultiFab mf_rho(ba,dm,1,IntVect(NGROW,NGROW,0));
    MultiFab mf_rhoS(ba,dm,1,IntVect(NGROW,NGROW,0));
    MultiFab mf_rhoA(ba,dm,1,IntVect(NGROW,NGROW,0));
    std::unique_ptr<MultiFab>& mf_ru = vec_ru[lev];
    std::unique_ptr<MultiFab>& mf_rv = vec_rv[lev];
    std::unique_ptr<MultiFab>& mf_rufrc = vec_rufrc[lev];
    std::unique_ptr<MultiFab>& mf_rvfrc = vec_rvfrc[lev];
    std::unique_ptr<MultiFab>& mf_sustr = vec_sustr[lev];
    std::unique_ptr<MultiFab>& mf_svstr = vec_svstr[lev];
    std::unique_ptr<MultiFab>& mf_rdrag = vec_rdrag[lev];
    std::unique_ptr<MultiFab>& mf_bustr = vec_bustr[lev];
    std::unique_ptr<MultiFab>& mf_bvstr = vec_bvstr[lev];
    std::unique_ptr<MultiFab>& mf_ubar = vec_ubar[lev];
    std::unique_ptr<MultiFab>& mf_vbar = vec_vbar[lev];
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
    /*
    mf_u.setVal(0.e34,IntVect(AMREX_D_DECL(NGROW-1,NGROW-1,0)));
    mf_v.setVal(0.e34,IntVect(AMREX_D_DECL(NGROW-1,NGROW-1,0)));
    mf_uold.setVal(0.e34,IntVect(AMREX_D_DECL(NGROW-1,NGROW-1,0)));
    mf_vold.setVal(0.e34,IntVect(AMREX_D_DECL(NGROW-1,NGROW-1,0)));*/
    mf_pden.setVal(0.e34,IntVect(AMREX_D_DECL(NGROW-1,NGROW-1,0)));
    mf_rho.setVal(0.e34,IntVect(AMREX_D_DECL(NGROW-1,NGROW-1,0)));
    mf_rhoS.setVal(0.e34,IntVect(AMREX_D_DECL(NGROW-1,NGROW-1,0)));
    mf_rhoA.setVal(0.e34,IntVect(AMREX_D_DECL(NGROW-1,NGROW-1,0)));
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
    if(iic==ntfirst&&false)
        MultiFab::Copy(S_new,S_old,0,0,S_new.nComp(),IntVect(AMREX_D_DECL(NGROW,NGROW,NGROW)));
    set_smflux(lev,t_old[lev]);

    auto N = Geom(lev).Domain().size()[2]-1; // Number of vertical "levs" aka, NZ

    const auto prob_lo          = Geom(lev).ProbLoArray();
    const auto dxi              = Geom(lev).InvCellSizeArray();
    const auto dx               = Geom(lev).CellSizeArray();
    const int Mm = Geom(lev).Domain().size()[1];
    auto geomdata = Geom(lev).data();

    //MFIter::allowMultipleMFIters(true);
    for ( MFIter mfi(mf_temp, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<Real> const& DC = mf_DC.array(mfi);
        Array4<Real> const& Akv = (vec_Akv[lev])->array(mfi);
        Array4<Real> const& h = (vec_hOfTheConfusingName[lev])->array(mfi);
        Array4<Real> const& Hz  = (vec_Hz[lev])->array(mfi);
        Array4<Real> const& Huon  = (vec_Huon[lev])->array(mfi);
        Array4<Real> const& Hvom  = (vec_Hvom[lev])->array(mfi);
        Array4<Real> const& z_r = (mf_z_r)->array(mfi);
        Array4<Real> const& z_w = (mf_z_w)->array(mfi);
        Array4<Real> const& uold = (mf_uold).array(mfi);
        Array4<Real> const& vold = (mf_vold).array(mfi);
        Array4<Real> const& u = (mf_u).array(mfi);
        Array4<Real> const& v = (mf_v).array(mfi);
        Array4<Real> const& pden = (mf_pden).array(mfi);
        Array4<Real> const& rho = (mf_rho).array(mfi);
        Array4<Real> const& rhoA = (mf_rhoA).array(mfi);
        Array4<Real> const& rhoS = (mf_rhoS).array(mfi);
        Array4<Real> const& tempold = (mf_tempold).array(mfi);
        Array4<Real> const& saltold = (mf_saltold).array(mfi);
        Array4<Real> const& temp = (mf_temp).array(mfi);
        Array4<Real> const& salt = (mf_salt).array(mfi);
        Array4<Real> const& tempstore = (vec_t3[lev])->array(mfi);
        Array4<Real> const& saltstore = (vec_s3[lev])->array(mfi);
        Array4<Real> const& ru = (mf_ru)->array(mfi);
        Array4<Real> const& rv = (mf_rv)->array(mfi);
        Array4<Real> const& rufrc = (mf_rufrc)->array(mfi);
        Array4<Real> const& rvfrc = (mf_rvfrc)->array(mfi);
        Array4<Real> const& W = (mf_W).array(mfi);
        Array4<Real> const& sustr = (mf_sustr)->array(mfi);
        Array4<Real> const& svstr = (mf_svstr)->array(mfi);
        Array4<Real> const& rdrag = (mf_rdrag)->array(mfi);
        Array4<Real> const& bustr = (mf_bustr)->array(mfi);
        Array4<Real> const& bvstr = (mf_bvstr)->array(mfi);
        Array4<Real> const& ubar = (mf_ubar)->array(mfi);
        Array4<Real> const& vbar = (mf_vbar)->array(mfi);
        Array4<Real> const& visc2_p = (mf_visc2_p)->array(mfi);
        Array4<Real> const& visc2_r = (mf_visc2_r)->array(mfi);
        Array4<Real> const& diff2_salt = (mf_diff2_salt)->array(mfi);
        Array4<Real> const& diff2_temp = (mf_diff2_temp)->array(mfi);

        Box bx = mfi.tilebox();
        Box gbx = mfi.growntilebox();
        Box gbx1 = mfi.growntilebox(IntVect(NGROW-1,NGROW-1,0));
        Box gbx2 = mfi.growntilebox(IntVect(NGROW,NGROW,0));
        Box gbx11 = mfi.growntilebox(IntVect(NGROW-1,NGROW-1,NGROW-1));
        //copy the tilebox
        //Box gbx1 = bx;
        //Box gbx11 = bx;
        //Box gbx2 = bx;
        //TODO: adjust for tiling
        //Box gbx3uneven(IntVect(AMREX_D_DECL(bx.smallEnd(0)-3,bx.smallEnd(1)-3,bx.smallEnd(2))),
        //               IntVect(AMREX_D_DECL(bx.bigEnd(0)+2,bx.bigEnd(1)+2,bx.bigEnd(2))));
        //Box gbx2uneven(IntVect(AMREX_D_DECL(bx.smallEnd(0)-2,bx.smallEnd(1)-2,bx.smallEnd(2))),
        //               IntVect(AMREX_D_DECL(bx.bigEnd(0)+1,bx.bigEnd(1)+1,bx.bigEnd(2))));
        Box ubx = surroundingNodes(bx,0);
        Box vbx = surroundingNodes(bx,1);
        //make only gbx be grown to match multifabs
        //gbx2.grow(IntVect(NGROW,NGROW,0));
        //gbx1.grow(IntVect(NGROW-1,NGROW-1,0));
        //gbx11.grow(IntVect(NGROW-1,NGROW-1,NGROW-1));

        Box bxD = bx;
        bxD.makeSlab(2,0);
        Box gbx1D = gbx1;
        gbx1D.makeSlab(2,0);
        Box gbx2D = gbx2;
        gbx2D.makeSlab(2,0);
        //gbx1D.grow(IntVect(NGROW-1,NGROW-1,0));
        //gbx2D.grow(IntVect(NGROW,NGROW,0));

        FArrayBox fab_FC(gbx2,1,amrex::The_Async_Arena()); //3D
        FArrayBox fab_FX(gbx2,1,amrex::The_Async_Arena()); //3D
        FArrayBox fab_FE(gbx2,1,amrex::The_Async_Arena()); //3D
        FArrayBox fab_BC(gbx2,1,amrex::The_Async_Arena());
        FArrayBox fab_CF(gbx2,1,amrex::The_Async_Arena());
        FArrayBox fab_pn(gbx2D,1,amrex::The_Async_Arena());
        FArrayBox fab_pm(gbx2D,1,amrex::The_Async_Arena());
        FArrayBox fab_on_u(gbx2D,1,amrex::The_Async_Arena());
        FArrayBox fab_om_v(gbx2D,1,amrex::The_Async_Arena());
        FArrayBox fab_om_u(gbx2D,1,amrex::The_Async_Arena());
        FArrayBox fab_on_v(gbx2D,1,amrex::The_Async_Arena());
        FArrayBox fab_om_r(gbx2D,1,amrex::The_Async_Arena());
        FArrayBox fab_on_r(gbx2D,1,amrex::The_Async_Arena());
        FArrayBox fab_om_p(gbx2D,1,amrex::The_Async_Arena());
        FArrayBox fab_on_p(gbx2D,1,amrex::The_Async_Arena());
        FArrayBox fab_pmon_u(gbx2D,1,amrex::The_Async_Arena());
        FArrayBox fab_pnom_u(gbx2D,1,amrex::The_Async_Arena());
        FArrayBox fab_pmon_v(gbx2D,1,amrex::The_Async_Arena());
        FArrayBox fab_pnom_v(gbx2D,1,amrex::The_Async_Arena());
        FArrayBox fab_fomn(gbx2D,1,amrex::The_Async_Arena());
        //FArrayBox fab_oHz(gbx11,1,amrex::The_Async_Arena());

        auto FC=fab_FC.array();
        auto FX=fab_FX.array();
        auto FE=fab_FE.array();
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

        //From ana_grid.h and metrics.F
        amrex::ParallelFor(gbx2D,
        [=] AMREX_GPU_DEVICE (int i, int j, int  )
            {
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

        amrex::ParallelFor(gbx2D,
        [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
          //Note: are the comment definitons right? Don't seem to match metrics.f90
          om_v(i,j,0)=1.0/dxi[0]; // 2/(pm(i,j-1)+pm(i,j))
          on_u(i,j,0)=1.0/dxi[1]; // 2/(pm(i,j-1)+pm(i,j))
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
            //DO j=Jstr,Jend
            //  DO i=IstrU,Iend
            bustr(i,j,0) = 0.5 * (rdrag(i-1,j,0)+rdrag(i,j,0))*(uold(i,j,0,nrhs));
            //DO j=JstrV,Jend
            //  DO i=Istr,Iend
            bvstr(i,j,0) = 0.5 * (rdrag(i,j-1,0)+rdrag(i,j,0))*(vold(i,j,0,nrhs));
        });

        // updates Huon/Hvom
        set_massflux_3d(gbx2,1,0,uold,Huon,Hz,on_u,nnew);
        set_massflux_3d(gbx2,0,1,vold,Hvom,Hz,om_v,nnew);

        rho_eos(gbx2,temp,salt,rho,rhoA,rhoS,pden,Hz,z_w,h,nrhs,N);
    }

    if(solverChoice.use_prestep) {
        prestep(lev, mf_uold, mf_vold, mf_u, mf_v, mf_ru, mf_rv, mf_tempold, mf_saltold,
                mf_temp, mf_salt, mf_Hz, vec_Akv[lev], vec_Huon[lev], vec_Hvom[lev], mf_W, mf_DC, vec_t3[lev],
                vec_s3[lev], mf_z_r, mf_z_w, mf_h, mf_sustr, mf_svstr, mf_bustr, mf_bvstr, iic, ntfirst, nnew,
                nstp, nrhs, N, dt_lev);
    }


    mf_W.FillBoundary(geom[lev].periodicity());

    for ( MFIter mfi(mf_temp, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<Real> const& DC = mf_DC.array(mfi);
        Array4<Real> const& Akv = (vec_Akv[lev])->array(mfi);
        Array4<Real> const& Hz  = (vec_Hz[lev])->array(mfi);
        Array4<Real> const& Huon  = (vec_Huon[lev])->array(mfi);
        Array4<Real> const& Hvom  = (vec_Hvom[lev])->array(mfi);
        Array4<Real> const& z_r = (mf_z_r)->array(mfi);
        Array4<Real> const& z_w = (mf_z_w)->array(mfi);
        Array4<Real> const& uold = (mf_uold).array(mfi);
        Array4<Real> const& vold = (mf_vold).array(mfi);
        Array4<Real> const& u = (mf_u).array(mfi);
        Array4<Real> const& v = (mf_v).array(mfi);
        Array4<Real> const& pden = (mf_pden).array(mfi);
        Array4<Real> const& rho = (mf_rho).array(mfi);
        Array4<Real> const& rhoA = (mf_rhoA).array(mfi);
        Array4<Real> const& rhoS = (mf_rhoS).array(mfi);
        Array4<Real> const& tempold = (mf_tempold).array(mfi);
        Array4<Real> const& saltold = (mf_saltold).array(mfi);
        Array4<Real> const& temp = (mf_temp).array(mfi);
        Array4<Real> const& salt = (mf_salt).array(mfi);
        Array4<Real> const& tempstore = (vec_t3[lev])->array(mfi);
        Array4<Real> const& saltstore = (vec_s3[lev])->array(mfi);
        Array4<Real> const& ru = (mf_ru)->array(mfi);
        Array4<Real> const& rv = (mf_rv)->array(mfi);
        Array4<Real> const& rufrc = (mf_rufrc)->array(mfi);
        Array4<Real> const& rvfrc = (mf_rvfrc)->array(mfi);
        Array4<Real> const& W = (mf_W).array(mfi);
        Array4<Real> const& sustr = (mf_sustr)->array(mfi);
        Array4<Real> const& svstr = (mf_svstr)->array(mfi);
        Array4<Real> const& rdrag = (mf_rdrag)->array(mfi);
        Array4<Real> const& bustr = (mf_bustr)->array(mfi);
        Array4<Real> const& bvstr = (mf_bvstr)->array(mfi);
        Array4<Real> const& ubar = (mf_ubar)->array(mfi);
        Array4<Real> const& vbar = (mf_vbar)->array(mfi);
        Array4<Real> const& visc2_p = (mf_visc2_p)->array(mfi);
        Array4<Real> const& visc2_r = (mf_visc2_r)->array(mfi);
        Array4<Real> const& diff2_salt = (mf_diff2_salt)->array(mfi);
        Array4<Real> const& diff2_temp = (mf_diff2_temp)->array(mfi);

        Box bx = mfi.tilebox();
        Box tbxp1 = bx;
        Box tbxp2 = bx;
        Box gbx = mfi.growntilebox();
        Box gbx1 = mfi.growntilebox(IntVect(NGROW-1,NGROW-1,0));
        Box gbx2 = mfi.growntilebox(IntVect(NGROW,NGROW,0));
        Box gbx11 = mfi.growntilebox(IntVect(NGROW-1,NGROW-1,NGROW-1));

        tbxp1.grow(IntVect(NGROW-1,NGROW-1,0));
        tbxp2.grow(IntVect(NGROW,NGROW,0));
        //copy the tilebox
        //Box gbx1 = bx;
        //Box gbx11 = bx;
        //Box gbx2 = bx;
        //TODO: adjust for tiling
        //Box gbx3uneven(IntVect(AMREX_D_DECL(bx.smallEnd(0)-3,bx.smallEnd(1)-3,bx.smallEnd(2))),
        //               IntVect(AMREX_D_DECL(bx.bigEnd(0)+2,bx.bigEnd(1)+2,bx.bigEnd(2))));
        //Box gbx2uneven(IntVect(AMREX_D_DECL(bx.smallEnd(0)-2,bx.smallEnd(1)-2,bx.smallEnd(2))),
        //               IntVect(AMREX_D_DECL(bx.bigEnd(0)+1,bx.bigEnd(1)+1,bx.bigEnd(2))));
        Box ubx = surroundingNodes(bx,0);
        Box vbx = surroundingNodes(bx,1);
        //make only gbx be grown to match multifabs
        //gbx2.grow(IntVect(NGROW,NGROW,0));
        //gbx1.grow(IntVect(NGROW-1,NGROW-1,0));
        //gbx11.grow(IntVect(NGROW-1,NGROW-1,NGROW-1));

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

        FArrayBox fab_FC(gbx2,1,amrex::The_Async_Arena()); //3D
        FArrayBox fab_FX(gbx2,1,amrex::The_Async_Arena()); //3D
        FArrayBox fab_FE(gbx2,1,amrex::The_Async_Arena()); //3D
        FArrayBox fab_BC(gbx2,1,amrex::The_Async_Arena());
        FArrayBox fab_CF(gbx2,1,amrex::The_Async_Arena());
        FArrayBox fab_pn(tbxp2D,1,amrex::The_Async_Arena());
        FArrayBox fab_pm(tbxp2D,1,amrex::The_Async_Arena());
        FArrayBox fab_on_u(tbxp2D,1,amrex::The_Async_Arena());
        FArrayBox fab_om_v(tbxp2D,1,amrex::The_Async_Arena());
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
        //FArrayBox fab_oHz(gbx11,1,amrex::The_Async_Arena());

        auto FC=fab_FC.array();
        auto FX=fab_FX.array();
        auto FE=fab_FE.array();
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

        //From ana_grid.h and metrics.F
        amrex::ParallelFor(tbxp2D,
        [=] AMREX_GPU_DEVICE (int i, int j, int  )
            {
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

        amrex::ParallelFor(tbxp2D,
        [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
          //Note: are the comment definitons right? Don't seem to match metrics.f90
          om_v(i,j,0)=1.0/dxi[0]; // 2/(pm(i,j-1)+pm(i,j))
          on_u(i,j,0)=1.0/dxi[1]; // 2/(pm(i,j-1)+pm(i,j))
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

        prsgrd(gbx1,ru,rv,on_u,om_v,rho,FC,Hz,z_r,z_w,nrhs,N);
        if (verbose > 2) {
           amrex::PrintToFile("ru_afterprsgrd").SetPrecision(18)<<FArrayBox(ru)<<std::endl;
           amrex::PrintToFile("rv_afterprsgrd").SetPrecision(18)<<FArrayBox(rv)<<std::endl;
        }

        if (verbose > 2) {
            amrex::PrintToFile("temp_pret3dmix").SetPrecision(18)<<FArrayBox(temp) << std::endl;
            amrex::PrintToFile("salt_pret3dmix").SetPrecision(18)<<FArrayBox(salt) << std::endl;
        }
        t3dmix(bx, temp, diff2_temp, Hz, pm, pn, pmon_u, pnom_v, nrhs, nnew, dt_lev);
        t3dmix(bx, salt, diff2_salt, Hz, pm, pn, pmon_u, pnom_v, nrhs, nnew, dt_lev);
        if (verbose > 2) {
            amrex::PrintToFile("temp_postt3dmix").SetPrecision(18)<<FArrayBox(temp) << std::endl;
            amrex::PrintToFile("salt_postt3dmix").SetPrecision(18)<<FArrayBox(salt) << std::endl;
        }

        if (solverChoice.use_coriolis) {
            //-----------------------------------------------------------------------
            // coriolis
            //-----------------------------------------------------------------------
            //
            // ru, rv updated
            // In ROMS, coriolis is the first (un-ifdefed) thing to happen in rhs3d_tile, which gets called after t3dmix
            coriolis(bx, gbx, uold, vold, ru, rv, Hz, fomn, nrhs, nrhs);
        }
        if (verbose > 2) {
            amrex::PrintToFile("ru_aftercor").SetPrecision(18)<<FArrayBox(ru)<<std::endl;
            amrex::PrintToFile("rv_aftercor").SetPrecision(18)<<FArrayBox(rv)<<std::endl;
            amrex::PrintToFile("u_aftercor").SetPrecision(18)<<FArrayBox(u)<<std::endl;
            amrex::PrintToFile("v_aftercor").SetPrecision(18)<<FArrayBox(v)<<std::endl;
            amrex::PrintToFile("temp_aftercor").SetPrecision(18)<<FArrayBox(temp)<<std::endl;
            amrex::PrintToFile("tempstore_aftercor").SetPrecision(18)<<FArrayBox(tempstore)<<std::endl;
            amrex::PrintToFile("salt_aftercor").SetPrecision(18)<<FArrayBox(salt)<<std::endl;
            amrex::PrintToFile("saltstore_aftercor").SetPrecision(18)<<FArrayBox(saltstore)<<std::endl;
        }

        //
        //-----------------------------------------------------------------------
        // rhs_3d
        //-----------------------------------------------------------------------
        //

        ////rufrc from 3d is set to ru, then the wind stress (and bottom stress) is added, then the mixing is added
        //rufrc=ru+sustr*om_u*on_u
        rhs_3d(bx, gbx, uold, vold, ru, rv, rufrc, rvfrc, sustr, svstr, bustr, bvstr, Huon, Hvom, on_u, om_v, om_u, on_v, W, FC, nrhs, N);
        if (verbose > 2) {
           amrex::PrintToFile("ru_afterrhs").SetPrecision(18)<<FArrayBox(ru)<<std::endl;
           amrex::PrintToFile("rv_afterrhs").SetPrecision(18)<<FArrayBox(rv)<<std::endl;
           amrex::PrintToFile("u_afterrhs").SetPrecision(18)<<FArrayBox(u)<<std::endl;
           amrex::PrintToFile("v_afterrhs").SetPrecision(18)<<FArrayBox(v)<<std::endl;
        }

        if(solverChoice.use_uv3dmix) {
        //u=u+(contributions from S-surfaces viscosity not scaled by dt)*dt*dx*dy
        //rufrc=rufrc + (contributions from S-surfaces viscosity not scaled by dt*dx*dy)
        //For debugging version: uv3dmix(bx, u, v, rufrc, rvfrc, visc2_p, visc2_r, Hz, om_r, on_r, om_p, on_p, pm, pn, nrhs, nnew, dt_lev,false);
        uv3dmix(bx, u, v, uold, vold, rufrc, rvfrc, visc2_p, visc2_r, Hz, om_r, on_r, om_p, on_p, pm, pn, nrhs, nnew, dt_lev);
        if (verbose > 2) {
            amrex::PrintToFile("ru_afteruvmix").SetPrecision(18)<<FArrayBox(ru)<<std::endl;
            amrex::PrintToFile("rv_afteruvmix").SetPrecision(18)<<FArrayBox(rv)<<std::endl;
            amrex::PrintToFile("u_afteruvmix").SetPrecision(18)<<FArrayBox(u)<<std::endl;
            amrex::PrintToFile("v_afteruvmix").SetPrecision(18)<<FArrayBox(v)<<std::endl;
            amrex::PrintToFile("Huon").SetPrecision(18)<<FArrayBox(Huon)<<std::endl;
            amrex::PrintToFile("temp_uv3dmix").SetPrecision(18)<<FArrayBox(temp)<<std::endl;
            amrex::PrintToFile("tempstore_uv3dmix").SetPrecision(18)<<FArrayBox(tempstore)<<std::endl;
            amrex::PrintToFile("salt_uv3dmix").SetPrecision(18)<<FArrayBox(salt)<<std::endl;
            amrex::PrintToFile("saltstore_uv3dmix").SetPrecision(18)<<FArrayBox(saltstore)<<std::endl;
        }
        }
    } // MFIter

    mf_temp.FillBoundary(geom[lev].periodicity());
    mf_salt.FillBoundary(geom[lev].periodicity());
    mf_tempold.FillBoundary(geom[lev].periodicity());
    mf_saltold.FillBoundary(geom[lev].periodicity());
    vec_t3[lev]->FillBoundary(geom[lev].periodicity());
    vec_s3[lev]->FillBoundary(geom[lev].periodicity());
    vec_Huon[lev]->FillBoundary(geom[lev].periodicity());
    vec_Hvom[lev]->FillBoundary(geom[lev].periodicity());

    if(solverChoice.use_barotropic) {
    bool predictor_2d_step=true;
    bool first_2d_step=(iic==ntfirst);
    int nfast_counter=predictor_2d_step ? nfast : nfast-1;
    nfast_counter += 1;
    //Compute fast timestep from dt_lev and ratio
    Real dtfast_lev=dt_lev/Real(fixed_ndtfast_ratio);
    int next_indx1 = 0;
    for(int my_iif = 0; my_iif < nfast_counter; my_iif++) {
        first_2d_step=(my_iif==0);
        //Predictor
        predictor_2d_step=true;
        advance_2d(lev, mf_u, mf_v, mf_rhoS, mf_rhoA, vec_ru[lev], vec_rv[lev],
                   vec_rufrc[lev], vec_rvfrc[lev],
                   vec_Zt_avg1[lev],
                   vec_DU_avg1[lev], vec_DU_avg2[lev],
                   vec_DV_avg1[lev], vec_DV_avg2[lev],
                   vec_rubar[lev], vec_rvbar[lev], vec_rzeta[lev],
                    vec_ubar[lev],  vec_vbar[lev],  vec_zeta[lev],
                   vec_hOfTheConfusingName[lev], vec_visc2_p[lev], vec_visc2_r[lev],
                   ncomp, dt_lev, dtfast_lev, predictor_2d_step, first_2d_step, my_iif, nfast, next_indx1);

       //Corrector. Skip it on last fast step
        predictor_2d_step=false;
        if (my_iif < nfast_counter - 1) {
            advance_2d(lev, mf_u, mf_v, mf_rhoS, mf_rhoA, vec_ru[lev], vec_rv[lev],
                   vec_rufrc[lev], vec_rvfrc[lev],
                   vec_Zt_avg1[lev],
                   vec_DU_avg1[lev], vec_DU_avg2[lev],
                   vec_DV_avg1[lev], vec_DV_avg2[lev],
                   vec_rubar[lev], vec_rvbar[lev], vec_rzeta[lev],
                    vec_ubar[lev],  vec_vbar[lev],  vec_zeta[lev],
                   vec_hOfTheConfusingName[lev], vec_visc2_p[lev], vec_visc2_r[lev],
                   ncomp, dt_lev, dtfast_lev, predictor_2d_step, first_2d_step, my_iif, nfast, next_indx1);
        }
    }
    }

    // Currently, bathymetry doesn't change, so we do not need to re-set h here. Just re-do stretch, etc
    // because zeta may have changed
    stretch_transform(lev);

    if(solverChoice.use_baroclinic) {
    advance_3d(lev, mf_u, mf_v, mf_tempold, mf_saltold, mf_temp, mf_salt, vec_t3[lev], vec_s3[lev], vec_ru[lev], vec_rv[lev],
               vec_DU_avg1[lev], vec_DU_avg2[lev],
               vec_DV_avg1[lev], vec_DV_avg2[lev],
               vec_ubar[lev],  vec_vbar[lev],
               mf_AK, mf_DC,
               mf_Hzk, vec_Akv[lev], vec_Hz[lev], vec_Huon[lev], vec_Hvom[lev], vec_z_w[lev], vec_hOfTheConfusingName[lev], ncomp, N, dt_lev);
    }

    U_new.FillBoundary(geom[lev].periodicity());
    V_new.FillBoundary(geom[lev].periodicity());

    U_old.FillBoundary(geom[lev].periodicity());
    V_old.FillBoundary(geom[lev].periodicity());

    mf_temp.FillBoundary(geom[lev].periodicity());
    mf_salt.FillBoundary(geom[lev].periodicity());

    mf_tempold.FillBoundary(geom[lev].periodicity());
    mf_saltold.FillBoundary(geom[lev].periodicity());

    vec_t3[lev]->FillBoundary(geom[lev].periodicity());
    vec_s3[lev]->FillBoundary(geom[lev].periodicity());

    for (int lev = 0; lev <= finest_level; ++lev) {
        FillPatch(lev, t_new[lev], vars_new[lev]);
    }

#ifdef ROMSX_USE_PARTICLES
    // Update tracer particles on level 0
    if (lev == 0 && use_tracer_particles) {
        MultiFab* Umac = &vars_new[lev][Vars::xvel];
        tracer_particles->AdvectWithUmac(Umac, lev, dt_lev, *vec_z_phys_nd[0]);
    }
#endif


    //We are not storing computed W aka Omega
    //    MultiFab::Copy(W_new,mf_w,0,0,W_new.nComp(),IntVect(AMREX_D_DECL(NGROW-1,NGROW-1,0)));
    //    W_new.FillBoundary();
    //    MultiFab::Copy(mf_W,S_old,Omega_comp,0,mf_W.nComp(),mf_w.nGrowVect());

//    set_drag(lev);

}
