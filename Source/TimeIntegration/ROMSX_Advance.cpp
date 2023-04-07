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

    U_old.FillBoundary();
    V_old.FillBoundary();
    W_old.FillBoundary();
    MultiFab::Copy(S_new,S_old,0,0,S_new.nComp(),S_new.nGrowVect());
    MultiFab::Copy(U_new,U_old,0,0,U_new.nComp(),U_new.nGrowVect());
    MultiFab::Copy(V_new,V_old,0,0,V_new.nComp(),V_new.nGrowVect());
    MultiFab::Copy(W_new,W_old,0,0,W_new.nComp(),W_new.nGrowVect());
    //////////    //pre_step3d corrections to boundaries

    auto& lev_old = vars_old[lev];
    // Moving terrain
    FillPatch(lev, time, lev_old);

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
    MultiFab mf_AK(ba,dm,1,IntVect(2,2,0)); //2d missing j coordinate
    MultiFab mf_DC(ba,dm,1,IntVect(2,2,1)); //2d missing j coordinate
    MultiFab mf_Hzk(ba,dm,1,IntVect(2,2,1)); //2d missing j coordinate
    std::unique_ptr<MultiFab>& mf_z_r = z_r[lev];
    //Consider passing these into the advance function or renaming relevant things
    MultiFab mf_u(ba,dm,1,IntVect(2,2,0));
    MultiFab mf_v(ba,dm,1,IntVect(2,2,0));
    MultiFab mf_w(ba,dm,1,IntVect(2,2,0));
    std::unique_ptr<MultiFab>& mf_ru = ru[lev];
    std::unique_ptr<MultiFab>& mf_rv = rv[lev];
    std::unique_ptr<MultiFab>& mf_sustr = sustr[lev];
    std::unique_ptr<MultiFab>& mf_svstr = svstr[lev];
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
    MultiFab mf_rw(ba,dm,1,IntVect(2,2,0));
    MultiFab mf_W(ba,dm,1,IntVect(3,3,0));
    // We need to set these because otherwise in the first call to romsx_advance we may
    //    read uninitialized data on ghost values in setting the bc's on the velocities
    mf_u.setVal(0.e34,IntVect(AMREX_D_DECL(1,1,0)));
    mf_v.setVal(0.e34,IntVect(AMREX_D_DECL(1,1,0)));
    mf_w.setVal(0);
    mf_DC.setVal(0);
    mf_w.setVal(0.e34,IntVect(AMREX_D_DECL(1,1,0)));
    MultiFab::Copy(mf_u,U_new,0,0,U_new.nComp(),IntVect(AMREX_D_DECL(2,2,0)));
    MultiFab::Copy(mf_v,V_new,0,0,V_new.nComp(),IntVect(AMREX_D_DECL(2,2,0)));
    MultiFab::Copy(mf_w,W_new,0,0,W_new.nComp(),IntVect(AMREX_D_DECL(2,2,0)));
    MultiFab::Copy(mf_W,S_old,Omega_comp,0,mf_W.nComp(),IntVect(AMREX_D_DECL(2,2,0)));
    //    MultiFab::Copy(S_new,S_old,0,0,S_new.nComp(),IntVect(AMREX_D_DECL(2,2,2)));
    mf_u.FillBoundary();
    mf_v.FillBoundary();
    mf_w.FillBoundary();
    mf_W.FillBoundary();

    mf_rw.setVal(0.0);
    mf_W.setVal(0.0);
    U_old.FillBoundary();
    V_old.FillBoundary();

    int iic = istep[lev];
    int ntfirst = 0;
    set_smflux(lev,time);
    /*
!
!  Set linear bottom stress.
!
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          bustr(i,j)=0.5_r8*(rdrag(i-1,j)+rdrag(i,j))*                  &
     &               ubar(i,j,krhs)
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
          bvstr(i,j)=0.5_r8*(rdrag(i,j-1)+rdrag(i,j))*                  &
     &               vbar(i,j,krhs)
        END DO
      END DO*/

    auto N = Geom(lev).Domain().size()[2]-1; // Number of vertical "levs" aka, NZ

    const auto dxi              = Geom(lev).InvCellSizeArray();
    const int Mm = Geom(lev).Domain().size()[1];
    auto geomdata = Geom(lev).data();

    for ( MFIter mfi(mf_u, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<Real> const& DC = mf_DC.array(mfi);
        Array4<Real> const& Akv_arr = (Akv[lev])->array(mfi);
        Array4<Real> const& Hz_arr  = (Hz[lev])->array(mfi);
        Array4<Real> const& z_r_arr = (mf_z_r)->array(mfi);
        Array4<Real> const& uold = (U_old).array(mfi);
        Array4<Real> const& vold = (V_old).array(mfi);
        Array4<Real> const& u = (mf_u).array(mfi);
        Array4<Real> const& v = (mf_v).array(mfi);
        Array4<Real> const& tempold = (mf_tempold).array(mfi);
        Array4<Real> const& saltold = (mf_saltold).array(mfi);
        Array4<Real> const& temp = (mf_temp).array(mfi);
        Array4<Real> const& salt = (mf_salt).array(mfi);
        Array4<Real> const& ru_arr = (mf_ru)->array(mfi);
        Array4<Real> const& rv_arr = (mf_rv)->array(mfi);
        Array4<Real> const& W = (mf_W).array(mfi);
        Array4<Real> const& sustr_arr = (mf_sustr)->array(mfi);
        Array4<Real> const& svstr_arr = (mf_svstr)->array(mfi);

        Box bx = mfi.tilebox();
        //copy the tilebox
        Box gbx1 = bx;
        Box gbx11 = bx;
        Box gbx2 = bx;
        Box gbx3uneven(IntVect(AMREX_D_DECL(bx.smallEnd(0)-3,bx.smallEnd(1)-3,bx.smallEnd(2))),
                       IntVect(AMREX_D_DECL(bx.bigEnd(0)+2,bx.bigEnd(1)+2,bx.bigEnd(2))));
        Box gbx2uneven(IntVect(AMREX_D_DECL(bx.smallEnd(0)-2,bx.smallEnd(1)-2,bx.smallEnd(2))),
                       IntVect(AMREX_D_DECL(bx.bigEnd(0)+1,bx.bigEnd(1)+1,bx.bigEnd(2))));
        //make only gbx be grown to match multifabs
        gbx2.grow(IntVect(2,2,0));
        gbx1.grow(IntVect(1,1,0));
        gbx11.grow(IntVect(1,1,1));

        FArrayBox fab_FC(gbx2,1,amrex::The_Async_Arena());
        FArrayBox fab_FX(gbx2,1,amrex::The_Async_Arena());
        FArrayBox fab_FE(gbx2,1,amrex::The_Async_Arena());
        FArrayBox fab_tempstore(gbx2,1,amrex::The_Async_Arena());
	FArrayBox fab_BC(gbx2,1,amrex::The_Async_Arena());
        FArrayBox fab_CF(gbx2,1,amrex::The_Async_Arena());
        FArrayBox fab_pn(gbx2,1,amrex::The_Async_Arena());
        FArrayBox fab_pm(gbx2,1,amrex::The_Async_Arena());
        FArrayBox fab_on_u(gbx2,1,amrex::The_Async_Arena());
        FArrayBox fab_om_v(gbx2,1,amrex::The_Async_Arena());
        FArrayBox fab_fomn(gbx2,1,amrex::The_Async_Arena());
        FArrayBox fab_Huon(gbx2,1,amrex::The_Async_Arena()); fab_Huon.setVal(0.);
        FArrayBox fab_Hvom(gbx2,1,amrex::The_Async_Arena()); fab_Hvom.setVal(0.);
        FArrayBox fab_oHz(gbx11,1,amrex::The_Async_Arena());

        auto FC=fab_FC.array();
        auto FX=fab_FX.array();
        auto FE=fab_FE.array();
        auto tempstore=fab_tempstore.array();
        auto pn=fab_pn.array();
        auto pm=fab_pm.array();
        auto on_u=fab_on_u.array();
        auto om_v=fab_om_v.array();
        auto fomn=fab_fomn.array();
        auto Huon=fab_Huon.array();
        auto Hvom=fab_Hvom.array();

        //From ana_grid.h and metrics.F
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

        Real lambda = 1.0;
        //
        //-----------------------------------------------------------------------
        // prestep_uv_3d
        //-----------------------------------------------------------------------
        //
        prestep_t_3d(bx, uold, vold, u, v, tempold, saltold, temp, salt, ru_arr, rv_arr, Hz_arr, Akv_arr, on_u, om_v, Huon, Hvom,
                     pm, pn, W, DC, FC, tempstore, FX, FE, iic, ntfirst, nnew, nstp, nrhs, N,
                          lambda, dt_lev);

        //
        //-----------------------------------------------------------------------
        // prestep_uv_3d
        //-----------------------------------------------------------------------
        //
        prestep_uv_3d(bx, uold, vold, u, v, ru_arr, rv_arr, Hz_arr, Akv_arr, on_u, om_v, Huon, Hvom,
                          pm, pn, W, DC, FC, z_r_arr, sustr_arr, svstr_arr, iic, ntfirst, nnew, nstp, nrhs, N,
                          lambda, dt_lev);

#ifdef UV_COR
        //
        //-----------------------------------------------------------------------
        // coriolis
        //-----------------------------------------------------------------------
        //
        coriolis(bx, uold, vold, ru_arr, rv_arr, Hz_arr, fomn, nrhs);
#endif

        //
        //-----------------------------------------------------------------------
        // rhs_3d
        //-----------------------------------------------------------------------
        //
        rhs_3d(bx, uold, vold, ru_arr, rv_arr, Huon, Hvom, W, FC, nrhs, N);

    } // MFIter

    advance_2d(lev, mf_u, mf_v, ru[lev], rv[lev],
               Zt_avg1[lev],
               DU_avg1[lev], DU_avg2[lev],
               DV_avg1[lev], DV_avg2[lev],
               rubar[lev], rvbar[lev], rzeta[lev],
                ubar[lev],  vbar[lev],  zeta[lev],
               hOfTheConfusingName[lev], ncomp, dt_lev);

    advance_3d(lev, mf_u, mf_v, mf_temp, mf_salt, ru[lev], rv[lev],
               DU_avg1[lev], DU_avg2[lev],
               DV_avg1[lev], DV_avg2[lev],
               ubar[lev],  vbar[lev],
               mf_AK, mf_DC,
               mf_Hzk, Akv[lev], Hz[lev], ncomp, N, dt_lev);

    MultiFab::Copy(U_new,mf_u,0,0,U_new.nComp(),IntVect(AMREX_D_DECL(1,1,0)));
    U_new.FillBoundary();

    MultiFab::Copy(V_new,mf_v,0,0,V_new.nComp(),IntVect(AMREX_D_DECL(1,1,0)));
    V_new.FillBoundary();

    //    MultiFab::Copy(W_new,mf_w,0,0,W_new.nComp(),IntVect(AMREX_D_DECL(1,1,0)));
    //    W_new.FillBoundary();
    //    MultiFab::Copy(mf_W,S_old,Omega_comp,0,mf_W.nComp(),mf_w.nGrowVect());

}
