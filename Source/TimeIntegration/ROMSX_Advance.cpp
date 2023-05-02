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
    //    MultiFab::Copy(S_new,S_old,0,0,S_new.nComp(),S_new.nGrowVect());
    //    MultiFab::Copy(U_new,U_old,0,0,U_new.nComp(),U_new.nGrowVect());
    //    MultiFab::Copy(V_new,V_old,0,0,V_new.nComp(),V_new.nGrowVect());
    //    MultiFab::Copy(W_new,W_old,0,0,W_new.nComp(),W_new.nGrowVect());
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
    std::unique_ptr<MultiFab>& mf_z_r = vec_z_r[lev];
    //Consider passing these into the advance function or renaming relevant things
    MultiFab mf_u(ba,dm,1,IntVect(2,2,0));
    MultiFab mf_v(ba,dm,1,IntVect(2,2,0));
    MultiFab mf_w(ba,dm,1,IntVect(2,2,0));
    std::unique_ptr<MultiFab>& mf_ru = vec_ru[lev];
    std::unique_ptr<MultiFab>& mf_rv = vec_rv[lev];
    std::unique_ptr<MultiFab>& mf_sustr = vec_sustr[lev];
    std::unique_ptr<MultiFab>& mf_svstr = vec_svstr[lev];
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
    mf_u.FillBoundary();
    mf_v.FillBoundary();
    mf_w.FillBoundary();
    mf_W.FillBoundary();
    mf_tempold.FillBoundary();
    mf_temp.FillBoundary();
    mf_saltold.FillBoundary();
    mf_salt.FillBoundary();

    mf_rw.setVal(0.0);
    mf_W.setVal(0.0);
    U_old.FillBoundary();
    V_old.FillBoundary();

    int iic = istep[lev];
    int ntfirst = 0;
    if(iic==ntfirst&&false)
        MultiFab::Copy(S_new,S_old,0,0,S_new.nComp(),IntVect(AMREX_D_DECL(2,2,2)));
    set_smflux(lev,t_old[lev]);
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

    const auto prob_lo          = Geom(lev).ProbLoArray();
    const auto dxi              = Geom(lev).InvCellSizeArray();
    const auto dx               = Geom(lev).CellSizeArray();
    const int Mm = Geom(lev).Domain().size()[1];
    auto geomdata = Geom(lev).data();

    for ( MFIter mfi(mf_u, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<Real> const& DC = mf_DC.array(mfi);
        Array4<Real> const& Akv = (vec_Akv[lev])->array(mfi);
        Array4<Real> const& Hz  = (vec_Hz[lev])->array(mfi);
        Array4<Real> const& Huon  = (vec_Huon[lev])->array(mfi);
        Array4<Real> const& Hvom  = (vec_Hvom[lev])->array(mfi);
        Array4<Real> const& z_r = (mf_z_r)->array(mfi);
        Array4<Real> const& uold = (U_old).array(mfi);
        Array4<Real> const& vold = (V_old).array(mfi);
        Array4<Real> const& u = (mf_u).array(mfi);
        Array4<Real> const& v = (mf_v).array(mfi);
        Array4<Real> const& tempold = (mf_tempold).array(mfi);
        Array4<Real> const& saltold = (mf_saltold).array(mfi);
        Array4<Real> const& temp = (mf_temp).array(mfi);
        Array4<Real> const& salt = (mf_salt).array(mfi);
        Array4<Real> const& tempstore = (vec_t3[lev])->array(mfi);
        Array4<Real> const& saltstore = (vec_s3[lev])->array(mfi);
        Array4<Real> const& ru = (mf_ru)->array(mfi);
        Array4<Real> const& rv = (mf_rv)->array(mfi);
        Array4<Real> const& W = (mf_W).array(mfi);
        Array4<Real> const& sustr = (mf_sustr)->array(mfi);
        Array4<Real> const& svstr = (mf_svstr)->array(mfi);

        Box bx = mfi.tilebox();
        //copy the tilebox
        Box gbx1 = bx;
        Box gbx11 = bx;
        Box gbx2 = bx;
        Box gbx3uneven(IntVect(AMREX_D_DECL(bx.smallEnd(0)-3,bx.smallEnd(1)-3,bx.smallEnd(2))),
                       IntVect(AMREX_D_DECL(bx.bigEnd(0)+2,bx.bigEnd(1)+2,bx.bigEnd(2))));
        Box gbx2uneven(IntVect(AMREX_D_DECL(bx.smallEnd(0)-2,bx.smallEnd(1)-2,bx.smallEnd(2))),
                       IntVect(AMREX_D_DECL(bx.bigEnd(0)+1,bx.bigEnd(1)+1,bx.bigEnd(2))));
        Box ubx = surroundingNodes(bx,0);
        Box vbx = surroundingNodes(bx,1);
        //make only gbx be grown to match multifabs
        gbx2.grow(IntVect(2,2,0));
        gbx1.grow(IntVect(1,1,0));
        gbx11.grow(IntVect(1,1,1));

        FArrayBox fab_FC(gbx2,1,amrex::The_Async_Arena());
        FArrayBox fab_FX(gbx2,1,amrex::The_Async_Arena());
        FArrayBox fab_FE(gbx2,1,amrex::The_Async_Arena());
        FArrayBox fab_BC(gbx2,1,amrex::The_Async_Arena());
        FArrayBox fab_CF(gbx2,1,amrex::The_Async_Arena());
        FArrayBox fab_pn(gbx2,1,amrex::The_Async_Arena());
        FArrayBox fab_pm(gbx2,1,amrex::The_Async_Arena());
        FArrayBox fab_on_u(gbx2,1,amrex::The_Async_Arena());
        FArrayBox fab_om_v(gbx2,1,amrex::The_Async_Arena());
        FArrayBox fab_fomn(gbx2,1,amrex::The_Async_Arena());
        FArrayBox fab_oHz(gbx11,1,amrex::The_Async_Arena());

        auto FC=fab_FC.array();
        auto FX=fab_FX.array();
        auto FE=fab_FE.array();
        auto pn=fab_pn.array();
        auto pm=fab_pm.array();
        auto on_u=fab_on_u.array();
        auto om_v=fab_om_v.array();
        auto fomn=fab_fomn.array();

        //From ana_grid.h and metrics.F
        amrex::LoopConcurrentOnCpu(gbx2,
        [=] (int i, int j, int  )
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

        amrex::LoopConcurrentOnCpu(gbx2,
        [=] (int i, int j, int k)
        {
          om_v(i,j,0)=1.0/dxi[0];
          on_u(i,j,0)=1.0/dxi[1];
          Huon(i,j,k)=0.0;
          Hvom(i,j,k)=0.0;
        });
        set_massflux_3d(Box(Huon),1,0,uold,Huon,Hz,on_u,nnew);
        set_massflux_3d(Box(Hvom),0,1,vold,Hvom,Hz,om_v,nnew);
        Real lambda = 1.0;
        //
        //-----------------------------------------------------------------------
        // prestep_t_3d
        //-----------------------------------------------------------------------
        //
        //Test this after advection included in 3d time, consider refactoring to call once per tracer
#if 1
        prestep_t_3d(bx, uold, vold, u, v, tempold, saltold, temp, salt, ru, rv, Hz, Akv, on_u, om_v, Huon, Hvom,
                     pm, pn, W, DC, FC, tempstore, saltstore, FX, FE, z_r, iic, ntfirst, nnew, nstp, nrhs, N,
                          lambda, dt_lev);
        prestep_t_3d(bx, uold, vold, u, v, saltold, saltold, salt, salt, ru, rv, Hz, Akv, on_u, om_v, Huon, Hvom,
                     pm, pn, W, DC, FC, saltstore, saltstore, FX, FE, z_r, iic, ntfirst, nnew, nstp, nrhs, N,
                          lambda, dt_lev);
#endif

        //
        //-----------------------------------------------------------------------
        // prestep_uv_3d
        //-----------------------------------------------------------------------
        //
        prestep_uv_3d(bx, uold, vold, u, v, ru, rv, Hz, Akv, on_u, om_v, Huon, Hvom,
                          pm, pn, W, DC, FC, z_r, sustr, svstr, iic, ntfirst, nnew, nstp, nrhs, N,
                          lambda, dt_lev);

#ifdef UV_COR
        //
        //-----------------------------------------------------------------------
        // coriolis
        //-----------------------------------------------------------------------
        //
        coriolis(bx, uold, vold, ru, rv, Hz, fomn, nrhs);
#endif

        //
        //-----------------------------------------------------------------------
        // rhs_3d
        //-----------------------------------------------------------------------
        //
        rhs_3d(bx, uold, vold, ru, rv, Huon, Hvom, W, FC, nrhs, N);

    } // MFIter
    mf_temp.FillBoundary();
    mf_salt.FillBoundary();
    mf_tempold.FillBoundary();
    mf_saltold.FillBoundary();
    vec_t3[lev]->FillBoundary();
    vec_s3[lev]->FillBoundary();
    vec_Huon[lev]->FillBoundary();
    vec_Hvom[lev]->FillBoundary();

    bool predictor_2d_step=true;
    bool first_2d_step=true;
    int nfast=2;
    int nfast_counter=predictor_2d_step ? nfast : nfast-1;
    int next_indx1 = 0;
    for(int my_iif = 0; my_iif < nfast_counter; my_iif++) {
        first_2d_step=(my_iif==0);
        predictor_2d_step=true;
        advance_2d(lev, mf_u, mf_v, vec_ru[lev], vec_rv[lev],
                   vec_Zt_avg1[lev],
                   vec_DU_avg1[lev], vec_DU_avg2[lev],
                   vec_DV_avg1[lev], vec_DV_avg2[lev],
                   vec_rubar[lev], vec_rvbar[lev], vec_rzeta[lev],
                    vec_ubar[lev],  vec_vbar[lev],  vec_zeta[lev],
                   vec_hOfTheConfusingName[lev], ncomp, dt_lev, predictor_2d_step, first_2d_step, my_iif, nfast, next_indx1);
#if 1
        predictor_2d_step=false;
        advance_2d(lev, mf_u, mf_v, vec_ru[lev], vec_rv[lev],
                   vec_Zt_avg1[lev],
                   vec_DU_avg1[lev], vec_DU_avg2[lev],
                   vec_DV_avg1[lev], vec_DV_avg2[lev],
                   vec_rubar[lev], vec_rvbar[lev], vec_rzeta[lev],
                    vec_ubar[lev],  vec_vbar[lev],  vec_zeta[lev],
                   vec_hOfTheConfusingName[lev], ncomp, dt_lev, predictor_2d_step, first_2d_step, my_iif, nfast, next_indx1);
#endif
    }

    advance_3d(lev, mf_u, mf_v, mf_tempold, mf_saltold, mf_temp, mf_salt, vec_t3[lev], vec_s3[lev], vec_ru[lev], vec_rv[lev],
               vec_DU_avg1[lev], vec_DU_avg2[lev],
               vec_DV_avg1[lev], vec_DV_avg2[lev],
               vec_ubar[lev],  vec_vbar[lev],
               mf_AK, mf_DC,
               mf_Hzk, vec_Akv[lev], vec_Hz[lev], vec_Huon[lev], vec_Hvom[lev], ncomp, N, dt_lev);

    MultiFab::Copy(U_new,mf_u,0,0,U_new.nComp(),IntVect(AMREX_D_DECL(1,1,0)));
    U_new.FillBoundary();

    MultiFab::Copy(V_new,mf_v,0,0,V_new.nComp(),IntVect(AMREX_D_DECL(1,1,0)));
    V_new.FillBoundary();

    mf_temp.FillBoundary();
    mf_salt.FillBoundary();
    mf_tempold.FillBoundary();
    mf_saltold.FillBoundary();
    vec_t3[lev]->FillBoundary();
    vec_s3[lev]->FillBoundary();
    //    MultiFab::Copy(W_new,mf_w,0,0,W_new.nComp(),IntVect(AMREX_D_DECL(1,1,0)));
    //    W_new.FillBoundary();
    //    MultiFab::Copy(mf_W,S_old,Omega_comp,0,mf_W.nComp(),mf_w.nGrowVect());

}
