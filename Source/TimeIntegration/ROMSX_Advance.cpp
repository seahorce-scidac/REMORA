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
    Real time_mt = t_new[lev] - 0.5*dt[lev];
    FillPatch(lev, time, time_mt, dt[lev], lev_old);

    const BoxArray&            ba = S_old.boxArray();
    const DistributionMapping& dm = S_old.DistributionMap();

    int nvars = S_old.nComp();

    // Place-holder for source array -- for now just set to 0
    MultiFab source(ba,dm,nvars,1);
    source.setVal(0.0);

    //This is primarily to make a constant "old" state
    // We don't need to call FillPatch on cons_mf because we have fillpatch'ed S_old above
    // MultiFab cons_mf(ba,dm,nvars,S_old.nGrowVect());
    // MultiFab::Copy(cons_mf,S_old,0,0,S_old.nComp(),S_old.nGrowVect());

    // *****************************************************************
    // Update the cell-centered state and face-based velocity using
    // a time integrator.
    // Inputs:
    //          S_old    (state on cell centers)
    //          U_old    (x-velocity on x-faces)
    //          V_old    (y-velocity on y-faces)
    //          W_old    (z-velocity on z-faces)
    //          source   (source term on cell centers)
    // Outputs:
    //          S_new    (state on cell centers)
    //          U_new    (x-velocity on x-faces)
    //          V_new    (y-velocity on y-faces)
    //          W_new    (z-velocity on z-faces)
    // *****************************************************************

    //-----------------------------------------------------------------------
    //  Time step momentum equation in the XI-direction.
    //-----------------------------------------------------------------------

    //Only used locally, probably should be rearranged into FArrayBox declaration
    MultiFab mf_AK(ba,dm,1,IntVect(2,2,0)); //2d missing j coordinate
    MultiFab mf_DC(ba,dm,1,IntVect(2,2,0)); //2d missing j coordinate
    MultiFab mf_Hzk(ba,dm,1,IntVect(2,2,1)); //2d missing j coordinate
    std::unique_ptr<MultiFab>& mf_Hz = Hz[lev];
    std::unique_ptr<MultiFab>& mf_z_r = z_r[lev];
    //Consider passing these into the advance function or renaming relevant things
    MultiFab mf_u(ba,dm,1,IntVect(2,2,0));
    MultiFab mf_v(ba,dm,1,IntVect(2,2,0));
    MultiFab mf_w(ba,dm,1,IntVect(2,2,0));
    std::unique_ptr<MultiFab>& mf_ru = ru[lev];
    std::unique_ptr<MultiFab>& mf_rv = rv[lev];
    std::unique_ptr<MultiFab>& mf_sustr = sustr[lev];
    std::unique_ptr<MultiFab>& mf_svstr = svstr[lev];
    std::unique_ptr<MultiFab>& mf_ubar = ubar[lev];
    std::unique_ptr<MultiFab>& mf_vbar = vbar[lev];

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

    mf_rw.setVal(0.0);
    mf_W.setVal(0.0);
    U_old.FillBoundary();
    V_old.FillBoundary();

    int ncomp = 1;
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
    //check this////////////
    const int nrhs = ncomp-1;
    const int nnew = ncomp-1;
    const int nstp = ncomp-1;
    const Real Gadv = -0.25;
    auto N = Geom(lev).Domain().size()[2]-1; // Number of vertical "levs" aka, NZ

    const auto dxi              = Geom(lev).InvCellSizeArray();
    //    const auto dx               = Geom(lev).CellSizeArray();
    const int Lm = Geom(lev).Domain().size()[0];
    const int Mm = Geom(lev).Domain().size()[1];
    auto geomdata = Geom(lev).data();

    //
    //-----------------------------------------------------------------------
    // prestep_uv_3d
    //-----------------------------------------------------------------------
    //
    for ( MFIter mfi(mf_u, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<Real> const& AK = mf_AK.array(mfi);
        Array4<Real> const& DC = mf_DC.array(mfi);
        Array4<Real> const& Hzk_arr = mf_Hzk.array(mfi);
        Array4<Real> const& Akv_arr = (Akv[lev])->array(mfi);
        Array4<Real> const& Hz_arr  = (Hz[lev])->array(mfi);
        Array4<Real> const& z_r = (mf_z_r)->array(mfi);
        Array4<Real> const& uold = (U_old).array(mfi);
        Array4<Real> const& vold = (V_old).array(mfi);
        Array4<Real> const& u = (mf_u).array(mfi);
        Array4<Real> const& v = (mf_v).array(mfi);
        Array4<Real> const& w = (mf_w).array(mfi);
        Array4<Real> const& ru_arr = (mf_ru)->array(mfi);
        Array4<Real> const& rv_arr = (mf_rv)->array(mfi);
        Array4<Real> const& rw = (mf_rw).array(mfi);
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
        //make only gbx be grown to match multifabs
        gbx2.grow(IntVect(2,2,0));
        gbx1.grow(IntVect(1,1,0));
        gbx11.grow(IntVect(1,1,1));
        Box gbx=gbx2;

        FArrayBox fab_FC(gbx2,1,amrex::The_Async_Arena);
        FArrayBox fab_BC(gbx2,1,amrex::The_Async_Arena);
        FArrayBox fab_CF(gbx2,1,amrex::The_Async_Arena);
        FArrayBox fab_pn(gbx2,1,amrex::The_Async_Arena);
        FArrayBox fab_pm(gbx2,1,amrex::The_Async_Arena);
        FArrayBox fab_on_u(gbx2,1,amrex::The_Async_Arena);
        FArrayBox fab_om_v(gbx2,1,amrex::The_Async_Arena);
        FArrayBox fab_fomn(gbx2,1,amrex::The_Async_Arena);
        FArrayBox fab_Huon(gbx2,1,amrex::The_Async_Arena);
        FArrayBox fab_Hvom(gbx2,1,amrex::The_Async_Arena);
        FArrayBox fab_oHz(gbx11,1,amrex::The_Async_Arena);

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

        auto FC=fab_FC.array();
        auto BC=fab_BC.array();
        auto CF=fab_CF.array();
        auto pn=fab_pn.array();
        auto pm=fab_pm.array();
        auto on_u=fab_on_u.array();
        auto om_v=fab_om_v.array();
        auto fomn=fab_fomn.array();
        auto Huon=fab_Huon.array();
        auto Hvom=fab_Hvom.array();
        auto oHz_arr=fab_oHz.array();
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
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
          om_v(i,j,0)=1.0/dxi[0];
          on_u(i,j,0)=1.0/dxi[1];
        });

        fab_Huon.setVal(0.0);
        fab_Hvom.setVal(0.0);
        fab_Huxx.setVal(0.0);
        fab_Huee.setVal(0.0);
        fab_Hvxx.setVal(0.0);
        fab_Hvee.setVal(0.0);
        fab_uxx.setVal(0.0);
        fab_uee.setVal(0.0);
        fab_UFx.setVal(0.0);
        fab_UFe.setVal(0.0);
        fab_vxx.setVal(0.0);
        fab_vee.setVal(0.0);
        fab_VFx.setVal(0.0);
        fab_VFe.setVal(0.0);

        Real lambda = 1.0;

        prestep_uv_3d(bx, uold, vold, u, v, ru_arr, rv_arr, Hz_arr, Akv_arr, on_u, om_v, Huon, Hvom,
                          pm, pn, W, DC, FC, z_r, sustr, svstr, iic, ntfirst, nnew, nstp, nrhs, N,
                          lambda, dt_lev);

#ifdef UV_COR
        coriolis(bx, uold, vold, ru_arr, rv_arr, Hz_arr, fomn, nrhs);
#endif

    //
    //-----------------------------------------------------------------------
    // rhs_3d
    //-----------------------------------------------------------------------
    //

        amrex::ParallelFor(gbx1, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
        //should not include grow cells
              uxx(i,j,k)=uold(i-1,j,k,nrhs)-2.0*uold(i,j,k,nrhs)+uold(i+1,j,k,nrhs);
              //neglecting terms about periodicity since testing only periodic for now
              Huxx(i,j,k)=Huon(i-1,j,k)-2.0*Huon(i,j,k)+Huon(i+1,j,k);
            });
        amrex::ParallelFor(gbx1, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
        {
              Real cff;
              Real cff1=uold(i  ,j,k,nrhs)+uold(i+1,j,k,nrhs);
              if (cff1 > 0.0)
                cff=uxx(i,j,k);
              else
                cff=uxx(i+1,j,k);
              UFx(i,j,k)=0.25*(cff1+Gadv*cff)*
                (Huon(i  ,j,k)+
                 Huon(i+1,j,k)+
                 Gadv*0.5*(Huxx(i  ,j,k)+
                           Huxx(i+1,j,k)));
                //should not include grow cells
              uee(i,j,k)=uold(i,j-1,k,nrhs)-2.0*uold(i,j,k,nrhs)+uold(i,j+1,k,nrhs);
            });
        amrex::ParallelFor(gbx1,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
              /////////////MIGHT NEED NEW LOOP HERE
              //neglecting terms about periodicity since testing only periodic for now
              Hvxx(i,j,k)=Hvom(i-1,j,k)-2.0*Hvom(i,j,k)+Hvom(i+1,j,k);
            });
        amrex::ParallelFor(gbx1,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
              Real cff;
              Real cff1=uold(i,j  ,k,nrhs)+uold(i,j-1,k,nrhs);
              Real cff2=Hvom(i,j,k)+Hvom(i-1,j,k);
              if (cff2>0.0)
                cff=uee(i,j-1,k);
              else
                cff=uee(i,j,k);
              UFe(i,j,k)=0.25*(cff1+Gadv*cff)*
                (cff2+Gadv*0.5*(Hvxx(i  ,j,k)+
                                Hvxx(i-1,j,k)));
              vxx(i,j,k)=vold(i-1,j,k,nrhs)-2.0*vold(i,j,k,nrhs)+
                vold(i+1,j,k,nrhs);
              //neglecting terms about periodicity since testing only periodic for now
              Huee(i,j,k)=Huon(i,j-1,k)-2.0*Huon(i,j,k)+Huon(i,j+1,k);
              cff1=vold(i  ,j,k,nrhs)+vold(i-1,j,k,nrhs);
              cff2=Huon(i,j,k)+Huon(i,j-1,k);
              if (cff2>0.0)
                cff=vxx(i-1,j,k);
              else
                cff=vxx(i,j,k);
              VFx(i,j,k)=0.25*(cff1+Gadv*cff)*
                (cff2+Gadv*0.5*(Huee(i,j  ,k)+
                                Huee(i,j-1,k)));
              vee(i,j,k)=vold(i,j-1,k,nrhs)-2.0*vold(i,j,k,nrhs)+
                vold(i,j+1,k,nrhs);
              Hvee(i,j,k)=Hvom(i,j-1,k)-2.0*Hvom(i,j,k)+Hvom(i,j+1,k);
            });
        amrex::ParallelFor(gbx1, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
              //neglecting terms about periodicity since testing only periodic for now
              Real cff;
              Real cff1=vold(i,j  ,k,nrhs)+vold(i,j+1,k,nrhs);
              if (cff1>0.0)
                cff=vee(i,j,k);
              else
                cff=vee(i,j+1,k);
              VFe(i,j,k)=0.25*(cff1+Gadv*cff)*
                    (Hvom(i,j  ,k)+
                     Hvom(i,j+1,k)+
                     Gadv*0.5*(Hvee(i,j  ,k)+
                               Hvee(i,j+1,k)));
            });
        amrex::ParallelFor(gbx1, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
              //
              //  Add in horizontal advection.
              //

              Real cff1=UFx(i,j  ,k)-UFx(i-1,j,k);
              Real cff2=UFe(i,j+1,k)-UFe(i  ,j,k);
              Real cff=cff1+cff2;

              ru_arr(i,j,k,nrhs) -= cff;

              cff1=VFx(i+1,j,k)-VFx(i  ,j,k);
              cff2=VFe(i  ,j,k)-VFe(i,j-1,k);
              cff=cff1+cff2;
              rv_arr(i,j,k,nrhs) -= cff;

              //-----------------------------------------------------------------------
              //  Add in vertical advection.
              //-----------------------------------------------------------------------
              cff1=9.0/16.0;
              cff2=1.0/16.0;
              //              if(i>=0)
              {
              if(k>=1&&k<=N-2)
              {
                      FC(i,j,k)=(cff1*(uold(i,j,k  ,nrhs)+
                             uold(i,j,k+1,nrhs))-
                       cff2*(uold(i,j,k-1,nrhs)+
                             uold(i,j,k+2,nrhs)))*
                      (cff1*(W(i  ,j,k)+
                             W(i-1,j,k))-
                       cff2*(W(i+1,j,k)+
                             W(i-2,j,k)));
              }
              else // this needs to be split up so that the following can be concurent
                {
                  FC(i,j,N)=0.0;
                  FC(i,j,N-1)=(cff1*(uold(i,j,N-1,nrhs)+
                                   uold(i,j,N  ,nrhs))-
                             cff2*(uold(i,j,N-2,nrhs)+
                                   uold(i,j,N  ,nrhs)))*
                            (cff1*(W(i  ,j,N-1)+
                                   W(i-1,j,N-1))-
                             cff2*(W(i+1,j,N-1)+
                                   W(i-2,j,N-1)));
                  FC(i,j,0)=(cff1*(uold(i,j,1,nrhs)+
                                 uold(i,j,2,nrhs))-
                           cff2*(uold(i,j,1,nrhs)+
                                 uold(i,j,3,nrhs)))*
                          (cff1*(W(i  ,j,1)+
                                 W(i-1,j,1))-
                           cff2*(W(i+1,j,1)+
                                 W(i-2,j,1)));
                  //              FC(i,0,-1)=0.0;
                }
              }

              if(k-1>=0) {
                  cff=FC(i,j,k)-FC(i,j,k-1);
              } else {
                  cff=FC(i,j,k);
              }

              ru_arr(i,j,k,nrhs) -= cff;

              //              if(j>=0)
              {
              if (k>=1 && k<=N-2)
              {
                  FC(i,j,k)=(cff1*(vold(i,j,k  ,nrhs)+
                             vold(i,j,k+1,nrhs))-
                       cff2*(vold(i,j,k-1,nrhs)+
                             vold(i,j,k+2,nrhs)))*
                      (cff1*(W(i,j  ,k)+
                             W(i,j-1,k))-
                       cff2*(W(i,j+1,k)+
                             W(i,j-2,k)));
              }
              else // this needs to be split up so that the following can be concurent
              {
                  FC(i,j,N)=0.0;
                  FC(i,j,N-1)=(cff1*(vold(i,j,N-1,nrhs)+
                                   vold(i,j,N  ,nrhs))-
                             cff2*(vold(i,j,N-2,nrhs)+
                                   vold(i,j,N  ,nrhs)))*
                            (cff1*(W(i,j  ,N-1)+
                                   W(i,j-1,N-1))-
                             cff2*(W(i,j+1,N-1)+
                                   W(i,j-2,N-1)));
                  FC(i,j,0)=(cff1*(vold(i,j,1,nrhs)+
                                 vold(i,j,2,nrhs))-
                           cff2*(vold(i,j,1,nrhs)+
                                 vold(i,j,3,nrhs)))*
                          (cff1*(W(i,j  ,1)+
                                 W(i,j-1,1))-
                           cff2*(W(i,j+1,1)+
                                 W(i,j-2,1)));
                  //              FC(i,0,-1)=0.0;
              }

              if(k-1>=0) {
                  cff=FC(i,j,k)-FC(i,j,k-1);
              } else {
                  cff=FC(i,j,k);
              }
              rv_arr(i,j,k,nrhs) -= cff;
              }

            });

        // End rhs3d_tile
    } // MFIter

    advance_2d(lev, mf_u, mf_v, ru[lev], rv[lev],
               DU_avg1[lev], DU_avg2[lev],
               DV_avg1[lev], DV_avg2[lev],
               rubar[lev], rvbar[lev], rzeta[lev],
                ubar[lev],  vbar[lev],  zeta[lev],
               dt_lev);

    advance_3d(lev, mf_u, mf_v, ru[lev], rv[lev],
               DU_avg1[lev], DU_avg2[lev],
               DV_avg1[lev], DV_avg2[lev],
               ubar[lev],  vbar[lev],
               mf_AK, mf_DC,
               mf_Hzk, Akv[lev], Hz[lev], dt_lev);

    MultiFab::Copy(U_new,mf_u,0,0,U_new.nComp(),IntVect(AMREX_D_DECL(1,1,0)));
    U_new.FillBoundary();

    MultiFab::Copy(V_new,mf_v,0,0,V_new.nComp(),IntVect(AMREX_D_DECL(1,1,0)));
    V_new.FillBoundary();

    //    MultiFab::Copy(W_new,mf_w,0,0,W_new.nComp(),IntVect(AMREX_D_DECL(1,1,0)));
    //    W_new.FillBoundary();
    //    MultiFab::Copy(mf_W,S_old,Omega_comp,0,mf_W.nComp(),mf_w.nGrowVect());

}
