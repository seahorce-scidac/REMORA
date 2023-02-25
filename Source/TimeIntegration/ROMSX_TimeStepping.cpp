#include <ROMSX.H>
#include <Utils.H>

using namespace amrex;

// Advance a level by dt
// includes a recursive call for finer levels
void
ROMSX::timeStep (int lev, Real time, int iteration)
{
    // Update what we call "old" and "new" time
    t_old[lev] = t_new[lev];
    t_new[lev] += dt[lev];

    if (Verbose()) {
        amrex::Print() << "[Level " << lev << " step " << istep[lev]+1 << "] ";
        amrex::Print() << "ADVANCE from time = " << t_old[lev] << " to " << t_new[lev]
                       << " with dt = " << dt[lev] << std::endl;
    }

    // Advance a single level for a single time step
    Advance(lev, time, dt[lev], iteration, nsubsteps[lev]);

    ++istep[lev];

    if (Verbose())
    {
        amrex::Print() << "[Level " << lev << " step " << istep[lev] << "] ";
        amrex::Print() << "Advanced " << CountCells(lev) << " cells" << std::endl;
    }

    if (lev < finest_level)
    {
        // recursive call for next-finer level
        for (int i = 1; i <= nsubsteps[lev+1]; ++i)
        {
            timeStep(lev+1, time+(i-1)*dt[lev+1], i);
        }

        AverageDownTo(lev); // average lev+1 down to lev
    }
}

// advance a single level for a single time step
void
ROMSX::Advance (int lev, Real time, Real dt_lev, int /*iteration*/, int /*ncycle*/)
{
    BL_PROFILE("ROMSX::Advance()");

    // We must swap the pointers so the previous step's "new" is now this step's "old"
    std::swap(vars_old[lev], vars_new[lev]);

    MultiFab& S_old = vars_old[lev][Vars::cons];
    MultiFab& S_new = vars_new[lev][Vars::cons];

    MultiFab& U_old = vars_old[lev][Vars::xvel];
    MultiFab& V_old = vars_old[lev][Vars::yvel];
    MultiFab& W_old = vars_old[lev][Vars::zvel];

    MultiFab& U_new = vars_new[lev][Vars::xvel];
    MultiFab& V_new = vars_new[lev][Vars::yvel];
    MultiFab& W_new = vars_new[lev][Vars::zvel];

    /*
    U_old.setVal(0.e34,U_new.nGrowVect());
    V_old.setVal(0.e34,V_new.nGrowVect());
    W_old.setVal(0.e34,W_new.nGrowVect());
    */
    U_old.FillBoundary();
    V_old.FillBoundary();
    W_old.FillBoundary();
    MultiFab::Copy(S_new,S_old,0,0,S_new.nComp(),S_new.nGrowVect());
    MultiFab::Copy(U_new,U_old,0,0,U_new.nComp(),U_new.nGrowVect());
    MultiFab::Copy(V_new,V_old,0,0,V_new.nComp(),V_new.nGrowVect());
    MultiFab::Copy(W_new,W_old,0,0,W_new.nComp(),W_new.nGrowVect());
    //////////    //pre_step3d corrections to boundaries
    /*
    // We need to set these because otherwise in the first call to romsx_advance we may
    //    read uninitialized data on ghost values in setting the bc's on the velocities
    U_new.setVal(0.e34,U_new.nGrowVect());
    V_new.setVal(0.e34,V_new.nGrowVect());
    W_new.setVal(0.e34,W_new.nGrowVect());
    */
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
    MultiFab cons_mf(ba,dm,nvars,S_old.nGrowVect());
    MultiFab::Copy(cons_mf,S_old,0,0,S_old.nComp(),S_old.nGrowVect());

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

    romsx_advance(lev,
                  cons_mf, S_new,
                  U_old, V_old, W_old,
                  U_new, V_new, W_new,
                  source,
                  Geom(lev), dt_lev, time
    );
}

    // Interface for advancing the data at one level by one "slow" timestep
void ROMSX::romsx_advance(int level,
                          amrex::MultiFab& cons_old,  amrex::MultiFab& cons_new,
                          amrex::MultiFab& xvel_old,  amrex::MultiFab& yvel_old,  amrex::MultiFab& zvel_old,
                          amrex::MultiFab& xvel_new,  amrex::MultiFab& yvel_new,  amrex::MultiFab& zvel_new,
                          amrex::MultiFab& source,
                          const amrex::Geometry fine_geom,
                          const amrex::Real dt, const amrex::Real time
                          )
{

    //-----------------------------------------------------------------------
    //  Time step momentum equation in the XI-direction.
    //-----------------------------------------------------------------------

    const BoxArray & ba = cons_old.boxArray();
    const DistributionMapping & dm = cons_old.DistributionMap();
    //Only used locally, probably should be rearranged into FArrayBox declaration
    MultiFab mf_AK(ba,dm,1,IntVect(2,2,0)); //2d missing j coordinate
    MultiFab mf_DC(ba,dm,1,IntVect(2,2,0)); //2d missing j coordinate
    MultiFab mf_Hzk(ba,dm,1,IntVect(2,2,1)); //2d missing j coordinate
    std::unique_ptr<MultiFab>& mf_Akv = Akv[level];
    std::unique_ptr<MultiFab>& mf_Hz = Hz[level];
    std::unique_ptr<MultiFab>& mf_z_r = z_r[level];
    //Consider passing these into the advance function or renaming relevant things
    MultiFab mf_u(ba,dm,1,IntVect(2,2,0));
    MultiFab mf_v(ba,dm,1,IntVect(2,2,0));
    MultiFab mf_w(ba,dm,1,IntVect(2,2,0));
    std::unique_ptr<MultiFab>& mf_ru = ru[level];
    std::unique_ptr<MultiFab>& mf_rv = rv[level];
    //    MultiFab mf_ru(ba,dm,1,IntVect(2,2,0));
    //    MultiFab mf_rv(ba,dm,1,IntVect(2,2,0));
    MultiFab mf_rw(ba,dm,1,IntVect(2,2,0));
    MultiFab mf_W(ba,dm,1,IntVect(2,2,0));
    // We need to set these because otherwise in the first call to romsx_advance we may
    //    read uninitialized data on ghost values in setting the bc's on the velocities
    mf_u.setVal(0.e34,IntVect(AMREX_D_DECL(1,1,0)));
    mf_v.setVal(0.e34,IntVect(AMREX_D_DECL(1,1,0)));
    mf_w.setVal(0);
    mf_DC.setVal(0);
    mf_w.setVal(0.e34,IntVect(AMREX_D_DECL(1,1,0)));
    MultiFab::Copy(mf_u,xvel_new,0,0,xvel_new.nComp(),IntVect(AMREX_D_DECL(2,2,0)));
    MultiFab::Copy(mf_v,yvel_new,0,0,yvel_new.nComp(),IntVect(AMREX_D_DECL(2,2,0)));
    MultiFab::Copy(mf_w,zvel_new,0,0,zvel_new.nComp(),IntVect(AMREX_D_DECL(2,2,0)));
    MultiFab::Copy(mf_W,cons_old,Omega_comp,0,mf_W.nComp(),IntVect(AMREX_D_DECL(2,2,0)));
    mf_u.FillBoundary();
    mf_v.FillBoundary();
    mf_w.FillBoundary();
    mf_W.FillBoundary();
    //    mf_ru.setVal(0.0);
    //    mf_rv.setVal(0.0);
    mf_rw.setVal(0.0);
    mf_W.setVal(0.0);
    xvel_old.FillBoundary();
    yvel_old.FillBoundary();

    int ncomp = 1;
    int iic = istep[level];
    int ntfirst = 0;
    //If we had wind stress and bottom stress we would need to set these:
    Real tdays;
    Real dstart;
    Real rho0;
    /*
       IF ((tdays(ng)-dstart).le.2.0_r8) THEN
          windamp=-0.1_r8*SIN(pi*(tdays(ng)-dstart)/4.0_r8)/rho0
       ELSE
          windamp=-0.1_r8/rho0
       END IF     */
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
    auto N = Geom(level).Domain().size()[2]-1; // Number of vertical "levels" aka, NZ

    const auto test_point=IntVect(AMREX_D_DECL(1-1,2-1,4-1));
    print_state(xvel_new,test_point);

    const auto dxi              = Geom(level).InvCellSizeArray();
    //    const auto dx               = Geom(level).CellSizeArray();
    const int Lm = Geom(level).Domain().size()[0]-1;
    const int Mm = Geom(level).Domain().size()[1]-1;
    auto geomdata = Geom(level).data();
    for ( MFIter mfi(mf_u, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<Real> const& AK = (mf_AK).array(mfi);
        Array4<Real> const& DC = (mf_DC).array(mfi);
        Array4<Real> const& Hzk = (mf_Hzk).array(mfi);
        Array4<Real> const& Akv = (mf_Akv)->array(mfi);
        Array4<Real> const& Hz = (mf_Hz)->array(mfi);
        Array4<Real> const& z_r = (mf_z_r)->array(mfi);
        Array4<Real> const& uold = (xvel_old).array(mfi);
        Array4<Real> const& vold = (yvel_old).array(mfi);
        //Array4<Real> const& uold = (mf_u).array(mfi);
        //Array4<Real> const& vold = (mf_v).array(mfi);
        Array4<Real> const& u = (mf_u).array(mfi);
        Array4<Real> const& v = (mf_v).array(mfi);
        Array4<Real> const& w = (mf_w).array(mfi);
        Array4<Real> const& ru = (mf_ru)->array(mfi);
        Array4<Real> const& rv = (mf_rv)->array(mfi);
        Array4<Real> const& rw = (mf_rw).array(mfi);
        Array4<Real> const& W = (mf_W).array(mfi);

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
        auto oHz=fab_oHz.array();
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
        amrex::ParallelFor(gbx2, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {

	      const auto prob_lo         = geomdata.ProbLo();
              const auto dx              = geomdata.CellSize();

              pm(i,j,0)=dxi[0];
              pn(i,j,0)=dxi[1];
	      //defined UPWELLING
	      Real f0=-8.26e-5;
	      Real beta=0.0;
	      Real Esize=1000*(Mm+1);
	      Real y = prob_lo[1] + (j + 0.5) * dx[1];
	      Real f=fomn(i,j,0)=f0+beta*(y-.5*Esize);
	      fomn(i,j,0)=f*(1.0/(pm(i,j,0)*pn(i,j,0)));
            });
        amrex::ParallelFor(gbx2, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
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
        amrex::ParallelFor(gbx2, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
              //-----------------------------------------------------------------------
              //  Compute horizontal mass fluxes, Hz*u/n and Hz*v/m.
              //-----------------------------------------------------------------------
                if(k+1<=N)
                if(i-1>=-2)
                    {
                  Huon(i,j,k)=0.5*(Hz(i,j,k)+Hz(i-1,j,k))*uold(i,j,k,nrhs)*
                on_u(i,j,0);
                    }
                else
                    Huon(i,j,k)=(Hz(i,j,k))*uold(i,j,k,nrhs)*
                on_u(i,j,0);
              if(k+1<=N)
                  if(j-1>=-2)
                      {
                      Hvom(i,j,k)=0.5*(Hz(i,j,k)+Hz(i,j-1,k))*vold(i,j,k,nrhs)*
                          om_v(i,j,0);
                      }
                  else
                      Hvom(i,j,k)=(Hz(i,j,k))*vold(i,j,k,nrhs)*
                          om_v(i,j,0);
                    });
        //Need to include pre_step3d.F terms

        //
        //  Weighting coefficient for the newest (implicit) time step derivatives
        //  using either a Crack-Nicolson implicit scheme (lambda=0.5) or a
        //  backward implicit scheme (lambda=1.0).
        //
#if 1
        //  Except the commented out part means its always 1.0
        Real lambda = 1.0;
        amrex::ParallelFor(gbx1, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                Real cff3=dt*(1.0-lambda);
                Real cff, cff1, cff2, cff4;

                if(k<=N&&k>=0)
                {
                    cff=1.0/(z_r(i,j,k+1)+z_r(i-1,j,k+1)-
                             z_r(i,j,k  )-z_r(i-1,j,k  ));
                    FC(i,j,k)=cff3*cff*(uold(i,j,k,nstp)-uold(i,j,k-1,nstp))*
                        (Akv(i,j,k)+Akv(i-1,j,k));
                }
                else
                {
                    //              FC(i,j,-1)=0.0;//dt*bustr(i,j,0);
                    //              FC(i,j,N)=0.0;//dt*sustr(i,j,0);
                }
                cff=dt*.25;
                DC(i,j,k)=cff*(pm(i,j,0)+pm(i-1,j,0))*(pn(i,j,0)+pn(i-1,j,0));
            });
        amrex::ParallelFor(gbx1, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                Real cff3=dt*(1.0-lambda);
                Real cff, cff1, cff2, cff4;

                int indx=0; //nrhs-3

                if(iic==ntfirst)
                {
                    //Hz still might need adjusting
                    if(k+1<=N&&k>=1)
                    {
                        cff1=uold(i,j,k,nstp)*0.5*(Hz(i,j,k)+Hz(i-1,j,k));
                        cff2=FC(i,j,k)-FC(i,j,k-1);
                        u(i,j,k,nnew)=cff1+cff2;
                    }
                    else if(k==0)
                    {
                        cff1=uold(i,j,k,nstp)*0.5*(Hz(i,j,k)+Hz(i-1,j,k));
                        cff2=FC(i,j,k);//-bustr(i,j,0);
                        u(i,j,k,nnew)=cff1+cff2;
                    }
                    else if(k==N)
                    {
                        cff1=uold(i,j,k,nstp)*0.5*(Hz(i,j,k)+Hz(i-1,j,k));
                        cff2=-FC(i,j,k);//+sustr(i,j,0);
                        u(i,j,k,nnew)=cff1+cff2;
                    }
		    if(IntVect(AMREX_D_DECL(i,j,k))==test_point)
			{
			    amrex::Print()<<test_point<<u(i,j,k,nnew)<<std::endl;
			    //amrex::Abort("word");
			}
                }
                else if(iic==ntfirst+1)
                {
                    if(k<N&&k>0) {
                        cff1=uold(i,j,k,nstp)*0.5*(Hz(i,j,k)+Hz(i-1,j,k));                     
                        cff2=FC(i,j,k)-FC(i,j,k-1);
                    }
                    else if(k==0) {
                        cff1=uold(i,j,k,nstp)*0.5*(Hz(i,j,k)+Hz(i-1,j,k));
                        cff2=FC(i,j,k);//-bustr(i,j,0);
                    }
                    else if(k==N) {
                        cff1=uold(i,j,k,nstp)*0.5*(Hz(i,j,k)+Hz(i-1,j,k));
                        cff2=-FC(i,j,k);//+sustr(i,j,0);
                    }
                    cff3=0.5*DC(i,j,k);
                    Real r_swap= ru(i,j,k,indx);
                    indx=nrhs ? 0 : 1;
                    ru(i,j,k,indx) = ru(i,j,k,nrhs);
                    ru(i,j,k,nrhs) = r_swap;
                    u(i,j,k,nnew)=cff1-
                                  cff3*ru(i,j,k,indx)+
                                  cff2;
		    if(IntVect(AMREX_D_DECL(i,j,k))==test_point)
			{
			    amrex::Print()<<test_point<<u(i,j,k,nnew)<<std::endl;
			    //	amrex::Abort("word");
			}
                }
                else
                {
                    cff1= 5.0/12.0;
                    cff2=16.0/12.0;
                    if(k<N&&k>0) {
                        cff3=uold(i,j,k,nstp)*0.5*(Hz(i,j,k)+Hz(i-1,j,k));
                        cff4=FC(i,j,k)-FC(i,j,k-1);
                    }
                    else if(k==0) {
                        cff3=uold(i,j,k,nstp)*0.5*(Hz(i,j,k)+Hz(i-1,j,k));
                        cff4=FC(i,j,k);//-bustr(i,j,0);
                    }
                    else if(k==N) {
                        cff3=uold(i,j,k,nstp)*0.5*(Hz(i,j,k)+Hz(i-1,j,k));
                        cff4=-FC(i,j,k);//+sustr(i,j,0);
                    }
                    Real r_swap= ru(i,j,k,indx);
                    indx=nrhs ? 0 : 1;
                    ru(i,j,k,indx) = ru(i,j,k,nrhs);
                    ru(i,j,k,nrhs) = r_swap;
                    u(i,j,k,nnew)=cff3+
                        DC(i,j,k)*(cff1*ru(i,j,k,indx)+
                                   cff2*ru(i,j,k,nrhs))+
                        cff4;
		    if(IntVect(AMREX_D_DECL(i,j,k))==test_point)
			{
			    amrex::Print()<<test_point<<u(i,j,k,nnew)<<std::endl;
			    //	amrex::Abort("word");
			}
		  }
            });
	lambda = 1.0;
        amrex::ParallelFor(gbx1, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                Real cff3=dt*(1.0-lambda);
                Real cff, cff1, cff2, cff4;

                if(k<=N&&k>=0)
                {
                    cff=1.0/(z_r(i,j,k+1)+z_r(i,j-1,k+1)-
                             z_r(i,j,k  )-z_r(i,j-1,k  ));
                    FC(i,j,k)=cff3*cff*(vold(i,j,k,nstp)-vold(i,j,k-1,nstp))*
                        (Akv(i,j,k)+Akv(i,j-1,k));
                }
                else
                {
                    //              FC(i,j,-1)=0.0;//dt*bustr(i,j,0);
                    //              FC(i,j,N)=0.0;//dt*sustr(i,j,0);
                }
                cff=dt*.25;
                DC(i,j,k)=cff*(pm(i,j,0)+pm(i,j-1,0))*(pn(i,j,0)+pn(i,j-1,0));
            });
        amrex::ParallelFor(gbx1, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                Real cff3=dt*(1.0-lambda);
                Real cff, cff1, cff2, cff4;

                int indx=0; //nrhs-3

                if(iic==ntfirst)
                {
                    //Hz still might need adjusting
                    if(k+1<=N&&k>=1)
                    {
                        cff1=vold(i,j,k,nstp)*0.5*(Hz(i,j,k)+Hz(i,j-1,k));
                        cff2=FC(i,j,k)-FC(i,j,k-1);
                        v(i,j,k,nnew)=cff1+cff2;
                    }
                    else if(k==0)
                    {
                        cff1=vold(i,j,k,nstp)*0.5*(Hz(i,j,k)+Hz(i,j-1,k));
                        cff2=FC(i,j,k);//-bustr(i,j,0);
                        v(i,j,k,nnew)=cff1+cff2;
                    }
                    else if(k==N)
                    {
                        cff1=vold(i,j,k,nstp)*0.5*(Hz(i,j,k)+Hz(i,j-1,k));
                        cff2=-FC(i,j,k);//+sustr(i,j,0);
                        v(i,j,k,nnew)=cff1+cff2;
                    }
		    if(IntVect(AMREX_D_DECL(i,j,k))==test_point)
			{
			    amrex::Print()<<test_point<<v(i,j,k,nnew)<<"v"<<std::endl;
			    //	amrex::Abort("word");
			}
                }
                else if(iic==ntfirst+1)
                {
                    if(k<N&&k>0) {
                        cff1=vold(i,j,k,nstp)*0.5*(Hz(i,j,k)+Hz(i,j-1,k));
                        cff2=FC(i,j,k)-FC(i,j,k-1);
                    }
                    else if(k==0) {
                        cff1=vold(i,j,k,nstp)*0.5*(Hz(i,j,k)+Hz(i,j-1,k));
                        cff2=FC(i,j,k);//-bustr(i,j,0);
                    }
                    else if(k==N) {
                        cff1=vold(i,j,k,nstp)*0.5*(Hz(i,j,k)+Hz(i,j-1,k));
                        cff2=-FC(i,j,k);//+sustr(i,j,0);
                    }
                    cff3=0.5*DC(i,j,k);
                    Real r_swap= rv(i,j,k,indx);
                    indx=nrhs ? 0 : 1;
                    rv(i,j,k,indx) = rv(i,j,k,nrhs);
                    rv(i,j,k,nrhs) = r_swap;
                    v(i,j,k,nnew)=cff1-
                                  cff3*rv(i,j,k,indx)+
                                  cff2;
		    		    if(IntVect(AMREX_D_DECL(i,j,k))==test_point)
			{
			    amrex::Print()<<test_point<<vold(i,j,k,nnew)<<"vold "<<std::endl;
			    amrex::Print()<<test_point<<v(i,j,k,nnew)<<"v "<<std::endl;
			    //	amrex::Abort("word");
			}
                }
                else
                {
                    cff1= 5.0/12.0;
                    cff2=16.0/12.0;
                    if(k<N&&k>0) {
                        cff3=vold(i,j,k,nstp)*0.5*(Hz(i,j,k)+Hz(i,j-1,k));
                        cff4=FC(i,j,k)-FC(i,j,k-1);
                    }
                    else if(k==0) {
                        cff3=vold(i,j,k,nstp)*0.5*(Hz(i,j,k)+Hz(i,j-1,k));
                        cff4=FC(i,j,k);//-bustr(i,j,0);
                    }
                    else if(k==N) {
                        cff3=vold(i,j,k,nstp)*0.5*(Hz(i,j,k)+Hz(i,j-1,k));
                        cff4=-FC(i,j,k);//+sustr(i,j,0);
                    }
                    Real r_swap= rv(i,j,k,indx);
                    indx=nrhs ? 0 : 1;
                    rv(i,j,k,indx) = rv(i,j,k,nrhs);
                    rv(i,j,k,nrhs) = r_swap;
                    v(i,j,k,nnew)=cff3+
                        DC(i,j,k)*(cff1*rv(i,j,k,indx)+
                                   cff2*rv(i,j,k,nrhs))+
                        cff4;
		    		    if(IntVect(AMREX_D_DECL(i,j,k))==test_point)
			{
			    amrex::Print()<<test_point<<v(i,j,k,nnew)<<"v"<<std::endl;
			    //	amrex::Abort("word");
			}
		}
            });
#endif  
#ifdef UV_COR
	//
	//-----------------------------------------------------------------------
	//  Add in Coriolis terms.
	//-----------------------------------------------------------------------
        //
        amrex::ParallelFor(gbx1, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
		Real cff=0.5*Hz(i,j,k)*fomn(i,j,0);
		UFx(i,j,k)=cff*(vold(i,j  ,k,nrhs)+
				vold(i,j+1,k,nrhs));
		VFe(i,j,k)=cff*(uold(i  ,j,k,nrhs)+
				uold(i+1,j,k,nrhs));
		
	    });
	amrex::ParallelFor(gbx1, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
		Real cff1=0.5*(UFx(i,j,k)-UFx(i-1,j,k));
              		    if(IntVect(AMREX_D_DECL(i,j,k))==test_point)
			{
			    ///None of the vs are correct in the second step, which seems to make the u grow incorrectly
			    amrex::Print()<<test_point<<vold(i,j,k,nnew)<<"cor "<<vold(i,j+1,k,nnew)<<" "<<v(i,j,k,nrhs)<<std::endl;
			    amrex::Print()<<test_point<<fomn(i,j,0)<<"cor "<<Hz(i,j,k)<<std::endl;
			    amrex::Print()<<test_point<<ru(i,j,k,nnew)<<"cor "<<std::endl;
			    //	amrex::Abort("word");
			}

		ru(i,j,k,nrhs)=ru(i,j,k,nrhs)+cff1;

		cff1=0.5*(VFe(i,j,k)+VFe(i,j-1,k));

		rv(i,j,k,nrhs)=rv(i,j,k,nrhs)-cff1;
	      if(IntVect(AMREX_D_DECL(i,j,k))==test_point)
		  {

			    amrex::Print()<<test_point<<VFx(i+1,j,k)<<"vVFx1 "<<std::endl;
			    amrex::Print()<<test_point<<VFx(i,j,k)<<"vVFx "<<std::endl;
			    amrex::Print()<<test_point<<VFe(i,j,k)<<"vVFe "<<std::endl;
			    amrex::Print()<<test_point<<VFe(i,j-1,k)<<"vVFe1  "<<std::endl;
			    amrex::Print()<<test_point<<rv(i,j,k,nrhs)<<"rv "<<std::endl;
			    //   amrex::Abort("rv");
			}		

            });
#endif
        //Need to include pre_step3d.F terms
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
        amrex::ParallelFor(gbx1, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
              /////////////MIGHT NEED NEW LOOP HERE
              //neglecting terms about periodicity since testing only periodic for now
              Hvxx(i,j,k)=Hvom(i-1,j,k)-2.0*Hvom(i,j,k)+Hvom(i+1,j,k);
            });
        amrex::ParallelFor(gbx1, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
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

              Real cff1=UFx(i,j,k)-UFx(i-1,j,k);
              Real cff2=UFe(i,j+1,k)-UFe(i,j,k);
              Real cff=cff1+cff2;

              ru(i,j,k,nrhs)=ru(i,j,k,nrhs)-cff;

              cff1=VFx(i+1,j,k)-VFx(i,j,k);
              cff2=VFe(i,j,k)-VFe(i,j-1,k);
              cff=cff1+cff2;
	      if(IntVect(AMREX_D_DECL(i,j,k))==test_point)
		  {

			    amrex::Print()<<test_point<<VFx(i+1,j,k)<<"vVFx1 "<<std::endl;
			    amrex::Print()<<test_point<<VFx(i,j,k)<<"vVFx "<<std::endl;
			    amrex::Print()<<test_point<<VFe(i,j,k)<<"vVFe "<<std::endl;
			    amrex::Print()<<test_point<<VFe(i,j-1,k)<<"vVFe1  "<<std::endl;
			    amrex::Print()<<test_point<<rv(i,j,k,nrhs)<<"rv "<<std::endl;
			    //amrex::Abort("rv");
			}
              rv(i,j,k,nrhs)=rv(i,j,k,nrhs)-cff;
	      if(IntVect(AMREX_D_DECL(i,j,k))==test_point)
		  {

			    amrex::Print()<<test_point<<VFx(i+1,j,k)<<"vVFx1 "<<std::endl;
			    amrex::Print()<<test_point<<VFx(i,j,k)<<"vVFx "<<std::endl;
			    amrex::Print()<<test_point<<VFe(i,j,k)<<"vVFe "<<std::endl;
			    amrex::Print()<<test_point<<VFe(i,j-1,k)<<"vVFe1  "<<std::endl;
			    amrex::Print()<<test_point<<rv(i,j,k,nrhs)<<"rv "<<std::endl;
			    //amrex::Abort("rv");
			}
              //-----------------------------------------------------------------------
              //  Add in vertical advection.
              //-----------------------------------------------------------------------
              cff1=9.0/16.0;
              cff2=1.0/16.0;
              if(i>=0)
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

              if(k-1>=0)
                  cff=FC(i,j,k)-FC(i,j,k-1);
              else
                  cff=FC(i,j,k);

              ru(i,j,k,nrhs)=ru(i,j,k,nrhs)-cff;

              if(j>=0)
              {
              if(k>=1&&k<=N-2)
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
              if(k-1>=0)
                  cff=FC(i,j,k)-FC(i,j,k-1);
              else
                  cff=FC(i,j,k);
              rv(i,j,k,nrhs)=rv(i,j,k,nrhs)-cff;
              }

            });

        // End rhs3d_tile
        // Need to include uv3dmix
        // Begin step3d_uv.F
        amrex::ParallelFor(gbx1, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                    AK(i,j,k)=0.5*(Akv(i-1,j,k)+
                                   Akv(i  ,j,k));
            });
        amrex::ParallelFor(gbx11, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                    Hzk(i,j,k)=0.5*(Hz(i-1,j,k)+
                                    Hz(i  ,j,k));
            });
        amrex::ParallelFor(gbx11, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                oHz(i,j,k) = 1.0/Hzk(i,j,k);
            });
        amrex::ParallelFor(gbx1, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                Real cff;
                if(iic==ntfirst)
                  cff=0.25*dt;
                else if(iic==ntfirst+1)
                  cff=0.25*dt*3.0/2.0;
                else
                  cff=0.25*dt*23.0/12.0;
                DC(i,j,k)=cff*(pm(i,j,0)+pm(i-1,j,0))*(pn(i,j,0)+pn(i-1,j,0));

                u(i,j,k)=u(i,j,k)+
                         DC(i,j,k)*ru(i,j,k,nrhs);
		    if(IntVect(AMREX_D_DECL(i,j,k))==test_point)
			{
			    amrex::Print()<<test_point<<u(i,j,k,nnew)<<std::endl;
			    //amrex::Abort("word");
			}
                v(i,j,k)=v(i,j,k)+
                         DC(i,j,k)*rv(i,j,k,nrhs);
		    		    if(IntVect(AMREX_D_DECL(i,j,k))==test_point)
			{
			    amrex::Print()<<test_point<<DC(i,j,k)<<"vDC "<<std::endl;
			    amrex::Print()<<test_point<<rv(i,j,k,nrhs)<<"rv "<<std::endl;
			    amrex::Print()<<test_point<<vold(i,j,k,nnew)<<"vold "<<std::endl;
			    amrex::Print()<<test_point<<v(i,j,k,nnew)<<"v "<<std::endl;
			    //	amrex::Abort("word");
			}
                //ifdef SPLINES_VVISC is true
                u(i,j,k)=u(i,j,k)*oHz(i,j,k);
		    if(IntVect(AMREX_D_DECL(i,j,k))==test_point)
			{
			    amrex::Print()<<test_point<<u(i,j,k,nnew)<<std::endl;
			    //amrex::Abort("word");
			}
                v(i,j,k)=v(i,j,k)*oHz(i,j,k);
		    		    if(IntVect(AMREX_D_DECL(i,j,k))==test_point)
			{
			    amrex::Print()<<test_point<<oHz(i,j,k)<<"oHz "<<std::endl;
			    amrex::Print()<<test_point<<vold(i,j,k,nnew)<<"vold "<<std::endl;
			    amrex::Print()<<test_point<<v(i,j,k,nnew)<<"v "<<std::endl;
			    //	amrex::Abort("word");
			}
            });
        // End previous
	#if 1
       // Begin vertical viscosity term
       //should be gbx1, but need to fix some bounds inside this loop:
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
               	   FC(i,j,k)=cff1*Hzk(i,j,k  )-dt*AK(i,j,k-1)*oHz(i,j,k  );
	       else
               	   FC(i,j,k)=cff1*Hzk(i,j,k  );
               if(k<=N-1)
                {
                    CF(i,j,k)=cff1*Hzk(i,j,k+1)-dt*AK(i,j,k+1)*oHz(i,j,k+1);
                }
             if(IntVect(AMREX_D_DECL(i,j,k))==test_point)
                 {
                     printf("%d %d %d %d %15.15g %15.15g %15.15g\n",i,j,k,n,cff1*Hzk(i,j,k),dt*AK(i,j,k-1)*oHz(i,j,k  ),FC(i,j,k));
                     amrex::Print()<<"splines"<<std::endl;
                     printf("%d %d %d %d %15.15g %15.15g %15.15g %15.5g\n",i,j,k,n,dt,AK(i,j,k-1),oHz(i,j,k  ),FC(i,j,k));
                     //              amrex::Abort("STOP");
               }

               //
               //  LU decomposition and forward substitution.
               //
               cff1=1.0/3.0;
               if(k==0)
               {
                   BC(i,j,k)=cff1*(Hzk(i,j,k)+Hzk(i,j,k+1))+
                       dt*AK(i,j,k)*(oHz(i,j,k)+oHz(i,j,k+1));
                   cff=1.0/(BC(i,j,k)-FC(i,j,k)*0.0);
                   CF(i,j,k)=cff*CF(i,j,k);
                   DC(i,j,k)=cff*(uold(i,j,k+1,nnew)-uold(i,j,k,nnew)-
                                  FC(i,j,k)*0.0);
               }
	       if(k+1<=N&&k>=1)
               {
                   BC(i,j,k)=cff1*(Hzk(i,j,k)+Hzk(i,j,k+1))+
                       dt*AK(i,j,k)*(oHz(i,j,k)+oHz(i,j,k+1));
                   cff=1.0/(BC(i,j,k)-FC(i,j,k)*CF(i,j,k-1));
                   CF(i,j,k)=cff*CF(i,j,k);
                   DC(i,j,k)=cff*(uold(i,j,k+1,nnew)-uold(i,j,k,nnew)-
                                  FC(i,j,k)*DC(i,j,k-1));
               }

               if(IntVect(AMREX_D_DECL(i,j,k))==test_point)
                  {
	       printf("%d %d %d %d %15.15g %15.15g %15.15g\n",i,j,k,n,cff1,Hzk(i,j,k),Hzk(i,j,k+1));
	       printf("%d %d %d %d %15.15g %15.15g %15.15g %15.15g\n",i,j,k,n,dt,AK(i,j,k),oHz(i,j,k),oHz(i,j,k+1));
                     printf("%d %d %d %d %15.15g %15.15g %15.15g\n",i,j,k,n,uold(i,j,k+1),uold(i,j,k),DC(i,j,k));
		     
                     printf("%d %d %d %d %15.15g %15.15g %15.15g\n",i,j,k,n,cff,FC(i,j,k),DC(i,j,k));
                     printf("%d %d %d %d %15.15g %15.15g %15.15g\n",i,j,k,n,BC(i,j,k),CF(i,j,k),DC(i,j,k));
               }
	    });
       amrex::ParallelFor(gbx1, ncomp,
       [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
               //
               //  Backward substitution.
               //
               DC(i,j,N)=0.0;
	       
               if(N-k+1<=N&&N-k>=0) //-N,1,-1 => kidx =N-k+1
               {
		   //		   amrex::Print()<<"index k: "<<k<<"corresponds to : "<<N-k<<"prev: "<<DC(i,j,N-k+1)<<std::endl;
                   if(N-k+1<0||N-k+2<0)
                       amrex::Abort("-1 here");
                   DC(i,j,N-k)=DC(i,j,N-k)-CF(i,j,N-k)*DC(i,j,N-k+1);
		   //amrex::Print()<<"index k: "<<k<<"corresponds to : "<<N-k<<"cur:  "<<DC(i,j,N-k)<<std::endl;
                   //              DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1);
               }
	    });

       amrex::ParallelFor(gbx1, ncomp,
       [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
               DC(i,j,k)=DC(i,j,k)*AK(i,j,k);
	    });
       amrex::ParallelFor(gbx1, ncomp,
       [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
		Real cff;
               if(k-1>=0)
                   cff=dt*oHz(i,j,k)*(DC(i,j,k)-DC(i,j,k-1));
               else
                   cff=dt*oHz(i,j,k)*(DC(i,j,k));
	       u(i,j,k)=u(i,j,k)+cff;
	       		    if(IntVect(AMREX_D_DECL(i,j,k))==test_point)
			{
			    amrex::Print()<<test_point<<u(i,j,k,nnew)<<std::endl;
			    //		amrex::Abort("backword");
			}
             if(i==3-1&&j==3-1&&k==1-1)
                 {
                     printf("%d %d %d %d %15.15g %15.15g %15.15g\n",i,j,k,n,DC(i,j,k),cff,u(i,j,k));
		     //amrex::Abort("STOP");
                     }
            });
#endif
	#if 1
       // Begin vertical viscosity term
       //should be gbx1, but need to fix some bounds inside this loop:
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
               	   FC(i,j,k)=cff1*Hzk(i,j,k  )-dt*AK(i,j,k-1)*oHz(i,j,k  );
	       else
               	   FC(i,j,k)=cff1*Hzk(i,j,k  );
               if(k<=N-1)
                {
                    CF(i,j,k)=cff1*Hzk(i,j,k+1)-dt*AK(i,j,k+1)*oHz(i,j,k+1);
                }
             if(IntVect(AMREX_D_DECL(i,j,k))==test_point)
                 {
                     printf("%d %d %d %d %15.15g %15.15g %15.15g\n",i,j,k,n,cff1*Hzk(i,j,k),dt*AK(i,j,k-1)*oHz(i,j,k  ),FC(i,j,k));
                     amrex::Print()<<"splines"<<std::endl;
                     printf("%d %d %d %d %15.15g %15.15g %15.15g %15.5g\n",i,j,k,n,dt,AK(i,j,k-1),oHz(i,j,k  ),FC(i,j,k));
                     //              amrex::Abort("STOP");
               }

               //
               //  LU decomposition and forward substitution.
               //
               cff1=1.0/3.0;
               if(k==0)
               {
                   BC(i,j,k)=cff1*(Hzk(i,j,k)+Hzk(i,j,k+1))+
                       dt*AK(i,j,k)*(oHz(i,j,k)+oHz(i,j,k+1));
                   cff=1.0/(BC(i,j,k)-FC(i,j,k)*0.0);
                   CF(i,j,k)=cff*CF(i,j,k);
                   DC(i,j,k)=cff*(vold(i,j,k+1,nnew)-vold(i,j,k,nnew)-
                                  FC(i,j,k)*0.0);
               }
	       if(k+1<=N&&k>=1)
               {
                   BC(i,j,k)=cff1*(Hzk(i,j,k)+Hzk(i,j,k+1))+
                       dt*AK(i,j,k)*(oHz(i,j,k)+oHz(i,j,k+1));
                   cff=1.0/(BC(i,j,k)-FC(i,j,k)*CF(i,j,k-1));
                   CF(i,j,k)=cff*CF(i,j,k);
                   DC(i,j,k)=cff*(vold(i,j,k+1,nnew)-vold(i,j,k,nnew)-
                                  FC(i,j,k)*DC(i,j,k-1));
               }

               if(IntVect(AMREX_D_DECL(i,j,k))==test_point)
                  {
	       printf("%d %d %d %d %15.15g %15.15g %15.15g\n",i,j,k,n,cff1,Hzk(i,j,k),Hzk(i,j,k+1));
	       printf("%d %d %d %d %15.15g %15.15g %15.15g %15.15g\n",i,j,k,n,dt,AK(i,j,k),oHz(i,j,k),oHz(i,j,k+1));
                     printf("%d %d %d %d %15.15g %15.15g %15.15g\n",i,j,k,n,uold(i,j,k+1),uold(i,j,k),DC(i,j,k));
		     
                     printf("%d %d %d %d %15.15g %15.15g %15.15g\n",i,j,k,n,cff,FC(i,j,k),DC(i,j,k));
                     printf("%d %d %d %d %15.15g %15.15g %15.15g\n",i,j,k,n,BC(i,j,k),CF(i,j,k),DC(i,j,k));
               }
	    });
       amrex::ParallelFor(gbx1, ncomp,
       [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
               //
               //  Backward substitution.
               //
               DC(i,j,N)=0.0;
	       
               if(N-k+1<=N&&N-k>=0) //-N,1,-1 => kidx =N-k+1
               {
		   //		   amrex::Print()<<"index k: "<<k<<"corresponds to : "<<N-k<<"prev: "<<DC(i,j,N-k+1)<<std::endl;
                   if(N-k+1<0||N-k+2<0)
                       amrex::Abort("-1 here");
                   DC(i,j,N-k)=DC(i,j,N-k)-CF(i,j,N-k)*DC(i,j,N-k+1);
		   //amrex::Print()<<"index k: "<<k<<"corresponds to : "<<N-k<<"cur:  "<<DC(i,j,N-k)<<std::endl;
                   //              DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1);
               }
	    });

       amrex::ParallelFor(gbx1, ncomp,
       [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
               DC(i,j,k)=DC(i,j,k)*AK(i,j,k);
	    });
       amrex::ParallelFor(gbx1, ncomp,
       [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
		Real cff;
               if(k-1>=0)
                   cff=dt*oHz(i,j,k)*(DC(i,j,k)-DC(i,j,k-1));
               else
                   cff=dt*oHz(i,j,k)*(DC(i,j,k));
	       v(i,j,k)=v(i,j,k)+cff;
	       		    if(IntVect(AMREX_D_DECL(i,j,k))==test_point)
			{
			    amrex::Print()<<test_point<<u(i,j,k,nnew)<<std::endl;
			    //		amrex::Abort("backword");
			}
             if(i==3-1&&j==3-1&&k==1-1)
                 {
                     printf("%d %d %d %d %15.15g %15.15g %15.15g\n",i,j,k,n,DC(i,j,k),cff,u(i,j,k));
		     //amrex::Abort("STOP");
                     }
            });
#endif
    }
    print_state(xvel_new,test_point);
    MultiFab::Copy(xvel_new,mf_u,0,0,xvel_new.nComp(),IntVect(AMREX_D_DECL(1,1,0)));
    xvel_new.FillBoundary();
    print_state(xvel_new,test_point);
    MultiFab::Copy(yvel_new,mf_v,0,0,yvel_new.nComp(),IntVect(AMREX_D_DECL(1,1,0)));
    yvel_new.FillBoundary();
    print_state(xvel_new,test_point);
    //    MultiFab::Copy(zvel_new,mf_w,0,0,zvel_new.nComp(),IntVect(AMREX_D_DECL(1,1,0)));
    //    zvel_new.FillBoundary();
    //    MultiFab::Copy(mf_W,cons_old,Omega_comp,0,mf_W.nComp(),mf_w.nGrowVect());

}
