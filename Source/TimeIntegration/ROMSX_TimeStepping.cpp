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

    MultiFab::Copy(S_new,S_old,0,0,S_new.nComp(),0);
    MultiFab::Copy(U_new,U_old,0,0,U_new.nComp(),U_new.nGrowVect());
    MultiFab::Copy(V_new,V_old,0,0,V_new.nComp(),V_new.nGrowVect());
    MultiFab::Copy(W_new,W_old,0,0,W_new.nComp(),W_new.nGrowVect());

    //////////    //pre_step3d corrections to boundaries
    /*
    // We need to set these because otherwise in the first call to romsx_advance we may
    //    read uninitialized data on ghost values in setting the bc's on the velocities
    U_new.setVal(1.e34,U_new.nGrowVect());
    V_new.setVal(1.e34,V_new.nGrowVect());
    W_new.setVal(1.e34,W_new.nGrowVect());
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

  print_state(xvel_old,amrex::IntVect(AMREX_D_DECL(3,3,3)));
    //-----------------------------------------------------------------------
    //  Time step momentum equation in the XI-direction.
    //-----------------------------------------------------------------------

    const BoxArray & ba = cons_old.boxArray();
    const DistributionMapping & dm = cons_old.DistributionMap();
    //Only used locally, probably should be rearranged into FArrayBox declaration
    MultiFab mf_AK(ba,dm,1,IntVect(1,1,0)); //2d missing j coordinate
    MultiFab mf_DC(ba,dm,1,IntVect(1,1,0)); //2d missing j coordinate
    MultiFab mf_Hzk(ba,dm,1,IntVect(1,1,0)); //2d missing j coordinate
    std::unique_ptr<MultiFab>& mf_Akv = Akv[level];
    std::unique_ptr<MultiFab>& mf_Hz = Hz[level];
    //Consider passing these into the advance function or renaming relevant things
    MultiFab mf_u(ba,dm,1,IntVect(1,1,0));
    MultiFab mf_v(ba,dm,1,IntVect(1,1,0));
    MultiFab mf_w(ba,dm,1,IntVect(1,1,0));
    MultiFab mf_ru(ba,dm,1,IntVect(1,1,0));
    MultiFab mf_rv(ba,dm,1,IntVect(1,1,0));
    MultiFab mf_rw(ba,dm,1,IntVect(1,1,0));
    MultiFab mf_W(ba,dm,1,IntVect(1,1,0));
    // We need to set these because otherwise in the first call to romsx_advance we may
    //    read uninitialized data on ghost values in setting the bc's on the velocities
    mf_u.setVal(1.e34,mf_u.nGrowVect());
    mf_v.setVal(1.e34,mf_v.nGrowVect());
    mf_w.setVal(1.e34,mf_w.nGrowVect());
    MultiFab::Copy(mf_u,xvel_new,0,0,xvel_new.nComp(),mf_u.nGrowVect());
    MultiFab::Copy(mf_v,yvel_new,0,0,yvel_new.nComp(),mf_v.nGrowVect());
    MultiFab::Copy(mf_w,zvel_new,0,0,zvel_new.nComp(),mf_w.nGrowVect());
    MultiFab::Copy(mf_W,cons_old,Omega_comp,0,mf_W.nComp(),mf_w.nGrowVect());

    int ncomp = 1;
    int iic = istep[level] - 1;
    int ntfirst = 1;
    //check this////////////
    const int nrhs = ncomp-1;
    const int nnew = ncomp-1;
    const Real Gadv = -0.25;
    auto N = Geom(level).ProbLength(2); // Number of vertical "levels" aka, NZ

    const auto dxi              = Geom(level).InvCellSizeArray();
    for ( MFIter mfi(*(mf_Akv), TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<Real> const& AK = (mf_AK).array(mfi);
	Array4<Real> const& DC = (mf_DC).array(mfi);
	Array4<Real> const& Hzk = (mf_Hzk).array(mfi);
	Array4<Real> const& Akv = (mf_Akv)->array(mfi);
	Array4<Real> const& Hz = (mf_Hz)->array(mfi);
	Array4<Real> const& u = (mf_u).array(mfi);
	Array4<Real> const& v = (mf_v).array(mfi);
	Array4<Real> const& w = (mf_w).array(mfi);
	Array4<Real> const& ru = (mf_ru).array(mfi);
	Array4<Real> const& rv = (mf_rv).array(mfi);
	Array4<Real> const& rw = (mf_rw).array(mfi);
	Array4<Real> const& W = (mf_W).array(mfi);

	Box bx = mfi.tilebox();
	//copy the tilebox
	Box gbx = bx;
	//make only gbx be grown to match multifabs
	gbx.grow(IntVect(1,1,0));
	FArrayBox fab_FC(bx,1,amrex::The_Async_Arena);
	FArrayBox fab_pn(gbx,1,amrex::The_Async_Arena);
	FArrayBox fab_pm(gbx,1,amrex::The_Async_Arena);
	FArrayBox fab_on_u(bx,1,amrex::The_Async_Arena);
	FArrayBox fab_om_v(bx,1,amrex::The_Async_Arena);
	FArrayBox fab_Huon(bx,1,amrex::The_Async_Arena);
	FArrayBox fab_Hvom(bx,1,amrex::The_Async_Arena);
	//rhs3d work arrays
	FArrayBox fab_Huxx(bx,1,amrex::The_Async_Arena);
	FArrayBox fab_Huee(bx,1,amrex::The_Async_Arena);
	FArrayBox fab_Hvxx(bx,1,amrex::The_Async_Arena);
	FArrayBox fab_Hvee(bx,1,amrex::The_Async_Arena);
	FArrayBox fab_uxx(bx,1,amrex::The_Async_Arena);
	FArrayBox fab_uee(bx,1,amrex::The_Async_Arena);
	FArrayBox fab_vxx(bx,1,amrex::The_Async_Arena);
	FArrayBox fab_vee(bx,1,amrex::The_Async_Arena);
	FArrayBox fab_UFx(bx,1,amrex::The_Async_Arena);
	FArrayBox fab_UFe(bx,1,amrex::The_Async_Arena);
	FArrayBox fab_VFx(bx,1,amrex::The_Async_Arena);
	FArrayBox fab_VFe(bx,1,amrex::The_Async_Arena);

	auto FC=fab_FC.array();
	auto pn=fab_pn.array();
	auto pm=fab_pm.array();
	auto on_u=fab_on_u.array();
	auto om_v=fab_om_v.array();
	auto Huon=fab_Huon.array();
	auto Hvom=fab_Hvom.array();
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

	amrex::ParallelFor(gbx, ncomp,
	[=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
	    {
	      pm(i,j,0)=dxi[0];
	      pn(i,j,0)=dxi[1];
	    });
	amrex::ParallelFor(bx, ncomp,
	[=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
	    {
	      om_v(i,j,0)=2.0/(pm(i,j-1,0)+pm(i,j,0));
	      //	      printf("%d %d %d %d %15.5g %15.5g %15.5g\n",i,j,k,n,om_v(i,j,0),pm(i,j-1,0),pm(i,j,0));
	      on_u(i,j,0)=2.0/(pn(i-1,j,0)+pn(i,j,0));
	      //	      printf("%d %d %d %d %15.5g %15.5g %15.5g\n",i,j,k,n,on_u(i,j,0),pn(i-1,j,0),pn(i,j,0));
	      //	      printf("%d %d %d %d %15.5g %15.5g %15.5g\n",i,j,k,n,Hz(i,j,k),Hz(i-1,j,k),u(i,j,k,nrhs));
	      //	      printf("%d %d %d %d %15.5g %15.5g %15.5g\n",i,j,k,n,Hz(i,j,k),Hz(i,j-1,k),v(i,j,k,nrhs));
	      //-----------------------------------------------------------------------
	      //  Compute horizontal mass fluxes, Hz*u/n and Hz*v/m.
	      //-----------------------------------------------------------------------
	      Huon(i,j,k)=0.5*(Hz(i,j,k)+Hz(i-1,j,k))*u(i,j,k,nrhs)*   
		on_u(i,j,0);
	      Hvom(i,j,k)=0.5*(Hz(i,j,k)+Hz(i,j-1,k))*v(i,j,k,nrhs)*   
		om_v(i,j,0);
	  	    });
	amrex::ParallelFor(bx, ncomp,
	[=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
	    {
	      //	      if(i==0&&j==1)
		{
		  //		  printf("%d %d %d %d %15.5g %15.5g %15.5g\n",i,j,k,n,u(i-1,j,k,nrhs),-2.0*u(i,j,k,nrhs),u(i+1,j,k,nrhs));
		}
	//should not include grow cells	      
	      uxx(i,j,0)=u(i-1,j,k,nrhs)-2.0*u(i,j,k,nrhs)+u(i+1,j,k,nrhs);
	      //neglecting terms about periodicity since testing only periodic for now
	      Huxx(i,j,0)=Huon(i-1,j,k)-2.0*Huon(i,j,k)+Huon(i+1,j,k);
	      //	      printf("%d %d %d %d %15.5g %15.5g %15.5g\n",i,j,k,n,Huon(i+1,j,k),uxx(i,j,0),Huxx(i,j,0));
	    });
	amrex::ParallelFor(bx, ncomp,
	[=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
	    {
      	      Real cff;
	      Real cff1=u(i  ,j,k,nrhs)+u(i+1,j,k,nrhs);
	      if (cff1 > 0.0)
		cff=uxx(i,j,0);
	      else
		cff=uxx(i+1,j,0);
	      UFx(i,j,0)=0.25*(cff1+Gadv*cff)*
		(Huon(i  ,j,k)+
		 Huon(i+1,j,k)+
		 Gadv*0.5*(Huxx(i  ,j,0)+
			   Huxx(i+1,j,0)));
		//should not include grow cells
	      uee(i,j,0)=u(i,j-1,k,nrhs)-2.0*u(i,j,k,nrhs)+u(i,j+1,k,nrhs);
	      //neglecting terms about periodicity since testing only periodic for now
	      Hvxx(i,j,0)=Hvom(i-1,j,k)-2.0*Hvom(i,j,k)+Hvom(i+1,j,k);
	      cff1=u(i,j  ,k,nrhs)+u(i,j-1,k,nrhs);
	      Real cff2=Hvom(i,j,k)+Hvom(i-1,j,k);
	      if (cff2>0.0)
		cff=uee(i,j-1,0);
	      else
		cff=uee(i,j,0);
	      UFe(i,j,0)=0.25*(cff1+Gadv*cff)*
		(cff2+Gadv*0.5*(Hvxx(i  ,j,0)+
				Hvxx(i-1,j,0)));
	      vxx(i,j,0)=v(i-1,j,k,nrhs)-2.0*v(i,j,k,nrhs)+
		v(i+1,j,k,nrhs);
	      //neglecting terms about periodicity since testing only periodic for now
	      Huee(i,j,0)=Huon(i,j-1,k)-2.0*Huon(i,j,k)+Huon(i,j+1,k);
	      cff1=v(i  ,j,k,nrhs)+v(i-1,j,k,nrhs);
	      cff2=Huon(i,j,k)+Huon(i,j-1,k);
	      if (cff2>0.0)
		cff=vxx(i-1,j,0);
	      else
		cff=vxx(i,j,0);
	      VFx(i,j,0)=0.25*(cff1+Gadv*cff)*
		(cff2+Gadv*0.5*(Huee(i,j  ,0)+
				Huee(i,j-1,0)));
	      vee(i,j,0)=v(i,j-1,k,nrhs)-2.0*v(i,j,k,nrhs)+
		v(i,j+1,k,nrhs);
	      Hvee(i,j,0)=Hvom(i,j-1,k)-2.0*Hvom(i,j,k)+Hvom(i,j+1,k);
	      //neglecting terms about periodicity since testing only periodic for now
	      cff1=v(i,j  ,k,nrhs)+v(i,j+1,k,nrhs);
	      if (cff1>0.0)
		cff=vee(i,j,0);
	      else
		cff=vee(i,j+1,0);
	      VFe(i,j,0)=0.25*(cff1+Gadv*cff)*
                    (Hvom(i,j  ,k)+
                     Hvom(i,j+1,k)+
                     Gadv*0.5*(Hvee(i,j  ,0)+
                               Hvee(i,j+1,0)));
	    });
	amrex::ParallelFor(bx, ncomp,
	[=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
	    {
	      //
	      //  Add in horizontal advection.
	      //

	      Real cff1=UFx(i,j,0)-UFx(i-1,j,0);
	      Real cff2=UFe(i,j+1,0)-UFe(i,j,0);
	      Real cff=cff1+cff2;
	      ru(i,j,k,nrhs)=ru(i,j,k,nrhs)-cff;
	      cff1=VFx(i+1,j,0)-VFx(i,j,0);
	      cff2=VFe(i,j,0)-VFe(i,j-1,0);
	      cff=cff1+cff2;
	      rv(i,j,k,nrhs)=rv(i,j,k,nrhs)-cff;

	      //-----------------------------------------------------------------------
	      //  Add in vertical advection.
	      //-----------------------------------------------------------------------
	      cff1=9.0/16.0;
	      cff2=1.0/16.0;
	      if(k>=2&&k<=N-2)
	      {
		FC(i,0,k)=(cff1*(u(i,j,k  ,nrhs)+
			     u(i,j,k+1,nrhs))-
		       cff2*(u(i,j,k-1,nrhs)+
			     u(i,j,k+2,nrhs)))*
		      (cff1*(W(i  ,j,k)+
			     W(i-1,j,k))-
		       cff2*(W(i+1,j,k)+
			     W(i-2,j,k)));
	      }
	      else if(k==0) // this needs to be split up so that the following can be concurent
		{
		  FC(i,0,N)=0.0;
		  FC(i,0,N-1)=(cff1*(u(i,j,N-1,nrhs)+
				   u(i,j,N  ,nrhs))-
			     cff2*(u(i,j,N-2,nrhs)+
				   u(i,j,N  ,nrhs)))*
		            (cff1*(W(i  ,j,N-1)+
				   W(i-1,j,N-1))-
			     cff2*(W(i+1,j,N-1)+
				   W(i-2,j,N-1)));
		  FC(i,0,1)=(cff1*(u(i,j,1,nrhs)+
				 u(i,j,2,nrhs))-
			   cff2*(u(i,j,1,nrhs)+
				 u(i,j,3,nrhs)))*
		          (cff1*(W(i  ,j,1)+
				 W(i-1,j,1))-
			   cff2*(W(i+1,j,1)+
				 W(i-2,j,1)));
		  FC(i,0,0)=0.0;
		}
	      cff=FC(i,0,k)-FC(i,0,k-1);
	      ru(i,j,k,nrhs)=ru(i,j,k,nrhs)-cff;
	      if(j>=2)
	      {
	      if(k>=2&&k<=N-2)
	      {
		FC(i,0,k)=(cff1*(v(i,j,k  ,nrhs)+
			     v(i,j,k+1,nrhs))-
		       cff2*(v(i,j,k-1,nrhs)+
			     v(i,j,k+2,nrhs)))*
		      (cff1*(W(i,j  ,k)+
			     W(i,j-1,k))-
		       cff2*(W(i,j+1,k)+
			     W(i,j-2,k)));
	      }
	      else if(k==0) // this needs to be split up so that the following can be concurent
		{
		  FC(i,0,N)=0.0;
		  FC(i,0,N-1)=(cff1*(v(i,j,N-1,nrhs)+
				   v(i,j,N  ,nrhs))-
			     cff2*(v(i,j,N-2,nrhs)+
				   v(i,j,N  ,nrhs)))*
		            (cff1*(W(i,j  ,N-1)+
				   W(i,j-1,N-1))-
			     cff2*(W(i,j+1,N-1)+
				   W(i,j-2,N-1)));
		  FC(i,0,1)=(cff1*(v(i,j,1,nrhs)+
				 v(i,j,2,nrhs))-
			   cff2*(v(i,j,1,nrhs)+
				 v(i,j,3,nrhs)))*
		          (cff1*(W(i,j  ,1)+
				 W(i,j-1,1))-
			   cff2*(W(i,j+1,1)+
				 W(i,j-2,1)));
		  FC(i,0,0)=0.0;
		}
	      cff=FC(i,0,k)-FC(i,0,k-1);
	      rv(i,j,k,nrhs)=rv(i,j,k,nrhs)-cff;
	      }
	    });
	amrex::ParallelFor(gbx, ncomp,
	[=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
	        Real cff;
	        AK(i,0,k)=0.5*(Akv(i-1,j,k)+
			       Akv(i  ,j,k));
		if(k!=0)
		  Hzk(i,0,k)=0.5*(Hz(i-1,j,k)+
				  Hz(i  ,j,k));

		if(iic==ntfirst)
		  cff=0.25*dt;
		else if(iic==ntfirst+1)
		  cff=0.25*dt*3.0/2.0;
		else
		  cff=0.25*dt*23.0/12.0;

		DC(i,0,0)=cff*(2*dxi[0])*(2*dxi[1]);
		
		//rhs contributions are in rhs3d.F and are from coriolis, horizontal advection, and vertical advection
		printf("%d %d %d %d %15.5g %15.5g %15.5g\n",i,j,k,n,u(i,j,k),DC(i,0,0),ru(i,j,k,nrhs));
		u(i,j,k)=u(i,j,k)+
		  DC(i,0,0)*ru(i,j,k,nrhs);
	      printf("%d %d %d %d %15.5g %15.5g %15.5g\n",i,j,k,n,u(i,j,k),DC(i,0,0),ru(i,j,k,nrhs));
    //  Time step right-hand-side terms.
    //            u(i,j,k,nnew)=u(i,j,k,nnew)+                                &
    //     &                    DC(i,0)*ru(i,j,k,nrhs)
    
		amrex::Abort("testing");

    //  Backward substitution.
    //            u(i,j,k,nnew)=u(i,j,k,nnew)+cff

    //  Compute new solution by back substitution.
    //            u(i,j,k,nnew)=DC(i,k)

    //  Compute new solution by back substitution.
    //            u(i,j,k,nnew)=DC(i,k)

    //  Couple and update new solution.
    //            u(i,j,k,nnew)=u(i,j,k,nnew)-DC(i,0)
    //# ifdef MASKING
    //            u(i,j,k,nnew)=u(i,j,k,nnew)*umask(i,j)
    //# endif
    //# ifdef WET_DRY
    //            u(i,j,k,nnew)=u(i,j,k,nnew)*umask_wet(i,j)
    //            ru(i,j,k,nrhs)=ru(i,j,k,nrhs)*umask_wet(i,j)
    //# endif

    //-----------------------------------------------------------------------
    //  Time step momentum equation in the ETA-direction.
    //-----------------------------------------------------------------------
    //  Time step right-hand-side terms.
    //                v(i,j,k,nnew)=v(i,j,k,nnew)+DC(i,0)*rv(i,j,k,nrhs)

    //  Backward substitution.
    //            v(i,j,k,nnew)=v(i,j,k,nnew)+cff

    //  Compute new solution by back substitution.
    //            v(i,j,N(ng),nnew)=DC(i,N(ng))

    //  Compute new solution by back substitution.
    //              v(i,j,k,nnew)=DC(i,k)
    //back-substition happens differently if OMEGA_IMPLICIT

    //  Couple and update new solution.
    //              v(i,j,k,nnew)=v(i,j,k,nnew)-DC(i,0)
    //# ifdef MASKING
    //              v(i,j,k,nnew)=v(i,j,k,nnew)*vmask(i,j)
    //# endif
    //# ifdef WET_DRY
    //              v(i,j,k,nnew)=v(i,j,k,nnew)*vmask_wet(i,j)
    //              rv(i,j,k,nrhs)=rv(i,j,k,nrhs)*vmask_wet(i,j)
    //# endif

    //    # ifdef UV_WAVEDRAG
    //                  u(i,j,k,nnew)=u(i,j,k,nstp)-dt(ng)*cff*               &
    // &                (u(i,j,k,nnew)-TIDES(ng)%filt_ubot(i,j))

    //-----------------------------------------------------------------------
    //  Set lateral boundary conditions.
    //-----------------------------------------------------------------------
    //      CALL u3dbc_tile (ng, tile,                                        &
    //      CALL v3dbc_tile (ng, tile,                                        &

    //-----------------------------------------------------------------------
    //  Apply momentum transport point sources (like river runoff), if any.
    //-----------------------------------------------------------------------
    //u(i,j,k,nnew)=SOURCES(ng)%Qsrc(is,k)*cff1
    //v(i,j,k,nnew)=SOURCES(ng)%Qsrc(is,k)*cff1

    //-----------------------------------------------------------------------
    //  Couple 2D and 3D momentum equations.
    //-----------------------------------------------------------------------

    //  Replace only BOUNDARY POINTS incorrect vertical mean with more
    //  accurate barotropic component, ubar=DU_avg1/(D*on_u). Recall that,
    //  D=CF(:,0).

    //  Replace only BOUNDARY POINTS incorrect vertical mean with more
    //  accurate barotropic component, vbar=DV_avg1/(D*om_v).  Recall that,
    //  D=CF(:,0).

    //-----------------------------------------------------------------------
    //  Exchange boundary data.
    //-----------------------------------------------------------------------
	    });
    }
}
