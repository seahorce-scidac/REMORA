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
    MultiFab::Copy(U_new,U_old,0,0,U_new.nComp(),0);
    MultiFab::Copy(V_new,V_old,0,0,V_new.nComp(),0);
    MultiFab::Copy(W_new,W_old,0,0,W_new.nComp(),0);

    //////////    //pre_step3d corrections to boundaries

    // We need to set these because otherwise in the first call to romsx_advance we may
    //    read uninitialized data on ghost values in setting the bc's on the velocities
    U_new.setVal(1.e34,U_new.nGrowVect());
    V_new.setVal(1.e34,V_new.nGrowVect());
    W_new.setVal(1.e34,W_new.nGrowVect());

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

    // Scratch space for time integrator
    amrex::Vector<amrex::MultiFab> rU_old;
    amrex::Vector<amrex::MultiFab> rU_new;
    amrex::Vector<amrex::MultiFab> rV_old;
    amrex::Vector<amrex::MultiFab> rV_new;
    amrex::Vector<amrex::MultiFab> rW_old;
    amrex::Vector<amrex::MultiFab> rW_new;

    const BoxArray & ba = cons_old.boxArray();
    const DistributionMapping & dm = cons_old.DistributionMap();
    //Only used locally, probably should be rearranged into FArrayBox declaration
    MultiFab mf_AK(ba,dm,1,IntVect(1,1,0)); //2d missing j coordinate
    MultiFab mf_DC(ba,dm,1,IntVect(1,1,0)); //2d missing j coordinate
    MultiFab mf_Hzk(ba,dm,1,IntVect(1,1,0)); //2d missing j coordinate
    std::unique_ptr<MultiFab>& mf_Akv = Akv[level];
    std::unique_ptr<MultiFab>& mf_Hz = Hz[level];

    int ncomp = 1;
    int iic = istep[level] - 1;
    int ntfirst = 1;

    const auto dxi              = Geom(level).InvCellSizeArray();
    for ( MFIter mfi(*(mf_Akv), TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<Real> const& AK = (mf_AK).array(mfi);
	Array4<Real> const& DC = (mf_DC).array(mfi);
	Array4<Real> const& Hzk = (mf_Hzk).array(mfi);
	Array4<Real> const& Akv = (mf_Akv)->array(mfi);
	Array4<Real> const& Hz = (mf_Hz)->array(mfi);
	Box bx = mfi.tilebox();
	bx.grow(IntVect(1,1,0));
	FArrayBox Huon(bx,1,amrex::The_Async_Arena);
	//rhs3d work arrays
	FArrayBox Huxx(bx,1,amrex::The_Async_Arena);
	FArrayBox Huee(bx,1,amrex::The_Async_Arena);
	FArrayBox Hvxx(bx,1,amrex::The_Async_Arena);
	FArrayBox Hvee(bx,1,amrex::The_Async_Arena);
	FArrayBox uxx(bx,1,amrex::The_Async_Arena);
	FArrayBox uee(bx,1,amrex::The_Async_Arena);
	FArrayBox vxx(bx,1,amrex::The_Async_Arena);
	FArrayBox vee(bx,1,amrex::The_Async_Arena);
	FArrayBox UFx(bx,1,amrex::The_Async_Arena);
	FArrayBox UFe(bx,1,amrex::The_Async_Arena);
	FArrayBox VFx(bx,1,amrex::The_Async_Arena);
	FArrayBox VFe(bx,1,amrex::The_Async_Arena);
	amrex::ParallelFor(bx, ncomp,
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
		//		xvel_new(i,j,k)=xvel_new(i,j,k)+
		//		  DC(i,0,0)*ru(i,j,k,nrhs);

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
