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

    const Vector<MultiFab*> mfs = {&vars_old[lev][Vars::cons], &vars_old[lev][Vars::xvel], &vars_old[lev][Vars::yvel], &vars_old[lev][Vars::zvel]};
    FillPatch(lev, time, mfs);

#if 0
    MultiFab* S_crse;
    MultiFab rU_crse, rV_crse, rW_crse;
    // Scratch space for time integrator
    amrex::Vector<amrex::MultiFab> rU_old;
    amrex::Vector<amrex::MultiFab> rU_new;
    amrex::Vector<amrex::MultiFab> rV_old;
    amrex::Vector<amrex::MultiFab> rV_new;
    amrex::Vector<amrex::MultiFab> rW_old;
    amrex::Vector<amrex::MultiFab> rW_new;

    if (lev > 0)
    {
        S_crse = &vars_old[lev-1][Vars::cons];

        MultiFab& U_crse = vars_old[lev-1][Vars::xvel];
        MultiFab& V_crse = vars_old[lev-1][Vars::yvel];
        MultiFab& W_crse = vars_old[lev-1][Vars::zvel];

        rU_crse.define(U_crse.boxArray(), U_crse.DistributionMap(), 1, U_crse.nGrow());
        rV_crse.define(V_crse.boxArray(), V_crse.DistributionMap(), 1, V_crse.nGrow());
        rW_crse.define(W_crse.boxArray(), W_crse.DistributionMap(), 1, W_crse.nGrow());

        VelocityToMomentum(U_crse, U_crse.nGrowVect(),
                           V_crse, V_crse.nGrowVect(),
                           W_crse, W_crse.nGrowVect(),
                          *S_crse,rU_crse,rV_crse,rW_crse);
    }

    const auto& local_ref_ratio = (lev > 0) ? ref_ratio[lev-1] : IntVect(1,1,1);

    InterpFaceRegister ifr;
    if (lev > 0)
    {
        ifr.define(S_old.boxArray(), S_old.DistributionMap(), Geom(lev), local_ref_ratio);
    }
#endif

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
#if 0
                  rU_old[lev], rV_old[lev], rW_old[lev],
                  rU_new[lev], rV_new[lev], rW_new[lev],
                  rU_crse, rV_crse, rW_crse,
#endif
                  source,
                  Geom(lev), dt_lev, time
#if 0
		  ,&ifr);
#else
    );
#endif
}

    // Interface for advancing the data at one level by one "slow" timestep
void ROMSX::romsx_advance(int level,
                          amrex::MultiFab& cons_old,  amrex::MultiFab& cons_new,
                          amrex::MultiFab& xvel_old,  amrex::MultiFab& yvel_old,  amrex::MultiFab& zvel_old,
                          amrex::MultiFab& xvel_new,  amrex::MultiFab& yvel_new,  amrex::MultiFab& zvel_new,
#if 0
                          amrex::MultiFab& xmom_old,  amrex::MultiFab& ymom_old,  amrex::MultiFab& zmom_old,
                          amrex::MultiFab& xmom_new,  amrex::MultiFab& ymom_new,  amrex::MultiFab& zmom_new,
                          amrex::MultiFab& xmom_crse, amrex::MultiFab& ymom_crse, amrex::MultiFab& zmom_crse,
#endif
		          amrex::MultiFab& source,
                          const amrex::Geometry fine_geom,
                          const amrex::Real dt, const amrex::Real time
#if 0
                          ,amrex::InterpFaceRegister* ifr);
#else
                          )
#endif
{

    //-----------------------------------------------------------------------
    //  Time step momentum equation in the XI-direction.
    //-----------------------------------------------------------------------

    //  Time step right-hand-side terms.
    //            u(i,j,k,nnew)=u(i,j,k,nnew)+                                &
    //     &                    DC(i,0)*ru(i,j,k,nrhs)

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
  }
