#include <REMORA.H>

using namespace amrex;

#ifdef REMORA_USE_FUNWAVE_FORT
#include <REMORA_funwave_Fortran_Interface.H>
#endif

/**
 * @param[in] lev            level of refinement
 * @param[in] time           simulation time at start of step
 * @param[in] dt_lev         baroclinic time step at level
 * @param[in] iteration      iteration in subcycling, if using
 * @param[in] ncycle         total number of subcycles, if using
 */
 void
REMORA::Advance (int lev, Real time, Real dt_lev, int /*iteration*/, int /*ncycle*/)
{
    BL_PROFILE("REMORA::Advance()");

    // Start the correction lev+1 will contribute to, before any flux is added to it.
    if (do_reflux && do_substep && lev < finest_level &&
        solverChoice.coupling_type == CouplingType::two_way) {
        getAdvFluxReg(lev+1)->reset();
    }

    setup_step(lev, time, dt_lev);

    int nfast_counter=nfast + 1;

    //***************************************************
    //Compute fast timestep from dt_lev and ratio
    //***************************************************
    Real dtfast_lev=dt_lev/Real(ndtfast);

    //***************************************************
    //Advance nfast_counter steps of the 2d integrator
    //***************************************************
    for (int my_iif = 0; my_iif < nfast_counter; my_iif++) {
        advance_2d_onestep(lev, dt_lev, dtfast_lev, my_iif, nfast_counter);
    }

#ifdef REMORA_USE_FUNWAVE_FORT
    MultiFab* mf_rhoS = vec_rhoS[lev].get();
    for ( MFIter mfi(*mf_rhoS, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Box bx = mfi.validbox();
        int ims = bx.smallEnd(0);
        int jms = bx.smallEnd(1);
        int kms = bx.smallEnd(2);
        int ime = bx.bigEnd(0);
        int jme = bx.bigEnd(1);
        int kme = bx.bigEnd(2);

        Array4<Real> const& rho_salt = mf_rhoS->array(mfi);

        funwave_advance_c(rho_salt.dataPtr(), ims, ime, jms, jme, kms, kme);
    }
#endif

    //***************************************************
    //Advance one step of the 3d integrator
    //***************************************************
    advance_3d_ml(lev, dt_lev);

    // vert_mean_3d has just replaced the depth mean of the fine 3D velocity with the
    // barotropic transport the parent imposed, while its shear is still the profile the
    // parent handed over at the start of the step. Nothing reconciles the two, so the
    // interface carries a depth-uniform offset that AverageDownTo then pushes onto the
    // parent. timeStepML re-imposes the parent's profile here; timeStep does not.
    if (cf_fill_vel_after && lev > 0) {
        FillPatch(lev, time + dt_lev, *xvel_new[lev], xvel_new, xvel_bc(),
                  BdyVars::u, 0, true, true);
        FillPatch(lev, time + dt_lev, *yvel_new[lev], yvel_new, yvel_bc(),
                  BdyVars::v, 0, true, true);
    }

    //***************************************************
    //Hand the completed step to a finer level. Only timeStep reaches this;
    //timeStepML registers its own coarse data inline.
    //***************************************************
    register_coarse_data(lev, time, dt_lev);
}

/**
 * Store this level's old and new state in the coarse/fine fill patchers so the next finer
 * level can interpolate its contact points to its own sub-times.
 *
 * @param[in] lev            level of refinement
 * @param[in] time           simulation time at start of the step just taken
 * @param[in] dt_lev         baroclinic time step at level
 */
void
REMORA::register_coarse_data (int lev, Real time, Real dt_lev)
{
    // At the end of Advance, once the new state is valid. The parent completes its whole
    // step before any child substep, so {t_old, t_new} brackets every time the child asks for.
    if (lev >= finest_level) { return; }

    if (cf_width > 0) {
        // The parallel copy inside RegisterCoarseData needs filled ghost cells
        cons_old[lev]->FillBoundary(geom[lev].periodicity());
        cons_new[lev]->FillBoundary(geom[lev].periodicity());
        FPr_c[lev].RegisterCoarseData({cons_old[lev], cons_new[lev]}, {time, time + dt_lev});
    }

    if (cf_width >= 0) {
        xvel_old[lev]->FillBoundary(geom[lev].periodicity());
        xvel_new[lev]->FillBoundary(geom[lev].periodicity());
        FPr_u[lev].RegisterCoarseData({xvel_old[lev], xvel_new[lev]}, {time, time + dt_lev});

        yvel_old[lev]->FillBoundary(geom[lev].periodicity());
        yvel_new[lev]->FillBoundary(geom[lev].periodicity());
        FPr_v[lev].RegisterCoarseData({yvel_old[lev], yvel_new[lev]}, {time, time + dt_lev});

        zvel_old[lev]->FillBoundary(geom[lev].periodicity());
        zvel_new[lev]->FillBoundary(geom[lev].periodicity());
        FPr_w[lev].RegisterCoarseData({zvel_old[lev], zvel_new[lev]}, {time, time + dt_lev});

        // ubar and vbar carry their time levels as components of one MultiFab, so there is
        // no old/new pair. These serve the tangential and interior contact points; the
        // normal ones on the interface come from the mass flux below.
        vec_ubar[lev]->FillBoundary(geom[lev].periodicity());
        FPr_ubar[lev].RegisterCoarseData({vec_ubar[lev].get(), vec_ubar[lev].get()},
                                         {time, time + dt_lev});

        vec_vbar[lev]->FillBoundary(geom[lev].periodicity());
        FPr_vbar[lev].RegisterCoarseData({vec_vbar[lev].get(), vec_vbar[lev].get()},
                                         {time, time + dt_lev});

        // ROMS uses only the newest flux record unless TIME_INTERP_FLUX is defined
        // (nesting.F): holding the step average constant over the child's substeps conserves
        // mass across the parent interval, interpolating between two averages does not.
        store_2d_flux(lev);

        MultiFab* Du_old = time_interp_flux ? vec_Dubar_old[lev].get() : vec_Dubar_new[lev].get();
        MultiFab* Dv_old = time_interp_flux ? vec_Dvbar_old[lev].get() : vec_Dvbar_new[lev].get();

        vec_Dubar_new[lev]->FillBoundary(geom[lev].periodicity());
        Du_old->FillBoundary(geom[lev].periodicity());
        FPr_Dubar[lev].RegisterCoarseData({Du_old, vec_Dubar_new[lev].get()},
                                          {time, time + dt_lev});

        vec_Dvbar_new[lev]->FillBoundary(geom[lev].periodicity());
        Dv_old->FillBoundary(geom[lev].periodicity());
        FPr_Dvbar[lev].RegisterCoarseData({Dv_old, vec_Dvbar_new[lev].get()},
                                          {time, time + dt_lev});

        // The free surface has no old/new pair above, so the generic FillPatch hands a
        // subcycled child this level's state frozen at the end of the step, for every substep.
        // ROMS's put_refine2d interpolates the donor between two stored snapshots onto the
        // child's own time, so keep the pair it needs. Zt_avg1 is the end-of-step surface and
        // set_zeta_to_Ztavg puts it in every leapfrog component, so broadcasting it to all
        // three matches both what ROMS stores and what the child reads, whichever component
        // its own knew names.
        if (cf_time_interp_zeta) {
            roll_2d_snapshot(vec_zeta_crse_old, vec_zeta_crse_new, lev,
                             *vec_Zt_avg1[lev], *vec_zeta[lev]);
        }
    }
}

/**
 * Roll a two-snapshot history of one level's end-of-step 2D state: the previous "new" becomes
 * "old" and src becomes the new one, so the pair brackets [t_old, t_new] of the step just
 * taken. src has one component and is broadcast to all of them. On the first call, and
 * whenever the level's layout changes, both snapshots are set to src, making the child's time
 * interpolation a no-op for that step rather than reading uninitialised data -- the same
 * special case ROMS takes when RollingIndex is still zero.
 *
 * @param[inout] old_v   the older snapshot, one entry per level
 * @param[inout] new_v   the newer snapshot, one entry per level
 * @param[in]    lev     level of refinement
 * @param[in]    src     end-of-step value to store, single component
 * @param[in]    like    MultiFab whose layout, component count and ghost width to match
 */
void
REMORA::roll_2d_snapshot (Vector<std::unique_ptr<MultiFab>>& old_v,
                          Vector<std::unique_ptr<MultiFab>>& new_v,
                          int lev, const MultiFab& src, const MultiFab& like)
{
    if (int(old_v.size()) <= lev) { old_v.resize(lev+1); }
    if (int(new_v.size()) <= lev) { new_v.resize(lev+1); }

    const int ncomp = like.nComp();
    const IntVect ng = like.nGrowVect();
    const bool fresh = !new_v[lev] || new_v[lev]->boxArray() != like.boxArray()
                       || new_v[lev]->DistributionMap() != like.DistributionMap();

    if (fresh) {
        old_v[lev].reset(new MultiFab(like.boxArray(), like.DistributionMap(), ncomp, ng));
        new_v[lev].reset(new MultiFab(like.boxArray(), like.DistributionMap(), ncomp, ng));
    } else {
        MultiFab::Copy(*old_v[lev], *new_v[lev], 0, 0, ncomp, ng);
    }

    // src carries one component; every component of the snapshot gets it
    for (int n = 0; n < ncomp; ++n) {
        MultiFab::Copy(*new_v[lev], src, 0, n, 1, amrex::min(ng, src.nGrowVect()));
    }
    new_v[lev]->FillBoundary(geom[lev].periodicity());

    if (fresh) {
        MultiFab::Copy(*old_v[lev], *new_v[lev], 0, 0, ncomp, ng);
    }
}
