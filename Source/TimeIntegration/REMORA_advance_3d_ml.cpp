#include <REMORA.H>

using namespace amrex;

/**
 * @param[in] lev     refinement level
 * @param[in] dt_lev  time step at refinement level
 */
void REMORA::advance_3d_ml (int lev, Real dt_lev)
{
    BL_PROFILE("REMORA::advance_3d_ml()");
    // Fill in three ways: 1) interpolate from coarse grid if lev > 0; 2) fill from physical boundaries;
    //                     3) fine-fine fill of ghost cells with FillBoundary call
    //FillPatch(lev, t_old[lev], *cons_new[lev], cons_new, BCVars::cons_bc, BdyVars::t);
    amrex::FillBoundary(Vector<MultiFab*>{xvel_new[lev], yvel_new[lev]}, geom[lev].periodicity());
    //FillPatch(lev, t_old[lev], *zvel_new[lev], zvel_new, BCVars::zvel_bc, BdyVars::null);

    //FillPatch(lev, t_old[lev], *vec_sstore[lev], GetVecOfPtrs(vec_sstore), BCVars::cons_bc, BdyVars::t);
    //FillPatch(lev, time, *vec_sstore[lev], GetVecOfPtrs(vec_sstore), BCVars::cons_bc, BdyVars::t,0,true,true,0,0,0.0,*cons_old[lev]);

    auto N = Geom(lev).Domain().size()[2]-1; // Number of vertical "levs" aka, NZ

    advance_3d(lev, *cons_new[lev], *xvel_new[lev], *yvel_new[lev],
               vec_sstore[lev].get(),
               vec_ru[lev].get(), vec_rv[lev].get(),
               vec_DU_avg1[lev], vec_DU_avg2[lev],
               vec_DV_avg1[lev], vec_DV_avg2[lev],
               vec_ubar[lev],  vec_vbar[lev],
               vec_Akv[lev], vec_Akt[lev], vec_Hz[lev], vec_Huon[lev], vec_Hvom[lev],
               vec_z_w[lev], vec_h[lev].get(),
               vec_pm[lev].get(), vec_pn[lev].get(), vec_mskr[lev].get(),
               vec_msku[lev].get(), vec_mskv[lev].get(),
               N, dt_lev);

    // All three time levels of ubar and vbar: one exchange at level 0, one two-level fill
    // per field above it. Filling all components hands the interpolator the BCRecs of the
    // next bccomps for components 1 and 2, which face_cons_linear_interp ignores.
    if (lev == 0) {
        amrex::FillBoundary(Vector<MultiFab*>{vec_ubar[lev].get(), vec_vbar[lev].get()},
                            geom[lev].periodicity());
    } else {
        FillPatchNoBC(lev, t_old[lev], *vec_ubar[lev], GetVecOfPtrs(vec_ubar), bdy_ubar(),0,true,false);
        FillPatchNoBC(lev, t_old[lev], *vec_vbar[lev], GetVecOfPtrs(vec_vbar), bdy_vbar(),0,true,false);
    }
    //FillPatch(lev, t_old[lev], *vec_sstore[lev], GetVecOfPtrs(vec_sstore), BCVars::cons_bc, BdyVars::t);


    // Fill in three ways: 1) interpolate from coarse grid if lev > 0; 2) fill from physical boundaries;
    //                     3) fine-fine fill of ghost cells with FillBoundary call
    // Note that we need the fine-fine and physical bc's in order to correctly move the particles
    amrex::FillBoundary(Vector<MultiFab*>{xvel_new[lev], yvel_new[lev]}, geom[lev].periodicity());
    // The true vertical velocity is a diagnostic: nothing in the dynamics reads it. A single
    // level leaves it for update_wvel to compute when a reader asks. With refinement, coarse w
    // is averaged down from the fine level and the fine ghost cells are
    // interpolated from the parent's w at both ends of its step, so w is formed every step
    // -- except on the finest level's substeps before the last (wvel_skip_substep). Nothing
    // reads that w: the parent averages down the last substep's, and the next substep's
    // zvel_old, which holds it, is read only by a checkpoint, before which every substep
    // forms w. A regrid or tagging that reads it first calls update_wvel, which forms it from
    // the same state the step would have.
    if (wvel_is_lazy() || wvel_skip_this_step) {
        wvel_current[lev] = 0;
    } else {
        compute_and_fill_wvel(lev);
    }

    // Apply land/sea mask to tracers
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(*cons_new[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<Real> const& cons = cons_new[lev]->array(mfi);
        Array4<Real> const& mskr = vec_mskr[lev]->array(mfi);

        Box bx = mfi.tilebox();

        ParallelFor(bx, ncons, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
        {
            cons(i,j,k,n) *= mskr(i,j,0);
        });
    }

#ifdef REMORA_USE_PARTICLES
    // **************************************************************************************
    // Update the particle positions
    // **************************************************************************************
   evolveTracers(lev, dt_lev);
#endif
}
