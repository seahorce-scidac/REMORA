#include <REMORA.H>

using namespace amrex;

/**
 * Advance a barotropic step, as in ROMS main2d.F and main3d.F.
 *
 * @param[in] lev            level of refinement
 * @param[in] dt_lev         model time step at level
 * @param[in] dtfast_lev     barotropic time step at level
 * @param[in] my_iif         current barotropic time step
 * @param[in] nfast_counter  number of barotropic calls in this model step
 */
void REMORA::advance_2d_onestep (int lev, Real /*dt_lev*/, Real dtfast_lev, int my_iif, int nfast_counter)
{
    bool first_2d_step=solverChoice.solve3d ? my_iif==0 : istep[lev]==0;

    //Predictor
    bool predictor_2d_step=true;
    int next_indx1 = 0;
    advance_2d(lev,
               vec_rhoS[lev].get() , vec_rhoA[lev].get(),
               vec_ru2d[lev].get()   , vec_rv2d[lev].get(),
               vec_rufrc[lev].get(), vec_rvfrc[lev].get(),
               vec_Zt_avg1[lev].get(),
               vec_DU_avg1[lev], vec_DU_avg2[lev],
               vec_DV_avg1[lev], vec_DV_avg2[lev],
               vec_rubar[lev], vec_rvbar[lev], vec_rzeta[lev],
                vec_ubar[lev],  vec_vbar[lev],  vec_zeta[lev].get(),
               vec_h[lev].get(),
               vec_pm[lev].get(), vec_pn[lev].get(),
               vec_fcor[lev].get(),
               vec_visc2_p[lev].get(), vec_visc2_r[lev].get(),
               vec_mskr[lev].get(), vec_msku[lev].get(), vec_mskv[lev].get(),
               vec_mskp[lev].get(),
               dtfast_lev, predictor_2d_step, first_2d_step, my_iif, next_indx1);

    //Corrector. Skip it on last fast step
    predictor_2d_step=false;
    if (!solverChoice.solve3d || my_iif < nfast_counter - 1) {

        // These are needed to pass ctests

        advance_2d(lev,
                   vec_rhoS[lev].get(), vec_rhoA[lev].get(),
                   vec_ru2d[lev].get(), vec_rv2d[lev].get(),
                   vec_rufrc[lev].get(), vec_rvfrc[lev].get(),
                   vec_Zt_avg1[lev].get(),
                   vec_DU_avg1[lev], vec_DU_avg2[lev],
                   vec_DV_avg1[lev], vec_DV_avg2[lev],
                   vec_rubar[lev], vec_rvbar[lev], vec_rzeta[lev],
                    vec_ubar[lev],  vec_vbar[lev],  vec_zeta[lev].get(),
                   vec_h[lev].get(),
                   vec_pm[lev].get(), vec_pn[lev].get(),
                    vec_fcor[lev].get(),
                   vec_visc2_p[lev].get(), vec_visc2_r[lev].get(),
                   vec_mskr[lev].get(), vec_msku[lev].get(), vec_mskv[lev].get(),
                   vec_mskp[lev].get(),
                   dtfast_lev, predictor_2d_step, first_2d_step, my_iif, next_indx1);
    }

    if (!solverChoice.solve3d) {
        // Keep the current state in component 0 and the previous state in component 1.
        // The predictor uses component 2; the corrector has just written component 1.
        for (MultiFab* mf : {vec_ubar[lev].get(), vec_vbar[lev].get(), vec_zeta[lev].get(),
                             vec_rubar[lev].get(), vec_rvbar[lev].get(), vec_rzeta[lev].get()}) {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(*mf, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.growntilebox();
                auto const& a = mf->array(mfi);
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    const Real previous = a(i,j,k,0);
                    a(i,j,k,0) = a(i,j,k,1);
                    a(i,j,k,1) = previous;
                });
            }
        }
        MultiFab::Copy(*vec_Zt_avg1[lev], *vec_zeta[lev], 0, 0, 1,
                       vec_Zt_avg1[lev]->nGrowVect());
        stretch_transform(lev);

        // The existing output fields show the depth-uniform horizontal velocity.
        for (int dir = 0; dir < 2; ++dir) {
            MultiFab& vel = dir == 0 ? *xvel_new[lev] : *yvel_new[lev];
            const MultiFab& bar = dir == 0 ? *vec_ubar[lev] : *vec_vbar[lev];
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(vel, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                auto const& a = vel.array(mfi);
                auto const& b = bar.const_array(mfi);
                ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    a(i,j,k) = b(i,j,0,0);
                });
            }
            FillPatch(lev, t_new[lev], vel, dir == 0 ? xvel_new : yvel_new,
                      foextrap_periodic_bc());
        }
        zvel_new[lev]->setVal(Real(0.0));
    }
}
