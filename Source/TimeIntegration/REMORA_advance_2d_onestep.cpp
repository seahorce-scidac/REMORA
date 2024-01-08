#include <REMORA.H>

using namespace amrex;

// advance a single level one 2D step
void REMORA::advance_2d_onestep (int lev, Real /*dt_lev*/, Real dtfast_lev, int my_iif, int nfast_counter)
{
    bool first_2d_step=(my_iif==0);

    // These are needed to pass ctests
    Real dummy_time = 0.0;
    FillPatch(lev, dummy_time, *vec_ubar[lev], GetVecOfPtrs(vec_ubar), BdyVars::ubar);
    FillPatch(lev, dummy_time, *vec_vbar[lev], GetVecOfPtrs(vec_vbar), BdyVars::vbar);

    //Predictor
    bool predictor_2d_step=true;
    int next_indx1 = 0;
    advance_2d(lev,
               vec_rhoS[lev].get() , vec_rhoA[lev].get(),
               vec_ru[lev].get()   , vec_rv[lev].get(),
               vec_rufrc[lev].get(), vec_rvfrc[lev].get(),
               vec_Zt_avg1[lev].get(),
               vec_DU_avg1[lev], vec_DU_avg2[lev],
               vec_DV_avg1[lev], vec_DV_avg2[lev],
               vec_rubar[lev], vec_rvbar[lev], vec_rzeta[lev],
                vec_ubar[lev],  vec_vbar[lev],  vec_zeta[lev].get(),
               vec_hOfTheConfusingName[lev].get(),
               vec_pm[lev].get(), vec_pn[lev].get(),
               vec_fcor[lev].get(),
               vec_visc2_p[lev].get(), vec_visc2_r[lev].get(),
               dtfast_lev, predictor_2d_step, first_2d_step, my_iif, next_indx1);

    //Corrector. Skip it on last fast step
    predictor_2d_step=false;
    if (my_iif < nfast_counter - 1) {

        // These are needed to pass ctests
        FillPatch(lev, dummy_time, *vec_ubar[lev], GetVecOfPtrs(vec_ubar), BdyVars::ubar);
        FillPatch(lev, dummy_time, *vec_vbar[lev], GetVecOfPtrs(vec_vbar), BdyVars::vbar);

        advance_2d(lev,
                   vec_rhoS[lev].get(), vec_rhoA[lev].get(),
                   vec_ru[lev].get(), vec_rv[lev].get(),
                   vec_rufrc[lev].get(), vec_rvfrc[lev].get(),
                   vec_Zt_avg1[lev].get(),
                   vec_DU_avg1[lev], vec_DU_avg2[lev],
                   vec_DV_avg1[lev], vec_DV_avg2[lev],
                   vec_rubar[lev], vec_rvbar[lev], vec_rzeta[lev],
                    vec_ubar[lev],  vec_vbar[lev],  vec_zeta[lev].get(),
                   vec_hOfTheConfusingName[lev].get(),
                   vec_pm[lev].get(), vec_pn[lev].get(),
                    vec_fcor[lev].get(),
                   vec_visc2_p[lev].get(), vec_visc2_r[lev].get(),
                   dtfast_lev, predictor_2d_step, first_2d_step, my_iif, next_indx1);
    }
}
