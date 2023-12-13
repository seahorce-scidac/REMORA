#include <ROMSX.H>

using namespace amrex;

// Advance a single 3D level for a single time step
void ROMSX::advance_3d_ml (int lev, Real dt_lev)
{
    // Fill in three ways: 1) interpolate from coarse grid if lev > 0; 2) fill from physical boundaries;
    //                     3) fine-fine fill of ghost cells with FillBoundary call
    FillPatch(lev, t_old[lev], cons_old[lev], cons_old);
    FillPatch(lev, t_old[lev], xvel_old[lev], xvel_old);
    FillPatch(lev, t_old[lev], yvel_old[lev], yvel_old);
    FillPatch(lev, t_old[lev], zvel_old[lev], zvel_old);

    MultiFab mf_temp(*cons_new[lev], amrex::make_alias, Temp_comp, 1);
    MultiFab mf_salt(*cons_new[lev], amrex::make_alias, Salt_comp, 1);

    auto N = Geom(lev).Domain().size()[2]-1; // Number of vertical "levs" aka, NZ

    advance_3d(lev, *xvel_new[lev], *yvel_new[lev],
               mf_temp, mf_salt,
               vec_t3[lev].get(), vec_s3[lev].get(),
               vec_ru[lev].get(), vec_rv[lev].get(),
               vec_DU_avg1[lev], vec_DU_avg2[lev],
               vec_DV_avg1[lev], vec_DV_avg2[lev],
               vec_ubar[lev],  vec_vbar[lev],
               vec_Akv[lev], vec_Akt[lev], vec_Hz[lev], vec_Huon[lev], vec_Hvom[lev],
               vec_z_w[lev], vec_hOfTheConfusingName[lev], N, dt_lev);

    vec_ubar[lev]->FillBoundary(geom[lev].periodicity());
    vec_vbar[lev]->FillBoundary(geom[lev].periodicity());

    vec_t3[lev]->FillBoundary(geom[lev].periodicity());
    vec_s3[lev]->FillBoundary(geom[lev].periodicity());

    // Fill in three ways: 1) interpolate from coarse grid if lev > 0; 2) fill from physical boundaries;
    //                     3) fine-fine fill of ghost cells with FillBoundary call
    // Note that we need the fine-fine and physical bc's in order to correctly move the particles
    FillPatch(lev, t_new[lev], cons_new[lev], cons_new);
    FillPatch(lev, t_new[lev], xvel_new[lev], xvel_new);
    FillPatch(lev, t_new[lev], yvel_new[lev], yvel_new);
    FillPatch(lev, t_new[lev], zvel_new[lev], zvel_new);

#ifdef ROMSX_USE_PARTICLES
    //***************************************************
    //Advance particles
    //***************************************************
    particleData.advance_particles(lev, dt_lev, xvel_new[lev], yvel_new[lev], zvel_new[lev], vec_z_phys_nd);
#endif
}
