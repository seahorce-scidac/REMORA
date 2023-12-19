#include <ROMSX.H>

using namespace amrex;

// Advance a single 3D level for a single time step
void ROMSX::advance_3d_ml (int lev, Real dt_lev)
{
    // Fill in three ways: 1) interpolate from coarse grid if lev > 0; 2) fill from physical boundaries;
    //                     3) fine-fine fill of ghost cells with FillBoundary call
    FillPatch(lev, t_new[lev], *cons_new[lev], cons_new);
    FillPatch(lev, t_new[lev], *xvel_new[lev], xvel_new);
    FillPatch(lev, t_new[lev], *yvel_new[lev], yvel_new);
    FillPatch(lev, t_new[lev], *zvel_new[lev], zvel_new);

    FillPatch(lev, t_new[lev], *vec_t3[lev], GetVecOfPtrs(vec_t3));
    FillPatch(lev, t_new[lev], *vec_s3[lev], GetVecOfPtrs(vec_s3));

    auto N = Geom(lev).Domain().size()[2]-1; // Number of vertical "levs" aka, NZ

    advance_3d(lev, *cons_new[lev], *xvel_new[lev], *yvel_new[lev],
               vec_sstore[lev].get(),
               vec_ru[lev].get(), vec_rv[lev].get(),
               vec_DU_avg1[lev], vec_DU_avg2[lev],
               vec_DV_avg1[lev], vec_DV_avg2[lev],
               vec_ubar[lev],  vec_vbar[lev],
               vec_Akv[lev], vec_Akt[lev], vec_Hz[lev], vec_Huon[lev], vec_Hvom[lev],
               vec_z_w[lev], vec_hOfTheConfusingName[lev], N, dt_lev);

    FillPatch(lev, t_new[lev], *vec_ubar[lev], GetVecOfPtrs(vec_ubar));
    FillPatch(lev, t_new[lev], *vec_vbar[lev], GetVecOfPtrs(vec_vbar));
    FillPatch(lev, t_new[lev], *vec_sstore[lev], GetVecOfPtrs(vec_sstore));

    // Fill in three ways: 1) interpolate from coarse grid if lev > 0; 2) fill from physical boundaries;
    //                     3) fine-fine fill of ghost cells with FillBoundary call
    // Note that we need the fine-fine and physical bc's in order to correctly move the particles
    FillPatch(lev, t_new[lev], *cons_new[lev], cons_new);
    FillPatch(lev, t_new[lev], *xvel_new[lev], xvel_new);
    FillPatch(lev, t_new[lev], *yvel_new[lev], yvel_new);
    FillPatch(lev, t_new[lev], *zvel_new[lev], zvel_new);

#ifdef ROMSX_USE_PARTICLES
    //***************************************************
    //Advance particles
    //***************************************************
    particleData.advance_particles(lev, dt_lev, xvel_new[lev], yvel_new[lev], zvel_new[lev], vec_z_phys_nd);
#endif
}
