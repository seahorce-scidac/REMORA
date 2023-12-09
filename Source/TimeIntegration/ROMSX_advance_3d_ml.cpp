#include <ROMSX.H>

using namespace amrex;

// advance a single 3D level for a single time step
void ROMSX::advance_3d_ml (int lev, Real dt_lev)
{
    MultiFab& S_old = *cons_old[lev];
    MultiFab& S_new = *cons_new[lev];

    MultiFab& U_old = *xvel_old[lev];
    MultiFab& V_old = *yvel_old[lev];

    MultiFab& U_new = *xvel_new[lev];
    MultiFab& V_new = *yvel_new[lev];

    MultiFab mf_u(U_new, amrex::make_alias, 0, 1);
    MultiFab mf_v(V_new, amrex::make_alias, 0, 1);

    MultiFab mf_temp(S_new, amrex::make_alias, Temp_comp, 1);
#ifdef ROMSX_USE_SALINITY
    MultiFab mf_salt(S_new, amrex::make_alias, Salt_comp, 1);
#else
    MultiFab mf_salt(S_new, amrex::make_alias, Temp_comp, 1);
#endif

    const BoxArray&            ba = S_old.boxArray();
    const DistributionMapping& dm = S_old.DistributionMap();

    //Only used locally, probably should be rearranged into FArrayBox declaration
    MultiFab mf_AK(ba,dm,1,IntVect(NGROW,NGROW,0)); //2d missing j coordinate
    MultiFab mf_DC(ba,dm,1,IntVect(NGROW,NGROW,NGROW-1)); //2d missing j coordinate
    MultiFab mf_Hzk(ba,dm,1,IntVect(NGROW,NGROW,NGROW-1)); //2d missing j coordinate

    auto N = Geom(lev).Domain().size()[2]-1; // Number of vertical "levs" aka, NZ

    const int ncomp = 1;

    advance_3d(lev, mf_u, mf_v, mf_temp, mf_salt, vec_t3[lev], vec_s3[lev], vec_ru[lev], vec_rv[lev],
               vec_DU_avg1[lev], vec_DU_avg2[lev],
               vec_DV_avg1[lev], vec_DV_avg2[lev],
               vec_ubar[lev],  vec_vbar[lev],
               mf_AK, mf_DC,
               mf_Hzk, vec_Akv[lev], vec_Akt[lev], vec_Hz[lev], vec_Huon[lev], vec_Hvom[lev],
               vec_z_w[lev], vec_hOfTheConfusingName[lev], ncomp, N, dt_lev);

    U_new.FillBoundary(geom[lev].periodicity());
    V_new.FillBoundary(geom[lev].periodicity());

    U_old.FillBoundary(geom[lev].periodicity());
    V_old.FillBoundary(geom[lev].periodicity());

    vec_ubar[lev]->FillBoundary(geom[lev].periodicity());
    vec_vbar[lev]->FillBoundary(geom[lev].periodicity());

    mf_temp.FillBoundary(geom[lev].periodicity());
    mf_salt.FillBoundary(geom[lev].periodicity());

    vec_t3[lev]->FillBoundary(geom[lev].periodicity());
    vec_s3[lev]->FillBoundary(geom[lev].periodicity());

    // Fill in three ways: 1) interpolate from coarse grid if lev > 0; 2) fill from physical boundaries;
    //                     3) fine-fine fill of ghost cells with FillBoundary call
    FillPatch(lev, t_new[lev], cons_new[lev], cons_new);
    FillPatch(lev, t_new[lev], xvel_new[lev], xvel_new);
    FillPatch(lev, t_new[lev], yvel_new[lev], yvel_new);
    FillPatch(lev, t_new[lev], zvel_new[lev], zvel_new);

#ifdef ROMSX_USE_PARTICLES
    particleData.advance_particles(lev, dt_lev, {cons_new, xvel_new, yvel_new, zvel_new}, vec_z_phys_nd);
#endif

}
