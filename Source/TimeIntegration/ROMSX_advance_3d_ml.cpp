#include <ROMSX.H>

using namespace amrex;

// advance a single 3D level for a single time step
void ROMSX::advance_3d_ml (int lev, Real dt_lev)
{
    MultiFab& S_old = vars_old[lev][Vars::cons];
    MultiFab& S_new = vars_new[lev][Vars::cons];

    MultiFab& U_old = vars_old[lev][Vars::xvel];
    MultiFab& V_old = vars_old[lev][Vars::yvel];
    MultiFab& W_old = vars_old[lev][Vars::zvel];

    MultiFab& U_new = vars_new[lev][Vars::xvel];
    MultiFab& V_new = vars_new[lev][Vars::yvel];
    MultiFab& W_new = vars_new[lev][Vars::zvel];

    MultiFab mf_u(U_new, amrex::make_alias, 0, 1);
    MultiFab mf_v(V_new, amrex::make_alias, 0, 1);

    MultiFab mf_temp(S_new, amrex::make_alias, Temp_comp, 1);
#ifdef ROMSX_USE_SALINITY
    MultiFab mf_salt(S_new, amrex::make_alias, Salt_comp, 1);
#else
    MultiFab mf_salt(S_new, amrex::make_alias, Temp_comp, 1);
#endif
    MultiFab mf_tempold(S_old, amrex::make_alias, Temp_comp, 1);
#ifdef ROMSX_USE_SALINITY
    MultiFab mf_saltold(S_old, amrex::make_alias, Salt_comp, 1);
#else
    MultiFab mf_saltold(S_old, amrex::make_alias, Temp_comp, 1);
#endif

    const BoxArray&            ba = S_old.boxArray();
    const DistributionMapping& dm = S_old.DistributionMap();

    //Only used locally, probably should be rearranged into FArrayBox declaration
    MultiFab mf_AK(ba,dm,1,IntVect(NGROW,NGROW,0)); //2d missing j coordinate
    MultiFab mf_DC(ba,dm,1,IntVect(NGROW,NGROW,NGROW-1)); //2d missing j coordinate
    MultiFab mf_Hzk(ba,dm,1,IntVect(NGROW,NGROW,NGROW-1)); //2d missing j coordinate

    auto N = Geom(lev).Domain().size()[2]-1; // Number of vertical "levs" aka, NZ

    const int ncomp = 1;

    advance_3d(lev, mf_u, mf_v, mf_tempold, mf_saltold, mf_temp, mf_salt, vec_t3[lev], vec_s3[lev], vec_ru[lev], vec_rv[lev],
               vec_DU_avg1[lev], vec_DU_avg2[lev],
               vec_DV_avg1[lev], vec_DV_avg2[lev],
               vec_ubar[lev],  vec_vbar[lev],
               mf_AK, mf_DC,
               mf_Hzk, vec_Akv[lev], vec_Akt[lev], vec_Hz[lev], vec_Huon[lev], vec_Hvom[lev], vec_z_w[lev], vec_hOfTheConfusingName[lev], ncomp, N, dt_lev);

    U_new.FillBoundary(geom[lev].periodicity());
    V_new.FillBoundary(geom[lev].periodicity());

    U_old.FillBoundary(geom[lev].periodicity());
    V_old.FillBoundary(geom[lev].periodicity());

    mf_temp.FillBoundary(geom[lev].periodicity());
    mf_salt.FillBoundary(geom[lev].periodicity());

    mf_tempold.FillBoundary(geom[lev].periodicity());
    mf_saltold.FillBoundary(geom[lev].periodicity());

    vec_t3[lev]->FillBoundary(geom[lev].periodicity());
    vec_s3[lev]->FillBoundary(geom[lev].periodicity());

    // Not sure why this FillPatch is here??
    FillPatch(lev, t_new[lev], vars_new[lev]);

#ifdef ROMSX_USE_PARTICLES
    particleData.advance_particles(lev, dt_lev, vars_new, vec_z_phys_nd);
#endif

}
