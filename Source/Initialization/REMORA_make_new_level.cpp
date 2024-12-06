/**
 * \file REMORA_make_new_level.cpp
 */

#include <REMORA.H>
#include <REMORA_prob_common.H>

#include <AMReX_buildInfo.H>

using namespace amrex;

// Make a new level using provided BoxArray and DistributionMapping and
// fill with interpolated coarse level data (overrides the pure virtual function in AmrCore)
// regrid  --> RemakeLevel            (if level already existed)
// regrid  --> MakeNewLevelFromCoarse (if adding new level)
void
REMORA::MakeNewLevelFromCoarse (int lev, Real time, const BoxArray& ba,
                               const DistributionMapping& dm)
{
    BoxList bl2d = ba.boxList();
    for (auto& b : bl2d) {
        b.setRange(2,0);
    }
    BoxArray ba2d(std::move(bl2d));

    cons_new[lev] = new MultiFab(ba, dm, NCONS, cons_new[lev-1]->nGrowVect());
    cons_old[lev] = new MultiFab(ba, dm, NCONS, cons_new[lev-1]->nGrowVect());

    xvel_new[lev] = new MultiFab(convert(ba, IntVect(1,0,0)), dm, 1, xvel_new[lev-1]->nGrowVect());
    xvel_old[lev] = new MultiFab(convert(ba, IntVect(1,0,0)), dm, 1, xvel_new[lev-1]->nGrowVect());

    yvel_new[lev] = new MultiFab(convert(ba, IntVect(0,1,0)), dm, 1, yvel_new[lev-1]->nGrowVect());
    yvel_old[lev] = new MultiFab(convert(ba, IntVect(0,1,0)), dm, 1, yvel_new[lev-1]->nGrowVect());

    zvel_new[lev] = new MultiFab(convert(ba, IntVect(0,0,1)), dm, 1, zvel_new[lev-1]->nGrowVect());
    zvel_old[lev] = new MultiFab(convert(ba, IntVect(0,0,1)), dm, 1, zvel_new[lev-1]->nGrowVect());

    resize_stuff(lev);

    vec_Zt_avg1[lev].reset(new MultiFab(ba2d ,dm,1,IntVect(NGROW+1,NGROW+1,0))); //2d, average of the free surface (zeta)
    vec_hOfTheConfusingName[lev].reset(new MultiFab(ba2d ,dm,2,IntVect(NGROW+1,NGROW+1,0))); //2d, average of the free surface (zeta)
    vec_ubar[lev].reset(new MultiFab(convert(ba2d,IntVect(1,0,0)),dm,3,IntVect(NGROW,NGROW,0)));
    vec_vbar[lev].reset(new MultiFab(convert(ba2d,IntVect(0,1,0)),dm,3,IntVect(NGROW,NGROW,0)));

    vec_ru[lev].reset(new MultiFab(convert(ba,IntVect(1,0,0)),dm,2,IntVect(NGROW,NGROW,0))); // RHS u (incl horizontal and vertical advection)
    vec_rv[lev].reset(new MultiFab(convert(ba,IntVect(0,1,0)),dm,2,IntVect(NGROW,NGROW,0))); // RHS v

    vec_ru2d[lev].reset(new MultiFab(convert(ba2d,IntVect(1,0,0)),dm,2,IntVect(NGROW,NGROW,0))); // RHS u for 2d
    vec_rv2d[lev].reset(new MultiFab(convert(ba2d,IntVect(0,1,0)),dm,2,IntVect(NGROW,NGROW,0))); // RHS v for 2d

    t_new[lev] = time;
    t_old[lev] = time - 1.e200_rt;

    init_masks(lev, ba, dm);

    FillCoarsePatch(lev, time, cons_new[lev], cons_new[lev-1]);
    FillCoarsePatch(lev, time, xvel_new[lev], xvel_new[lev-1]);
    FillCoarsePatch(lev, time, yvel_new[lev], yvel_new[lev-1]);
    FillCoarsePatch(lev, time, zvel_new[lev], zvel_new[lev-1]);

    init_stuff(lev, ba, dm);

    FillCoarsePatch(lev, time, vec_hOfTheConfusingName[lev].get(), vec_hOfTheConfusingName[lev-1].get());
    FillCoarsePatch(lev, time, vec_Zt_avg1[lev].get(), vec_Zt_avg1[lev-1].get());
    for (int icomp=0; icomp<3; icomp++) {
        FillCoarsePatch(lev, time, vec_ubar[lev].get(), vec_ubar[lev-1].get(),icomp,false);
        FillCoarsePatch(lev, time, vec_vbar[lev].get(), vec_vbar[lev-1].get(),icomp,false);
    }
    for (int icomp=0; icomp<2; icomp++) {
        FillCoarsePatch(lev, time, vec_ru[lev].get(), vec_ru[lev-1].get(),icomp,false);
        FillCoarsePatch(lev, time, vec_rv[lev].get(), vec_rv[lev-1].get(),icomp,false);
        FillCoarsePatch(lev, time, vec_ru2d[lev].get(), vec_ru2d[lev-1].get(),icomp,false);
        FillCoarsePatch(lev, time, vec_rv2d[lev].get(), vec_rv2d[lev-1].get(),icomp,false);
    }


    set_pm_pn(lev);
    stretch_transform(lev);

    init_set_vmix(lev);
    set_hmixcoef(lev);
    set_coriolis(lev);
    set_zeta_to_Ztavg(lev);
    init_custom_smflux(geom[lev], time, *vec_sustr[lev], *vec_svstr[lev], solverChoice);

    // ********************************************************************************************
    // If we are making a new level then the FillPatcher for this level hasn't been allocated yet
    // ********************************************************************************************
    if (cf_width >= 0) {
        Construct_REMORAFillPatchers(lev);
           Define_REMORAFillPatchers(lev);
    }

#ifdef REMORA_USE_PARTICLES
    // particleData.Redistribute();
#endif
}

// Remake an existing level using provided BoxArray and DistributionMapping and
// fill with existing fine and coarse data.
// overrides the pure virtual function in AmrCore
void
REMORA::RemakeLevel (int lev, Real time, const BoxArray& ba, const DistributionMapping& dm)
{
    BoxArray            ba_old(cons_new[lev]->boxArray());
    DistributionMapping dm_old(cons_new[lev]->DistributionMap());

    BoxList bl2d = ba.boxList();
    for (auto& b : bl2d) {
        b.setRange(2,0);
    }
    BoxArray ba2d(std::move(bl2d));

#if (NGROW==2)
    int ngrow_state   = ComputeGhostCells(solverChoice.spatial_order)+1;
    int ngrow_vels    = ComputeGhostCells(solverChoice.spatial_order)+1;
    int ngrow_zeta    = ComputeGhostCells(solverChoice.spatial_order)+1;
    int ngrow_h       = ComputeGhostCells(solverChoice.spatial_order)+1;
    int ngrow_velbar  = ComputeGhostCells(solverChoice.spatial_order);
#else
    int ngrow_state   = ComputeGhostCells(solverChoice.spatial_order)+2;
    int ngrow_vels    = ComputeGhostCells(solverChoice.spatial_order)+2;
    int ngrow_zeta    = ComputeGhostCells(solverChoice.spatial_order)+2;
    int ngrow_h       = ComputeGhostCells(solverChoice.spatial_order)+2;
    int ngrow_velbar  = ComputeGhostCells(solverChoice.spatial_order)+1;
#endif

    MultiFab tmp_cons_new(ba, dm, NCONS, ngrow_state);
    MultiFab tmp_cons_old(ba, dm, NCONS, ngrow_state);

    MultiFab tmp_xvel_new(convert(ba, IntVect(1,0,0)), dm, 1, ngrow_vels);
    MultiFab tmp_xvel_old(convert(ba, IntVect(1,0,0)), dm, 1, ngrow_vels);

    MultiFab tmp_yvel_new(convert(ba, IntVect(0,1,0)), dm, 1, ngrow_vels);
    MultiFab tmp_yvel_old(convert(ba, IntVect(0,1,0)), dm, 1, ngrow_vels);

    MultiFab tmp_zvel_new(convert(ba, IntVect(0,0,1)), dm, 1, IntVect(ngrow_vels,ngrow_vels,0));
    MultiFab tmp_zvel_old(convert(ba, IntVect(0,0,1)), dm, 1, IntVect(ngrow_vels,ngrow_vels,0));

    MultiFab tmp_Zt_avg1_new(ba2d, dm, 1, IntVect(ngrow_zeta,ngrow_zeta,0));
    MultiFab tmp_Zt_avg1_old(ba2d, dm, 1, IntVect(ngrow_zeta,ngrow_zeta,0));
    MultiFab tmp_h(ba2d, dm, 2, IntVect(ngrow_h,ngrow_h,0));

    MultiFab tmp_ubar_new(convert(ba2d, IntVect(1,0,0)), dm, 3, IntVect(ngrow_velbar,ngrow_velbar,0));
    MultiFab tmp_ubar_old(convert(ba2d, IntVect(1,0,0)), dm, 3, IntVect(ngrow_velbar,ngrow_velbar,0));

    MultiFab tmp_vbar_new(convert(ba2d, IntVect(0,1,0)), dm, 3, IntVect(ngrow_velbar,ngrow_velbar,0));
    MultiFab tmp_vbar_old(convert(ba2d, IntVect(0,1,0)), dm, 3, IntVect(ngrow_velbar,ngrow_velbar,0));

    MultiFab tmp_ru_new(convert(ba, IntVect(1,0,0)),dm,2,IntVect(NGROW,NGROW,0));
    MultiFab tmp_rv_new(convert(ba, IntVect(0,1,0)),dm,2,IntVect(NGROW,NGROW,0));

    MultiFab tmp_ru2d_new(convert(ba2d, IntVect(1,0,0)),dm,2,IntVect(NGROW,NGROW,0));
    MultiFab tmp_rv2d_new(convert(ba2d, IntVect(0,1,0)),dm,2,IntVect(NGROW,NGROW,0));

    init_masks(lev, ba, dm);

    // This will fill the temporary MultiFabs with data from previous fine data as well as coarse where needed
    FillPatch(lev, time, tmp_cons_new, cons_new, BCVars::cons_bc, BdyVars::t,0,true,false);
    FillPatch(lev, time, tmp_xvel_new, xvel_new, BCVars::xvel_bc, BdyVars::u,0,true,false,0,0,0.0,tmp_xvel_new);
    FillPatch(lev, time, tmp_yvel_new, yvel_new, BCVars::yvel_bc, BdyVars::v,0,true,false,0,0,0.0,tmp_yvel_new);
    FillPatch(lev, time, tmp_zvel_new, zvel_new, BCVars::zvel_bc, BdyVars::null,0,true,false);

    FillPatch(lev, time, tmp_h, GetVecOfPtrs(vec_hOfTheConfusingName), BCVars::cons_bc, BdyVars::null,0,false,false);
    FillPatch(lev, time, tmp_h, GetVecOfPtrs(vec_hOfTheConfusingName), BCVars::cons_bc, BdyVars::null,1,false,false);
    FillPatch(lev, time, tmp_Zt_avg1_new, GetVecOfPtrs(vec_Zt_avg1), BCVars::zeta_bc, BdyVars::null,0,true,false);
    for (int icomp=0; icomp<3; icomp++) {
        FillPatch(lev, time, tmp_ubar_new, GetVecOfPtrs(vec_ubar), BCVars::ubar_bc, BdyVars::ubar, icomp,false,false);
        FillPatch(lev, time, tmp_vbar_new, GetVecOfPtrs(vec_vbar), BCVars::vbar_bc, BdyVars::vbar, icomp,false,false);
    }
    for (int icomp=0; icomp<2; icomp++) {
        FillPatch(lev, time, tmp_ru_new, GetVecOfPtrs(vec_ru),BCVars::xvel_bc, BdyVars::null, icomp,false,false);
        FillPatch(lev, time, tmp_rv_new, GetVecOfPtrs(vec_rv),BCVars::yvel_bc, BdyVars::null, icomp,false,false);
        // These might want to have BCVars::ubar_bc and vbar_bc
        FillPatch(lev, time, tmp_ru2d_new, GetVecOfPtrs(vec_ru2d),BCVars::xvel_bc, BdyVars::null, icomp,false,false);
        FillPatch(lev, time, tmp_rv2d_new, GetVecOfPtrs(vec_rv2d),BCVars::yvel_bc, BdyVars::null, icomp,false,false);
    }

    MultiFab::Copy(tmp_cons_old,tmp_cons_new,0,0,NCONS,tmp_cons_new.nGrowVect());
    MultiFab::Copy(tmp_xvel_old,tmp_xvel_new,0,0,    1,tmp_xvel_new.nGrowVect());
    MultiFab::Copy(tmp_yvel_old,tmp_yvel_new,0,0,    1,tmp_yvel_new.nGrowVect());
    MultiFab::Copy(tmp_zvel_old,tmp_zvel_new,0,0,    1,tmp_zvel_new.nGrowVect());

    std::swap(tmp_cons_new, *cons_new[lev]);
    std::swap(tmp_cons_old, *cons_old[lev]);
    std::swap(tmp_xvel_new, *xvel_new[lev]);
    std::swap(tmp_xvel_old, *xvel_old[lev]);
    std::swap(tmp_yvel_new, *yvel_new[lev]);
    std::swap(tmp_yvel_old, *yvel_old[lev]);
    std::swap(tmp_zvel_new, *zvel_new[lev]);
    std::swap(tmp_zvel_old, *zvel_old[lev]);
    std::swap(tmp_Zt_avg1_new, *vec_Zt_avg1[lev]);
    std::swap(tmp_h,           *vec_hOfTheConfusingName[lev]);
    std::swap(tmp_ubar_new,    *vec_ubar[lev]);
    std::swap(tmp_vbar_new,    *vec_vbar[lev]);
    std::swap(tmp_ru_new,    *vec_ru[lev]);
    std::swap(tmp_rv_new,    *vec_rv[lev]);
    std::swap(tmp_ru2d_new,    *vec_ru2d[lev]);
    std::swap(tmp_rv2d_new,    *vec_rv2d[lev]);

    t_new[lev] = time;
    t_old[lev] = time - 1.e200_rt;

    init_stuff(lev, ba, dm);

    set_pm_pn(lev);
    stretch_transform(lev);

    init_set_vmix(lev);
    set_hmixcoef(lev);
    set_coriolis(lev);
    set_zeta_to_Ztavg(lev);
    init_custom_smflux(geom[lev], time, *vec_sustr[lev], *vec_svstr[lev], solverChoice);

    // We need to re-define the FillPatcher if the grids have changed
    if (lev > 0 && cf_width >= 0) {
        bool ba_changed = (ba != ba_old);
        bool dm_changed = (dm != dm_old);
        if (ba_changed || dm_changed) {
          Define_REMORAFillPatchers(lev);
        }
    }

#ifdef REMORA_USE_PARTICLES
    particleData.Redistribute();
#endif
}

// Make a new level from scratch using provided BoxArray and DistributionMapping.
// This is called both for initialization and for restart
// (overrides the pure virtual function in AmrCore)
// main.cpp --> REMORA::InitData --> InitFromScratch --> MakeNewGrids --> MakeNewLevelFromScratch
//                                         restart  --> MakeNewGrids --> MakeNewLevelFromScratch
void REMORA::MakeNewLevelFromScratch (int lev, Real time, const BoxArray& ba,
                                     const DistributionMapping& dm)
{
    // Set BoxArray grids and DistributionMapping dmap in AMReX_AmrMesh.H class
    SetBoxArray(lev, ba);
    SetDistributionMap(lev, dm);

    BoxList bl2d = ba.boxList();
    for (auto& b : bl2d) {
        b.setRange(2,0);
    }
    BoxArray ba2d(std::move(bl2d));

    amrex::Print() << "GRIDS AT LEVEL " << lev << " ARE " << ba << std::endl;

    // The number of ghost cells for density must be 1 greater than that for velocity
    //     so that we can go back in forth between velocity and momentum on all faces
#if NGROW==2
    int ngrow_state = ComputeGhostCells(solverChoice.spatial_order)+1;
    int ngrow_vels  = ComputeGhostCells(solverChoice.spatial_order)+1;
#else
    int ngrow_state = ComputeGhostCells(solverChoice.spatial_order)+2;
    int ngrow_vels  = ComputeGhostCells(solverChoice.spatial_order)+2;
#endif

    cons_old[lev] = new MultiFab(ba, dm, NCONS, ngrow_state);
    cons_new[lev] = new MultiFab(ba, dm, NCONS, ngrow_state);

    xvel_new[lev] = new MultiFab(convert(ba, IntVect(1,0,0)), dm, 1, ngrow_vels);
    xvel_old[lev] = new MultiFab(convert(ba, IntVect(1,0,0)), dm, 1, ngrow_vels);

    yvel_new[lev] = new MultiFab(convert(ba, IntVect(0,1,0)), dm, 1, ngrow_vels);
    yvel_old[lev] = new MultiFab(convert(ba, IntVect(0,1,0)), dm, 1, ngrow_vels);

    zvel_new[lev] = new MultiFab(convert(ba, IntVect(0,0,1)), dm, 1, IntVect(ngrow_vels,ngrow_vels,0));
    zvel_old[lev] = new MultiFab(convert(ba, IntVect(0,0,1)), dm, 1, IntVect(ngrow_vels,ngrow_vels,0));

    resize_stuff(lev);

    vec_Zt_avg1[lev].reset(new MultiFab(ba2d ,dm,1,IntVect(NGROW+1,NGROW+1,0))); //2d, average of the free surface (zeta)
    vec_hOfTheConfusingName[lev].reset(new MultiFab(ba2d ,dm,2,IntVect(NGROW+1,NGROW+1,0))); //2d, bathymetry
    vec_ubar[lev].reset(new MultiFab(convert(ba2d,IntVect(1,0,0)),dm,3,IntVect(NGROW,NGROW,0)));
    vec_vbar[lev].reset(new MultiFab(convert(ba2d,IntVect(0,1,0)),dm,3,IntVect(NGROW,NGROW,0)));

    vec_ru[lev].reset(new MultiFab(convert(ba,IntVect(1,0,0)),dm,2,IntVect(NGROW,NGROW,0))); // RHS u (incl horizontal and vertical advection)
    vec_rv[lev].reset(new MultiFab(convert(ba,IntVect(0,1,0)),dm,2,IntVect(NGROW,NGROW,0))); // RHS v

    vec_ru2d[lev].reset(new MultiFab(convert(ba2d,IntVect(1,0,0)),dm,2,IntVect(NGROW,NGROW,0))); // RHS u (incl horizontal and vertical advection)
    vec_rv2d[lev].reset(new MultiFab(convert(ba2d,IntVect(0,1,0)),dm,2,IntVect(NGROW,NGROW,0))); // RHS v

    init_masks(lev, ba, dm);
    init_stuff(lev, ba, dm);

    init_only(lev, time);

#ifdef REMORA_USE_PARTICLES
    if (restart_chkfile.empty()) {
        if (lev == 0) {
            initializeTracers((ParGDBBase*)GetParGDB(),vec_z_phys_nd);
        } else {
            particleData.Redistribute();
        }
    }
#endif
}

void REMORA::resize_stuff(int lev)
{
    vec_z_phys_nd.resize(lev+1);

    vec_hOfTheConfusingName.resize(lev+1);
    vec_Zt_avg1.resize(lev+1);
    vec_s_r.resize(lev+1);
    vec_s_w.resize(lev+1);
    vec_z_w.resize(lev+1);
    vec_z_r.resize(lev+1);
    vec_Hz.resize(lev+1);
    vec_Huon.resize(lev+1);
    vec_Hvom.resize(lev+1);
    vec_Akv.resize(lev+1);
    vec_Akt.resize(lev+1);
    vec_visc2_p.resize(lev+1);
    vec_visc2_r.resize(lev+1);
    vec_diff2.resize(lev+1);
    vec_ru.resize(lev+1);
    vec_rv.resize(lev+1);
    vec_ru2d.resize(lev+1);
    vec_rv2d.resize(lev+1);
    vec_rufrc.resize(lev+1);
    vec_rvfrc.resize(lev+1);
    vec_sustr.resize(lev+1);
    vec_svstr.resize(lev+1);
    vec_rdrag.resize(lev+1);
    vec_bustr.resize(lev+1);
    vec_bvstr.resize(lev+1);

    vec_DU_avg1.resize(lev+1);
    vec_DU_avg2.resize(lev+1);
    vec_DV_avg1.resize(lev+1);
    vec_DV_avg2.resize(lev+1);
    vec_rubar.resize(lev+1);
    vec_rvbar.resize(lev+1);
    vec_rzeta.resize(lev+1);
    vec_ubar.resize(lev+1);
    vec_vbar.resize(lev+1);
    vec_zeta.resize(lev+1);
    vec_mskr.resize(lev+1);
    vec_msku.resize(lev+1);
    vec_mskv.resize(lev+1);
    vec_mskp.resize(lev+1);
    vec_sstore.resize(lev+1);

    vec_pm.resize(lev+1);
    vec_pn.resize(lev+1);
    vec_fcor.resize(lev+1);

    vec_xr.resize(lev+1);
    vec_yr.resize(lev+1);
    vec_xu.resize(lev+1);
    vec_yu.resize(lev+1);
    vec_xv.resize(lev+1);
    vec_yv.resize(lev+1);
    vec_xp.resize(lev+1);
    vec_yp.resize(lev+1);

    vec_rhoS.resize(lev+1);
    vec_rhoA.resize(lev+1);
    vec_bvf.resize(lev+1);

    mapfac_m.resize(lev+1);
    mapfac_u.resize(lev+1);
    mapfac_v.resize(lev+1);

    vec_tke.resize(lev+1);
    vec_gls.resize(lev+1);
    vec_Lscale.resize(lev+1);
    vec_Akk.resize(lev+1);
    vec_Akp.resize(lev+1);
}
void REMORA::init_masks (int lev, const BoxArray& ba, const DistributionMapping& dm)
{
    BoxList bl2d = ba.boxList();
    for (auto& b : bl2d) {
        b.setRange(2,0);
    }

    BoxArray ba2d(std::move(bl2d));
    vec_mskr[lev].reset(new MultiFab(ba2d,dm,1,IntVect(NGROW+1,NGROW+1,0)));
    vec_msku[lev].reset(new MultiFab(convert(ba2d,IntVect(1,0,0)),dm,1,IntVect(NGROW+1,NGROW+1,0)));
    vec_mskv[lev].reset(new MultiFab(convert(ba2d,IntVect(0,1,0)),dm,1,IntVect(NGROW+1,NGROW+1,0)));
    vec_mskp[lev].reset(new MultiFab(convert(ba2d,IntVect(1,1,0)),dm,1,IntVect(NGROW+1,NGROW+1,0)));

    vec_mskr[lev]->setVal(1.0_rt);
    vec_msku[lev]->setVal(1.0_rt);
    vec_mskv[lev]->setVal(1.0_rt);
    vec_mskp[lev]->setVal(1.0_rt);
}

void REMORA::init_stuff (int lev, const BoxArray& ba, const DistributionMapping& dm)
{
    // ********************************************************************************************
    // Initialize the boundary conditions
    // ********************************************************************************************
    physbcs[lev] = std::make_unique<REMORAPhysBCFunct> (lev, geom[lev], domain_bcs_type, domain_bcs_type_d,
                                                       m_bc_extdir_vals);

    BoxList bl2d = ba.boxList();
    for (auto& b : bl2d) {
        b.setRange(2,0);
    }
    BoxArray ba2d(std::move(bl2d));

    BoxList bl1d = ba.boxList();
    for (auto& b : bl1d) {
        b.setRange(0,0);
        b.setRange(1,0);
    }
    BoxArray ba1d(std::move(bl1d));

    // Map factors
    mapfac_m[lev].reset(new MultiFab(ba2d,dm,1,0));
    mapfac_u[lev].reset(new MultiFab(convert(ba2d,IntVect(1,0,0)),dm,1,0));
    mapfac_v[lev].reset(new MultiFab(convert(ba2d,IntVect(0,1,0)),dm,1,0));
    mapfac_m[lev]->setVal(1.);
    mapfac_u[lev]->setVal(1.);
    mapfac_v[lev]->setVal(1.);

    BoxArray ba_nd(ba);
    ba_nd.surroundingNodes();
    BoxArray ba_w(ba);
    ba_w.surroundingNodes(2);

    vec_z_phys_nd[lev].reset          (new MultiFab(ba_nd,dm,1,IntVect(NGROW,NGROW,1))); // z at psi points (nodes) MIGHT NEED NGROW+1

    vec_s_r[lev].reset                (new MultiFab(ba1d,dm,1,IntVect(    0,    0,0))); // scaled vertical coordinate [0,1] , transforms to z

    vec_s_w[lev].reset                (new MultiFab(convert(ba1d,IntVect(0,0,1)),dm,1,IntVect(    0,    0,0))); // scaled vertical coordinate at w-points [0,1] , transforms to z

    vec_z_w[lev].reset                (new MultiFab(convert(ba,IntVect(0,0,1)),dm,1,IntVect(NGROW+1,NGROW+1,0))); // z at w points (cell faces)
    vec_z_r[lev].reset                (new MultiFab(ba,dm,1,IntVect(NGROW+1,NGROW+1,0))); // z at r points (cell center)
    vec_Hz[lev].reset                 (new MultiFab(ba,dm,1,IntVect(NGROW+1,NGROW+1,NGROW+1))); // like in ROMS, thickness of cell in z

    vec_Huon[lev].reset               (new MultiFab(convert(ba,IntVect(1,0,0)),dm,1,IntVect(NGROW,NGROW,0))); // mass flux for u component
    vec_Hvom[lev].reset               (new MultiFab(convert(ba,IntVect(0,1,0)),dm,1,IntVect(NGROW,NGROW,0))); // mass flux for v component

    vec_Akv[lev].reset                (new MultiFab(convert(ba,IntVect(0,0,1)),dm,1,IntVect(NGROW,NGROW,0))); // vertical mixing coefficient (.in)
    vec_Akt[lev].reset                (new MultiFab(convert(ba,IntVect(0,0,1)),dm,NCONS,IntVect(NGROW,NGROW,0))); // vertical mixing coefficient (.in)

    // check dimensionality
    vec_visc2_p[lev].reset(new MultiFab(ba,dm,1,IntVect(NGROW,NGROW,0))); // harmonic viscosity at psi points -- difference to 3d?
    vec_visc2_r[lev].reset(new MultiFab(ba,dm,1,IntVect(NGROW,NGROW,0))); // harmonic viscosity at rho points
    vec_diff2[lev].reset(new MultiFab(ba,dm,NCONS,IntVect(NGROW,NGROW,0))); // harmonic diffusivity temperature/salt

    //2d, (incl advection terms and surface/bottom stresses, integral over the whole column, k=0)
    vec_rufrc[lev].reset(new MultiFab(convert(ba2d,IntVect(1,0,0)),dm,2,IntVect(NGROW,NGROW,0)));
    vec_rvfrc[lev].reset(new MultiFab(convert(ba2d,IntVect(0,1,0)),dm,2,IntVect(NGROW,NGROW,0))); //2d, same as above but v

    vec_sustr[lev].reset(new MultiFab(convert(ba2d,IntVect(1,0,0)),dm,1,IntVect(NGROW,NGROW,0))); //2d, surface stress
    vec_svstr[lev].reset(new MultiFab(convert(ba2d,IntVect(0,1,0)),dm,1,IntVect(NGROW,NGROW,0))); //2d

    //2d, linear drag coefficient [m/s], defined at rho, somehow related to rdrg
    vec_rdrag[lev].reset(new MultiFab(ba2d,dm,1,IntVect(NGROW,NGROW,0)));

    vec_bustr[lev].reset(new MultiFab(convert(ba2d,IntVect(1,0,0)),dm,1,IntVect(NGROW,NGROW,0))); //2d, bottom stress
    vec_bvstr[lev].reset(new MultiFab(convert(ba2d,IntVect(0,1,0)),dm,1,IntVect(NGROW,NGROW,0)));

    //all 2d -- all associated with the 2D advance
    //2d DU: sum(height[incl free surface?] * u)
    vec_DU_avg1[lev].reset(new MultiFab(convert(ba2d,IntVect(1,0,0)),dm,1,IntVect(NGROW,NGROW,0)));

    //2d like above, but correct(or)?
    vec_DU_avg2[lev].reset(new MultiFab(convert(ba2d,IntVect(1,0,0)),dm,1,IntVect(NGROW,NGROW,0)));

    vec_DV_avg1[lev].reset(new MultiFab(convert(ba2d,IntVect(0,1,0)),dm,1,IntVect(NGROW,NGROW,0)));
    vec_DV_avg2[lev].reset(new MultiFab(convert(ba2d,IntVect(0,1,0)),dm,1,IntVect(NGROW,NGROW,0)));

    vec_rubar[lev].reset(new MultiFab(convert(ba2d,IntVect(1,0,0)),dm,4,IntVect(NGROW,NGROW,0))); // 2d RHS ubar
    vec_rvbar[lev].reset(new MultiFab(convert(ba2d,IntVect(0,1,0)),dm,4,IntVect(NGROW,NGROW,0)));
    vec_rzeta[lev].reset(new MultiFab(ba2d,dm,4,IntVect(NGROW,NGROW,0))); // 2d RHS zeta

    // starts off kind of like a depth-averaged u, but exists at more points and more timesteps (b/c fast 2D update) than full u
    vec_zeta[lev].reset(new MultiFab(ba2d,dm,3,IntVect(NGROW+1,NGROW+1,0)));  // 2d free surface

    vec_pm[lev].reset(new MultiFab(ba2d,dm,1,IntVect(NGROW+1,NGROW+2,0)));
    vec_pn[lev].reset(new MultiFab(ba2d,dm,1,IntVect(NGROW+2,NGROW+1,0)));
    vec_fcor[lev].reset(new MultiFab(ba2d,dm,1,IntVect(NGROW+1,NGROW+1,0)));

    vec_xr[lev].reset(new MultiFab(ba2d,dm,1,IntVect(NGROW+1,NGROW+1,0)));
    vec_yr[lev].reset(new MultiFab(ba2d,dm,1,IntVect(NGROW+1,NGROW+1,0)));

    vec_xu[lev].reset(new MultiFab(convert(ba2d,IntVect(1,0,0)),dm,1,IntVect(NGROW,NGROW,0)));
    vec_yu[lev].reset(new MultiFab(convert(ba2d,IntVect(1,0,0)),dm,1,IntVect(NGROW,NGROW,0)));
    vec_xv[lev].reset(new MultiFab(convert(ba2d,IntVect(0,1,0)),dm,1,IntVect(NGROW,NGROW,0)));
    vec_yv[lev].reset(new MultiFab(convert(ba2d,IntVect(0,1,0)),dm,1,IntVect(NGROW,NGROW,0)));
    vec_xp[lev].reset(new MultiFab(convert(ba2d,IntVect(1,1,0)),dm,1,IntVect(NGROW,NGROW,0)));
    vec_yp[lev].reset(new MultiFab(convert(ba2d,IntVect(1,1,0)),dm,1,IntVect(NGROW,NGROW,0)));


    // tempstore, saltstore, etc
    vec_sstore[lev].reset(new MultiFab(ba,dm,NCONS,IntVect(NGROW,NGROW,0)));

    vec_rhoS[lev].reset(new MultiFab(ba,dm,1,IntVect(NGROW,NGROW,0)));
    vec_rhoA[lev].reset(new MultiFab(ba,dm,1,IntVect(NGROW,NGROW,0)));
    vec_bvf[lev].reset(new MultiFab(convert(ba,IntVect(0,0,1)),dm,1,IntVect(NGROW,NGROW,0)));

    vec_tke[lev].reset(new MultiFab(convert(ba,IntVect(0,0,1)),dm,3,IntVect(NGROW,NGROW,0)));
    vec_gls[lev].reset(new MultiFab(convert(ba,IntVect(0,0,1)),dm,3,IntVect(NGROW,NGROW,0)));
    vec_Lscale[lev].reset(new MultiFab(convert(ba,IntVect(0,0,1)),dm,1,IntVect(NGROW,NGROW,0)));
    vec_Akk[lev].reset(new MultiFab(convert(ba,IntVect(0,0,1)),dm,1,IntVect(NGROW,NGROW,0)));
    vec_Akp[lev].reset(new MultiFab(convert(ba,IntVect(0,0,1)),dm,1,IntVect(NGROW,NGROW,0)));

    set_weights(lev);

    vec_DU_avg1[lev]->setVal(0.0_rt);
    vec_DU_avg2[lev]->setVal(0.0_rt);
    vec_DV_avg1[lev]->setVal(0.0_rt);
    vec_DV_avg2[lev]->setVal(0.0_rt);
    vec_rubar[lev]->setVal(0.0_rt);
    vec_rvbar[lev]->setVal(0.0_rt);
    vec_rzeta[lev]->setVal(0.0_rt);

    // Initialize these vars even if we aren't using GLS to
    // avoid issues on e.g. checkpoint
    vec_tke[lev]->setVal(solverChoice.gls_Kmin);
    vec_gls[lev]->setVal(solverChoice.gls_Pmin);
    vec_Lscale[lev]->setVal(0.0_rt);
    vec_Akk[lev]->setVal(solverChoice.Akk_bak);
    vec_Akp[lev]->setVal(solverChoice.Akp_bak);

    // NOTE: Used to set vec_pm and vec_pn to 1e34 here to make foextrap work
    // when init_type = real. However, this does not appear to be necessary so removing

    // Set initial linear drag coefficient
    vec_rdrag[lev]->setVal(solverChoice.rdrag);

    // ********************************************************************************************
    // Create the REMORAFillPatcher object
    // ********************************************************************************************
    if (lev > 0 && cf_width >= 0) {
        Construct_REMORAFillPatchers(lev);
           Define_REMORAFillPatchers(lev);
    }
}

// Delete level data
// overrides the pure virtual function in AmrCore
void
REMORA::ClearLevel (int lev)
{
    delete cons_new[lev]; delete xvel_new[lev];  delete yvel_new[lev];  delete zvel_new[lev];
    delete cons_old[lev]; delete xvel_old[lev];  delete yvel_old[lev];  delete zvel_old[lev];
}

void
REMORA::set_pm_pn (int lev)
{
    AMREX_ASSERT(solverChoice.ic_bc_type == IC_BC_Type::Custom);
    const auto dxi = Geom(lev).InvCellSize();
    vec_pm[lev]->setVal(dxi[0]); vec_pm[lev]->FillBoundary(geom[lev].periodicity());
    vec_pn[lev]->setVal(dxi[1]); vec_pn[lev]->FillBoundary(geom[lev].periodicity());

    for ( MFIter mfi(*vec_xr[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<Real> const& xr = vec_xr[lev]->array(mfi);
        Array4<Real> const& yr = vec_yr[lev]->array(mfi);
        Array4<Real> const& xu = vec_xu[lev]->array(mfi);
        Array4<Real> const& yu = vec_yu[lev]->array(mfi);
        Array4<Real> const& xv = vec_xv[lev]->array(mfi);
        Array4<Real> const& yv = vec_yv[lev]->array(mfi);
        Array4<Real> const& xp = vec_xp[lev]->array(mfi);
        Array4<Real> const& yp = vec_yp[lev]->array(mfi);

        Box bx = mfi.growntilebox(IntVect(NGROW-1,NGROW-1,0));
        Real dx = 1.0_rt / dxi[0];
        Real dy = 1.0_rt / dxi[1];

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            xr(i,j,0) = (i + 0.5_rt) * dx;
            yr(i,j,0) = (j + 0.5_rt) * dy;
        });

        ParallelFor(convert(bx,IntVect(1,0,0)), [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            xu(i,j,0) = i * dx;
            yu(i,j,0) = (j + 0.5_rt) * dy;
        });

        ParallelFor(convert(bx,IntVect(0,1,0)), [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            xv(i,j,0) = (i + 0.5_rt) * dx;
            yv(i,j,0) = j * dy;
        });

        ParallelFor(convert(bx,IntVect(1,1,0)), [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            xp(i,j,0) = i * dx;
            yp(i,j,0) = j * dy;
        });
    }
}

void
REMORA::set_zeta_to_Ztavg (int lev)
{
    std::unique_ptr<MultiFab>& mf_zeta = vec_zeta[lev];
    std::unique_ptr<MultiFab>& mf_Zt_avg1  = vec_Zt_avg1[lev];
    for ( MFIter mfi(*vec_zeta[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<const Real> const& Zt_avg1 = (mf_Zt_avg1)->const_array(mfi);
        Array4<Real> const& zeta = mf_zeta->array(mfi);

        Box  bx3 = mfi.tilebox(); bx3.grow(IntVect(NGROW+1,NGROW+1,0));

        ParallelFor(bx3, 3, [=] AMREX_GPU_DEVICE (int i, int j, int , int n)
        {
            zeta(i,j,0,n) = Zt_avg1(i,j,0);
        });

    }
}

void
REMORA::update_mskp (int lev)
{
    for ( MFIter mfi(*vec_mskr[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<const Real> const& mskr = vec_mskr[lev]->const_array(mfi);
        Array4<      Real> const& mskp = vec_mskp[lev]->array(mfi);

        Box bx = mfi.tilebox(); bx.grow(IntVect(1,1,0)); bx.makeSlab(2,0);

        Real cff1 = 1.0_rt;
        Real cff2 = 2.0_rt;

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            if ((mskr(i-1,j,0) > 0.5) and (mskr(i,j,0) > 0.5) and (mskr(i-1,j-1,0) > 0.5) and (mskr(i,j-1,0) > 0.5)) {
                mskp(i,j,0) = 1.0_rt;
            } else if ((mskr(i-1,j,0) < 0.5) and (mskr(i,j,0) > 0.5) and (mskr(i-1,j-1,0) > 0.5) and (mskr(i,j-1,0) > 0.5)) {
                mskp(i,j,0) = cff1;
            } else if ((mskr(i-1,j,0) > 0.5) and (mskr(i,j,0) < 0.5) and (mskr(i-1,j-1,0) > 0.5) and (mskr(i,j-1,0) > 0.5)) {
                mskp(i,j,0) = cff1;
            } else if ((mskr(i-1,j,0) > 0.5) and (mskr(i,j,0) > 0.5) and (mskr(i-1,j-1,0) < 0.5) and (mskr(i,j-1,0) > 0.5)) {
                mskp(i,j,0) = cff1;
            } else if ((mskr(i-1,j,0) > 0.5) and (mskr(i,j,0) > 0.5) and (mskr(i-1,j-1,0) > 0.5) and (mskr(i,j-1,0) < 0.5)) {
                mskp(i,j,0) = cff1;
            } else if ((mskr(i-1,j,0) > 0.5) and (mskr(i,j,0) < 0.5) and (mskr(i-1,j-1,0) > 0.5) and (mskr(i,j-1,0) < 0.5)) {
                mskp(i,j,0) = cff2;
            } else if ((mskr(i-1,j,0) < 0.5) and (mskr(i,j,0) > 0.5) and (mskr(i-1,j-1,0) < 0.5) and (mskr(i,j-1,0) > 0.5)) {
                mskp(i,j,0) = cff2;
            } else if ((mskr(i-1,j,0) > 0.5) and (mskr(i,j,0) > 0.5) and (mskr(i-1,j-1,0) < 0.5) and (mskr(i,j-1,0) < 0.5)) {
                mskp(i,j,0) = cff2;
            } else if ((mskr(i-1,j,0) < 0.5) and (mskr(i,j,0) < 0.5) and (mskr(i-1,j-1,0) > 0.5) and (mskr(i,j-1,0) > 0.5)) {
                mskp(i,j,0) = cff2;
            } else {
                mskp(i,j,0) = 0.0_rt;
            }

        });
    }
}
