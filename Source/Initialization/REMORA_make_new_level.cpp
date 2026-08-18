/**
 * \file REMORA_make_new_level.cpp
 */

#include <REMORA.H>
#include <REMORA_prob_common.H>

#include <AMReX_buildInfo.H>

using namespace amrex;

/**
 * Make a new level using provided BoxArray and DistributionMapping and
 * fill with interpolated coarse level data (overrides the pure virtual function in AmrCore)
 * regrid  --> RemakeLevel            (if level already existed)
 * regrid  --> MakeNewLevelFromCoarse (if adding new level)
 *
 * @param[in   ] lev     level to make
 * @param[in   ] time    current time
 * @param[in   ] ba      BoxArray for the level
 * @param[in   ] dm      DistributionMapping for the level
 */
void
REMORA::MakeNewLevelFromCoarse (int lev, Real time, const BoxArray& ba,
                               const DistributionMapping& dm)
{
    BoxList bl2d = ba.boxList();
    for (auto& b : bl2d) {
        b.setRange(2,0);
    }
    BoxArray ba2d(std::move(bl2d));

    amrex::Print() << "Making level " << lev << " from coarse" << std::endl;
    amrex::Print() << "GRIDS AT LEVEL " << lev << " ARE " << ba << std::endl;

    cons_new[lev] = new MultiFab(ba, dm, ncons, cons_new[lev-1]->nGrowVect());
    cons_old[lev] = new MultiFab(ba, dm, ncons, cons_new[lev-1]->nGrowVect());

    xvel_new[lev] = new MultiFab(convert(ba, IntVect(1,0,0)), dm, 1, xvel_new[lev-1]->nGrowVect());
    xvel_old[lev] = new MultiFab(convert(ba, IntVect(1,0,0)), dm, 1, xvel_new[lev-1]->nGrowVect());

    yvel_new[lev] = new MultiFab(convert(ba, IntVect(0,1,0)), dm, 1, yvel_new[lev-1]->nGrowVect());
    yvel_old[lev] = new MultiFab(convert(ba, IntVect(0,1,0)), dm, 1, yvel_new[lev-1]->nGrowVect());

    zvel_new[lev] = new MultiFab(convert(ba, IntVect(0,0,1)), dm, 1, zvel_new[lev-1]->nGrowVect());
    zvel_old[lev] = new MultiFab(convert(ba, IntVect(0,0,1)), dm, 1, zvel_new[lev-1]->nGrowVect());

    resize_stuff(lev);

    vec_Zt_avg1[lev].reset(new MultiFab(ba2d ,dm,1,IntVect(NGROW+1,NGROW+1,0))); //2d, average of the free surface (zeta)
    vec_h[lev].reset(new MultiFab(ba2d ,dm,2,IntVect(NGROW+1,NGROW+1,0))); //2d, average of the free surface (zeta)
    vec_ubar[lev].reset(new MultiFab(convert(ba2d,IntVect(1,0,0)),dm,3,IntVect(NGROW,NGROW,0)));
    vec_vbar[lev].reset(new MultiFab(convert(ba2d,IntVect(0,1,0)),dm,3,IntVect(NGROW,NGROW,0)));

    vec_ru[lev].reset(new MultiFab(convert(ba,IntVect(1,0,0)),dm,2,IntVect(NGROW,NGROW,0))); // RHS u (incl horizontal and vertical advection)
    vec_rv[lev].reset(new MultiFab(convert(ba,IntVect(0,1,0)),dm,2,IntVect(NGROW,NGROW,0))); // RHS v

    vec_ru2d[lev].reset(new MultiFab(convert(ba2d,IntVect(1,0,0)),dm,2,IntVect(NGROW,NGROW,0))); // RHS u for 2d
    vec_rv2d[lev].reset(new MultiFab(convert(ba2d,IntVect(0,1,0)),dm,2,IntVect(NGROW,NGROW,0))); // RHS v for 2d

    t_new[lev] = time;
    t_old[lev] = time - bogus_large_value;

    init_masks(lev, ba, dm);

    init_stuff(lev, ba, dm);

    cons_new[lev]->setVal(zero);
    xvel_new[lev]->setVal(zero);
    yvel_new[lev]->setVal(zero);
    zvel_new[lev]->setVal(zero);

    cons_old[lev]->setVal(zero);
    xvel_old[lev]->setVal(zero);
    yvel_old[lev]->setVal(zero);
    zvel_old[lev]->setVal(zero);

    vec_ru[lev]->setVal(zero);
    vec_rv[lev]->setVal(zero);

    vec_ru2d[lev]->setVal(zero);
    vec_rv2d[lev]->setVal(zero);

    vec_ubar[lev]->setVal(zero);
    vec_vbar[lev]->setVal(zero);


    FillCoarsePatch(lev, time, cons_new[lev], cons_new[lev-1],BCVars::Temp_bc_comp,BdyVars::t);
    FillCoarsePatch(lev, time, xvel_new[lev], xvel_new[lev-1], xvel_bc(), BdyVars::u);
    FillCoarsePatch(lev, time, yvel_new[lev], yvel_new[lev-1], yvel_bc(), BdyVars::v);
    FillCoarsePatch(lev, time, zvel_new[lev], zvel_new[lev-1], zvel_bc(), BdyVars::null);

    if (lev > hires_grid_level) {
        FillCoarsePatch(lev, time, vec_h[lev].get(), vec_h[lev-1].get(),
                        foextrap_periodic_bc(),BdyVars::null,0,false);
        FillCoarsePatch(lev, time, vec_h[lev].get(), vec_h[lev-1].get(),
                        foextrap_periodic_bc(),BdyVars::null,1,false);
    } else {
        set_bathymetry_averaged_down(lev);
    }

    FillCoarsePatch(lev, time, vec_Zt_avg1[lev].get(), vec_Zt_avg1[lev-1].get(),BCVars::cons_bc);
    for (int icomp=0; icomp<3; icomp++) {
        FillCoarsePatch(lev, time, vec_ubar[lev].get(), vec_ubar[lev-1].get(), ubar_bc(),
                bdy_ubar(),icomp,false);
        FillCoarsePatch(lev, time, vec_vbar[lev].get(), vec_vbar[lev-1].get(), vbar_bc(),
                bdy_vbar(),icomp,false);
    }
    for (int icomp=0; icomp<2; icomp++) {
        FillCoarsePatch(lev, time, vec_ru[lev].get(), vec_ru[lev-1].get(), xvel_bc(),
                BdyVars::null,icomp,false);
        FillCoarsePatch(lev, time, vec_rv[lev].get(), vec_rv[lev-1].get(), yvel_bc(),
                BdyVars::null,icomp,false);
        FillCoarsePatch(lev, time, vec_ru2d[lev].get(), vec_ru2d[lev-1].get(), xvel_bc(),
                BdyVars::null,icomp,false);
        FillCoarsePatch(lev, time, vec_rv2d[lev].get(), vec_rv2d[lev-1].get(), yvel_bc(),
                BdyVars::null,icomp,false);
    }

    // Not totally sure foextrap is right here
    FillCoarsePatchPC(lev, time, vec_mskr[lev].get(), vec_mskr[lev-1].get(),
            foextrap_bc());

    calculate_nodal_masks(lev);


    set_grid_scale(lev);
    stretch_transform(lev);

    init_set_vmix(lev);
    set_hmixcoef(lev);
    set_coriolis(lev);
    bool apply_eminusp = false;
    set_zeta_to_Ztavg(lev,apply_eminusp);
    // Previously set smflux

#ifdef REMORA_USE_NETCDF
    if (solverChoice.do_rivers) {
        init_riv_pos_from_netcdf(lev);
    }
#endif

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

/**
 * Remake an existing level using provided BoxArray and DistributionMapping and
 * fill with existing fine and coarse data.
 * overrides the pure virtual function in AmrCore
 * @param[in   ] lev     level to make
 * @param[in   ] time    current time
 * @param[in   ] ba      BoxArray for the level
 * @param[in   ] dm      DistributionMapping for the level
 */
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

    amrex::Print() << "Remaking level " << lev << std::endl;
    amrex::Print() << "GRIDS AT LEVEL " << lev << " ARE " << ba << std::endl;

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

    MultiFab tmp_cons_new(ba, dm, ncons, ngrow_state);
    MultiFab tmp_cons_old(ba, dm, ncons, ngrow_state);

    MultiFab tmp_xvel_new(convert(ba, IntVect(1,0,0)), dm, 1, ngrow_vels);
    MultiFab tmp_xvel_old(convert(ba, IntVect(1,0,0)), dm, 1, ngrow_vels);

    MultiFab tmp_yvel_new(convert(ba, IntVect(0,1,0)), dm, 1, ngrow_vels);
    MultiFab tmp_yvel_old(convert(ba, IntVect(0,1,0)), dm, 1, ngrow_vels);

    MultiFab tmp_zvel_new(convert(ba, IntVect(0,0,1)), dm, 1, IntVect(ngrow_vels,ngrow_vels,0));
    MultiFab tmp_zvel_old(convert(ba, IntVect(0,0,1)), dm, 1, IntVect(ngrow_vels,ngrow_vels,0));

    MultiFab tmp_Zt_avg1_new(ba2d, dm, 1, IntVect(ngrow_zeta,ngrow_zeta,0));
    MultiFab tmp_h(ba2d, dm, 2, IntVect(ngrow_h,ngrow_h,0));

    MultiFab tmp_ubar_new(convert(ba2d, IntVect(1,0,0)), dm, 3, IntVect(ngrow_velbar,ngrow_velbar,0));

    MultiFab tmp_vbar_new(convert(ba2d, IntVect(0,1,0)), dm, 3, IntVect(ngrow_velbar,ngrow_velbar,0));

    MultiFab tmp_ru_new(convert(ba, IntVect(1,0,0)),dm,2,IntVect(NGROW,NGROW,0));
    MultiFab tmp_rv_new(convert(ba, IntVect(0,1,0)),dm,2,IntVect(NGROW,NGROW,0));

    MultiFab tmp_ru2d_new(convert(ba2d, IntVect(1,0,0)),dm,2,IntVect(NGROW,NGROW,0));
    MultiFab tmp_rv2d_new(convert(ba2d, IntVect(0,1,0)),dm,2,IntVect(NGROW,NGROW,0));

    init_masks(lev, ba, dm);

    tmp_cons_new.setVal(zero);
    tmp_xvel_new.setVal(zero);
    tmp_yvel_new.setVal(zero);
    tmp_zvel_new.setVal(zero);

    tmp_cons_old.setVal(zero);
    tmp_xvel_old.setVal(zero);
    tmp_yvel_old.setVal(zero);
    tmp_zvel_old.setVal(zero);

    tmp_ru_new.setVal(zero);
    tmp_rv_new.setVal(zero);

    tmp_ru2d_new.setVal(zero);
    tmp_rv2d_new.setVal(zero);

    tmp_ubar_new.setVal(zero);
    tmp_vbar_new.setVal(zero);


    // This will fill the temporary MultiFabs with data from previous fine data as well as coarse where needed
    FillPatch(lev, time, tmp_cons_new, cons_new, BCVars::cons_bc, BdyVars::t,0,true,false);
    FillPatch(lev, time, tmp_xvel_new, xvel_new, xvel_bc(), BdyVars::u,0,true,false,0,0,zero,tmp_xvel_new);
    FillPatch(lev, time, tmp_yvel_new, yvel_new, yvel_bc(), BdyVars::v,0,true,false,0,0,zero,tmp_yvel_new);
    FillPatch(lev, time, tmp_zvel_new, zvel_new, zvel_bc(), BdyVars::null,0,true,false);
    FillPatch(lev, time, tmp_Zt_avg1_new, GetVecOfPtrs(vec_Zt_avg1), zeta_bc(), BdyVars::null,0,true,false);

    for (int icomp=0; icomp<3; icomp++) {
        FillPatch(lev, time, tmp_ubar_new, GetVecOfPtrs(vec_ubar), ubar_bc(), bdy_ubar(), icomp,false,false);
        FillPatch(lev, time, tmp_vbar_new, GetVecOfPtrs(vec_vbar), vbar_bc(), bdy_vbar(), icomp,false,false);
    }
    for (int icomp=0; icomp<2; icomp++) {
        FillPatch(lev, time, tmp_ru_new, GetVecOfPtrs(vec_ru), xvel_bc(), BdyVars::null, icomp,false,false);
        FillPatch(lev, time, tmp_rv_new, GetVecOfPtrs(vec_rv), yvel_bc(), BdyVars::null, icomp,false,false);
        // These might want to have BCVars::ubar_bc and vbar_bc
        FillPatch(lev, time, tmp_ru2d_new, GetVecOfPtrs(vec_ru2d), xvel_bc(), BdyVars::null, icomp,false,false);
        FillPatch(lev, time, tmp_rv2d_new, GetVecOfPtrs(vec_rv2d), yvel_bc(), BdyVars::null, icomp,false,false);
    }

    MultiFab::Copy(tmp_cons_old,tmp_cons_new,0,0,ncons,tmp_cons_new.nGrowVect());
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
    std::swap(tmp_ubar_new,    *vec_ubar[lev]);
    std::swap(tmp_vbar_new,    *vec_vbar[lev]);
    std::swap(tmp_ru_new,    *vec_ru[lev]);
    std::swap(tmp_rv_new,    *vec_rv[lev]);
    std::swap(tmp_ru2d_new,    *vec_ru2d[lev]);
    std::swap(tmp_rv2d_new,    *vec_rv2d[lev]);

    // Handle bathymetry separately
    if (lev > hires_grid_level) {
        FillPatch(lev, time, tmp_h, GetVecOfPtrs(vec_h), foextrap_periodic_bc(), BdyVars::null,0,false,false);
        FillPatch(lev, time, tmp_h, GetVecOfPtrs(vec_h), foextrap_periodic_bc(), BdyVars::null,1,false,false);
        std::swap(tmp_h,           *vec_h[lev]);
    } else {
        set_bathymetry_averaged_down(lev);
    }

    t_new[lev] = time;
    t_old[lev] = time - bogus_large_value;

    init_masks(lev, ba, dm);
    FillCoarsePatchPC(lev, time, vec_mskr[lev].get(), vec_mskr[lev-1].get(),
            foextrap_bc());
    calculate_nodal_masks(lev);

    init_stuff(lev, ba, dm);

    set_grid_scale(lev);
    stretch_transform(lev);

    init_set_vmix(lev);
    set_hmixcoef(lev);
    set_coriolis(lev);
    bool apply_eminusp = false;
    set_zeta_to_Ztavg(lev,apply_eminusp);
    // Previously set smflux here

#ifdef REMORA_USE_NETCDF
    if (solverChoice.do_rivers) {
        init_riv_pos_from_netcdf(lev);
    }
#endif

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

/**
 * Make a new level from scratch using provided BoxArray and DistributionMapping.
 * This is called both for initialization and for restart
 * (overrides the pure virtual function in AmrCore)
 * main.cpp --> REMORA::InitData --> InitFromScratch --> MakeNewGrids --> MakeNewLevelFromScratch
 *                                        restart  --> MakeNewGrids --> MakeNewLevelFromScratch
 *
 * @param[in   ] lev     level to make
 * @param[in   ] time    current time
 * @param[in   ] ba      BoxArray for the level
 * @param[in   ] dm      DistributionMapping for the level
 */
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

    amrex::Print() << "Making level " << lev << " from scratch" << std::endl;
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

    cons_old[lev] = new MultiFab(ba, dm, ncons, ngrow_state);
    cons_new[lev] = new MultiFab(ba, dm, ncons, ngrow_state);

    xvel_new[lev] = new MultiFab(convert(ba, IntVect(1,0,0)), dm, 1, ngrow_vels);
    xvel_old[lev] = new MultiFab(convert(ba, IntVect(1,0,0)), dm, 1, ngrow_vels);

    yvel_new[lev] = new MultiFab(convert(ba, IntVect(0,1,0)), dm, 1, ngrow_vels);
    yvel_old[lev] = new MultiFab(convert(ba, IntVect(0,1,0)), dm, 1, ngrow_vels);

    zvel_new[lev] = new MultiFab(convert(ba, IntVect(0,0,1)), dm, 1, IntVect(ngrow_vels,ngrow_vels,0));
    zvel_old[lev] = new MultiFab(convert(ba, IntVect(0,0,1)), dm, 1, IntVect(ngrow_vels,ngrow_vels,0));

    resize_stuff(lev);

    vec_Zt_avg1[lev].reset(new MultiFab(ba2d ,dm,1,IntVect(NGROW+1,NGROW+1,0))); //2d, average of the free surface (zeta)
    vec_h[lev].reset(new MultiFab(ba2d ,dm,2,IntVect(NGROW+1,NGROW+1,0))); //2d, bathymetry
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

/**
 * @param[in   ] lev    level do operate on
 */
void REMORA::resize_stuff(int lev)
{
    vec_z_phys_nd.resize(lev+1);

    vec_h_full_domain.resize(hires_grid_level+1);

    vec_h.resize(lev+1);
    vec_Zt_avg1.resize(lev+1);
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
    vec_btflx.resize(lev+1);
    vec_stflx.resize(lev+1);
    vec_btflux.resize(lev+1);
    vec_stflux.resize(lev+1);
    vec_lrflx.resize(lev+1);
    vec_longwave_down.resize(lev+1);
    vec_lhflx.resize(lev+1);
    vec_shflx.resize(lev+1);
    vec_rain.resize(lev+1);
    vec_evap.resize(lev+1);
    vec_rdrag.resize(lev+1);
    vec_rdrag2.resize(lev+1);
    vec_ZoBot.resize(lev+1);
    vec_bustr.resize(lev+1);
    vec_bvstr.resize(lev+1);
    vec_uwind.resize(lev+1);
    vec_vwind.resize(lev+1);
    vec_Tair.resize(lev+1);
    vec_qair.resize(lev+1);
    vec_Pair.resize(lev+1);
    vec_srflx.resize(lev+1);
    vec_cloud.resize(lev+1);
    vec_EminusP.resize(lev+1);
    vec_alpha.resize(lev+1);
    vec_beta.resize(lev+1);

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
    vec_zeta_full_domain.resize(hires_init_level+1);
    vec_mskr.resize(lev+1);
    vec_msku.resize(lev+1);
    vec_mskv.resize(lev+1);
    vec_mskp.resize(lev+1);
    vec_mskr3d.resize(lev+1);
    vec_sstore.resize(lev+1);

    vec_cons_full_domain.resize(hires_init_level+1);
    vec_xvel_full_domain.resize(hires_init_level+1);
    vec_yvel_full_domain.resize(hires_init_level+1);

    vec_pm.resize(lev+1);
    vec_pn.resize(lev+1);
    vec_fcor.resize(lev+1);
    vec_pm_full_domain.resize(hires_grid_level+1);
    vec_pn_full_domain.resize(hires_grid_level+1);

    vec_xr.resize(lev+1);
    vec_yr.resize(lev+1);
    vec_xu.resize(lev+1);
    vec_yu.resize(lev+1);
    vec_xv.resize(lev+1);
    vec_yv.resize(lev+1);
    vec_xp.resize(lev+1);
    vec_yp.resize(lev+1);
    vec_lonp.resize(lev+1);
    vec_latp.resize(lev+1);

    vec_dndx.resize(lev+1);
    vec_dmde.resize(lev+1);

    vec_rhoS.resize(lev+1);
    vec_rhoA.resize(lev+1);
    vec_bvf.resize(lev+1);

    vec_tke.resize(lev+1);
    vec_gls.resize(lev+1);
    vec_Lscale.resize(lev+1);
    vec_Akk.resize(lev+1);
    vec_Akp.resize(lev+1);

    vec_river_position.resize(lev+1);

    if (lev==0) vec_nudg_coeff.resize(num_bdy_vars());

    vec_nudg_coeff[BdyVars::u].resize(lev+1);
    vec_nudg_coeff[BdyVars::v].resize(lev+1);
    for (int icomp = 0; icomp < ncons; ++icomp) {
        vec_nudg_coeff[BdyVars::cons(icomp)].resize(lev+1);
    }
    vec_nudg_coeff[bdy_ubar()].resize(lev+1);
    vec_nudg_coeff[bdy_vbar()].resize(lev+1);
    vec_nudg_coeff[bdy_zeta()].resize(lev+1);
}

/**
 * @param[in   ] lev    level to operate on
 * @param[in   ] ba     BoxArray for the level
 * @param[in   ] dm     DistributionMapping for the level
 */
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

    vec_mskr3d[lev].reset(new MultiFab(ba,dm,1,IntVect(NGROW+1,NGROW+1,0)));

    vec_mskr[lev]->setVal(one);
    vec_msku[lev]->setVal(one);
    vec_mskv[lev]->setVal(one);
    vec_mskp[lev]->setVal(one);

    vec_mskr3d[lev]->setVal(one);
}

/**
 * @param[in   ] lev    level to operate on
 * @param[in   ] ba     BoxArray for the level
 * @param[in   ] dm     DistributionMapping for the level
 */
void REMORA::init_stuff (int lev, const BoxArray& ba, const DistributionMapping& dm)
{
    // ********************************************************************************************
    // Initialize the boundary conditions
    // ********************************************************************************************
    physbcs[lev] = std::make_unique<REMORAPhysBCFunct> (lev, geom[lev], domain_bcs_type, domain_bcs_type_d,
                                                       m_bc_extdir_vals, ncons);

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

    BoxArray ba_nd(ba);
    ba_nd.surroundingNodes();
    BoxArray ba_w(ba);
    ba_w.surroundingNodes(2);

    vec_z_phys_nd[lev].reset          (new MultiFab(ba_nd,dm,1,IntVect(NGROW,NGROW,1))); // z at psi points (nodes) MIGHT NEED NGROW+1
    vec_z_w[lev].reset                (new MultiFab(convert(ba,IntVect(0,0,1)),dm,1,IntVect(NGROW+1,NGROW+1,0))); // z at w points (cell faces)
    vec_z_r[lev].reset                (new MultiFab(ba,dm,1,IntVect(NGROW+1,NGROW+1,0))); // z at r points (cell center)
    vec_Hz[lev].reset                 (new MultiFab(ba,dm,1,IntVect(NGROW+1,NGROW+1,NGROW+1))); // like in ROMS, thickness of cell in z

    vec_Huon[lev].reset               (new MultiFab(convert(ba,IntVect(1,0,0)),dm,1,IntVect(NGROW,NGROW,0))); // mass flux for u component
    vec_Hvom[lev].reset               (new MultiFab(convert(ba,IntVect(0,1,0)),dm,1,IntVect(NGROW,NGROW,0))); // mass flux for v component

    vec_Akv[lev].reset                (new MultiFab(convert(ba,IntVect(0,0,1)),dm,1,IntVect(NGROW,NGROW,0))); // vertical mixing coefficient (.in)
    // NAT components, not ncons: vertical mixing coefficients exist for the active tracers
    // alone, and passive tracers mix with the salinity one. See akt_comp().
    vec_Akt[lev].reset                (new MultiFab(convert(ba,IntVect(0,0,1)),dm,NAT,IntVect(NGROW,NGROW,0))); // vertical mixing coefficient (.in)

    // check dimensionality
    vec_visc2_p[lev].reset(new MultiFab(ba2d,dm,1,IntVect(NGROW,NGROW,0))); // harmonic viscosity at psi points -- difference to 3d?
    vec_visc2_r[lev].reset(new MultiFab(ba2d,dm,1,IntVect(NGROW,NGROW,0))); // harmonic viscosity at rho points
    vec_diff2[lev].reset(new MultiFab(ba2d,dm,ncons,IntVect(NGROW,NGROW,0))); // harmonic diffusivity temperature/salt

    //2d, (incl advection terms and surface/bottom stresses, integral over the whole column, k=0)
    vec_rufrc[lev].reset(new MultiFab(convert(ba2d,IntVect(1,0,0)),dm,2,IntVect(NGROW,NGROW,0)));
    vec_rvfrc[lev].reset(new MultiFab(convert(ba2d,IntVect(0,1,0)),dm,2,IntVect(NGROW,NGROW,0))); //2d, same as above but v

    vec_sustr[lev].reset(new MultiFab(convert(ba2d,IntVect(1,0,0)),dm,1,IntVect(NGROW,NGROW,0))); //2d, surface stress
    vec_svstr[lev].reset(new MultiFab(convert(ba2d,IntVect(0,1,0)),dm,1,IntVect(NGROW,NGROW,0))); //2d

    if (solverChoice.bottom_stress_type == BottomStressType::linear) {
        //2d, linear drag coefficient [m/s], defined at rho, somehow related to rdrg
        vec_rdrag[lev].reset(new MultiFab(ba2d,dm,1,IntVect(NGROW,NGROW,0)));
    } else if (solverChoice.bottom_stress_type == BottomStressType::quadratic) {
        vec_rdrag2[lev].reset(new MultiFab(ba2d,dm,1,IntVect(NGROW,NGROW,0)));
    }

    if (solverChoice.bottom_stress_type == BottomStressType::logarithmic ||
        solverChoice.vert_mixing_type == VertMixingType::GLS) {
        vec_ZoBot[lev].reset(new MultiFab(ba2d,dm,1,IntVect(NGROW,NGROW,0)));
    }

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

    if (solverChoice.use_curvilinear_grid) {
        vec_dndx[lev].reset(new MultiFab(convert(ba2d,IntVect(1,0,0)),dm,1,IntVect(NGROW,NGROW,0)));
        vec_dmde[lev].reset(new MultiFab(convert(ba2d,IntVect(1,0,0)),dm,1,IntVect(NGROW,NGROW,0)));
    }

    // tempstore, saltstore, etc
    vec_sstore[lev].reset(new MultiFab(ba,dm,ncons,IntVect(NGROW,NGROW,0)));

    vec_rhoS[lev].reset(new MultiFab(ba,dm,1,IntVect(NGROW,NGROW,0)));
    vec_rhoA[lev].reset(new MultiFab(ba,dm,1,IntVect(NGROW,NGROW,0)));
    vec_bvf[lev].reset(new MultiFab(convert(ba,IntVect(0,0,1)),dm,1,IntVect(NGROW,NGROW,0)));

    vec_tke[lev].reset(new MultiFab(convert(ba,IntVect(0,0,1)),dm,3,IntVect(NGROW,NGROW,0)));
    vec_gls[lev].reset(new MultiFab(convert(ba,IntVect(0,0,1)),dm,3,IntVect(NGROW,NGROW,0)));
    vec_Lscale[lev].reset(new MultiFab(convert(ba,IntVect(0,0,1)),dm,1,IntVect(NGROW,NGROW,0)));
    vec_Akk[lev].reset(new MultiFab(convert(ba,IntVect(0,0,1)),dm,1,IntVect(NGROW,NGROW,0)));
    vec_Akp[lev].reset(new MultiFab(convert(ba,IntVect(0,0,1)),dm,1,IntVect(NGROW,NGROW,0)));

    // surface/bottom tracer fluxes for update
    vec_stflx[lev].reset(new MultiFab(ba2d,dm,ncons,IntVect(NGROW,NGROW,0)));
    vec_btflx[lev].reset(new MultiFab(ba2d,dm,ncons,IntVect(NGROW,NGROW,0)));
    // surface/bottom tracer fluxes to be filled by inputs
    vec_stflux[lev].reset(new MultiFab(ba2d,dm,ncons,IntVect(NGROW,NGROW,0)));
    vec_btflux[lev].reset(new MultiFab(ba2d,dm,ncons,IntVect(NGROW,NGROW,0)));

    if (solverChoice.bulk_fluxes || solverChoice.atm2ocn_flux_mode) {
        vec_uwind[lev].reset(new MultiFab(ba2d,dm,1,IntVect(NGROW,NGROW,0))); //2d, surface wind u
        vec_vwind[lev].reset(new MultiFab(ba2d,dm,1,IntVect(NGROW,NGROW,0))); //2d, surface wind v
        vec_Tair[lev].reset(new MultiFab(ba2d,dm,1,IntVect(NGROW,NGROW,0)));  //2d, air temperature
        vec_qair[lev].reset(new MultiFab(ba2d,dm,1,IntVect(NGROW,NGROW,0)));  //2d, specific humidity
        vec_Pair[lev].reset(new MultiFab(ba2d,dm,1,IntVect(NGROW,NGROW,0)));  //2d, air pressure
        vec_srflx[lev].reset(new MultiFab(ba2d,dm,1,IntVect(NGROW,NGROW,0))); //2d, shortwave radiation flux
        vec_longwave_down[lev].reset(new MultiFab(ba2d,dm,1,IntVect(NGROW,NGROW,0))); //2d, external longwave radiation flux
        vec_cloud[lev].reset(new MultiFab(ba2d,dm,1,IntVect(NGROW,NGROW,0))); //2d, cloud cover fraction
        vec_EminusP[lev].reset(new MultiFab(ba2d,dm,1,IntVect(NGROW,NGROW,0))); //2d, evaporation minus precipitation
        vec_alpha[lev].reset(new MultiFab(ba2d,dm,1,IntVect(NGROW,NGROW,0)));
        vec_beta[lev].reset(new MultiFab(ba2d,dm,1,IntVect(NGROW,NGROW,0)));
        vec_lrflx[lev].reset(new MultiFab(ba2d,dm,1,IntVect(NGROW,NGROW,0)));
        vec_lhflx[lev].reset(new MultiFab(ba2d,dm,1,IntVect(NGROW,NGROW,0)));
        vec_shflx[lev].reset(new MultiFab(ba2d,dm,1,IntVect(NGROW,NGROW,0)));
        vec_rain[lev].reset(new MultiFab(ba2d,dm,1,IntVect(NGROW,NGROW,0)));
        vec_evap[lev].reset(new MultiFab(ba2d,dm,1,IntVect(NGROW,NGROW,0)));

        vec_uwind[lev]->setVal(solverChoice.Uwind);
        vec_vwind[lev]->setVal(solverChoice.Vwind);
        vec_Tair[lev]->setVal(solverChoice.Tair);
        vec_qair[lev]->setVal(solverChoice.Hair); // Hair can be specific humidity or RH
        vec_Pair[lev]->setVal(solverChoice.Pair);
        vec_srflx[lev]->setVal(solverChoice.srflux);
        vec_longwave_down[lev]->setVal(solverChoice.longwave_rad);
        vec_cloud[lev]->setVal(solverChoice.cloud);
        vec_EminusP[lev]->setVal(solverChoice.EminusP);
        vec_rain[lev]->setVal(solverChoice.rain);

        // Set flux vars that will be computed in bulk_fluxes to zero so initial plotting works
        vec_stflx[lev]->setVal(zero);
        vec_sustr[lev]->setVal(zero);
        vec_svstr[lev]->setVal(zero);
        vec_lhflx[lev]->setVal(zero);
        vec_shflx[lev]->setVal(zero);
        vec_lrflx[lev]->setVal(zero); // possibly this should be set to longwave_rad like longwave_down
    }

    if (solverChoice.do_rivers) {
        vec_river_position[lev].reset(new iMultiFab(ba2d,dm,1,IntVect(NGROW,NGROW,0)));
        vec_river_position[lev]->setVal(-1);
    }

    vec_nudg_coeff[BdyVars::u][lev].reset(new MultiFab(ba,dm,1,IntVect(NGROW+1,NGROW+1,0)));
    vec_nudg_coeff[BdyVars::v][lev].reset(new MultiFab(ba,dm,1,IntVect(NGROW+1,NGROW+1,0)));
    for (int icomp = 0; icomp < ncons; ++icomp) {
        vec_nudg_coeff[BdyVars::cons(icomp)][lev].reset(new MultiFab(ba,dm,1,IntVect(NGROW+1,NGROW+1,0)));
    }
    vec_nudg_coeff[bdy_ubar()][lev].reset(new MultiFab(ba2d,dm,1,IntVect(NGROW+1,NGROW+1,0)));
    vec_nudg_coeff[bdy_vbar()][lev].reset(new MultiFab(ba2d,dm,1,IntVect(NGROW+1,NGROW+1,0)));
    vec_nudg_coeff[bdy_zeta()][lev].reset(new MultiFab(ba2d,dm,1,IntVect(NGROW+1,NGROW+1,0)));

    set_weights(lev);

    vec_DU_avg1[lev]->setVal(zero);
    vec_DU_avg2[lev]->setVal(zero);
    vec_DV_avg1[lev]->setVal(zero);
    vec_DV_avg2[lev]->setVal(zero);
    vec_rubar[lev]->setVal(zero);
    vec_rvbar[lev]->setVal(zero);
    vec_rzeta[lev]->setVal(zero);

    // Initialize these vars even if we aren't using GLS to
    // avoid issues on e.g. checkpoint
    vec_tke[lev]->setVal(solverChoice.gls_Kmin);
    vec_gls[lev]->setVal(solverChoice.gls_Pmin);
    vec_Lscale[lev]->setVal(zero);
    vec_Akk[lev]->setVal(solverChoice.Akk_bak);
    vec_Akp[lev]->setVal(solverChoice.Akp_bak);

    vec_stflx[lev]->setVal(zero);
    vec_btflx[lev]->setVal(zero);
    vec_stflux[lev]->setVal(zero);
    vec_btflux[lev]->setVal(zero);

    // NOTE: Used to set vec_pm and vec_pn to 1e34 here to make foextrap work
    // when init_type = real. However, this does not appear to be necessary so removing

    // Set initial linear drag coefficient
    if (solverChoice.bottom_stress_type == BottomStressType::linear) {
        vec_rdrag[lev]->setVal(solverChoice.rdrag);
    } else if (solverChoice.bottom_stress_type == BottomStressType::quadratic) {
        vec_rdrag2[lev]->setVal(solverChoice.rdrag2);
    }

    if (solverChoice.bottom_stress_type == BottomStressType::logarithmic ||
        solverChoice.vert_mixing_type == VertMixingType::GLS) {
        vec_ZoBot[lev]->setVal(solverChoice.Zob);
    }


    // ********************************************************************************************
    // Create the REMORAFillPatcher object
    // ********************************************************************************************
    if (lev > 0 && cf_width >= 0) {
        Construct_REMORAFillPatchers(lev);
           Define_REMORAFillPatchers(lev);
    }
}

/**
 * Delete level data. Overrides the pure virtual function in AmrCore
 *
 * @param[in   ] lev    level to operate on
 */
void
REMORA::ClearLevel (int lev)
{
    delete cons_new[lev]; delete xvel_new[lev];  delete yvel_new[lev];  delete zvel_new[lev];
    delete cons_old[lev]; delete xvel_old[lev];  delete yvel_old[lev];  delete zvel_old[lev];
}

/**
 * @param[in   ] lev    level to operate on
 */
void
REMORA::set_grid_scale (int lev)
{
    // Even if we're using high-resolution grid initialization, don't set it up with average-down
    if (solverChoice.ic_type == IC_Type::analytic) {
        if (solverChoice.grid_scale_type == GridScaleType::constant) {
            const auto dxi = Geom(lev).InvCellSize();
            vec_pm[lev]->setVal(dxi[0]); vec_pm[lev]->FillBoundary(geom[lev].periodicity());
            vec_pn[lev]->setVal(dxi[1]); vec_pn[lev]->FillBoundary(geom[lev].periodicity());
        } else if (solverChoice.grid_scale_type == GridScaleType::analytic) {
            prob->init_analytic_grid_scale(lev, Geom(lev), solverChoice, *this, *vec_pm[lev].get(), *vec_pn[lev].get());
            vec_pm[lev]->FillBoundary(geom[lev].periodicity());
            vec_pn[lev]->FillBoundary(geom[lev].periodicity());
        }
        set_grid_coords_from_grid_scale(lev);
#ifdef REMORA_USE_NETCDF
    } else if (solverChoice.ic_type == IC_Type::netcdf) {
        if (lev == 0 && hires_grid_level < 0) {
            init_grid_vars_from_netcdf(lev);
        } else if (lev > hires_grid_level) {
            Real dummy_time = zero;
            FillCoarsePatch(lev,dummy_time,vec_pm[lev].get(), vec_pm[lev-1].get(), foextrap_bc());
            FillCoarsePatch(lev,dummy_time,vec_pn[lev].get(), vec_pn[lev-1].get(), foextrap_bc());

            int rrx = ref_ratio[lev-1][0];
            int rry = ref_ratio[lev-1][1];
            // pm and pn need to be rescaled by the refinement ratio
            for ( MFIter mfi(*cons_new[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi )
            {
                Array4<Real> const& pm   = vec_pm[lev]->array(mfi);
                Array4<Real> const& pn   = vec_pn[lev]->array(mfi);
                Box ubx = mfi.growntilebox(IntVect(NGROW+1,NGROW+2,0));
                Box vbx = mfi.growntilebox(IntVect(NGROW+2,NGROW+1,0));
                ParallelFor(makeSlab(ubx,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int ) {
                    pm(i,j,0) = pm(i,j,0) * (rrx);
                });
                ParallelFor(makeSlab(vbx,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int ) {
                    pn(i,j,0) = pn(i,j,0) * (rry);
                });
            }
            set_grid_coords_from_grid_scale(lev);
        } else {
            set_grid_vars_averaged_down(lev);
            set_grid_coords_from_grid_scale(lev);
        }
#endif
    }
    if (solverChoice.use_curvilinear_grid) {
        set_curvilinear_terms_from_grid_scale(lev);
    }
}

/**
 * @param[in   ]lev   level to operate on
 */
void
REMORA::set_grid_coords_from_grid_scale (int lev) {
    for ( MFIter mfi(*vec_xr[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<const Real> const& pm = vec_pm[lev]->const_array(mfi);
        Array4<const Real> const& pn = vec_pn[lev]->const_array(mfi);
        Array4<Real> const& xr = vec_xr[lev]->array(mfi);
        Array4<Real> const& yr = vec_yr[lev]->array(mfi);
        Array4<Real> const& xu = vec_xu[lev]->array(mfi);
        Array4<Real> const& yu = vec_yu[lev]->array(mfi);
        Array4<Real> const& xv = vec_xv[lev]->array(mfi);
        Array4<Real> const& yv = vec_yv[lev]->array(mfi);
        Array4<Real> const& xp = vec_xp[lev]->array(mfi);
        Array4<Real> const& yp = vec_yp[lev]->array(mfi);

        Box bx = mfi.growntilebox(IntVect(NGROW,NGROW,0));
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            xr(i,j,0) = (i + Real(0.5)) / pm(i,j,0);
            yr(i,j,0) = (j + Real(0.5)) / pn(i,j,0);
        });

        ParallelFor(grow(convert(bx,IntVect(1,0,0)),IntVect(-1,0,0)), [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            xu(i,j,0) = i / pm(i,j,0);
            yu(i,j,0) = (j + Real(0.5)) / pn(i,j,0);
        });

        ParallelFor(grow(convert(bx,IntVect(0,1,0)),IntVect(0,-1,0)), [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            xv(i,j,0) = (i + Real(0.5)) / pm(i,j,0);
            yv(i,j,0) = j / pn(i,j,0);
        });

        ParallelFor(grow(convert(bx,IntVect(1,1,0)),IntVect(-1,-1,0)), [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            xp(i,j,0) = i / pm(i,j,0);
            yp(i,j,0) = j / pn(i,j,0);
        });
    }
}

/**
 * @param[in   ]lev   level to operate on
 */
void
REMORA::set_curvilinear_terms_from_grid_scale (int lev) {
    for ( MFIter mfi(*vec_dndx[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<const Real> const& pm = vec_pm[lev]->const_array(mfi);
        Array4<const Real> const& pn = vec_pn[lev]->const_array(mfi);
        Array4<Real> const& dndx = vec_dndx[lev]->array(mfi);
        Array4<Real> const& dmde = vec_dmde[lev]->array(mfi);

        Box bx = mfi.growntilebox(IntVect(NGROW,NGROW,0));
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            dndx(i,j,0) = Real(0.5) * (one / pn(i+1,j  ,0) - one / pn(i-1,j  ,0));
            dmde(i,j,0) = Real(0.5) * (one / pm(i  ,j+1,0) - one / pn(i  ,j-1,0));
        });
    }
}

/**
 * @param[in   ] lev    level to operate on
 * @param[in   ] apply_eminusp    whether to apply the evaporation-minus-precipitation correction to Zt_avg1
 */
void
REMORA::set_zeta_to_Ztavg (int lev, bool apply_eminusp)
{
    BL_PROFILE("REMORA::set_zeta_to_Ztavg()");
    std::unique_ptr<MultiFab>& mf_zeta = vec_zeta[lev];
    std::unique_ptr<MultiFab>& mf_Zt_avg1  = vec_Zt_avg1[lev];
    if (solverChoice.eminusp_correct_ssh && apply_eminusp) {
        for ( MFIter mfi(*vec_zeta[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi )
        {
            Array4<Real> const& Zt_avg1 = (mf_Zt_avg1)->array(mfi);
            Array4<const Real> const& evap = vec_evap[lev]->const_array(mfi);
            Array4<const Real> const& rain = vec_rain[lev]->const_array(mfi);
            Array4<const Real> const& EminusP = vec_EminusP[lev]->const_array(mfi);
            bool use_EminusP_from_input = solverChoice.eminusp &&
                solverChoice.bulk_flux_type[BulkFlux::EminusP] != BulkForcingType::computed;

            Box  bx2 = mfi.growntilebox(IntVect(NGROW,NGROW,0));// bx2.grow(IntVect(NGROW,NGROW,0));

            Real cff = dt[lev] / rhow;
            Real dt_lev = dt[lev];

            ParallelFor(bx2, [=] AMREX_GPU_DEVICE (int i, int j, int )
            {
                if (use_EminusP_from_input) {
                    // EminusP is treated as a kinematic freshwater flux (m/s).
                    Zt_avg1(i,j,0) = Zt_avg1(i,j,0) - EminusP(i,j,0) * dt_lev;
                } else {
                    Zt_avg1(i,j,0) = Zt_avg1(i,j,0) - (evap(i,j,0) - rain(i,j,0)) * cff;
                }
            });
        }
    }
    Gpu::streamSynchronize();

    vec_Zt_avg1[lev]->FillBoundary(geom[lev].periodicity());

    for ( MFIter mfi(*vec_zeta[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Box  bx3 = mfi.tilebox(); bx3.grow(IntVect(NGROW+1,NGROW+1,0));
        Array4<Real> const& zeta = mf_zeta->array(mfi);
        Array4<Real> const& Zt_avg1 = (mf_Zt_avg1)->array(mfi);

        ParallelFor(bx3, 3, [=] AMREX_GPU_DEVICE (int i, int j, int , int n)
        {
            zeta(i,j,0,n) = Zt_avg1(i,j,0);
        });
    }
}

/**
 * @param[in   ] lev    level to operate on
 */
void
REMORA::update_mskp (int lev)
{
    for ( MFIter mfi(*vec_mskr[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<const Real> const& mskr = vec_mskr[lev]->const_array(mfi);
        Array4<      Real> const& mskp = vec_mskp[lev]->array(mfi);

        Box bx = mfi.tilebox(); bx.grow(IntVect(1,1,0)); bx.makeSlab(2,0);

        Real cff1 = one;
        Real cff2 = two;

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            if ((mskr(i-1,j,0) > Real(0.5)) and (mskr(i,j,0) > Real(0.5)) and (mskr(i-1,j-1,0) > Real(0.5)) and (mskr(i,j-1,0) > Real(0.5))) {
                mskp(i,j,0) = one;
            } else if ((mskr(i-1,j,0) < Real(0.5)) and (mskr(i,j,0) > Real(0.5)) and (mskr(i-1,j-1,0) > Real(0.5)) and (mskr(i,j-1,0) > Real(0.5))) {
                mskp(i,j,0) = cff1;
            } else if ((mskr(i-1,j,0) > Real(0.5)) and (mskr(i,j,0) < Real(0.5)) and (mskr(i-1,j-1,0) > Real(0.5)) and (mskr(i,j-1,0) > Real(0.5))) {
                mskp(i,j,0) = cff1;
            } else if ((mskr(i-1,j,0) > Real(0.5)) and (mskr(i,j,0) > Real(0.5)) and (mskr(i-1,j-1,0) < Real(0.5)) and (mskr(i,j-1,0) > Real(0.5))) {
                mskp(i,j,0) = cff1;
            } else if ((mskr(i-1,j,0) > Real(0.5)) and (mskr(i,j,0) > Real(0.5)) and (mskr(i-1,j-1,0) > Real(0.5)) and (mskr(i,j-1,0) < Real(0.5))) {
                mskp(i,j,0) = cff1;
            } else if ((mskr(i-1,j,0) > Real(0.5)) and (mskr(i,j,0) < Real(0.5)) and (mskr(i-1,j-1,0) > Real(0.5)) and (mskr(i,j-1,0) < Real(0.5))) {
                mskp(i,j,0) = cff2;
            } else if ((mskr(i-1,j,0) < Real(0.5)) and (mskr(i,j,0) > Real(0.5)) and (mskr(i-1,j-1,0) < Real(0.5)) and (mskr(i,j-1,0) > Real(0.5))) {
                mskp(i,j,0) = cff2;
            } else if ((mskr(i-1,j,0) > Real(0.5)) and (mskr(i,j,0) > Real(0.5)) and (mskr(i-1,j-1,0) < Real(0.5)) and (mskr(i,j-1,0) < Real(0.5))) {
                mskp(i,j,0) = cff2;
            } else if ((mskr(i-1,j,0) < Real(0.5)) and (mskr(i,j,0) < Real(0.5)) and (mskr(i-1,j-1,0) > Real(0.5)) and (mskr(i,j-1,0) > Real(0.5))) {
                mskp(i,j,0) = cff2;
            } else {
                mskp(i,j,0) = zero;
            }

        });
    }
}

/**
 * @param[in   ] lev    level to operate on
 */
void
REMORA::calculate_nodal_masks (int lev)
{
    for ( MFIter mfi(*vec_mskr[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<const Real> const& mskr = vec_mskr[lev]->const_array(mfi);
        Array4<      Real> const& msku = vec_msku[lev]->array(mfi);
        Array4<      Real> const& mskv = vec_mskv[lev]->array(mfi);
        Array4<      Real> const& mskp = vec_mskp[lev]->array(mfi);

        Box bx = mfi.tilebox(); bx.grow(IntVect(1,1,0)); bx.makeSlab(2,0);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            msku(i,j,0) = mskr(i-1,j  ,0) * mskr(i,j,0);
            mskv(i,j,0) = mskr(i  ,j-1,0) * mskr(i,j,0);
            mskp(i,j,0) = mskr(i-1,j-1,0) * mskr(i,j,0) * mskr(i-1,j,0) * mskr(i,j-1,0);
        });
    }
}

/**
 * @param[in   ] lev    level to operate on
 */
void
REMORA::fill_3d_masks (int lev)
{
    for ( MFIter mfi(*vec_mskr3d[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<const Real> const& mskr = vec_mskr[lev]->const_array(mfi);
        Array4<      Real> const& mskr3d = vec_mskr3d[lev]->array(mfi);

        Box bx = mfi.tilebox(); bx.grow(IntVect(1,1,0));

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            mskr3d(i,j,k) = mskr(i,j,0);
        });
    }
}
