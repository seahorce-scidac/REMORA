/**
 * \file REMORA_make_new_level.cpp
 */

#include <REMORA.H>
#include <REMORA_prob_common.H>

#include <AMReX_buildInfo.H>
#include <AMReX_OpenMP.H>

using namespace amrex;

namespace {
/**
 * By default OpenMP threads share the boxes of a level, not tiles (see the tile-size
 * default in main.cpp), so a rank with fewer boxes than threads leaves threads idle. Say
 * so once per level, naming the knobs that fix it. With tiling on the threads share the
 * tiles of a box instead, so there is nothing to warn about.
 */
void warn_if_fewer_boxes_than_threads (int lev, const BoxArray& ba, const DistributionMapping& dm)
{
#ifdef AMREX_USE_OMP
    // Tiling is on when a tile is smaller than some box in x or y; inputs often write
    // "1024 1024 1024" to mean no tiling, so compare with the boxes, not a fixed value.
    const int nthreads = OpenMP::get_max_threads();
    const IntVect& tile_size = FabArrayBase::mfiter_tile_size;
    int nlocal = 0;
    bool tiled = false;
    const int myproc = ParallelDescriptor::MyProc();
    for (int i = 0; i < ba.size(); ++i) {
        if (dm[i] == myproc) { ++nlocal; }
        tiled = tiled || ba[i].length(0) > tile_size[0] || ba[i].length(1) > tile_size[1];
    }
    // Collective: every rank reaches it, whatever its thread count
    ParallelDescriptor::ReduceIntMin(nlocal);
    if (nthreads > 1 && !tiled && nlocal < nthreads) {
        amrex::Print() << "WARNING: level " << lev << " has only " << nlocal
                       << " boxes on some rank for " << nthreads << " OpenMP threads;"
                       << " threads beyond the box count idle. Lower amr.max_grid_size"
                       << " (or amr.max_grid_size_x/y) to give every rank at least "
                       << nthreads << " boxes, or set remora.omp_tile_size so the"
                       << " threads share the tiles of a box.\n";
    }
#else
    amrex::ignore_unused(lev, ba, dm);
#endif
}
}

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
    warn_if_fewer_boxes_than_threads(lev, ba, dm);

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


    // Before the fills: FillPatch masks domain-boundary ghosts through physbcs, so these must
    // hold the real coastline rather than init_masks' all-water placeholder. Order matters:
    // init_stuff reallocates the coordinates set_grid_scale fills, which set_masks may read.
    set_grid_scale(lev);
    set_masks(lev);

    // init_stuff reallocated the nudging coefficients; init_only fills them only
    // for a level made from scratch. The radiation conditions below and the
    // climatology nudging read them.
    init_clim_nudg_coeff(lev);

#ifdef REMORA_USE_NETCDF
    // The per-level NetCDF boundary series are otherwise built only by init_only,
    // for a level made from scratch. The FillCoarsePatch calls below apply the
    // file's boundary data to this level and index boundary_series[lev], so a
    // level first created at a regrid needs its series built here.
    if (solverChoice.boundary_from_netcdf) {
        init_bdry_from_netcdf(lev);
    }
#endif

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

    // The free surface before the barotropic velocities: a Flather condition on
    // ubar or vbar at a domain edge reads zeta, and h, on this level.
    FillCoarsePatch(lev, time, vec_Zt_avg1[lev].get(), vec_Zt_avg1[lev-1].get(),BCVars::cons_bc);
    bool apply_eminusp = false;
    set_zeta_to_Ztavg(lev,apply_eminusp);
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

    stretch_transform(lev);

    init_set_vmix(lev);
    set_hmixcoef(lev);
    set_coriolis(lev);

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
    define_flux_register(lev);

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

    BoxList bl2d = ba.boxList();
    for (auto& b : bl2d) {
        b.setRange(2,0);
    }
    BoxArray ba2d(std::move(bl2d));

    amrex::Print() << "Remaking level " << lev << std::endl;
    amrex::Print() << "GRIDS AT LEVEL " << lev << " ARE " << ba << std::endl;
    warn_if_fewer_boxes_than_threads(lev, ba, dm);

    int ngrow_state   = 3;
    int ngrow_vels    = 3;
    int ngrow_zeta    = 3;
    int ngrow_h       = 3;
    int ngrow_velbar  = 2;

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


    // As in MakeNewLevelFromCoarse. init_stuff comes along because set_grid_scale follows it.
    init_stuff(lev, ba, dm);
    set_grid_scale(lev);
    set_masks(lev);
    // init_stuff reallocated the nudging coefficients the fills below may read
    init_clim_nudg_coeff(lev);

    // This will fill the temporary MultiFabs with data from previous fine data as well as coarse where needed
    FillPatch(lev, time, tmp_cons_new, cons_new, BCVars::cons_bc, BdyVars::t,0,true,false,0,0,zero,tmp_cons_new);
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
        // Swap first: vec_h[lev] is still on the pre-regrid BoxArray at this point, and the
        // branch above is what moves it onto the new one. set_bathymetry_averaged_down takes
        // its data from vec_h_full_domain rather than from the old grids, so the old contents
        // are not needed -- but its FillPatch uses the land masks, which init_masks has
        // already rebuilt on the new BoxArray, so leaving vec_h[lev] on the old one mixes two
        // BoxArrays inside REMORAPhysBCFunct.
        std::swap(tmp_h,           *vec_h[lev]);
        set_bathymetry_averaged_down(lev);
    }

    t_new[lev] = time;
    t_old[lev] = time - bogus_large_value;

    stretch_transform(lev);

    init_set_vmix(lev);
    set_hmixcoef(lev);
    set_coriolis(lev);
    bool apply_eminusp = false;
    set_zeta_to_Ztavg(lev,apply_eminusp);

#ifdef REMORA_USE_NETCDF
    if (solverChoice.do_rivers) {
        init_riv_pos_from_netcdf(lev);
    }
#endif

    // Both of these hold this level's layout and the one below it, so both have to be rebuilt
    // whenever either changes. Unconditionally: AmrCore::regrid calls RemakeLevel only when
    // this level's BoxArray changed or the coarser one's did, and in the latter case it
    // passes this level's own ba and dm unchanged -- so testing them here would skip exactly
    // the case where only the coarse side moved.
    if (lev > 0 && cf_width >= 0) {
        Define_REMORAFillPatchers(lev);
        define_flux_register(lev);
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
    warn_if_fewer_boxes_than_threads(lev, ba, dm);

    int ngrow_state = 3;
    int ngrow_vels  = 3;

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
    // Sized for both hires knobs: the initial-state cascade needs a mask at
    // hires_init_level, which may sit above hires_grid_level.
    vec_mskr_full_domain.resize(std::max(hires_grid_level, hires_init_level)+1);

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

    vec_DUon.resize(lev+1);
    vec_DVom.resize(lev+1);
    vec_scratch_2d.resize(lev+1);
    vec_DU_avg1.resize(lev+1);
    vec_DU_avg2.resize(lev+1);
    vec_DV_avg1.resize(lev+1);
    vec_DV_avg2.resize(lev+1);
    vec_Dubar_old.resize(lev+1);
    vec_Dubar_new.resize(lev+1);
    vec_Dvbar_old.resize(lev+1);
    vec_Dvbar_new.resize(lev+1);
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
    // Indexed by the coarse level of a pair, so lev entries would do; sized like the rest to
    // keep the indexing uniform.
    vec_mskr_crse_on_fine.resize(lev+1);
    vec_msku_crse_on_fine.resize(lev+1);
    vec_mskv_crse_on_fine.resize(lev+1);
    vec_sstore.resize(lev+1);

    vec_cons_full_domain.resize(hires_init_level+1);
    vec_xvel_full_domain.resize(hires_init_level+1);
    vec_yvel_full_domain.resize(hires_init_level+1);

    vec_pm.resize(lev+1);
    vec_pn.resize(lev+1);
    vec_metric_r.resize(lev+1);
    vec_metric_u.resize(lev+1);
    vec_metric_v.resize(lev+1);
    vec_metric_p.resize(lev+1);
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
    // Scratch of advance_2d, rewritten on every fast step
    vec_DUon[lev].reset(new MultiFab(convert(ba2d,IntVect(1,0,0)),dm,1,IntVect(NGROW,NGROW,0)));
    vec_DVom[lev].reset(new MultiFab(convert(ba2d,IntVect(0,1,0)),dm,1,IntVect(NGROW,NGROW,0)));
    vec_scratch_2d[lev].reset(new MultiFab(ba2d,dm,n_scratch_2d,IntVect(NGROW+1,NGROW+1,0)));

    //2d DU: sum(height[incl free surface?] * u)
    vec_DU_avg1[lev].reset(new MultiFab(convert(ba2d,IntVect(1,0,0)),dm,1,IntVect(NGROW,NGROW,0)));

    //2d like above, but correct(or)?
    vec_DU_avg2[lev].reset(new MultiFab(convert(ba2d,IntVect(1,0,0)),dm,1,IntVect(NGROW,NGROW,0)));

    vec_DV_avg1[lev].reset(new MultiFab(convert(ba2d,IntVect(0,1,0)),dm,1,IntVect(NGROW,NGROW,0)));
    vec_DV_avg2[lev].reset(new MultiFab(convert(ba2d,IntVect(0,1,0)),dm,1,IntVect(NGROW,NGROW,0)));

    // D*ubar and D*vbar: DU_avg2/on_u and DV_avg2/om_v, at the two ends of this level's step
    vec_Dubar_old[lev].reset(new MultiFab(convert(ba2d,IntVect(1,0,0)),dm,1,IntVect(NGROW,NGROW,0)));
    vec_Dubar_new[lev].reset(new MultiFab(convert(ba2d,IntVect(1,0,0)),dm,1,IntVect(NGROW,NGROW,0)));
    vec_Dvbar_old[lev].reset(new MultiFab(convert(ba2d,IntVect(0,1,0)),dm,1,IntVect(NGROW,NGROW,0)));
    vec_Dvbar_new[lev].reset(new MultiFab(convert(ba2d,IntVect(0,1,0)),dm,1,IntVect(NGROW,NGROW,0)));

    vec_rubar[lev].reset(new MultiFab(convert(ba2d,IntVect(1,0,0)),dm,4,IntVect(NGROW,NGROW,0))); // 2d RHS ubar
    vec_rvbar[lev].reset(new MultiFab(convert(ba2d,IntVect(0,1,0)),dm,4,IntVect(NGROW,NGROW,0)));
    vec_rzeta[lev].reset(new MultiFab(ba2d,dm,4,IntVect(NGROW,NGROW,0))); // 2d RHS zeta

    // starts off kind of like a depth-averaged u, but exists at more points and more timesteps (b/c fast 2D update) than full u
    vec_zeta[lev].reset(new MultiFab(ba2d,dm,3,IntVect(NGROW+1,NGROW+1,0)));  // 2d free surface

    vec_pm[lev].reset(new MultiFab(ba2d,dm,1,IntVect(NGROW+1,NGROW+2,0)));
    vec_pn[lev].reset(new MultiFab(ba2d,dm,1,IntVect(NGROW+2,NGROW+1,0)));
    // As wide as pm and pn allow: each product reads the cells either side of its point,
    // and pm/pn are filled one cell further than every reader of these needs.
    vec_metric_r[lev].reset(new MultiFab(ba2d,dm,MetricR::ncomp,IntVect(NGROW+1,NGROW+1,0)));
    vec_metric_u[lev].reset(new MultiFab(convert(ba2d,IntVect(1,0,0)),dm,MetricU::ncomp,IntVect(NGROW,NGROW,0)));
    vec_metric_v[lev].reset(new MultiFab(convert(ba2d,IntVect(0,1,0)),dm,MetricV::ncomp,IntVect(NGROW,NGROW,0)));
    vec_metric_p[lev].reset(new MultiFab(convert(ba2d,IntVect(1,1,0)),dm,MetricP::ncomp,IntVect(NGROW,NGROW,0)));
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
        vec_dndx[lev].reset(new MultiFab(ba2d,dm,1,IntVect(NGROW,NGROW,0)));
        vec_dmde[lev].reset(new MultiFab(ba2d,dm,1,IntVect(NGROW,NGROW,0)));
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
        vec_lhflx[lev]->setVal(zero);
        vec_shflx[lev]->setVal(zero);
        vec_lrflx[lev]->setVal(zero); // possibly this should be set to longwave_rad like longwave_down
        vec_evap[lev]->setVal(zero);
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
    vec_Dubar_old[lev]->setVal(zero);
    vec_Dubar_new[lev]->setVal(zero);
    vec_Dvbar_old[lev]->setVal(zero);
    vec_Dvbar_new[lev]->setVal(zero);
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
    vec_sustr[lev]->setVal(zero);
    vec_svstr[lev]->setVal(zero);

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
            for ( MFIter mfi(*vec_pm[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi )
            {
                Array4<Real> const& pm   = vec_pm[lev]->array(mfi);
                Array4<Real> const& pn   = vec_pn[lev]->array(mfi);
                Box ubx = mfi.growntilebox(IntVect(NGROW+1,NGROW+2,0));
                Box vbx = mfi.growntilebox(IntVect(NGROW+2,NGROW+1,0));
                ParallelFor(ubx, [=] AMREX_GPU_DEVICE (int i, int j, int ) {
                    pm(i,j,0) = pm(i,j,0) * (rrx);
                });
                ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int ) {
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
    // Every path that sets pm and pn on a level (analytic, NetCDF, averaged down from a
    // high-resolution grid, interpolated from the parent, regrid) ends here.
    set_derived_metrics(lev);
}

/**
 * The metric products the time step needs, computed once from pm and pn rather than
 * with divides inside every kernel (ROMS metrics.F). Each value uses the expression, and
 * operand order, of the kernels it replaces, so they read identical bits.
 *
 * @param[in   ]lev   level to operate on
 */
void
REMORA::set_derived_metrics (int lev)
{
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(*vec_metric_r[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<const Real> const& pm = vec_pm[lev]->const_array(mfi);
        Array4<const Real> const& pn = vec_pn[lev]->const_array(mfi);
        Array4<Real> const& mr = vec_metric_r[lev]->array(mfi);
        Array4<Real> const& mu = vec_metric_u[lev]->array(mfi);
        Array4<Real> const& mv = vec_metric_v[lev]->array(mfi);
        Array4<Real> const& mp = vec_metric_p[lev]->array(mfi);

        ParallelFor(mfi.growntilebox(vec_metric_r[lev]->nGrowVect()),
        [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            mr(i,j,0,MetricR::pmn)    = pm(i,j,0)*pn(i,j,0);
            mr(i,j,0,MetricR::omn)    = one/(pm(i,j,0)*pn(i,j,0));
            mr(i,j,0,MetricR::om_r)   = one/pm(i,j,0);
            mr(i,j,0,MetricR::on_r)   = one/pn(i,j,0);
            mr(i,j,0,MetricR::pmon_r) = pm(i,j,0)/pn(i,j,0);
            mr(i,j,0,MetricR::pnom_r) = pn(i,j,0)/pm(i,j,0);
        });

        ParallelFor(mfi.grownnodaltilebox(0, vec_metric_u[lev]->nGrowVect()),
        [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            const Real om_u = two/(pm(i,j,0)+pm(i-1,j,0));
            const Real on_u = two/(pn(i,j,0)+pn(i-1,j,0));
            mu(i,j,0,MetricU::om_u)   = om_u;
            mu(i,j,0,MetricU::on_u)   = on_u;
            mu(i,j,0,MetricU::omn_u)  = om_u*on_u;
            mu(i,j,0,MetricU::pmon_u) = (pm(i-1,j,0)+pm(i,j,0))/(pn(i-1,j,0)+pn(i,j,0));
            mu(i,j,0,MetricU::pmn4_u) = (pm(i,j,0)+pm(i-1,j,0))*(pn(i,j,0)+pn(i-1,j,0));
        });

        ParallelFor(mfi.grownnodaltilebox(1, vec_metric_v[lev]->nGrowVect()),
        [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            const Real om_v = two/(pm(i,j,0)+pm(i,j-1,0));
            const Real on_v = two/(pn(i,j,0)+pn(i,j-1,0));
            mv(i,j,0,MetricV::om_v)   = om_v;
            mv(i,j,0,MetricV::on_v)   = on_v;
            mv(i,j,0,MetricV::omn_v)  = om_v*on_v;
            mv(i,j,0,MetricV::pnom_v) = (pn(i,j-1,0)+pn(i,j,0))/(pm(i,j-1,0)+pm(i,j,0));
            mv(i,j,0,MetricV::pmn4_v) = (pm(i,j,0)+pm(i,j-1,0))*(pn(i,j,0)+pn(i,j-1,0));
        });

        Box pbx = mfi.tilebox(IntVect(1,1,0), vec_metric_p[lev]->nGrowVect());
        ParallelFor(pbx, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            const Real sum_pm = pm(i-1,j-1,0)+pm(i-1,j,0)+pm(i,j-1,0)+pm(i,j,0);
            const Real sum_pn = pn(i-1,j-1,0)+pn(i-1,j,0)+pn(i,j-1,0)+pn(i,j,0);
            mp(i,j,0,MetricP::om_p)   = Real(4.0)/sum_pm;
            mp(i,j,0,MetricP::on_p)   = Real(4.0)/sum_pn;
            mp(i,j,0,MetricP::pmon_p) = sum_pm/sum_pn;
            mp(i,j,0,MetricP::pnom_p) = sum_pn/sum_pm;
        });
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

        ParallelFor(mfi.tilebox(IntVect(1,0,0),IntVect(NGROW,NGROW,0)), [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            xu(i,j,0) = i / pm(i,j,0);
            yu(i,j,0) = (j + Real(0.5)) / pn(i,j,0);
        });

        ParallelFor(mfi.tilebox(IntVect(0,1,0),IntVect(NGROW,NGROW,0)), [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            xv(i,j,0) = (i + Real(0.5)) / pm(i,j,0);
            yv(i,j,0) = j / pn(i,j,0);
        });

        ParallelFor(mfi.tilebox(IntVect(1,1,0),IntVect(NGROW,NGROW,0)), [=] AMREX_GPU_DEVICE (int i, int j, int)
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
            dmde(i,j,0) = Real(0.5) * (one / pm(i  ,j+1,0) - one / pm(i  ,j-1,0));
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
REMORA::calculate_nodal_masks (int lev)
{
    for ( MFIter mfi(*vec_mskr[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<const Real> const& mskr = vec_mskr[lev]->const_array(mfi);
        Array4<      Real> const& msku = vec_msku[lev]->array(mfi);
        Array4<      Real> const& mskv = vec_mskv[lev]->array(mfi);
        Array4<      Real> const& mskp = vec_mskp[lev]->array(mfi);

        // NGROW rings, not one: the stencil reaches mskr(i-1,j-1), and mskr carries NGROW+1.
        Box bx = mfi.tilebox(); bx.grow(IntVect(NGROW,NGROW,0)); bx.makeSlab(2,0);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            msku(i,j,0) = mskr(i-1,j  ,0) * mskr(i,j,0);
            mskv(i,j,0) = mskr(i  ,j-1,0) * mskr(i,j,0);
            // mskp is 2 along a straight land-sea boundary
            if ((mskr(i-1,j,0) > Real(0.5)) and (mskr(i,j,0) > Real(0.5)) and (mskr(i-1,j-1,0) > Real(0.5)) and (mskr(i,j-1,0) > Real(0.5))) {
                mskp(i,j,0) = one;
            } else if ((mskr(i-1,j,0) < Real(0.5)) and (mskr(i,j,0) > Real(0.5)) and (mskr(i-1,j-1,0) > Real(0.5)) and (mskr(i,j-1,0) > Real(0.5))) {
                mskp(i,j,0) = one;
            } else if ((mskr(i-1,j,0) > Real(0.5)) and (mskr(i,j,0) < Real(0.5)) and (mskr(i-1,j-1,0) > Real(0.5)) and (mskr(i,j-1,0) > Real(0.5))) {
                mskp(i,j,0) = one;
            } else if ((mskr(i-1,j,0) > Real(0.5)) and (mskr(i,j,0) > Real(0.5)) and (mskr(i-1,j-1,0) < Real(0.5)) and (mskr(i,j-1,0) > Real(0.5))) {
                mskp(i,j,0) = one;
            } else if ((mskr(i-1,j,0) > Real(0.5)) and (mskr(i,j,0) > Real(0.5)) and (mskr(i-1,j-1,0) > Real(0.5)) and (mskr(i,j-1,0) < Real(0.5))) {
                mskp(i,j,0) = one;
            } else if ((mskr(i-1,j,0) > Real(0.5)) and (mskr(i,j,0) < Real(0.5)) and (mskr(i-1,j-1,0) > Real(0.5)) and (mskr(i,j-1,0) < Real(0.5))) {
                mskp(i,j,0) = two;
            } else if ((mskr(i-1,j,0) < Real(0.5)) and (mskr(i,j,0) > Real(0.5)) and (mskr(i-1,j-1,0) < Real(0.5)) and (mskr(i,j-1,0) > Real(0.5))) {
                mskp(i,j,0) = two;
            } else if ((mskr(i-1,j,0) > Real(0.5)) and (mskr(i,j,0) > Real(0.5)) and (mskr(i-1,j-1,0) < Real(0.5)) and (mskr(i,j-1,0) < Real(0.5))) {
                mskp(i,j,0) = two;
            } else if ((mskr(i-1,j,0) < Real(0.5)) and (mskr(i,j,0) < Real(0.5)) and (mskr(i-1,j-1,0) > Real(0.5)) and (mskr(i,j-1,0) > Real(0.5))) {
                mskp(i,j,0) = two;
            } else {
                mskp(i,j,0) = zero;
            }
        });
    }
}

/**
 * Rebuild the u-, v- and psi-point masks from vec_mskr and fill their ghost cells.
 *
 * Every path goes through here, so the levels cannot end up with different definitions.
 *
 * @param[in   ] lev    level to operate on
 */
void
REMORA::update_nodal_masks (int lev)
{
    calculate_nodal_masks(lev);

    vec_msku[lev]->FillBoundary(geom[lev].periodicity());
    vec_mskv[lev]->FillBoundary(geom[lev].periodicity());
    vec_mskp[lev]->FillBoundary(geom[lev].periodicity());
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

        // No stencil, so this can fill every ring mskr has.
        Box bx = mfi.tilebox(); bx.grow(IntVect(NGROW+1,NGROW+1,0));

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            mskr3d(i,j,k) = mskr(i,j,0);
        });
    }
}
