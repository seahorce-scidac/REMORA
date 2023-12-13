/**
 * \file ROMSX_make_new_level.cpp
 */

#include <prob_common.H>
#include <EOS.H>
#include <ROMSX.H>

#include <AMReX_buildInfo.H>

using namespace amrex;

// Make a new level using provided BoxArray and DistributionMapping and
// fill with interpolated coarse level data (overrides the pure virtual function in AmrCore)
// regrid  --> RemakeLevel            (if level already existed)
// regrid  --> MakeNewLevelFromCoarse (if adding new level)
void
ROMSX::MakeNewLevelFromCoarse (int lev, Real time, const BoxArray& ba,
                               const DistributionMapping& dm)
{
    cons_new[lev]->define(ba, dm, NCONS, cons_new[lev-1]->nGrowVect());
    cons_old[lev]->define(ba, dm, NCONS, cons_new[lev-1]->nGrowVect());

    xvel_new[lev]->define(convert(ba, IntVect(1,0,0)), dm, 1, xvel_new[lev-1]->nGrowVect());
    xvel_old[lev]->define(convert(ba, IntVect(1,0,0)), dm, 1, xvel_new[lev-1]->nGrowVect());

    yvel_new[lev]->define(convert(ba, IntVect(0,1,0)), dm, 1, yvel_new[lev-1]->nGrowVect());
    yvel_old[lev]->define(convert(ba, IntVect(0,1,0)), dm, 1, yvel_new[lev-1]->nGrowVect());

    zvel_new[lev]->define(convert(ba, IntVect(0,0,1)), dm, 1, zvel_new[lev-1]->nGrowVect());
    zvel_old[lev]->define(convert(ba, IntVect(0,0,1)), dm, 1, zvel_new[lev-1]->nGrowVect());

    resize_stuff(lev);
      init_stuff(lev, ba, dm);

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    FillCoarsePatch(lev, time, cons_new[lev], cons_new[lev-1]);
    FillCoarsePatch(lev, time, xvel_new[lev], xvel_new[lev-1]);
    FillCoarsePatch(lev, time, yvel_new[lev], yvel_new[lev-1]);
    FillCoarsePatch(lev, time, zvel_new[lev], zvel_new[lev-1]);
}

// Remake an existing level using provided BoxArray and DistributionMapping and
// fill with existing fine and coarse data.
// overrides the pure virtual function in AmrCore
void
ROMSX::RemakeLevel (int lev, Real time, const BoxArray& ba, const DistributionMapping& dm)
{
#if (NGROW==2)
    int ngrow_state = ComputeGhostCells(solverChoice.spatial_order)+1;
    int ngrow_vels  = ComputeGhostCells(solverChoice.spatial_order)+1;
#else
    int ngrow_state = ComputeGhostCells(solverChoice.spatial_order)+2;
    int ngrow_vels  = ComputeGhostCells(solverChoice.spatial_order)+2;
#endif

    MultiFab* tmp_cons_new; tmp_cons_new->define(ba, dm, NCONS, ngrow_state);
    MultiFab* tmp_cons_old; tmp_cons_old->define(ba, dm, NCONS, ngrow_state);

    MultiFab* tmp_xvel_new; tmp_xvel_new->define(convert(ba, IntVect(1,0,0)), dm, 1, ngrow_vels);
    MultiFab* tmp_xvel_old; tmp_xvel_old->define(convert(ba, IntVect(1,0,0)), dm, 1, ngrow_vels);

    MultiFab* tmp_yvel_new; tmp_yvel_new->define(convert(ba, IntVect(0,1,0)), dm, 1, ngrow_vels);
    MultiFab* tmp_yvel_old; tmp_yvel_old->define(convert(ba, IntVect(0,1,0)), dm, 1, ngrow_vels);

    MultiFab* tmp_zvel_new; tmp_zvel_new->define(convert(ba, IntVect(0,0,1)), dm, 1, IntVect(ngrow_vels,ngrow_vels,0));
    MultiFab* tmp_zvel_old; tmp_zvel_old->define(convert(ba, IntVect(0,0,1)), dm, 1, IntVect(ngrow_vels,ngrow_vels,0));

    // This will fill the temporary MultiFabs with data from previous fine data as well as coarse where needed
    FillPatch(lev, time, tmp_cons_new, cons_new);
    FillPatch(lev, time, tmp_xvel_new, xvel_new);
    FillPatch(lev, time, tmp_yvel_new, yvel_new);
    FillPatch(lev, time, tmp_zvel_new, zvel_new);

    MultiFab::Copy(*tmp_cons_old,*tmp_cons_new,0,0,NCONS,tmp_cons_new[lev].nGrowVect());
    MultiFab::Copy(*tmp_xvel_old,*tmp_xvel_new,0,0,   1,tmp_xvel_new[lev].nGrowVect());
    MultiFab::Copy(*tmp_yvel_old,*tmp_yvel_new,0,0,   1,tmp_yvel_new[lev].nGrowVect());
    MultiFab::Copy(*tmp_zvel_old,*tmp_zvel_new,0,0,   1,tmp_zvel_new[lev].nGrowVect());

    std::swap(tmp_cons_new, cons_new[lev]);
    std::swap(tmp_cons_old, cons_old[lev]);
    std::swap(tmp_xvel_new, xvel_new[lev]);
    std::swap(tmp_xvel_old, xvel_old[lev]);
    std::swap(tmp_yvel_new, yvel_new[lev]);
    std::swap(tmp_yvel_old, yvel_old[lev]);
    std::swap(tmp_zvel_new, zvel_new[lev]);
    std::swap(tmp_zvel_old, zvel_old[lev]);

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    init_stuff(lev, ba, dm);
}

// Make a new level from scratch using provided BoxArray and DistributionMapping.
// This is called both for initialization and for restart
// (overrides the pure virtual function in AmrCore)
// main.cpp --> ROMSX::InitData --> InitFromScratch --> MakeNewGrids --> MakeNewLevelFromScratch
//                                         restart  --> MakeNewGrids --> MakeNewLevelFromScratch
void ROMSX::MakeNewLevelFromScratch (int lev, Real /*time*/, const BoxArray& ba,
                                     const DistributionMapping& dm)
{
    // Set BoxArray grids and DistributionMapping dmap in AMReX_AmrMesh.H class
    SetBoxArray(lev, ba);
    SetDistributionMap(lev, dm);

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
      init_stuff(lev, ba, dm);
}

void ROMSX::resize_stuff(int lev)
{
    vec_z_phys_nd.resize(lev+1);

    vec_hOfTheConfusingName.resize(lev+1);
    vec_Zt_avg1.resize(lev+1);
    vec_s_r.resize(lev+1);
    vec_z_w.resize(lev+1);
    vec_z_r.resize(lev+1);
    vec_y_r.resize(lev+1);
    vec_x_r.resize(lev+1);
    vec_Hz.resize(lev+1);
    vec_Huon.resize(lev+1);
    vec_Hvom.resize(lev+1);
    vec_Akv.resize(lev+1);
    vec_Akt.resize(lev+1);
    vec_visc3d_r.resize(lev+1);
    vec_visc2_p.resize(lev+1);
    vec_visc2_r.resize(lev+1);
    vec_diff2.resize(lev+1);
    vec_ru.resize(lev+1);
    vec_rv.resize(lev+1);
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
    vec_t3.resize(lev+1);
    vec_s3.resize(lev+1);

    vec_rhoS.resize(lev+1);
    vec_rhoA.resize(lev+1);

    mapfac_m.resize(lev+1);
    mapfac_u.resize(lev+1);
    mapfac_v.resize(lev+1);
}

void ROMSX::init_stuff(int lev, const BoxArray& ba, const DistributionMapping& dm)
{
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

    vec_z_phys_nd[lev].reset          (new MultiFab(ba_nd,dm,1,IntVect(NGROW,NGROW,0))); // z at psi points (nodes) MIGHT NEED NGROW+1

    vec_hOfTheConfusingName[lev].reset(new MultiFab(ba2d ,dm,2,IntVect(NGROW+1,NGROW+1,0))); //2d, depth (double check if negative)
    vec_Zt_avg1[lev].reset            (new MultiFab(ba2d ,dm,1,IntVect(NGROW+1,NGROW+1,0))); //2d, average of the free surface (zeta)

    vec_x_r[lev].reset                (new MultiFab(ba2d,dm,1,IntVect(NGROW+1,NGROW+1,0))); // x at r points (cell center)
    vec_y_r[lev].reset                (new MultiFab(ba2d,dm,1,IntVect(NGROW+1,NGROW+1,0))); // y at r points (cell center)

    vec_s_r[lev].reset                (new MultiFab(ba1d,dm,1,IntVect(    0,    0,0))); // scaled vertical coordinate [0,1] , transforms to z

    vec_z_w[lev].reset                (new MultiFab(ba_w,dm,1,IntVect(NGROW+1,NGROW+1,0))); // z at w points (cell faces)
    vec_z_r[lev].reset                (new MultiFab(ba  ,dm,1,IntVect(NGROW+1,NGROW+1,0))); // z at r points (cell center)
    vec_Hz[lev].reset                 (new MultiFab(ba  ,dm,1,IntVect(NGROW+1,NGROW+1,NGROW+1))); // like in ROMS, thickness of cell in z

    vec_Huon[lev].reset               (new MultiFab(convert(ba,IntVect(1,0,0)),dm,1,IntVect(NGROW,NGROW,0))); // mass flux for u component
    vec_Hvom[lev].reset               (new MultiFab(convert(ba,IntVect(0,1,0)),dm,1,IntVect(NGROW,NGROW,0))); // mass flux for v component

    vec_Akv[lev].reset                (new MultiFab(ba  ,dm,1,IntVect(NGROW,NGROW,0))); // vertical mixing coefficient (.in)
    vec_Akt[lev].reset                (new MultiFab(ba  ,dm,1,IntVect(NGROW,NGROW,0))); // vertical mixing coefficient (.in)
    vec_visc3d_r[lev].reset           (new MultiFab(ba  ,dm,1,IntVect(NGROW,NGROW,0))); // not used


    // check dimensionality
    vec_visc2_p[lev].reset(new MultiFab(ba,dm,1,IntVect(NGROW,NGROW,0))); // harmonic viscosity at psi points -- difference to 3d?
    vec_visc2_r[lev].reset(new MultiFab(ba,dm,1,IntVect(NGROW,NGROW,0))); // harmonic viscosity at rho points
    vec_diff2[lev].reset(new MultiFab(ba,dm,2,IntVect(NGROW,NGROW,0))); // harmonic diffusivity temperature/salt

    // maybe TODO: clean up component indexing in prestep?
    vec_ru[lev].reset(new MultiFab(convert(ba,IntVect(1,0,0)),dm,2,IntVect(NGROW,NGROW,NGROW))); // RHS u (incl horizontal and vertical advection)
    vec_rv[lev].reset(new MultiFab(convert(ba,IntVect(0,1,0)),dm,2,IntVect(NGROW,NGROW,NGROW))); // RHS v
    vec_rufrc[lev].reset(new MultiFab(convert(ba2d,IntVect(1,0,0)),dm,2,IntVect(NGROW,NGROW,0))); //2d, (incl advection terms and surface/bottom stresses, integral over the whole columnn, k=0)
    vec_rvfrc[lev].reset(new MultiFab(convert(ba2d,IntVect(0,1,0)),dm,2,IntVect(NGROW,NGROW,0))); //2d, same as above but v
    vec_sustr[lev].reset(new MultiFab(convert(ba2d,IntVect(1,0,0)),dm,1,IntVect(NGROW,NGROW,0))); //2d, surface stress
    vec_svstr[lev].reset(new MultiFab(convert(ba2d,IntVect(0,1,0)),dm,1,IntVect(NGROW,NGROW,0))); //2d
    vec_rdrag[lev].reset(new MultiFab(ba2d,dm,1,IntVect(NGROW,NGROW,0))); //2d, linear drag coefficient [m/s], defined at rho, somehow related to rdrg
    vec_bustr[lev].reset(new MultiFab(convert(ba2d,IntVect(1,0,0)),dm,1,IntVect(NGROW,NGROW,0))); //2d, bottom stress
    vec_bvstr[lev].reset(new MultiFab(convert(ba2d,IntVect(0,1,0)),dm,1,IntVect(NGROW,NGROW,0)));

    //all 2d -- all associated with the 2D advance
    vec_DU_avg1[lev].reset(new MultiFab(convert(ba2d,IntVect(1,0,0)),dm,1,IntVect(NGROW,NGROW,0))); //2d DU: sum(height[incl free surface?] * u)
    vec_DU_avg2[lev].reset(new MultiFab(convert(ba2d,IntVect(1,0,0)),dm,1,IntVect(NGROW,NGROW,0))); //2d like above, but correct(or)?
    vec_DV_avg1[lev].reset(new MultiFab(convert(ba2d,IntVect(0,1,0)),dm,1,IntVect(NGROW,NGROW,0)));
    vec_DV_avg2[lev].reset(new MultiFab(convert(ba2d,IntVect(0,1,0)),dm,1,IntVect(NGROW,NGROW,0)));

    vec_rubar[lev].reset(new MultiFab(convert(ba2d,IntVect(1,0,0)),dm,4,IntVect(NGROW,NGROW,0))); // 2d RHS ubar
    vec_rvbar[lev].reset(new MultiFab(convert(ba2d,IntVect(0,1,0)),dm,4,IntVect(NGROW,NGROW,0)));
    vec_rzeta[lev].reset(new MultiFab(ba2d,dm,4,IntVect(NGROW,NGROW,0))); // 2d RHS zeta

    // starts off kind of like a depth-averaged u, but exists at more points and more timesteps (b/c fast 2D update) than full u
    vec_ubar[lev].reset(new MultiFab(convert(ba2d,IntVect(1,0,0)),dm,3,IntVect(NGROW,NGROW,0)));
    vec_vbar[lev].reset(new MultiFab(convert(ba2d,IntVect(0,1,0)),dm,3,IntVect(NGROW,NGROW,0)));
    vec_zeta[lev].reset(new MultiFab(ba2d,dm,3,IntVect(NGROW+1,NGROW+1,0)));  // 2d free surface

    vec_t3[lev].reset(new MultiFab(ba,dm,1,IntVect(NGROW,NGROW,0))); //tempstore
    vec_s3[lev].reset(new MultiFab(ba,dm,1,IntVect(NGROW,NGROW,0))); //saltstore

    vec_rhoS[lev].reset(new MultiFab(ba,dm,1,IntVect(NGROW,NGROW,0)));
    vec_rhoA[lev].reset(new MultiFab(ba,dm,1,IntVect(NGROW,NGROW,0)));

    set_bathymetry(lev);
    stretch_transform(lev);
    set_vmix(lev);
    set_hmixcoef(lev);
    set_weights(lev);

    //consider tracking ru and rv indexes more specifically or more similarly to indx
    vec_ru[lev]->setVal(0.0);
    vec_rv[lev]->setVal(0.0);

    vec_DU_avg1[lev]->setVal(0.0);
    vec_DU_avg2[lev]->setVal(0.0);
    vec_DV_avg1[lev]->setVal(0.0);
    vec_DV_avg2[lev]->setVal(0.0);
    vec_rubar[lev]->setVal(0.0);
    vec_rvbar[lev]->setVal(0.0);
    vec_rzeta[lev]->setVal(0.0);

    vec_ubar[lev]->setVal(0.0);
    vec_vbar[lev]->setVal(0.0);
    vec_zeta[lev]->setVal(0.0);

    // Set initial linear drag coefficient
    vec_rdrag[lev]->setVal(solverChoice.rdrag);
}

// Delete level data
// overrides the pure virtual function in AmrCore
void
ROMSX::ClearLevel (int lev)
{
    delete cons_new[lev]; delete xvel_new[lev];  delete yvel_new[lev];  delete zvel_new[lev];
    delete cons_old[lev]; delete xvel_old[lev];  delete yvel_old[lev];  delete zvel_old[lev];
}
