/**
 * \file ROMSX.cpp
 */

#include <prob_common.H>
#include <EOS.H>
#include <ROMSX.H>

#include <AMReX_buildInfo.H>

#include <Utils.H>
#include <TerrainMetrics.H>

#ifdef ROMSX_USE_MULTIBLOCK
#include <MultiBlockContainer.H>
#endif

using namespace amrex;

amrex::Real ROMSX::startCPUTime        = 0.0;
amrex::Real ROMSX::previousCPUTimeUsed = 0.0;

Vector<AMRErrorTag> ROMSX::ref_tags;

SolverChoice ROMSX::solverChoice;

// Time step control
amrex::Real ROMSX::cfl           =  0.8;
amrex::Real ROMSX::fixed_dt      = -1.0;
amrex::Real ROMSX::init_shrink   =  1.0;
amrex::Real ROMSX::change_max    =  1.1;

// Type of mesh refinement algorithm
int         ROMSX::do_reflux     = 0;
int         ROMSX::do_avg_down   = 0;

// Dictate verbosity in screen output
int         ROMSX::verbose       = 0;

// Frequency of diagnostic output
int         ROMSX::sum_interval  = -1;
amrex::Real ROMSX::sum_per       = -1.0;

// Native AMReX vs NetCDF
std::string ROMSX::plotfile_type    = "amrex";

// init_type:  "custom", "ideal", "real"
std::string ROMSX::init_type        = "custom";

// NetCDF wrfinput (initialization) file(s)
amrex::Vector<amrex::Vector<std::string>> ROMSX::nc_init_file = {{""}}; // Must provide via input
std::string ROMSX::nc_bdy_file;

amrex::Vector<std::string> BCNames = {"xlo", "ylo", "zlo", "xhi", "yhi", "zhi"};

// constructor - reads in parameters from inputs file
//             - sizes multilevel arrays and data structures
//             - initializes BCRec boundary condition object
ROMSX::ROMSX ()
{
    if (amrex::ParallelDescriptor::IOProcessor()) {
        const char* romsx_hash = amrex::buildInfoGetGitHash(1);
        const char* amrex_hash = amrex::buildInfoGetGitHash(2);
        const char* buildgithash = amrex::buildInfoGetBuildGitHash();
        const char* buildgitname = amrex::buildInfoGetBuildGitName();

        if (strlen(romsx_hash) > 0) {
          amrex::Print() << "\n"
                         << "ROMSX git hash: " << romsx_hash << "\n";
        }
        if (strlen(amrex_hash) > 0) {
          amrex::Print() << "AMReX git hash: " << amrex_hash << "\n";
        }
        if (strlen(buildgithash) > 0) {
          amrex::Print() << buildgitname << " git hash: " << buildgithash << "\n";
        }

        amrex::Print() << "\n";
    }

    ReadParameters();
    const std::string& pv1 = "plot_vars_1"; setPlotVariables(pv1,plot_var_names_1);
    const std::string& pv2 = "plot_vars_2"; setPlotVariables(pv2,plot_var_names_2);

    amrex_probinit(geom[0].ProbLo(),geom[0].ProbHi());

    // Geometry on all levels has been defined already.

    // No valid BoxArray and DistributionMapping have been defined.
    // But the arrays for them have been resized.

    int nlevs_max = max_level + 1;

    istep.resize(nlevs_max, 0);
    nsubsteps.resize(nlevs_max, 1);
    for (int lev = 1; lev <= max_level; ++lev) {
        nsubsteps[lev] = MaxRefRatio(lev-1);
    }

    t_new.resize(nlevs_max, 0.0);
    t_old.resize(nlevs_max, -1.e100);
    dt.resize(nlevs_max, 1.e100);

    vars_new.resize(nlevs_max);
    vars_old.resize(nlevs_max);

    for (int lev = 0; lev < nlevs_max; ++lev) {
        vars_new[lev].resize(Vars::NumTypes);
        vars_old[lev].resize(Vars::NumTypes);
    }

    flux_registers.resize(nlevs_max);

    // Initialize tagging criteria for mesh refinement
    refinement_criteria_setup();

    // We have already read in the ref_Ratio (via amr.ref_ratio =) but we need to enforce
    //     that there is no refinement in the vertical so we test on that here.
    for (int lev = 0; lev < max_level; ++lev)
    {
       amrex::Print() << "Refinement ratio at level " << lev << " set to be " <<
          ref_ratio[lev][0]  << " " << ref_ratio[lev][1]  <<  " " << ref_ratio[lev][2] << std::endl;

       if (ref_ratio[lev][2] != 1)
       {
           amrex::Error("We don't allow refinement in the vertical -- make sure to set ref_ratio = 1 in z");
       }
    }
}

ROMSX::~ROMSX ()
{
}

// advance solution to final time
void
ROMSX::Evolve ()
{
    Real cur_time = t_new[0];

    // Take one coarse timestep by calling timeStep -- which recursively calls timeStep
    //      for finer levels (with or without subcycling)
    for (int step = istep[0]; step < max_step && cur_time < stop_time; ++step)
    {
        amrex::Print() << "\nCoarse STEP " << step+1 << " starts ..." << std::endl;

        ComputeDt();

        int lev = 0;
        int iteration = 1;
        timeStep(lev, cur_time, iteration);

        cur_time  += dt[0];

        amrex::Print() << "Coarse STEP " << step+1 << " ends." << " TIME = " << cur_time
                       << " DT = " << dt[0]  << std::endl;
        print_state(vars_new[lev][Vars::xvel],IntVect(AMREX_D_DECL(2,2,2)));
        if (plot_int_1 > 0 && (step+1) % plot_int_1 == 0) {
            last_plot_file_step_1 = step+1;
            WritePlotFile(1,plot_var_names_1);
        }
        if (plot_int_2 > 0 && (step+1) % plot_int_2 == 0) {
            last_plot_file_step_2 = step+1;
            WritePlotFile(2,plot_var_names_2);
        }

        if (check_int > 0 && (step+1) % check_int == 0) {
            last_check_file_step = step+1;
#ifdef ROMSX_USE_NETCDF
            if (check_type == "netcdf") {
               WriteNCCheckpointFile();
            }
#endif
            if (check_type == "native") {
               WriteCheckpointFile();
            }
        }

        post_timestep(step, cur_time, dt[0]);

#ifdef AMREX_MEM_PROFILING
        {
            std::ostringstream ss;
            ss << "[STEP " << step+1 << "]";
            MemProfiler::report(ss.str());
        }
#endif

        if (cur_time >= stop_time - 1.e-6*dt[0]) break;
    }

    if (plot_int_1 > 0 && istep[0] > last_plot_file_step_1) {
        WritePlotFile(1,plot_var_names_1);
    }
    if (plot_int_2 > 0 && istep[0] > last_plot_file_step_2) {
        WritePlotFile(2,plot_var_names_2);
    }

    if (check_int > 0 && istep[0] > last_check_file_step) {
#ifdef ROMSX_USE_NETCDF
        if (check_type == "netcdf") {
           WriteNCCheckpointFile();
        }
#endif
        if (check_type == "native") {
           WriteCheckpointFile();
        }
    }

}

// Called after every coarse timestep
void
ROMSX::post_timestep (int nstep, Real time, Real dt_lev0)
{
    BL_PROFILE("ROMSX::post_timestep()");

    if (do_reflux)
    {
        for (int lev = finest_level-1; lev >= 0; lev--)
        {
            // This call refluxes from the lev/lev+1 interface onto lev
            get_flux_reg(lev+1).Reflux(vars_new[lev][Vars::cons],1.0, 0, 0, NVAR, geom[lev]);

            // We need to do this before anything else because refluxing changes the
            // values of coarse cells underneath fine grids with the assumption they'll
            // be over-written by averaging down
            AverageDownTo(lev);
        }
    }

    if (is_it_time_for_action(nstep, time, dt_lev0, sum_interval, sum_per)) {
        sum_integrated_quantities(time);
    }
}

// This is called from main.cpp and handles all initialization, whether from start or restart
void
ROMSX::InitData ()
{
    // Initialize the start time for our CPU-time tracker
    startCPUTime = amrex::ParallelDescriptor::second();

    // Map the words in the inputs file to BC types, then translate
    //     those types into what they mean for each variable
    init_bcs();

    last_plot_file_step_1 = -1;
    last_plot_file_step_2 = -1;
    last_check_file_step = -1;

    if (restart_chkfile == "") {
        // start simulation from the beginning

        const Real time = 0.0;
        InitFromScratch(time);

#ifdef ROMSX_USE_MULTIBLOCK
        // Multiblock: hook to set BL & comms once ba/dm are known
        if(domain_p[0].bigEnd(0) < 500 ) {
            m_mbc->SetBoxLists();
            m_mbc->SetBlockCommMetaData();
        }
#endif

        // For now we initialize rho_KE to 0
        for (int lev = finest_level-1; lev >= 0; --lev)
            vars_new[lev][Vars::cons].setVal(0.0,RhoKE_comp,1,0);

        if (solverChoice.use_terrain) {
            if (init_type != "real") {
                for (int lev = 0; lev <= finest_level; lev++)
                {
                    init_custom_terrain(geom[lev],*z_phys_nd[lev],time);
                    init_terrain_grid(geom[lev],*z_phys_nd[lev]);
                    make_metrics(geom[lev],*z_phys_nd[lev],*z_phys_cc[lev],*detJ_cc[lev]);
                }
            }
        }

        for (int lev = 0; lev <= finest_level; lev++)
            init_only(lev, time);

        AverageDown();

        if(solverChoice.use_terrain) {
            if (init_type == "real") {
                for (int lev = 0; lev <= finest_level; lev++)
                {
                    make_metrics(geom[lev],*z_phys_nd[lev],*z_phys_cc[lev],*detJ_cc[lev]);
                }
            }
        }

    } else { // Restart from a checkpoint

        restart();

        if(solverChoice.use_terrain) {
            // This must come after the call to restart because that
            //      is where we read in the mesh data
            for (int lev = finest_level-1; lev >= 0; --lev)
                make_metrics(geom[lev],*z_phys_nd[lev],*z_phys_cc[lev],*detJ_cc[lev]);
        }
    }

    // Initialize flux registers (whether we start from scratch or restart)
    if (do_reflux) {
        flux_registers[0] = 0;
        for (int lev = 1; lev <= finest_level; lev++)
        {
            flux_registers[lev] = new FluxRegister(grids[lev], dmap[lev], ref_ratio[lev-1], lev, NVAR);
        }
    }

    if (restart_chkfile == "" && check_int > 0)
    {
#ifdef ROMSX_USE_NETCDF
        if (check_type == "netcdf") {
           WriteNCCheckpointFile();
        }
#endif
        if (check_type == "native") {
           WriteCheckpointFile();
        }
        last_check_file_step = 0;
    }

    if ( (restart_chkfile == "") ||
         (restart_chkfile != "" && plot_file_on_restart) )
    {
        if (plot_int_1 > 0)
        {
            WritePlotFile(1,plot_var_names_1);
            last_plot_file_step_1 = istep[0];
        }
        if (plot_int_2 > 0)
        {
            WritePlotFile(2,plot_var_names_2);
            last_plot_file_step_2 = istep[0];
        }
    }

    if (is_it_time_for_action(istep[0], t_new[0], dt[0], sum_interval, sum_per)) {
        sum_integrated_quantities(t_new[0]);
    }

    ComputeDt();

    // Fill ghost cells/faces
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        auto& lev_new = vars_new[lev];
        auto& lev_old = vars_old[lev];

        FillPatch(lev, t_new[lev], lev_new);

        // Copy from new into old just in case
        int ngs   = lev_new[Vars::cons].nGrow();
        int ngvel = lev_new[Vars::xvel].nGrow();
        MultiFab::Copy(lev_old[Vars::cons],lev_new[Vars::cons],0,0,NVAR,ngs);
        MultiFab::Copy(lev_old[Vars::xvel],lev_new[Vars::xvel],0,0,1,ngvel);
        MultiFab::Copy(lev_old[Vars::yvel],lev_new[Vars::yvel],0,0,1,ngvel);
        MultiFab::Copy(lev_old[Vars::zvel],lev_new[Vars::zvel],0,0,1,IntVect(ngvel,ngvel,0));
    }
}

void
ROMSX::restart()
{
#ifdef ROMSX_USE_NETCDF
    if (restart_type == "netcdf") {
       ReadNCCheckpointFile();
    }
#endif
    if (restart_type == "native") {
       ReadCheckpointFile();
    }

    // We set this here so that we don't over-write the checkpoint file we just started from
    last_check_file_step = istep[0];
}

// Make a new level using provided BoxArray and DistributionMapping and
// fill with interpolated coarse level data.
// overrides the pure virtual function in AmrCore
void
ROMSX::MakeNewLevelFromCoarse (int lev, Real time, const BoxArray& ba,
                             const DistributionMapping& dm)
{
    const auto& crse_new = vars_new[lev-1];
    auto& lev_new = vars_new[lev];
    auto& lev_old = vars_old[lev];

    lev_new[Vars::cons].define(ba, dm, crse_new[Vars::cons].nComp(), crse_new[Vars::cons].nGrowVect());
    lev_old[Vars::cons].define(ba, dm, crse_new[Vars::cons].nComp(), crse_new[Vars::cons].nGrowVect());

    lev_new[Vars::xvel].define(convert(ba, IntVect(1,0,0)), dm, 1, crse_new[Vars::xvel].nGrowVect());
    lev_old[Vars::xvel].define(convert(ba, IntVect(1,0,0)), dm, 1, crse_new[Vars::xvel].nGrowVect());

    lev_new[Vars::yvel].define(convert(ba, IntVect(0,1,0)), dm, 1, crse_new[Vars::yvel].nGrowVect());
    lev_old[Vars::yvel].define(convert(ba, IntVect(0,1,0)), dm, 1, crse_new[Vars::yvel].nGrowVect());

    lev_new[Vars::zvel].define(convert(ba, IntVect(0,0,1)), dm, 1, crse_new[Vars::zvel].nGrowVect());
    lev_old[Vars::zvel].define(convert(ba, IntVect(0,0,1)), dm, 1, crse_new[Vars::zvel].nGrowVect());

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    FillCoarsePatchAllVars(lev, time, vars_new[lev]);
}

// Remake an existing level using provided BoxArray and DistributionMapping and
// fill with existing fine and coarse data.
// overrides the pure virtual function in AmrCore
void
ROMSX::RemakeLevel (int lev, Real time, const BoxArray& ba, const DistributionMapping& dm)
{
    Vector<MultiFab> tmp_lev_new(Vars::NumTypes);
    Vector<MultiFab> tmp_lev_old(Vars::NumTypes);

    int ngrow_state = ComputeGhostCells(solverChoice.spatial_order)+1;
    int ngrow_vels  = ComputeGhostCells(solverChoice.spatial_order);

    tmp_lev_new[Vars::cons].define(ba, dm, Cons::NumVars, ngrow_state);
    tmp_lev_old[Vars::cons].define(ba, dm, Cons::NumVars, ngrow_state);

    tmp_lev_new[Vars::xvel].define(convert(ba, IntVect(1,0,0)), dm, 1, ngrow_vels);
    tmp_lev_old[Vars::xvel].define(convert(ba, IntVect(1,0,0)), dm, 1, ngrow_vels);

    tmp_lev_new[Vars::yvel].define(convert(ba, IntVect(0,1,0)), dm, 1, ngrow_vels);
    tmp_lev_old[Vars::yvel].define(convert(ba, IntVect(0,1,0)), dm, 1, ngrow_vels);

    tmp_lev_new[Vars::zvel].define(convert(ba, IntVect(0,0,1)), dm, 1, IntVect(ngrow_vels,ngrow_vels,0));
    tmp_lev_old[Vars::zvel].define(convert(ba, IntVect(0,0,1)), dm, 1, IntVect(ngrow_vels,ngrow_vels,0));

    // This will fill the temporary MultiFabs with data from vars_new
    FillPatch(lev, time, tmp_lev_new);
    FillPatch(lev, time, tmp_lev_old);

    for (int var_idx = 0; var_idx < Vars::NumTypes; ++var_idx) {
        std::swap(tmp_lev_new[var_idx], vars_new[lev][var_idx]);
        std::swap(tmp_lev_old[var_idx], vars_old[lev][var_idx]);
    }

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;
}

// Delete level data
// overrides the pure virtual function in AmrCore
void
ROMSX::ClearLevel (int lev)
{
    for (int var_idx = 0; var_idx < Vars::NumTypes; ++var_idx) {
        vars_new[lev][var_idx].clear();
        vars_old[lev][var_idx].clear();
    }
}

// Make a new level from scratch using provided BoxArray and DistributionMapping.
// This is called both for initialization and for restart
// (overrides the pure virtual function in AmrCore)
// main.cpp --> ROMSX::InitData --> InitFromScratch --> MakeNewGrids --> MakeNewLevelFromScratch
//                                       restart  --> MakeNewGrids --> MakeNewLevelFromScratch
void ROMSX::MakeNewLevelFromScratch (int lev, Real /*time*/, const BoxArray& ba,
                                   const DistributionMapping& dm)
{
    // Set BoxArray grids and DistributionMapping dmap in AMReX_AmrMesh.H class
    SetBoxArray(lev, ba);
    SetDistributionMap(lev, dm);

    // The number of ghost cells for density must be 1 greater than that for velocity
    //     so that we can go back in forth between velocity and momentum on all faces
    int ngrow_state = ComputeGhostCells(solverChoice.spatial_order)+1;
    int ngrow_vels  = ComputeGhostCells(solverChoice.spatial_order);

    auto& lev_new = vars_new[lev];
    auto& lev_old = vars_old[lev];

    lev_new[Vars::cons].define(ba, dm, Cons::NumVars, ngrow_state);
    lev_old[Vars::cons].define(ba, dm, Cons::NumVars, ngrow_state);

    lev_new[Vars::xvel].define(convert(ba, IntVect(1,0,0)), dm, 1, ngrow_vels);
    lev_old[Vars::xvel].define(convert(ba, IntVect(1,0,0)), dm, 1, ngrow_vels);

    lev_new[Vars::yvel].define(convert(ba, IntVect(0,1,0)), dm, 1, ngrow_vels);
    lev_old[Vars::yvel].define(convert(ba, IntVect(0,1,0)), dm, 1, ngrow_vels);

    lev_new[Vars::zvel].define(convert(ba, IntVect(0,0,1)), dm, 1, IntVect(ngrow_vels,ngrow_vels,0));
    lev_old[Vars::zvel].define(convert(ba, IntVect(0,0,1)), dm, 1, IntVect(ngrow_vels,ngrow_vels,0));

    z_phys_nd.resize(lev+1);
    z_phys_cc.resize(lev+1);
    detJ_cc.resize(lev+1);

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

    mapfac_m.resize(lev+1);
    mapfac_u.resize(lev+1);
    mapfac_v.resize(lev+1);
    mapfac_m[lev].reset(new MultiFab(ba2d,dm,1,0));
    mapfac_u[lev].reset(new MultiFab(convert(ba2d,IntVect(1,0,0)),dm,1,0));
    mapfac_v[lev].reset(new MultiFab(convert(ba2d,IntVect(0,1,0)),dm,1,0));
    mapfac_m[lev]->setVal(1.);
    mapfac_u[lev]->setVal(1.);
    mapfac_v[lev]->setVal(1.);

    // Base state holds r_0, pres_0, pi_0 (in that order)
    base_state.resize(lev+1);
    base_state[lev].define(ba,dm,3,1);
    base_state[lev].setVal(0.);

    if (solverChoice.use_terrain) {
        z_phys_cc[lev].reset(new MultiFab(ba,dm,1,1));
          detJ_cc[lev].reset(new MultiFab(ba,dm,1,1));

        BoxArray ba_nd(ba);
        ba_nd.surroundingNodes();

        // We need this to be one greater than the ghost cells to handle levels > 0
        int ngrow = ComputeGhostCells(solverChoice.spatial_order)+2;
        z_phys_nd[lev].reset(new MultiFab(ba_nd,dm,1,IntVect(ngrow,ngrow,1)));

    } else {
        z_phys_nd[lev] = nullptr;
        z_phys_cc[lev] = nullptr;
          detJ_cc[lev] = nullptr;
    }

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
    vec_visc3d_r.resize(lev+1);
    vec_visc2_p.resize(lev+1);
    vec_visc2_r.resize(lev+1);
    vec_ru.resize(lev+1);
    vec_rv.resize(lev+1);
    vec_rufrc.resize(lev+1);
    vec_rvfrc.resize(lev+1);
    vec_sustr.resize(lev+1);
    vec_svstr.resize(lev+1);

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

    vec_hOfTheConfusingName[lev].reset(new MultiFab(ba,dm,1,IntVect(2,2,0)));
    vec_Zt_avg1[lev].reset(new MultiFab(ba,dm,1,IntVect(2,2,0)));
    vec_s_r[lev].reset(new MultiFab(ba1d,dm,1,IntVect(0,0,0)));
    vec_z_w[lev].reset(new MultiFab(ba,dm,1,IntVect(2,2,0)));
    vec_z_r[lev].reset(new MultiFab(ba,dm,1,IntVect(2,2,0)));
    vec_y_r[lev].reset(new MultiFab(ba2d,dm,1,IntVect(2,2,0)));
    vec_x_r[lev].reset(new MultiFab(ba2d,dm,1,IntVect(2,2,0)));
    vec_Hz[lev].reset(new MultiFab(ba,dm,1,IntVect(3,3,3)));
    vec_Huon[lev].reset(new MultiFab(ba,dm,1,IntVect(2,2,0)));
    vec_Hvom[lev].reset(new MultiFab(ba,dm,1,IntVect(2,2,0)));
    vec_Akv[lev].reset(new MultiFab(ba,dm,1,IntVect(2,2,0)));
    vec_visc3d_r[lev].reset(new MultiFab(ba,dm,1,IntVect(2,2,0)));
    vec_visc2_p[lev].reset(new MultiFab(ba,dm,1,IntVect(2,2,0)));
    vec_visc2_r[lev].reset(new MultiFab(ba,dm,1,IntVect(2,2,0)));
    vec_ru[lev].reset(new MultiFab(ba,dm,2,IntVect(2,2,0)));
    vec_rv[lev].reset(new MultiFab(ba,dm,2,IntVect(2,2,0)));
    vec_rufrc[lev].reset(new MultiFab(ba,dm,2,IntVect(2,2,0)));
    vec_rvfrc[lev].reset(new MultiFab(ba,dm,2,IntVect(2,2,0)));
    vec_sustr[lev].reset(new MultiFab(ba,dm,1,IntVect(2,2,0)));
    vec_svstr[lev].reset(new MultiFab(ba,dm,1,IntVect(2,2,0)));

    vec_DU_avg1[lev].reset(new MultiFab(ba2d,dm,1,IntVect(2,2,0)));
    vec_DU_avg2[lev].reset(new MultiFab(ba2d,dm,1,IntVect(2,2,0)));
    vec_DV_avg1[lev].reset(new MultiFab(ba2d,dm,1,IntVect(2,2,0)));
    vec_DV_avg2[lev].reset(new MultiFab(ba2d,dm,1,IntVect(2,2,0)));
    vec_rubar[lev].reset(new MultiFab(ba,dm,4,IntVect(2,2,0)));
    vec_rvbar[lev].reset(new MultiFab(ba,dm,4,IntVect(2,2,0)));
    vec_rzeta[lev].reset(new MultiFab(ba,dm,4,IntVect(2,2,0)));
    vec_ubar[lev].reset(new MultiFab(ba2d,dm,3,IntVect(2,2,0)));
    vec_vbar[lev].reset(new MultiFab(ba2d,dm,3,IntVect(2,2,0)));
    vec_zeta[lev].reset(new MultiFab(ba,dm,3,IntVect(2,2,0)));
    vec_t3[lev].reset(new MultiFab(ba,dm,1,IntVect(2,2,0)));
    vec_s3[lev].reset(new MultiFab(ba,dm,1,IntVect(2,2,0)));

    set_depth(lev);
    set_vmix(lev);
    set_hmixcoef(lev);

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

}

void
ROMSX::init_only(int lev, Real time)
{
    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    auto& lev_new = vars_new[lev];

    // Loop over grids at this level to initialize our grid data
    lev_new[Vars::cons].setVal(0.0);
    lev_new[Vars::xvel].setVal(0.0);
    lev_new[Vars::yvel].setVal(0.0);
    lev_new[Vars::zvel].setVal(0.0);

    if (init_type == "custom") {
        init_custom(lev);
#ifdef ROMSX_USE_NETCDF
    } else if (init_type == "ideal" || init_type == "real") {
        init_from_wrfinput(lev);
#endif
    }

    // Ensure that the face-based data are the same on both sides of a periodic domain.
    // The data associated with the lower grid ID is considered the correct value.
    lev_new[Vars::xvel].OverrideSync(geom[lev].periodicity());
    lev_new[Vars::yvel].OverrideSync(geom[lev].periodicity());
    lev_new[Vars::zvel].OverrideSync(geom[lev].periodicity());
}

// read in some parameters from inputs file
void
ROMSX::ReadParameters ()
{
    {
        ParmParse pp;  // Traditionally, max_step and stop_time do not have prefix.
        pp.query("max_step", max_step);
        pp.query("stop_time", stop_time);
    }

    ParmParse pp(pp_prefix);
    {
        // The type of the file we restart from
        pp.query("restart_type", restart_type);

        pp.query("regrid_int", regrid_int);
        pp.query("check_file", check_file);
        pp.query("check_type", check_type);
        pp.query("check_int", check_int);

        pp.query("restart", restart_chkfile);

        if (pp.contains("data_log"))
        {
            int num_datalogs = pp.countval("data_log");
            datalog.resize(num_datalogs);
            datalogname.resize(num_datalogs);
            pp.queryarr("data_log",datalogname,0,num_datalogs);
            for (int i = 0; i < num_datalogs; i++)
                setRecordDataInfo(i,datalogname[i]);
        }

        // Verbosity
        pp.query("v", verbose);

        // Frequency of diagnostic output
        pp.query("sum_interval", sum_interval);
        pp.query("sum_period"  , sum_per);

        // Time step controls
        pp.query("cfl", cfl);
        pp.query("init_shrink", init_shrink);
        pp.query("change_max", change_max);

        pp.query("fixed_dt", fixed_dt);

        AMREX_ASSERT(cfl > 0. || fixed_dt > 0.);

        // Mesh refinement
        do_reflux = 0;
        do_avg_down = 0;

        // How to initialize
        pp.query("init_type",init_type);
        if (init_type != "custom" &&
            init_type != "ideal" &&
            init_type != "real")
        {
            amrex::Error("init_type must be custom, ideal, or real ");
        }

        // We use this to keep track of how many boxes we read in from WRF initialization
        num_files_at_level.resize(max_level+1,0);

        // We use this to keep track of how many boxes are specified thru the refinement indicators
        num_boxes_at_level.resize(max_level+1,0);
            boxes_at_level.resize(max_level+1);

        // We always have exactly one file at level 0
        num_boxes_at_level[0] = 1;

#ifdef ROMSX_USE_NETCDF
        nc_init_file.resize(max_level+1);

        // NetCDF wrfinput initialization files -- possibly multiple files at each of multiple levels
        //        but we always have exactly one file at level 0
        for (int lev = 0; lev <= max_level; lev++)
        {
            const std::string nc_file_names = amrex::Concatenate("nc_init_file_",lev,1);
            if (pp.contains(nc_file_names.c_str()))
            {
                int num_files = pp.countval(nc_file_names.c_str());
                num_files_at_level[lev] = num_files;
                nc_init_file[lev].resize(num_files);
                pp.queryarr(nc_file_names.c_str(), nc_init_file[lev],0,num_files);
                for (int j = 0; j < num_files; j++)
                    amrex::Print() << "Reading NC init file names at level " << lev << " and index " << j << " : " << nc_init_file[lev][j] << std::endl;
            }
        }

        // NetCDF wrfbdy lateral boundary file
        pp.query("nc_bdy_file", nc_bdy_file);
#endif

        // Output format
        pp.query("plotfile_type", plotfile_type);
        if (plotfile_type != "amrex" &&
            plotfile_type != "netcdf" && plotfile_type != "NetCDF" &&
            plotfile_type != "hdf5"   && plotfile_type != "HDF5" )
        {
            amrex::Print() << "User selected plotfile_type = " << plotfile_type << std::endl;
            amrex::Abort("Dont know this plotfile_type");
        }
        pp.query("plot_file_1", plot_file_1);
        pp.query("plot_file_2", plot_file_2);
        pp.query("plot_int_1", plot_int_1);
        pp.query("plot_int_2", plot_int_2);
    }

#ifdef ROMSX_USE_MULTIBLOCK
    solverChoice.pp_prefix = pp_prefix;
#endif

    solverChoice.init_params();
}

// Create horizontal average quantities
void
ROMSX::MakeHorizontalAverages ()
{
    // We need to create horizontal averages for 5 variables:
    // density, temperature, pressure, qc, qv (if present)

    // First, average down all levels
    AverageDown();

    // get the number of cells in z at level 0
    int dir_z = AMREX_SPACEDIM-1;
    auto domain = geom[0].Domain();
    int size_z = domain.length(dir_z);
    int start_z = domain.smallEnd()[dir_z];
    Real area_z = static_cast<Real>(domain.length(0));
    area_z *= domain.length(1);

    // resize the level 0 horizontal average vectors
    h_havg_density.resize(size_z, 0.0_rt);
    h_havg_temperature.resize(size_z, 0.0_rt);
#ifdef ROMSX_USE_MOISTURE
    h_havg_qv.resize(size_z, 0.0_rt);
    h_havg_qc.resize(size_z, 0.0_rt);
#endif

    // get the cell centered data and construct sums
    auto& mf_cons = vars_new[0][Vars::cons];
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(mf_cons); mfi.isValid(); ++mfi) {
        const Box& box = mfi.validbox();
        const IntVect& se = box.smallEnd();
        const IntVect& be = box.bigEnd();
        auto  arr_cons = mf_cons[mfi].array();

#ifdef ROMSX_USE_MOISTURE
        FArrayBox fab_reduce(box, 5);
#else
        FArrayBox fab_reduce(box, 3);
#endif
        Elixir elx_reduce = fab_reduce.elixir();
        auto arr_reduce = fab_reduce.array();

        ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            Real dens = arr_cons(i, j, k, Cons::Rho);
            arr_reduce(i, j, k, 0) = dens;
            arr_reduce(i, j, k, 1) = arr_cons(i, j, k, Cons::Temp) / dens;
#ifdef ROMSX_USE_MOISTURE
            arr_reduce(i, j, k, 3) = arr_cons(i, j, k, Cons::RhoQv) / dens;
            arr_reduce(i, j, k, 4) = arr_cons(i, j, k, Cons::RhoQc) / dens;
#endif
        });

        for (int k=se[dir_z]; k <= be[dir_z]; ++k) {
            Box kbox(box); kbox.setSmall(dir_z,k); kbox.setBig(dir_z,k);
            h_havg_density     [k-start_z] += fab_reduce.sum<RunOn::Device>(kbox,0);
            h_havg_temperature [k-start_z] += fab_reduce.sum<RunOn::Device>(kbox,1);
#ifdef ROMSX_USE_MOISTURE
            h_havg_qv          [k-start_z] += fab_reduce.sum<RunOn::Device>(kbox,3);
            h_havg_qc          [k-start_z] += fab_reduce.sum<RunOn::Device>(kbox,4);
#endif
        }
    }

    // combine sums from different MPI ranks
    ParallelDescriptor::ReduceRealSum(h_havg_density.dataPtr(), h_havg_density.size());
    ParallelDescriptor::ReduceRealSum(h_havg_temperature.dataPtr(), h_havg_temperature.size());
#ifdef ROMSX_USE_MOISTURE
    ParallelDescriptor::ReduceRealSum(h_havg_qv.dataPtr(), h_havg_qv.size());
    ParallelDescriptor::ReduceRealSum(h_havg_qc.dataPtr(), h_havg_qc.size());
#endif

    // divide by the total number of cells we are averaging over
    for (int k = 0; k < size_z; ++k) {
        h_havg_density[k]     /= area_z;
        h_havg_temperature[k] /= area_z;
#ifdef ROMSX_USE_MOISTURE
        h_havg_qv[k]          /= area_z;
        h_havg_qc[k]          /= area_z;
#endif
    }

    // resize device vectors
    d_havg_density.resize(size_z, 0.0_rt);
    d_havg_temperature.resize(size_z, 0.0_rt);
    d_havg_pressure.resize(size_z, 0.0_rt);
#ifdef ROMSX_USE_MOISTURE
    d_havg_qv.resize(size_z, 0.0_rt);
    d_havg_qc.resize(size_z, 0.0_rt);
#endif

    // copy host vectors to device vectors
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_havg_density.begin(), h_havg_density.end(), d_havg_density.begin());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_havg_temperature.begin(), h_havg_temperature.end(), d_havg_temperature.begin());
#ifdef ROMSX_USE_MOISTURE
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_havg_qv.begin(), h_havg_qv.end(), d_havg_qv.begin());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_havg_qc.begin(), h_havg_qc.end(), d_havg_qc.begin());
#endif
}

// Set covered coarse cells to be the average of overlying fine cells for all levels
void
ROMSX::AverageDown ()
{
    for (int lev = finest_level-1; lev >= 0; --lev)
    {
        AverageDownTo(lev);
    }
}

// Set covered coarse cells to be the average of overlying fine cells at level crse_lev
void
ROMSX::AverageDownTo (int crse_lev)
{
    for (int var_idx = 0; var_idx < Vars::NumTypes; ++var_idx) {
        const BoxArray& ba(vars_new[crse_lev][var_idx].boxArray());
        if (ba[0].type() == IntVect::TheZeroVector())
            amrex::average_down(vars_new[crse_lev+1][var_idx], vars_new[crse_lev][var_idx],
                                0, vars_new[crse_lev][var_idx].nComp(), refRatio(crse_lev));
        else // We assume the arrays are face-centered if not cell-centered
            amrex::average_down_faces(vars_new[crse_lev+1][var_idx], vars_new[crse_lev][var_idx],
                                      refRatio(crse_lev),geom[crse_lev]);
    }
}
