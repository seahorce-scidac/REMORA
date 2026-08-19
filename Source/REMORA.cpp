/**
 * \file REMORA.cpp
 */

#include <REMORA_prob_common.H>
#include <REMORA.H>

#ifdef REMORA_USE_NETCDF
#include "REMORA_NCFile.H"
#endif

#include <AMReX_buildInfo.H>

using namespace amrex;

amrex::Real REMORA::startCPUTime        = zero;
amrex::Real REMORA::previousCPUTimeUsed = zero;

Vector<AMRErrorTag> REMORA::ref_tags;

SolverChoice REMORA::solverChoice;

// Time step control
amrex::Real REMORA::cfl           =  Real(0.8);
amrex::Real REMORA::fixed_dt      = -one;
amrex::Real REMORA::fixed_fast_dt = -one;
amrex::Real REMORA::change_max    =  Real(1.1);

int   REMORA::fixed_ndtfast_ratio = 0;

// Dictate verbosity in screen output
int         REMORA::verbose       = 0;

// Frequency of diagnostic output
int         REMORA::sum_interval  = -1;
amrex::Real REMORA::sum_per       = -one;

// Minimum number of digits in plotfile name
int         REMORA::file_min_digits = 5;

// Do we include staggered velocities in the plotfile?
bool        REMORA::plot_staggered_vels = false;

// Do we include nodal data (Nu_nd) in the plotfile?
bool        REMORA::plot_nodal_data = true;

// Native AMReX vs NetCDF
PlotfileType REMORA::plotfile_type    = PlotfileType::amrex;

#ifdef REMORA_USE_NETCDF

int   REMORA::total_nc_plot_file_step = 1;

// Do we write one file per timestep (false) or one file for all timesteps (true)
bool  REMORA::write_history_file      = true;

// NetCDF initialization file
amrex::Vector<std::string> REMORA::nc_bdry_file = {""}; // Must provide via input
amrex::Vector<amrex::Vector<std::string>> REMORA::nc_init_file = {{""}}; // Must provide via input
amrex::Vector<amrex::Vector<std::string>> REMORA::nc_grid_file = {{""}}; // Must provide via input
#endif

/**
 * constructor:
 * - reads in parameters from inputs file
 * - sizes multilevel arrays and data structures
 * - initializes BCRec boundary condition object
 */
REMORA::REMORA ()
{
    BL_PROFILE("REMORA::REMORA()");

    if (ParallelDescriptor::IOProcessor()) {
        const char* remora_hash = amrex::buildInfoGetGitHash(1);
        const char* amrex_hash = amrex::buildInfoGetGitHash(2);
        const char* buildgithash = amrex::buildInfoGetBuildGitHash();
        const char* buildgitname = amrex::buildInfoGetBuildGitName();

        if (strlen(remora_hash) > 0) {
          amrex::Print() << "\n"
                         << "REMORA git hash: " << remora_hash << "\n";
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

    // Blocking factor in z set to very large value to be > nz
    // This guarantees that there will be no domain decomposition in the z-direction
    // We have to set this by hand here because setting it in the input file will
    // cause checks in the AmrCore constructor to fail.
    Vector<IntVect> blocking_factor_vec = Vector<IntVect>();
    blocking_factor_vec.resize(max_level+1);
    for (int lev = 0; lev <= max_level; ++lev) {
        blocking_factor_vec[lev] = blockingFactor(lev);
        blocking_factor_vec[lev][2] = 4096;
    }
    SetBlockingFactor(blocking_factor_vec);

    const std::string& pv3d = "plot_vars_3d"; set3DPlotVariables(pv3d);
    const std::string& pv2d = "plot_vars_2d"; set2DPlotVariables(pv2d);

    prob = amrex_probinit(geom[0].ProbLo(),geom[0].ProbHi());

    // Geometry on all levels has been defined already.

    // No valid BoxArray and DistributionMapping have been defined.
    // But the arrays for them have been resized.

    int nlevs_max = max_level + 1;

    istep.resize(nlevs_max, 0);
    nsubsteps.resize(nlevs_max, 1);
    for (int lev = 1; lev <= max_level; ++lev) {
        nsubsteps[lev] = do_substep ? MaxRefRatio(lev-1) : 1;
    }

    physbcs.resize(nlevs_max);

    t_new.resize(nlevs_max, zero);
    t_old.resize(nlevs_max, -bogus_large_value);
    dt.resize(nlevs_max, bogus_large_value);

    cons_new.resize(nlevs_max);
    cons_old.resize(nlevs_max);
    xvel_new.resize(nlevs_max);
    xvel_old.resize(nlevs_max);
    yvel_new.resize(nlevs_max);
    yvel_old.resize(nlevs_max);
    zvel_new.resize(nlevs_max);
    zvel_old.resize(nlevs_max);

    advflux_reg.resize(nlevs_max);

    // Initialize tagging criteria for mesh refinement
    refinement_criteria_setup();

    IntVect cum_ref_ratio = IntVect(1,1,0);
    cum_ref_ratios.push_back(cum_ref_ratio);
    // We have already read in the ref_Ratio (via amr.ref_ratio =) but we need to enforce
    //     that there is no refinement in the vertical so we test on that here.
    for (int lev = 0; lev < max_level; ++lev)
    {
       amrex::Print() << "Refinement ratio at level " << lev << " set to be " <<
          ref_ratio[lev][0]  << " " << ref_ratio[lev][1]  <<  " " << ref_ratio[lev][2] << std::endl;

       if (ref_ratio[lev][2] != 1)
       {
           amrex::Print() << "********************************************************************************" << std::endl;
           amrex::Print() << "We don't allow refinement in the vertical -- make sure to set ref_ratio = 1 in z" << std::endl;
           amrex::Print() << "It's possible you set amr.ref_ratio when you meant to set amr.ref_ratio_vect    " << std::endl;
           amrex::Print() << "********************************************************************************" << std::endl;
           amrex::Abort();
       }

       cum_ref_ratio[0] *= ref_ratio[lev][0];
       cum_ref_ratio[1] *= ref_ratio[lev][1];
       cum_ref_ratios.push_back(cum_ref_ratio);
    }
}

REMORA::REMORA (const amrex::RealBox& rb, int max_level_in, const amrex::Vector<int>& n_cell_in, int coord, const amrex::Vector<amrex::IntVect>& ref_ratio_in, const amrex::Array<int,AMREX_SPACEDIM>& is_per, std::string prefix)
    : amrex::AmrCore (rb, max_level_in, n_cell_in, coord, ref_ratio_in, is_per)
{
    BL_PROFILE("REMORA::REMORA(explicit)");
    pp_prefix = prefix;

    if (ParallelDescriptor::IOProcessor()) {
        const char* remora_hash = amrex::buildInfoGetGitHash(1);
        const char* amrex_hash = amrex::buildInfoGetGitHash(2);
        const char* buildgithash = amrex::buildInfoGetBuildGitHash();
        const char* buildgitname = amrex::buildInfoGetBuildGitName();

        if (strlen(remora_hash) > 0) {
          amrex::Print() << "\n"
                         << "REMORA git hash: " << remora_hash << "\n";
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

    const std::string& pv3d = "plot_vars_3d"; set3DPlotVariables(pv3d);
    const std::string& pv2d = "plot_vars_2d"; set2DPlotVariables(pv2d);

    prob = amrex_probinit(geom[0].ProbLo(),geom[0].ProbHi());

    int nlevs_max = max_level + 1;

    istep.resize(nlevs_max, 0);
    nsubsteps.resize(nlevs_max, 1);
    for (int lev = 1; lev <= max_level; ++lev) {
        nsubsteps[lev] = do_substep ? MaxRefRatio(lev-1) : 1;
    }

    physbcs.resize(nlevs_max);

    t_new.resize(nlevs_max, zero);
    t_old.resize(nlevs_max, -bogus_large_value);
    dt.resize(nlevs_max, bogus_large_value);

    cons_new.resize(nlevs_max);
    cons_old.resize(nlevs_max);
    xvel_new.resize(nlevs_max);
    xvel_old.resize(nlevs_max);
    yvel_new.resize(nlevs_max);
    yvel_old.resize(nlevs_max);
    zvel_new.resize(nlevs_max);
    zvel_old.resize(nlevs_max);

    advflux_reg.resize(nlevs_max);

    refinement_criteria_setup();

    for (int lev = 0; lev < max_level; ++lev)
    {
       amrex::Print() << "Refinement ratio at level " << lev << " set to be " <<
          ref_ratio[lev][0]  << " " << ref_ratio[lev][1]  <<  " " << ref_ratio[lev][2] << std::endl;

       if (ref_ratio[lev][2] != 1)
       {
           amrex::Print() << "********************************************************************************" << std::endl;
           amrex::Print() << "We don't allow refinement in the vertical -- make sure to set ref_ratio = 1 in z" << std::endl;
           amrex::Print() << "It's possible you set amr.ref_ratio when you meant to set amr.ref_ratio_vect    " << std::endl;
           amrex::Print() << "********************************************************************************" << std::endl;
           amrex::Abort();
       }
    }
}

REMORA::~REMORA ()
{
}

void
REMORA::init_scalar_metadata ()
{
    cons_names.clear();
    cons_names.reserve(ncons);
    cons_names.emplace_back("temp");
    cons_names.emplace_back("salt");

    // Passive (dye) scalars come first, then the biology block, matching the component
    // layout: temp, salt, tracer, tracer_1, ..., NO3, NH4, ...
    if (nscalar > 0) {
        cons_names.emplace_back("tracer");
        for (int i = 1; i < nscalar; ++i) {
            cons_names.emplace_back("tracer_" + std::to_string(i));
        }
    }

    if (REMORABiology::has_biology(biology_model)) {
        const auto bio_names = REMORABiology::tracer_names(biology_model, fennel_params);
        for (const auto& name : bio_names) {
            cons_names.emplace_back(name);
        }
    }

    AMREX_ALWAYS_ASSERT(static_cast<int>(cons_names.size()) == ncons);
}

void
REMORA::Evolve ()
{
    BL_PROFILE_VAR("REMORA::Evolve()",evolve);
    Real cur_time = t_new[0];

    // Take one coarse timestep by calling timeStep -- which recursively calls timeStep
    //      for finer levels (with or without subcycling)
    for (int step = istep[0]; step < max_step && cur_time < stop_time; ++step)
    {
        amrex::Print() << "\nCoarse STEP " << step+1 << " starts ..." << std::endl;

        ComputeDt();

        int lev = 0;
        int iteration = 1;
        auto dEvolveTime0 = amrex::second();

        if (max_level == 0) {
            timeStep(lev, cur_time, iteration);
        }
        else {
            timeStepML(cur_time, iteration);
        }

        cur_time  += dt[0];

        amrex::Print() << "Coarse STEP " << step+1 << " ends." << " TIME = " << cur_time
                       << " DT = " << dt[0]  << std::endl;

        if (verbose > 0)
        {
            auto dEvolveTime = amrex::second() - dEvolveTime0;
            ParallelDescriptor::ReduceRealMax(dEvolveTime,ParallelDescriptor::IOProcessorNumber());
            amrex::Print() << "Timestep time = " << dEvolveTime << " seconds." << '\n';
        }

        WriteAtIntermediateTime(step, cur_time);

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

    BL_PROFILE_VAR_STOP(evolve);

    WriteAtFinalTime();
}

void
REMORA::WriteAtFinalTime()
{

    if ( (plot_int > 0 || plot_int_time > zero) && istep[0] > last_plot_file_step)
    {
        WritePlotFile(istep[0]);
        history_count++;
    }

    if ((check_int > 0 || check_int_time > zero) && istep[0] > last_check_file_step) {
        WriteCheckpointFile();
    }
}

void
REMORA::WriteAtIntermediateTime(int step, amrex::Real cur_time)
{
    if ( (plot_int > 0      && (step+1 - last_plot_file_step) == plot_int         ) ||
         (plot_int_time > 0 && (cur_time >= (last_plot_file_time + plot_int_time))) )
    {
        last_plot_file_step = step+1;
        last_plot_file_time = cur_time;
        WritePlotFile(step+1);
        history_count++;
    }

    if ((check_int > 0 && (step+1 - last_check_file_step) == check_int)
            || (check_int_time > 0 && cur_time >= (last_check_file_time + check_int_time))) {
        last_check_file_step = step+1;
        last_check_file_time = cur_time;
        WriteCheckpointFile();
    }
}

/**
 * @param[in   ] nstep    which step we're on
 * @param[in   ] time     current time
 * @param[in   ] dt_lev0  time step on level 0
 */
void
REMORA::post_timestep (int nstep, Real time, Real dt_lev0)
{
    BL_PROFILE("REMORA::post_timestep()");

#ifdef REMORA_USE_PARTICLES
    particleData.Redistribute();
#endif

    if (solverChoice.coupling_type == CouplingType::two_way)
    {
        for (int lev = finest_level-1; lev >= 0; lev--)
        {
            // This call refluxes from the lev/lev+1 interface onto lev
            //getAdvFluxReg(lev+1)->Reflux(*cons_new[lev], 0, 0, NCONS);

            // We need to do this before anything else because refluxing changes the
            // values of coarse cells underneath fine grids with the assumption they'll
            // be over-written by averaging down
            //
            AverageDownTo(lev);
        }
    }

    if (is_it_time_for_action(nstep, time, dt_lev0, sum_interval, sum_per)) {
        sum_integrated_quantities(time);
    }
}

/**
 * This is called from main.cpp and handles all initialization, whether from start or restart
 */
void
REMORA::InitData ()
{
    BL_PROFILE("REMORA::InitData()");
    if (running_with_coupling_driver) {
        amrex::Print() << "REMORA InitData: driver-managed atm2ocn coupling enabled"
                       << " two_way=" << (driver_uses_two_way_coupling ? 1 : 0)
                       << " active_contract="
                       << (DriverUsesStateForcing(driver_atmos_forcing_mode) ? "state" : "flux")
                       << "\n";
    }
    // Initialize the start time for our CPU-time tracker
    startCPUTime = Real(ParallelDescriptor::second());

    // Map the words in the inputs file to BC types, then translate
    //     those types into what they mean for each variable
    init_bcs();

    // Init vertical stretching coeffs
    init_stretch_coeffs();

    last_plot_file_step = -1;
    last_check_file_step = -1;
    last_plot_file_time = -one;
    last_check_file_time = -one;

    if (restart_chkfile == "") {
        // start simulation from the beginning

        InitFromScratch(start_time);

        if (solverChoice.coupling_type == CouplingType::two_way) {
            AverageDown();
        }

    } else { // Restart from a checkpoint

        restart();

    }
#ifdef REMORA_USE_MOAB
    InitMOABMesh();
#endif
    // Initialize flux registers (whether we start from scratch or restart)
    if (solverChoice.coupling_type == CouplingType::two_way) {
        advflux_reg[0] = nullptr;
        for (int lev = 1; lev <= finest_level; lev++)
        {
            advflux_reg[lev].reset( new YAFluxRegister(grids[lev], grids[lev-1],
                                                   dmap[lev],  dmap[lev-1],
                                                   geom[lev],  geom[lev-1],
                                              ref_ratio[lev-1], lev, ncons));
        }
    }

    // Fill ghost cells/faces
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        if (lev > 0 && cf_width >= 0) {
            Construct_REMORAFillPatchers(lev);
        }

        if (restart_chkfile == "") {
            FillPatch(lev, t_new[lev], *cons_new[lev], cons_new, BCVars::cons_bc, BdyVars::t, 0, true, false,0,0,zero,*cons_new[lev]);
            FillPatch(lev, t_new[lev], *xvel_new[lev], xvel_new, xvel_bc(), BdyVars::u, 0, true, false,0,0,zero,*xvel_new[lev]);
            FillPatch(lev, t_new[lev], *yvel_new[lev], yvel_new, yvel_bc(), BdyVars::v, 0, true, false,0,0,zero,*yvel_new[lev]);
            FillPatch(lev, t_new[lev], *zvel_new[lev], zvel_new, zvel_bc(), BdyVars::null, 0, true, false);

            // Copy from new into old just in case when initializing from scratch
            int ngs   = cons_new[lev]->nGrow();
            int ngvel = xvel_new[lev]->nGrow();
            MultiFab::Copy(*cons_old[lev],*cons_new[lev],0,0,ncons,ngs);
            MultiFab::Copy(*xvel_old[lev],*xvel_new[lev],0,0,1,ngvel);
            MultiFab::Copy(*yvel_old[lev],*yvel_new[lev],0,0,1,ngvel);
            MultiFab::Copy(*zvel_old[lev],*zvel_new[lev],0,0,1,IntVect(ngvel,ngvel,0));
        }
    } // lev

    // Check for additional plotting variables that are available after
    // particle containers are setup.
    const std::string& pv3d = "plot_vars_3d"; append3DPlotVariables(pv3d);
    const std::string& pv2d = "plot_vars_2d"; append2DPlotVariables(pv2d);

    if (restart_chkfile == "" && (check_int > 0 || check_int_time > zero))
    {
        WriteCheckpointFile();
        last_check_file_step = 0;
    }

    if ( (restart_chkfile == "") ||
         (restart_chkfile != "" && plot_file_on_restart) )
    {
        if (plot_int > 0 || plot_int_time > zero)
        {
            int step0 = 0;
            WritePlotFile(step0);
            history_count++;
            last_plot_file_step = istep[0];
        }
    }

    if (is_it_time_for_action(istep[0], t_new[0], dt[0], sum_interval, sum_per)) {
        sum_integrated_quantities(t_new[0]);
    }

    ComputeDt();

}

/**
 * @param[in   ] lev    level to operate on
 */
void
REMORA::Construct_REMORAFillPatchers (int lev)
{
    BL_PROFILE("REMORA::Construct_REMORAFillPatchers()");
    amrex::Print() << ":::Construct_REMORAFillPatchers " << lev << std::endl;

    auto& ba_fine  = cons_new[lev  ]->boxArray();
    auto& ba_crse  = cons_new[lev-1]->boxArray();
    auto& dm_fine  = cons_new[lev  ]->DistributionMap();
    auto& dm_crse  = cons_new[lev-1]->DistributionMap();

    BoxList bl2d_fine = ba_fine.boxList();
    for (auto& b : bl2d_fine) {
        b.setRange(2,0);
    }
    BoxArray ba2d_fine(std::move(bl2d_fine));

    BoxList bl2d_crse = ba_crse.boxList();
    for (auto& b : bl2d_crse) {
        b.setRange(2,0);
    }
    BoxArray ba2d_crse(std::move(bl2d_crse));

    int ncomp = cons_new[lev]->nComp();

    FPr_c.emplace_back(ba_fine, dm_fine, geom[lev]  ,
                       ba_crse, dm_crse, geom[lev-1],
                       -cf_width, -cf_set_width, ncomp, &cell_cons_interp);
    FPr_u.emplace_back(convert(ba_fine, IntVect(1,0,0)), dm_fine, geom[lev]  ,
                       convert(ba_crse, IntVect(1,0,0)), dm_crse, geom[lev-1],
                       -cf_width, -cf_set_width, 1, &face_cons_linear_interp);
    FPr_v.emplace_back(convert(ba_fine, IntVect(0,1,0)), dm_fine, geom[lev]  ,
                       convert(ba_crse, IntVect(0,1,0)), dm_crse, geom[lev-1],
                       -cf_width, -cf_set_width, 1, &face_cons_linear_interp);
    FPr_w.emplace_back(convert(ba_fine, IntVect(0,0,1)), dm_fine, geom[lev]  ,
                       convert(ba_crse, IntVect(0,0,1)), dm_crse, geom[lev-1],
                       -cf_width, -cf_set_width, 1, &face_cons_linear_interp);

    FPr_ubar.emplace_back(convert(ba2d_fine, IntVect(1,0,0)), dm_fine, geom[lev]  ,
                       convert(ba2d_crse, IntVect(1,0,0)), dm_crse, geom[lev-1],
                       -cf_width, -cf_set_width, 3, &face_cons_linear_interp);
    FPr_vbar.emplace_back(convert(ba2d_fine, IntVect(0,1,0)), dm_fine, geom[lev]  ,
                       convert(ba2d_crse, IntVect(0,1,0)), dm_crse, geom[lev-1],
                       -cf_width, -cf_set_width, 3, &face_cons_linear_interp);
}

/**
 * @param[in   ] lev    level to operate on
 */
void
REMORA::Define_REMORAFillPatchers (int lev)
{
    BL_PROFILE("REMORA::Define_REMORAFillPatchers()");
    amrex::Print() << ":::Define_REMORAFillPatchers " << lev << std::endl;

    auto& ba_fine  = cons_new[lev  ]->boxArray();
    auto& ba_crse  = cons_new[lev-1]->boxArray();
    auto& dm_fine  = cons_new[lev  ]->DistributionMap();
    auto& dm_crse  = cons_new[lev-1]->DistributionMap();

    BoxList bl2d_fine = ba_fine.boxList();
    for (auto& b : bl2d_fine) {
        b.setRange(2,0);
    }
    BoxArray ba2d_fine(std::move(bl2d_fine));

    BoxList bl2d_crse = ba_crse.boxList();
    for (auto& b : bl2d_crse) {
        b.setRange(2,0);
    }
    BoxArray ba2d_crse(std::move(bl2d_crse));


    int ncomp = cons_new[lev]->nComp();

    FPr_c[lev-1].Define(ba_fine, dm_fine, geom[lev]  ,
                        ba_crse, dm_crse, geom[lev-1],
                        -cf_width, -cf_set_width, ncomp, &cell_cons_interp);
    FPr_u[lev-1].Define(convert(ba_fine, IntVect(1,0,0)), dm_fine, geom[lev]  ,
                        convert(ba_crse, IntVect(1,0,0)), dm_crse, geom[lev-1],
                        -cf_width, -cf_set_width, 1, &face_cons_linear_interp);
    FPr_v[lev-1].Define(convert(ba_fine, IntVect(0,1,0)), dm_fine, geom[lev]  ,
                        convert(ba_crse, IntVect(0,1,0)), dm_crse, geom[lev-1],
                        -cf_width, -cf_set_width, 1, &face_cons_linear_interp);
    FPr_w[lev-1].Define(convert(ba_fine, IntVect(0,0,1)), dm_fine, geom[lev]  ,
                        convert(ba_crse, IntVect(0,0,1)), dm_crse, geom[lev-1],
                        -cf_width, -cf_set_width, 1, &face_cons_linear_interp);

    FPr_ubar[lev-1].Define(convert(ba2d_fine, IntVect(1,0,0)), dm_fine, geom[lev]  ,
                        convert(ba2d_crse, IntVect(1,0,0)), dm_crse, geom[lev-1],
                        -cf_width, -cf_set_width, 3, &face_cons_linear_interp);
    FPr_vbar[lev-1].Define(convert(ba2d_fine, IntVect(0,1,0)), dm_fine, geom[lev]  ,
                        convert(ba2d_crse, IntVect(0,1,0)), dm_crse, geom[lev-1],
                        -cf_width, -cf_set_width, 3, &face_cons_linear_interp);
}

void
REMORA::restart ()
{
    BL_PROFILE("REMORA::restart()");
    ReadCheckpointFile();

    // We set this here so that we don't over-write the checkpoint file we just started from
    last_check_file_step = istep[0];
}

/**
 * @param[in   ] lev    level to operate on
 */
void
REMORA::set_zeta (int lev)
{
    BL_PROFILE("REMORA::set_zeta()");
    if (lev==0) {
        if (hires_init_level < 0) {
            if (solverChoice.ic_type == IC_Type::analytic) {
                prob->init_analytic_zeta(lev, geom[lev], solverChoice, *this, *vec_zeta[lev]);
            } else if (solverChoice.ic_type == IC_Type::netcdf) {
#ifdef REMORA_USE_NETCDF
                amrex::Print() << "Calling init_zeta_from_netcdf on level " << lev << std::endl;
                init_zeta_from_netcdf(lev);
                amrex::Print() << "Sea surface height loaded from netcdf file \n " << std::endl;
#endif
            } else {
                amrex::Abort("Unknown IC_Type");
            }
        } else {
            set_zeta_averaged_down(lev);
        }
        vec_zeta[lev]->FillBoundary(geom[lev].periodicity());
    } else {
        // If our level is higher than the high resolution grid or initialization
        // is analytic, interpolate from level below. Otherwise, copy over the bathymetry
        // data that has been averaged down
        if (lev > hires_init_level) {
            Real dummy_time = zero;
            FillCoarsePatch(lev,dummy_time,vec_zeta[lev].get(), vec_zeta[lev-1].get(),BCVars::cons_bc);
        } else {
            set_zeta_averaged_down(lev);
            vec_zeta[lev]->FillBoundary(geom[lev].periodicity());
        }
    }
    set_zeta_average(lev);
}

/**
 * @param[in   ] lev    level to operate on
 */
void
REMORA::set_bathymetry (int lev)
{
    BL_PROFILE("REMORA::bathymetry()");
    // Only set bathymetry on level 0, and interpolate for finer levels
    if (lev==0) {
        if (solverChoice.flat_bathymetry) {
            init_flat_bathymetry(lev);
        // If grid data is not defined on a level > 0 (negative level) then
        // initialize from low-resolution grid normally. Otherwise use high-resolution
        // grid data averaged down to level 0
        } else if (hires_grid_level < 0) {
            if (solverChoice.ic_type == IC_Type::analytic) {
                prob->init_analytic_bathymetry(lev, geom[lev], solverChoice, *this, *vec_h[lev]);
            } else if (solverChoice.ic_type == IC_Type::netcdf) {
#ifdef REMORA_USE_NETCDF
                amrex::Print() << "Calling init_bathymetry_from_netcdf " << std::endl;
                init_bathymetry_from_netcdf(lev);
                amrex::Print() << "Bathymetry loaded from netcdf file \n " << std::endl;
                amrex::Print() << "Calling init_grid_vars_from_netcdf " << std::endl;
                init_grid_vars_from_netcdf(lev);
                amrex::Print() << "Grid variables loaded from netcdf file \n " << std::endl;
#endif
            } else {
                amrex::Abort("Unknown IC_Type");
            }
        } else {
            set_bathymetry_averaged_down(lev);
            // Only the netcdf path fills vec_pm/pn_full_domain; with analytic initialization
            // init_bathymetry_full_domain_from_analytic fills h alone, and set_grid_scale
            // below derives pm/pn from the geometry.
            if (solverChoice.ic_type == IC_Type::netcdf) {
                set_grid_vars_averaged_down(lev);
            }
        }
        // Need FillBoundary to fill at grid-grid boundaries, and EnforcePeriodicity
        // to make sure ghost cells in the domain corners are consistent.
        vec_h[lev]->FillBoundary(geom[lev].periodicity());
        vec_h[lev]->EnforcePeriodicity(geom[lev].periodicity());
    } else {
        // If our level is higher than the high resolution grid or initialization
        // is analytic, interpolate from level below. Otherwise, copy over the bathymetry
        // data that has been averaged down
        if (lev > hires_grid_level) {
            Real dummy_time = zero;
            FillCoarsePatch(lev,dummy_time,vec_h[lev].get(), vec_h[lev-1].get(),
                    foextrap_periodic_bc(),BdyVars::null,0,false);
            FillCoarsePatch(lev,dummy_time,vec_h[lev].get(), vec_h[lev-1].get(),
                    foextrap_periodic_bc(),BdyVars::null,1,false);
        } else {
            set_bathymetry_averaged_down(lev);
            vec_h[lev]->FillBoundary(geom[lev].periodicity());
            vec_h[lev]->EnforcePeriodicity(geom[lev].periodicity());
        }
    }
    set_grid_scale(lev);
}

/**
 * @param[in   ] lev   level to operate on
 */
void
REMORA::set_bathymetry_averaged_down (int lev) {
    Real dummy_time = zero;
    // Note: don't understand why the grow vector args aren't vec_h and then vec_h_full_domain
    ParallelCopy(*vec_h[lev].get(), *vec_h_full_domain[lev].get(), 0, 0, 1,vec_h_full_domain[lev]->nGrowVect(),vec_h[lev]->nGrowVect());
    ParallelCopy(*vec_h[lev].get(), *vec_h_full_domain[lev].get(), 0, 1, 1,vec_h_full_domain[lev]->nGrowVect(),vec_h[lev]->nGrowVect());
    FillPatch(lev,dummy_time,*vec_h[lev],GetVecOfPtrs(vec_h),
            foextrap_periodic_bc(),
            BdyVars::null,0,false,false,1);
    FillPatch(lev,dummy_time,*vec_h[lev],GetVecOfPtrs(vec_h),
            foextrap_periodic_bc(),
            BdyVars::null,1,false,false,1);
}

/**
 * @param[in   ] lev   level to operate on
 */
void
REMORA::set_grid_vars_averaged_down (int lev) {
    Real dummy_time = zero;
    ParallelCopy(*vec_pm[lev].get(), *vec_pm_full_domain[lev].get(), 0, 0, 1,
            vec_pm_full_domain[lev]->nGrowVect(),vec_pm[lev]->nGrowVect());
    ParallelCopy(*vec_pn[lev].get(), *vec_pn_full_domain[lev].get(), 0, 0, 1,
            vec_pn_full_domain[lev]->nGrowVect(),vec_pn[lev]->nGrowVect());
    FillPatch(lev,dummy_time,*vec_pm[lev],GetVecOfPtrs(vec_pm),
            foextrap_periodic_bc(),
            BdyVars::null,0,false);
    FillPatch(lev,dummy_time,*vec_pn[lev],GetVecOfPtrs(vec_pn),
            foextrap_periodic_bc(),
            BdyVars::null,0,false);
}

/**
 * @param[in   ] lev   level to operate on
 */
void
REMORA::set_zeta_averaged_down (int lev) {
    ParallelCopy(*vec_zeta[lev].get(), *vec_zeta_full_domain[lev].get(), 0, 0, 1,
            vec_zeta_full_domain[lev]->nGrowVect(),vec_zeta[lev]->nGrowVect());
    FillPatch(lev, t_new[lev], *vec_zeta[lev], GetVecOfPtrs(vec_zeta), zeta_bc(), bdy_zeta(),
                  0, false,false,0,0,zero,*vec_zeta[lev]);
}

/**
 * @param[in   ] lev   level to operate on
 */
void
REMORA::set_init_data_averaged_down (int lev) {
    ParallelCopy(*cons_new[lev], *vec_cons_full_domain[lev], 0, 0, ncons,
            vec_cons_full_domain[lev]->nGrowVect(),cons_new[lev]->nGrowVect());
    ParallelCopy(*xvel_new[lev], *vec_xvel_full_domain[lev], 0, 0, 1,
            vec_xvel_full_domain[lev]->nGrowVect(),xvel_new[lev]->nGrowVect());
    ParallelCopy(*yvel_new[lev], *vec_yvel_full_domain[lev], 0, 0, 1,
            vec_yvel_full_domain[lev]->nGrowVect(),yvel_new[lev]->nGrowVect());

    FillPatch(lev, t_new[lev], *cons_new[lev], cons_new, BCVars::cons_bc, BdyVars::t, 0, true, false,0,0,zero,*cons_new[lev]);
    FillPatch(lev, t_new[lev], *xvel_new[lev], xvel_new, xvel_bc(), BdyVars::u, 0, true, false,0,0,zero,*xvel_new[lev]);
    FillPatch(lev, t_new[lev], *yvel_new[lev], yvel_new, yvel_bc(), BdyVars::v, 0, true, false,0,0,zero,*yvel_new[lev]);
}

/**
 * @param[in   ] lev    level to operate on
 */
void
REMORA::set_coriolis(int lev) {
    BL_PROFILE("REMORA::set_coriolis()");
    if (solverChoice.use_coriolis) {
        if (solverChoice.coriolis_type == Cor_Type::analytic) {
            prob->init_analytic_coriolis(lev, geom[lev], solverChoice, *this, *vec_fcor[lev]);
        } else if (solverChoice.coriolis_type == Cor_Type::beta_plane) {
            init_beta_plane_coriolis(lev);
#ifdef REMORA_USE_NETCDF
        } else if (solverChoice.coriolis_type == Cor_Type::netcdf) {
            if (lev == 0) {
                amrex::Print() << "Calling init_coriolis_from_netcdf " << std::endl;
                init_coriolis_from_netcdf(lev);
                amrex::Print() << "Coriolis loaded from netcdf file \n" << std::endl;
            } else {
                Real dummy_time = zero;
                FillCoarsePatch(lev,dummy_time,vec_fcor[lev].get(), vec_fcor[lev-1].get(),BCVars::cons_bc);
            }
#endif
        } else {
            Abort("Don't know this coriolis_type!");
        }

        Real time = zero;
        FillPatch(lev, time, *vec_fcor[lev], GetVecOfPtrs(vec_fcor), foextrap_bc());
        vec_fcor[lev]->EnforcePeriodicity(geom[lev].periodicity());
    }
}

void
REMORA::init_set_vmix(int lev) {
    BL_PROFILE("REMORA::init_set_vmix()");
    if (solverChoice.vert_mixing_type == VertMixingType::analytic) {
        set_analytic_vmix(lev);
    } else if (solverChoice.vert_mixing_type == VertMixingType::GLS) {
        init_gls_vmix(lev, solverChoice);
        // The GLS initialization just sets the multifab to a value, so there's
        // no need to call FillPatch here
    } else {
        Abort("Don't know this vertical mixing type");
    }
}

/**
 * @param[in   ] lev    level to operate on
 */
void
REMORA::set_analytic_vmix(int lev) {
    BL_PROFILE("REMORA::set_analytic_vmix()");
    Real time = zero;
    vec_Akv[lev]->setVal(solverChoice.Akv_bak);
    for (int n = 0; n < NAT; n++) {
        vec_Akt[lev]->setVal(solverChoice.Akt_bak[n], n, 1);
    }
    prob->init_analytic_vmix(lev, geom[lev], solverChoice, *this,*vec_Akv[lev], *vec_Akt[lev]);
    FillPatch(lev, time, *vec_Akv[lev], GetVecOfPtrs(vec_Akv), zvel_bc(), BdyVars::null,0,true,false);
    for (int n = 0; n < NAT; n++) {
        FillPatch(lev, time, *vec_Akt[lev], GetVecOfPtrs(vec_Akt), zvel_bc(), BdyVars::null,n,false,false);
    }
}

/**
 * @param[in   ] lev    level to operate on
 */
void
REMORA::set_masks(int lev)
{
    if (solverChoice.mask_type == MaskType::analytic) {
        prob->init_analytic_masks(lev,geom[lev], solverChoice, *this, *vec_mskr[lev]);
        calculate_nodal_masks(lev);
    } else if (solverChoice.mask_type == MaskType::netcdf) {
#ifdef REMORA_USE_NETCDF
        if (lev == 0) {
            amrex::Print() << "Calling init_masks_from_netcdf level " << lev << std::endl;
            init_masks_from_netcdf(lev);
            amrex::Print() << "Masks loaded from netcdf file \n " << std::endl;
        } else {
            Real dummy_time = zero;
            FillCoarsePatchPC(lev, dummy_time, vec_mskr[lev].get(), vec_mskr[lev-1].get(),
                    foextrap_bc());
            calculate_nodal_masks(lev);
        }
#endif
    }
    fill_3d_masks(lev);
}

/**
 * @param[in   ] lev    level to operate on
 */
void
REMORA::set_hmixcoef(int lev)
{
    BL_PROFILE("REMORA::set_hmixcoef()");

    // Optional AMR scaling: decrease coefficients on refined levels linearly
    // with grid size (i.e., proportional to sqrt(cell area)). For a horizontal
    // refinement ratio rx x ry, the effective scale factor is 1/sqrt(rx*ry).
    Real lev_scale = one;
    if ((solverChoice.scaled_to_grid_amr_scaling == ScaledToGridAMRScaling::linear) && (lev > 0)) {
        Real rf = one;
        for (int l = 0; l < lev; ++l) {
            rf *= std::sqrt(static_cast<Real>(ref_ratio[l][0]) * static_cast<Real>(ref_ratio[l][1]));
        }
        lev_scale = one / rf;
    }

    if (solverChoice.horiz_mixing_type == HorizMixingType::analytic) {
        prob->init_analytic_hmix(lev, geom[lev], solverChoice,
                                 *this, *vec_visc2_p[lev], *vec_visc2_r[lev], *vec_diff2[lev]);

    } else if (solverChoice.horiz_mixing_type == HorizMixingType::constant) {
        vec_visc2_p[lev]->setVal(solverChoice.visc2 * lev_scale);
        vec_visc2_r[lev]->setVal(solverChoice.visc2 * lev_scale);
        for (int n = 0; n < ncons; n++) {
            vec_diff2[lev]->setVal(solverChoice.tnu2[n] * lev_scale, n, 1);
        }

    // Scale harmonic viscosity and diffusivity by the grid size as ROMS
    // does in Utility/ini_hmixcoef.F. Intended for curvilinear grids.
    //
    // Define the ROMS grid factor (grdscl):
    //     G(i,j) = sqrt( 1 / (pm(i,j) * pn(i,j)) )
    //            = sqrt(cell area)
    //     Gmax   = max over grid of G(i,j)
    //
    // Then horizontal harmonic mixing coefficients are scaled as:
    //     nu(i,j)       = nu0    * G(i,j) / Gmax
    //     kappa_n(i,j)  = kappa0 * G(i,j) / Gmax
    //
    // where:
    //     nu0     = solverChoice.visc2
    //     kappa0  = solverChoice.tnu2[n]
    //
    // This makes mixing strongest where grid spacing is largest.
    //
    // NOTE: The normalization (Gmax) is computed over the entire grid (ignoring masks).
    // Therefore, if the largest cell area occurs over land, the maximum over *wet* cells
    // (or in masked output files) may be smaller than the user-specified value.

    } else if (solverChoice.horiz_mixing_type == HorizMixingType::scaled_to_grid) {

        // ------------------------------------------------------------
        // Step 1: Compute grdmax over entire grid
        // ------------------------------------------------------------
        vec_visc2_r[lev]->setVal(solverChoice.visc2);
        vec_visc2_p[lev]->setVal(solverChoice.visc2);
        for (int n = 0; n < ncons; n++) {
            vec_diff2[lev]->setVal(solverChoice.tnu2[n], n, 1);
        }

        // NOTE: This must be GPU-safe. Do not dereference MultiFab data on host.
        // Force the reduction to run in the GPU launch region if GPUs are enabled.
        // (If the launch region is disabled at runtime, ReduceMax may fall back to
        // a host path that can try to read device-only data.)
        amrex::Gpu::LaunchSafeGuard lsg(true);
        Real denom_min = amrex::ReduceMin(*vec_pm[lev], *vec_pn[lev], 0,
            [=] AMREX_GPU_HOST_DEVICE (Box const& bx,
                                      Array4<Real const> const& pm,
                                      Array4<Real const> const& pn) -> Real
            {
                Real local_min = bogus_large_value;
                amrex::Loop(bx, [=,&local_min] (int i, int j, int) noexcept
                {
                    local_min = amrex::min(local_min, pm(i,j,0) * pn(i,j,0));
                });
                return local_min;
            });

        ParallelDescriptor::ReduceRealMin(denom_min);
        if (denom_min <= zero) {
            Abort("scaled_to_grid: found non-positive pm*pn (grid metrics must be > 0)");
        }

        Real grdmax = amrex::ReduceMax(*vec_pm[lev], *vec_pn[lev], 0,
            [=] AMREX_GPU_HOST_DEVICE (Box const& bx,
                                      Array4<Real const> const& pm,
                                      Array4<Real const> const& pn) -> Real
            {
                Real local_max = zero;
                amrex::Loop(bx, [=,&local_max] (int i, int j, int) noexcept
                {
                    Real denom = pm(i,j,0) * pn(i,j,0);
                    if (denom > zero) {
                        Real G = std::sqrt(one / denom);
                        local_max = amrex::max(local_max, G);
                    }
                });
                return local_max;
            });

        ParallelDescriptor::ReduceRealMax(grdmax);
        if (grdmax <= zero) {
            Abort("scaled_to_grid: grdmax <= 0");
        }

        // Optional AMR scaling: decrease coefficients on refined levels linearly
        // with grid size (i.e., proportional to sqrt(cell area)). For a horizontal
        // refinement ratio rx x ry, the effective scale factor is 1/sqrt(rx*ry).
        lev_scale = one;
        if ((solverChoice.scaled_to_grid_amr_scaling == ScaledToGridAMRScaling::linear) && (lev > 0)) {
            Real rf = one;
            for (int l = 0; l < lev; ++l) {
                rf *= std::sqrt(static_cast<Real>(ref_ratio[l][0]) * static_cast<Real>(ref_ratio[l][1]));
            }
            lev_scale = one / rf;
        }

        Real visc0 = solverChoice.visc2 * lev_scale;
        Real cff   = visc0 / grdmax;

        // ------------------------------------------------------------
        // Step 2: Set rho coefficients everywhere
        // ------------------------------------------------------------
        amrex::Gpu::DeviceVector<Real> diff0_d(ncons);
        amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                         solverChoice.tnu2.begin(), solverChoice.tnu2.begin() + ncons,
                         diff0_d.begin());
        Real const* diff0_ptr = diff0_d.data();

        for (MFIter mfi(*vec_visc2_r[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.validbox();
            auto pm    = vec_pm[lev]->const_array(mfi);
            auto pn    = vec_pn[lev]->const_array(mfi);
            auto visc2_r = vec_visc2_r[lev]->array(mfi);
            auto diff2   = vec_diff2[lev]->array(mfi);

            int ncons_local = ncons;
            ParallelFor(makeSlab(bx,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
            {
                Real denom  = pm(i,j,0) * pn(i,j,0);
                Real grdscl = (denom > zero) ? std::sqrt(one / denom) : zero;
                visc2_r(i,j,0) = cff * grdscl;

                for (int n = 0; n < ncons_local; n++) {
                    diff2(i,j,0,n) = ((diff0_ptr[n] * lev_scale) / grdmax) * grdscl;
                }
            });
        }

        // Fill ghost cells for rho coefficients BEFORE psi averaging
        Real time = zero;
        FillPatch(lev, time, *vec_visc2_r[lev], GetVecOfPtrs(vec_visc2_r), foextrap_periodic_bc());

        // ------------------------------------------------------------
        // Step 3: Psi coefficients = average of 4 surrounding rho
        // ------------------------------------------------------------
        for (MFIter mfi(*vec_visc2_p[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.validbox();
            auto visc2_p = vec_visc2_p[lev]->array(mfi);
            auto visc2_r = vec_visc2_r[lev]->const_array(mfi);

            ParallelFor(makeSlab(bx,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
            {
                visc2_p(i,j,0) = fourth * (
                    visc2_r(i-1,j-1,0) +
                    visc2_r(i  ,j-1,0) +
                    visc2_r(i-1,j  ,0) +
                    visc2_r(i  ,j  ,0)
                );
            });
        }

        FillPatch(lev, time, *vec_visc2_p[lev], GetVecOfPtrs(vec_visc2_p), foextrap_periodic_bc());

        // Diagnostics
        // NOTE: coefficients are computed everywhere (including land). Output routines may later
        // mask land points (e.g., to FillValue in NetCDF/plotfiles), and analysis tools may
        // additionally apply mask_rho (setting land to 0). Report both conventions.
        //
        // Global (MPI-reduced) extrema over all valid cells (no ghost).
        Real visc_min_all = vec_visc2_r[lev]->min(0,0,false);
        Real visc_max_all = vec_visc2_r[lev]->max(0,0,false);

        // Global extrema over *wet* rho points only, k=0.
        amrex::Gpu::LaunchSafeGuard lsg_diag(true);
        Real visc_min_wet = amrex::ReduceMin(*vec_visc2_r[lev], *vec_mskr[lev], 0,
            [=] AMREX_GPU_HOST_DEVICE (Box const& bx,
                                      Array4<Real const> const& visc2,
                                      Array4<Real const> const& mskr) -> Real
            {
                Real local_min = bogus_large_value;
                amrex::Loop(bx, [=,&local_min] (int i, int j, int) noexcept
                {
                    if (mskr(i,j,0) > zero) {
                        local_min = amrex::min(local_min, visc2(i,j,0));
                    }
                });
                return local_min;
            });
        ParallelDescriptor::ReduceRealMin(visc_min_wet);

        Real visc_max_wet = amrex::ReduceMax(*vec_visc2_r[lev], *vec_mskr[lev], 0,
            [=] AMREX_GPU_HOST_DEVICE (Box const& bx,
                                      Array4<Real const> const& visc2,
                                      Array4<Real const> const& mskr) -> Real
            {
                Real local_max = -bogus_large_value;
                amrex::Loop(bx, [=,&local_max] (int i, int j, int) noexcept
                {
                    if (mskr(i,j,0) > zero) {
                        local_max = amrex::max(local_max, visc2(i,j,0));
                    }
                });
                return local_max;
            });
        ParallelDescriptor::ReduceRealMax(visc_max_wet);

        // Mimic "apply mask_rho" convention (dry -> 0).
        Real visc_min_mask0 = amrex::ReduceMin(*vec_visc2_r[lev], *vec_mskr[lev], 0,
            [=] AMREX_GPU_HOST_DEVICE (Box const& bx,
                                      Array4<Real const> const& visc2,
                                      Array4<Real const> const& mskr) -> Real
            {
                Real local_min = bogus_large_value;
                amrex::Loop(bx, [=,&local_min] (int i, int j, int) noexcept
                {
                    const Real v = (mskr(i,j,0) > zero) ? visc2(i,j,0) : zero;
                    local_min = amrex::min(local_min, v);
                });
                return local_min;
            });
        ParallelDescriptor::ReduceRealMin(visc_min_mask0);

        Real visc_max_mask0 = amrex::ReduceMax(*vec_visc2_r[lev], *vec_mskr[lev], 0,
            [=] AMREX_GPU_HOST_DEVICE (Box const& bx,
                                      Array4<Real const> const& visc2,
                                      Array4<Real const> const& mskr) -> Real
            {
                Real local_max = -bogus_large_value;
                amrex::Loop(bx, [=,&local_max] (int i, int j, int) noexcept
                {
                    const Real v = (mskr(i,j,0) > zero) ? visc2(i,j,0) : zero;
                    local_max = amrex::max(local_max, v);
                });
                return local_max;
            });
        ParallelDescriptor::ReduceRealMax(visc_max_mask0);
        if (ParallelDescriptor::IOProcessor() && lev == 0)
        {
            Print() << "\nHorizontal mixing scaled by grid metric\n";
            Print() << "grdmax = " << grdmax << "\n";
            if (solverChoice.scaled_to_grid_amr_scaling == ScaledToGridAMRScaling::linear) {
                Print() << "AMR scaling (linear) lev_scale = " << lev_scale << "\n";
            }
            Print() << "visc2(all)      min/max = "
                    << visc_min_all << " / "
                    << visc_max_all << "\n";
            Print() << "visc2(wet,k=0)  min/max = "
                    << visc_min_wet << " / "
                    << visc_max_wet << "\n";
            Print() << "visc2(mask->0)  min/max = "
                    << visc_min_mask0 << " / "
                    << visc_max_mask0 << "\n";
        }

    } else {
        Abort("Don't know this horizontal mixing type");
    }

    // Final FillPatch for all fields
    Real time = zero;
    FillPatch(lev, time, *vec_visc2_p[lev], GetVecOfPtrs(vec_visc2_p), foextrap_periodic_bc());
    FillPatch(lev, time, *vec_visc2_r[lev], GetVecOfPtrs(vec_visc2_r), foextrap_periodic_bc());
    for (int n = 0; n < ncons; n++) {
        FillPatch(lev, time, *vec_diff2[lev], GetVecOfPtrs(vec_diff2),
                  foextrap_periodic_bc(), BdyVars::null, n, false);
    }
}

/**
 * @param[in   ] lev    level to operate on
 */
void
REMORA::init_flat_bathymetry(int lev)
{
    BL_PROFILE("REMORA::init_flat_bathymetry()");
    vec_h[lev]->setVal(-geom[0].ProbLo()[2]);
}

/**
 * @param[in   ] lev    level to operate on
 */
void
REMORA::set_smflux(int lev)
{
    BL_PROFILE("REMORA::set_smflux()");
    if (solverChoice.smflux_type == SMFluxType::analytic) {
        prob->init_analytic_smflux(lev, geom[lev], solverChoice, *this,*vec_sustr[lev], *vec_svstr[lev]);
    } else if (solverChoice.smflux_type == SMFluxType::netcdf) {
#ifdef REMORA_USE_NETCDF
        sustr_data_from_file->update_interpolated_to_time(t_old[lev], lev, vec_sustr[lev].get(), geom, ref_ratio);
        svstr_data_from_file->update_interpolated_to_time(t_old[lev], lev, vec_svstr[lev].get(), geom, ref_ratio);
        FillPatch(lev, t_old[lev], *vec_sustr[lev], GetVecOfPtrs(vec_sustr), foextrap_periodic_bc(), BdyVars::null,0,false,false);
        FillPatch(lev, t_old[lev], *vec_svstr[lev], GetVecOfPtrs(vec_svstr), foextrap_periodic_bc(), BdyVars::null,0,false,false);
#endif
    }
}

/**
 * @param[in   ] lev    level to operate on
 */
void
REMORA::set_surface_state (int lev)
{
    BL_PROFILE("REMORA::set_surface_state()");

    auto& bulk_flux_type = solverChoice.bulk_flux_type;

    for (int n=0; n < AtmosState::NumTypes; n++) {
        if (driver_atmos_state_from_driver[n]) {
            amrex::Abort("Reached set_surface_state() but variables have already been specified from driver!");
        }
    }
//    const bool driver_has_uwind = driver_atmos_state_from_driver[AtmosState::Uwind];
//    const bool driver_has_vwind = driver_atmos_state_from_driver[AtmosState::Vwind];
//
//    const bool use_analytic_uwind = bulk_flux_type[BulkFlux::Uwind] == BulkForcingType::analytic &&
//                                    !driver_has_uwind;
//    const bool use_analytic_vwind = bulk_flux_type[BulkFlux::Vwind] == BulkForcingType::analytic &&
//                                    !driver_has_vwind;
//
//    if (use_analytic_uwind || use_analytic_vwind) {
//        std::unique_ptr<MultiFab> tmp_uwind;
//        std::unique_ptr<MultiFab> tmp_vwind;
//        MultiFab* analytic_uwind = vec_uwind[lev].get();
//        MultiFab* analytic_vwind = vec_vwind[lev].get();
//
//        if (!use_analytic_uwind) {
//            tmp_uwind.reset(new MultiFab(vec_uwind[lev]->boxArray(), vec_uwind[lev]->DistributionMap(),
//                                         1, vec_uwind[lev]->nGrowVect()));
//            analytic_uwind = tmp_uwind.get();
//        }
//        if (!use_analytic_vwind) {
//            tmp_vwind.reset(new MultiFab(vec_vwind[lev]->boxArray(), vec_vwind[lev]->DistributionMap(),
//                                         1, vec_vwind[lev]->nGrowVect()));
//            analytic_vwind = tmp_vwind.get();
//        }
//
//        prob->init_analytic_wind(lev, geom[lev], solverChoice, *this, *analytic_uwind, *analytic_vwind);
//    }

#ifdef REMORA_USE_NETCDF
    auto update_from_netcdf = [&](std::unique_ptr<NCTimeSeries>& data_from_file,
                                  Vector<std::unique_ptr<MultiFab>>& mf_vec) {
        data_from_file->update_interpolated_to_time(t_old[lev], lev, mf_vec[lev].get(), geom, ref_ratio);
        FillPatch(lev, t_old[lev], *mf_vec[lev], GetVecOfPtrs(mf_vec),
                  foextrap_periodic_bc(), BdyVars::null, 0, false);
    };

    if (bulk_flux_type[BulkFlux::Uwind] == BulkForcingType::netcdf && !driver_atmos_state_from_driver[AtmosState::Uwind]) {
        update_from_netcdf(Uwind_data_from_file, vec_uwind);
    }
    if (bulk_flux_type[BulkFlux::Vwind] == BulkForcingType::netcdf && !driver_atmos_state_from_driver[AtmosState::Vwind]) {
        update_from_netcdf(Vwind_data_from_file, vec_vwind);
    }

    if (bulk_flux_type[BulkFlux::Tair] == BulkForcingType::netcdf && !driver_atmos_state_from_driver[AtmosState::Tair]) {
        update_from_netcdf(Tair_data_from_file, vec_Tair);
    }
    if (bulk_flux_type[BulkFlux::Qair] == BulkForcingType::netcdf && !driver_atmos_state_from_driver[AtmosState::Qair]) {
        update_from_netcdf(qair_data_from_file, vec_qair);
        if (solverChoice.qair_is_percent) {
            vec_qair[lev]->mult(amrex::Real(0.01));
        }
    }
    if (bulk_flux_type[BulkFlux::Pair] == BulkForcingType::netcdf && !driver_atmos_state_from_driver[AtmosState::Pair]) {
        update_from_netcdf(Pair_data_from_file, vec_Pair);
    }
    if (bulk_flux_type[BulkFlux::SWrad] == BulkForcingType::netcdf && !driver_atmos_state_from_driver[AtmosState::SWrad]) {
        update_from_netcdf(srflx_data_from_file, vec_srflx);
    }
    if (bulk_flux_type[BulkFlux::LWrad] == BulkForcingType::netcdf && !driver_atmos_state_from_driver[AtmosState::LWrad]) {
        update_from_netcdf(longwave_down_data_from_file, vec_longwave_down);
    }
    if (bulk_flux_type[BulkFlux::Rain] == BulkForcingType::netcdf && !driver_atmos_state_from_driver[AtmosState::Rain]) {
        update_from_netcdf(rain_data_from_file, vec_rain);
    }
    if (bulk_flux_type[BulkFlux::Cloud] == BulkForcingType::netcdf && !driver_atmos_state_from_driver[AtmosState::Cloud]) {
        update_from_netcdf(cloud_data_from_file, vec_cloud);
    }
    if (bulk_flux_type[BulkFlux::EminusP] == BulkForcingType::netcdf) {
        update_from_netcdf(EminusP_data_from_file, vec_EminusP);
    }
#else
    for (int idx = 0; idx < BulkFlux::NumTypes; ++idx) {
        if (bulk_flux_type[idx] == BulkForcingType::netcdf) {
            amrex::Abort("NetCDF bulk-flux forcing requires building with NetCDF");
        }
    }
#endif

    MultiFab* analytic_uwind = (bulk_flux_type[BulkFlux::Uwind] == BulkForcingType::analytic &&
                               !driver_atmos_state_from_driver[AtmosState::Uwind]) ? vec_uwind[lev].get() : nullptr;
    MultiFab* analytic_vwind = (bulk_flux_type[BulkFlux::Vwind] == BulkForcingType::analytic &&
                               !driver_atmos_state_from_driver[AtmosState::Vwind]) ? vec_vwind[lev].get() : nullptr;
    MultiFab* analytic_Tair = (bulk_flux_type[BulkFlux::Tair] == BulkForcingType::analytic &&
                               !driver_atmos_state_from_driver[AtmosState::Tair]) ? vec_Tair[lev].get() : nullptr;
    MultiFab* analytic_qair = (bulk_flux_type[BulkFlux::Qair] == BulkForcingType::analytic &&
                               !driver_atmos_state_from_driver[AtmosState::Qair]) ? vec_qair[lev].get() : nullptr;
    MultiFab* analytic_Pair = (bulk_flux_type[BulkFlux::Pair] == BulkForcingType::analytic &&
                               !driver_atmos_state_from_driver[AtmosState::Pair]) ? vec_Pair[lev].get() : nullptr;
    MultiFab* analytic_srflx = (bulk_flux_type[BulkFlux::SWrad] == BulkForcingType::analytic &&
                                !driver_atmos_state_from_driver[AtmosState::SWrad]) ? vec_srflx[lev].get() : nullptr;
    MultiFab* analytic_lwrad = (bulk_flux_type[BulkFlux::LWrad] == BulkForcingType::analytic &&
                                !driver_atmos_state_from_driver[AtmosState::LWrad]) ? vec_longwave_down[lev].get() : nullptr;
    MultiFab* analytic_rain = (bulk_flux_type[BulkFlux::Rain] == BulkForcingType::analytic &&
                               !driver_atmos_state_from_driver[AtmosState::Rain]) ? vec_rain[lev].get() : nullptr;
    MultiFab* analytic_cloud = (bulk_flux_type[BulkFlux::Cloud] == BulkForcingType::analytic &&
                                !driver_atmos_state_from_driver[AtmosState::Cloud]) ? vec_cloud[lev].get() : nullptr;
    MultiFab* analytic_EminusP = bulk_flux_type[BulkFlux::EminusP] == BulkForcingType::analytic ? vec_EminusP[lev].get() : nullptr;

    if (analytic_uwind != nullptr || analytic_vwind != nullptr ||
        analytic_Tair != nullptr || analytic_qair != nullptr || analytic_Pair != nullptr ||
        analytic_srflx != nullptr || analytic_lwrad != nullptr || analytic_rain != nullptr ||
        analytic_cloud != nullptr || analytic_EminusP != nullptr) {
        prob->init_analytic_surface_var(lev, geom[lev], solverChoice, *this,
                                        *analytic_uwind, *analytic_vwind,
                                        *analytic_Tair, *analytic_qair, *analytic_Pair,
                                        *analytic_srflx, *analytic_lwrad, *analytic_rain,
                                        *analytic_cloud, *analytic_EminusP);
    }

    if (vec_uwind[lev] != nullptr) { vec_uwind[lev]->FillBoundary(geom[lev].periodicity()); }
    if (vec_vwind[lev] != nullptr) { vec_vwind[lev]->FillBoundary(geom[lev].periodicity()); }
    if (vec_Tair[lev] != nullptr) { vec_Tair[lev]->FillBoundary(geom[lev].periodicity()); }
    if (vec_qair[lev] != nullptr) { vec_qair[lev]->FillBoundary(geom[lev].periodicity()); }
    if (vec_Pair[lev] != nullptr) { vec_Pair[lev]->FillBoundary(geom[lev].periodicity()); }
    if (vec_srflx[lev] != nullptr) { vec_srflx[lev]->FillBoundary(geom[lev].periodicity()); }
    if (vec_longwave_down[lev] != nullptr) { vec_longwave_down[lev]->FillBoundary(geom[lev].periodicity()); }
    if (vec_rain[lev] != nullptr) { vec_rain[lev]->FillBoundary(geom[lev].periodicity()); }
    if (vec_cloud[lev] != nullptr) { vec_cloud[lev]->FillBoundary(geom[lev].periodicity()); }
    if (vec_EminusP[lev] != nullptr) { vec_EminusP[lev]->FillBoundary(geom[lev].periodicity()); }
}

/**
 * @param[in   ] lev    level to operate on
 * @param[in   ] time   current time for initialization
 */
void
REMORA::init_only (int lev, Real time)
{
    BL_PROFILE("REMORA::init_only()");
    t_new[lev] = time;
    t_old[lev] = time - bogus_large_value;

    cons_new[lev]->setVal(zero);
    xvel_new[lev]->setVal(zero);
    yvel_new[lev]->setVal(zero);
    zvel_new[lev]->setVal(zero);

    xvel_old[lev]->setVal(zero);
    yvel_old[lev]->setVal(zero);
    zvel_old[lev]->setVal(zero);

    vec_ru[lev]->setVal(zero);
    vec_rv[lev]->setVal(zero);

    vec_ru2d[lev]->setVal(zero);
    vec_rv2d[lev]->setVal(zero);

    if (solverChoice.ic_type == IC_Type::analytic) {
        set_grid_scale(lev);
    }
    set_masks(lev);

#ifdef REMORA_USE_NETCDF
    if (solverChoice.ic_type == IC_Type::netcdf) {
        init_clim_nudg_coeff(lev);

        if (solverChoice.do_any_clim_nudg && lev == 0) {
            if (nc_clim_his_file.empty() || nc_clim_his_file[0].empty()) {
                amrex::Error("NetCDF climatology file name must be provided via input");
            }
            if (solverChoice.do_m2_clim_nudg) {
                ubar_clim_data_from_file.reset(new NCTimeSeries(nc_clim_his_file, "ubar",
                            clim_ubar_time_varname, geom[lev].Domain(),vec_ubar[lev].get(),true,true));
                vbar_clim_data_from_file.reset(new NCTimeSeries(nc_clim_his_file, "vbar",
                            clim_ubar_time_varname, geom[lev].Domain(),vec_vbar[lev].get(),true,true));
                ubar_clim_data_from_file->Initialize();
                vbar_clim_data_from_file->Initialize();
            }
            if (solverChoice.do_m3_clim_nudg) {
                u_clim_data_from_file.reset(new NCTimeSeries(nc_clim_his_file, "u", clim_u_time_varname, geom[lev].Domain(),xvel_new[lev],false,true));
                v_clim_data_from_file.reset(new NCTimeSeries(nc_clim_his_file, "v", clim_v_time_varname, geom[lev].Domain(),yvel_new[lev],false,true));
                u_clim_data_from_file->Initialize();
                v_clim_data_from_file->Initialize();
            }
            // Since the NCTimeSeries object isn't filling the cons_new MultiFab directly, we don't have to specify a component.
            // It just needs to know the shape of the MultiFab
            cons_clim_data_from_file.resize(ncons);
            for (int icomp = 0; icomp < ncons; ++icomp) {
                if (!solverChoice.do_cons_clim_nudg[icomp]) { continue; }
                // A tracer's climatology is stored in the file under the tracer's own
                // name, following the same convention ROMS uses. Check up front rather
                // than letting the read fail deep inside NCTimeSeries.
                for (const auto& fname : nc_clim_his_file) {
                    if (!QueryNetCDFHasVars(fname, {cons_names[icomp]})) {
                        amrex::Abort("Climatology file " + fname + " does not contain '" +
                                     cons_names[icomp] + "', which is required by remora.do_" +
                                     cons_names[icomp] + "_clim_nudg. Either add it to the file "
                                     "or turn that flag off.");
                    }
                }
                cons_clim_data_from_file[icomp].reset(new NCTimeSeries(nc_clim_his_file, cons_names[icomp],
                            clim_cons_time_varname[icomp],geom[lev].Domain(),cons_new[lev],false,true));
                cons_clim_data_from_file[icomp]->Initialize();
            }
        }
    }

    if (solverChoice.boundary_from_netcdf) {
        amrex::Print() << "Calling init_bdry_from_netcdf at level " << lev << std::endl;
        init_bdry_from_netcdf(lev);
        amrex::Print() << "Boundary data loaded from netcdf file \n " << std::endl;
    }

    // This will be a non-op if forcings specified analytically
    if (solverChoice.smflux_type == SMFluxType::netcdf) {
        if (lev==0) {
            if (nc_frc_file.empty() || nc_frc_file[0].empty()) {
                amrex::Error("NetCDF forcing file name must be provided via input for surface momentum fluxes");
            }
            sustr_data_from_file.reset(new NCTimeSeries(nc_frc_file, "sustr", frc_time_varname, geom[lev].Domain(),vec_sustr[lev].get(), true, false));
            svstr_data_from_file.reset(new NCTimeSeries(nc_frc_file, "svstr", frc_time_varname, geom[lev].Domain(),vec_svstr[lev].get(), true, false));
            sustr_data_from_file->Initialize();
            svstr_data_from_file->Initialize();
        } else {
            FillCoarsePatch(lev, time, vec_sustr[lev].get(), vec_sustr[lev-1].get(), foextrap_bc());
            FillCoarsePatch(lev, time, vec_svstr[lev].get(), vec_svstr[lev-1].get(), foextrap_bc());
        }
    }

    // Conditionally load atmospheric forcing fields from NetCDF based on source type.
    const auto& bulk_flux_type = solverChoice.bulk_flux_type;
    bool any_bulk_netcdf = false;
    for (int idx = 0; idx < BulkFlux::NumTypes; ++idx) {
        any_bulk_netcdf = any_bulk_netcdf || bulk_flux_type[idx] == BulkForcingType::netcdf;
    }
    if (lev == 0 && any_bulk_netcdf && (nc_frc_file.empty() || nc_frc_file[0].empty())) {
        amrex::Error("NetCDF forcing file name must be provided via input for bulk-flux atmospheric forcing");
    }

    if (lev==0) {
        if (bulk_flux_type[BulkFlux::Uwind] == BulkForcingType::netcdf) {
            Uwind_data_from_file.reset(new NCTimeSeries(nc_frc_file, "Uwind", frc_time_varname, geom[lev].Domain(),vec_uwind[lev].get(), true, false));
            Uwind_data_from_file->Initialize();
        }
        if (bulk_flux_type[BulkFlux::Vwind] == BulkForcingType::netcdf) {
            Vwind_data_from_file.reset(new NCTimeSeries(nc_frc_file, "Vwind", frc_time_varname, geom[lev].Domain(),vec_vwind[lev].get(), true, false));
            Vwind_data_from_file->Initialize();
        }
        if (bulk_flux_type[BulkFlux::Tair] == BulkForcingType::netcdf) {
            Tair_data_from_file.reset(new NCTimeSeries(nc_frc_file, "Tair", frc_time_varname, geom[lev].Domain(),vec_Tair[lev].get(), true, false));
            Tair_data_from_file->Initialize();
        }
        if (bulk_flux_type[BulkFlux::Qair] == BulkForcingType::netcdf) {
            qair_data_from_file.reset(new NCTimeSeries(nc_frc_file, "qair", frc_time_varname, geom[lev].Domain(),vec_qair[lev].get(), true, false));
            qair_data_from_file->Initialize();
        }
        if (bulk_flux_type[BulkFlux::Pair] == BulkForcingType::netcdf) {
            Pair_data_from_file.reset(new NCTimeSeries(nc_frc_file, "Pair", frc_time_varname, geom[lev].Domain(),vec_Pair[lev].get(), true, false));
            Pair_data_from_file->Initialize();
        }
        if (bulk_flux_type[BulkFlux::SWrad] == BulkForcingType::netcdf) {
            srflx_data_from_file.reset(new NCTimeSeries(nc_frc_file, "swrad", frc_time_varname, geom[lev].Domain(),vec_srflx[lev].get(), true, false));
            srflx_data_from_file->Initialize();
        }
        if (bulk_flux_type[BulkFlux::Rain] == BulkForcingType::netcdf) {
            rain_data_from_file.reset(new NCTimeSeries(nc_frc_file, "rain", frc_time_varname, geom[lev].Domain(),vec_rain[lev].get(), true, false));
            rain_data_from_file->Initialize();
        }
        if (bulk_flux_type[BulkFlux::Cloud] == BulkForcingType::netcdf) {
            cloud_data_from_file.reset(new NCTimeSeries(nc_frc_file, "cloud", frc_time_varname, geom[lev].Domain(),vec_cloud[lev].get(), true, false));
            cloud_data_from_file->Initialize();
        }
        if (bulk_flux_type[BulkFlux::EminusP] == BulkForcingType::netcdf) {
            EminusP_data_from_file.reset(new NCTimeSeries(nc_frc_file, "EminusP", frc_time_varname, geom[lev].Domain(),vec_EminusP[lev].get(), true, false));
            EminusP_data_from_file->Initialize();
        }
        if (bulk_flux_type[BulkFlux::LWrad] == BulkForcingType::netcdf) {
            longwave_down_data_from_file.reset(new NCTimeSeries(nc_frc_file, solverChoice.longwave_netcdf_varname, frc_time_varname,
                                                                geom[lev].Domain(), vec_longwave_down[lev].get(), true, false));
            longwave_down_data_from_file->Initialize();
        }
    } else {
        if (bulk_flux_type[BulkFlux::Uwind] == BulkForcingType::netcdf) {
            FillCoarsePatch(lev, time, vec_uwind[lev].get(), vec_uwind[lev-1].get(), foextrap_bc());
        }
        if (bulk_flux_type[BulkFlux::Vwind] == BulkForcingType::netcdf) {
            FillCoarsePatch(lev, time, vec_vwind[lev].get(), vec_vwind[lev-1].get(), foextrap_bc());
        }
        if (bulk_flux_type[BulkFlux::Tair] == BulkForcingType::netcdf) {
            FillCoarsePatch(lev, time, vec_Tair[lev].get(), vec_Tair[lev-1].get(), foextrap_bc());
        }
        if (bulk_flux_type[BulkFlux::Qair] == BulkForcingType::netcdf) {
            FillCoarsePatch(lev, time, vec_qair[lev].get(), vec_qair[lev-1].get(), foextrap_bc());
        }
        if (bulk_flux_type[BulkFlux::Pair] == BulkForcingType::netcdf) {
            FillCoarsePatch(lev, time, vec_Pair[lev].get(), vec_Pair[lev-1].get(), foextrap_bc());
        }
        if (bulk_flux_type[BulkFlux::SWrad] == BulkForcingType::netcdf) {
            FillCoarsePatch(lev, time, vec_srflx[lev].get(), vec_srflx[lev-1].get(), foextrap_bc());
        }
        if (bulk_flux_type[BulkFlux::Rain] == BulkForcingType::netcdf) {
            FillCoarsePatch(lev, time, vec_rain[lev].get(), vec_rain[lev-1].get(), foextrap_bc());
        }
        if (bulk_flux_type[BulkFlux::Cloud] == BulkForcingType::netcdf) {
            FillCoarsePatch(lev, time, vec_cloud[lev].get(), vec_cloud[lev-1].get(), foextrap_bc());
        }
        if (bulk_flux_type[BulkFlux::EminusP] == BulkForcingType::netcdf) {
            FillCoarsePatch(lev, time, vec_EminusP[lev].get(), vec_EminusP[lev-1].get(), foextrap_bc());
        }
        if (bulk_flux_type[BulkFlux::LWrad] == BulkForcingType::netcdf) {
            FillCoarsePatch(lev, time, vec_longwave_down[lev].get(), vec_longwave_down[lev-1].get(), foextrap_bc());
        }
    }

    // Only need to read in rivers on level 0
    // Will need to be on higher levels eventually
    if (solverChoice.do_rivers) {
        if (nc_riv_file.empty() || nc_riv_file[0].empty()) {
            amrex::Error("NetCDF river file name must be provided via input for rivers");
        }
        auto dom = geom[0].Domain();
        int nz = dom.length(2);
        // Every cell-centered tracer can take river input. The field is named for the
        // tracer, as in ROMS: river_temp, river_salt, river_tracer, river_NO3, ...
        river_source_cons.resize(ncons);
        for (int icomp = 0; icomp < ncons; ++icomp) {
            if (!solverChoice.do_rivers_cons[icomp]) { continue; }

            const std::string field = "river_" + cons_names[icomp];
            for (const auto& fname : nc_riv_file) {
                if (!QueryNetCDFHasVars(fname, {field})) {
                    // The flag may have come from remora.do_rivers_scalar rather than the
                    // per-tracer key, so name the tracer and the key that switches it off.
                    amrex::Abort("River file " + fname + " does not contain '" + field +
                                 "', but river input is enabled for tracer '" +
                                 cons_names[icomp] + "'. Either add that variable to the "
                                 "file, or set remora.do_rivers_" + cons_names[icomp] +
                                 " = false.");
                }
            }

            river_source_cons[icomp].reset(new NCTimeSeriesRiver(nc_riv_file, field, riv_time_varname, nz));
            river_source_cons[icomp]->Initialize();
        }
        river_source_transport.reset(new NCTimeSeriesRiver(nc_riv_file, "river_transport", riv_time_varname, nz));
        river_source_transport->Initialize();
        river_source_transportbar.reset(new NCTimeSeriesRiver(nc_riv_file, "river_transport", riv_time_varname, nz, 1));
        river_source_transportbar->Initialize();
        init_riv_pos_from_netcdf(lev);
    }

    if (lev==0 and hires_grid_level > 0 and solverChoice.ic_type == IC_Type::netcdf) {
        amrex::Print() << "Reading high resolution bathymetry and grid data" << std::endl;
        allocate_bathymetry_grid_vars_full_domain();
        init_bathymetry_full_domain_from_netcdf();
        init_grid_vars_full_domain_from_netcdf();
        amrex::Print() << "Done reading in high resolution bathymetry and grid data" << std::endl;
    }
    if (lev==0 and hires_init_level > 0 and solverChoice.ic_type == IC_Type::netcdf) {
        amrex::Print() << "Reading high resolution initial data" << std::endl;
        allocate_init_full_domain();
        init_data_full_domain_from_netcdf();
        // Biology source is chosen by remora.biology_ic_type, not by ic_type,
        // so this goes through the same dispatcher as the per-level path.
        // Must follow init_data_full_domain_from_netcdf: the analytic biology
        // profiles read temperature.
        init_biology_ic_full_domain();
        init_zeta_full_domain_from_netcdf();
        amrex::Print() << "Done reading in high resolution initial data" << std::endl;
    }
#else
    if (solverChoice.ic_type == IC_Type::netcdf) {
        Abort("Not compiled with NetCDF, but remora.ic_type = netcdf reads initial and grid data from file");
    }
    // No guard on hires_grid_level here: with analytic initialization it needs no NetCDF at all --
    // the bathymetry comes from prob->init_analytic_bathymetry evaluated at the fine level and
    // averaged down. hires_init_level needs no guard either: it is rejected for analytic
    // initialization in ReadParameters, and the netcdf case is caught just above.
    if (solverChoice.boundary_from_netcdf) {
        Abort("Not compiled with NetCDF, but selected boundary conditions require NetCDF");
    }
    if (solverChoice.do_rivers) {
        Abort("Not compiled with NetCDF, but using river sources requires NetCDF");
    }
#endif

    if (lev==0 and hires_grid_level > 0 and solverChoice.ic_type == IC_Type::analytic) {
        allocate_bathymetry_grid_vars_full_domain();
        init_bathymetry_full_domain_from_analytic();
    }

    if (lev==0 and hires_init_level > 0 and solverChoice.ic_type == IC_Type::analytic) {
        allocate_init_full_domain();
        init_full_domain_zeta_from_analytic();
    }

    set_bathymetry(lev);
    set_zeta(lev);
    stretch_transform(lev);

    if (lev==0 and hires_init_level > 0 and solverChoice.ic_type == IC_Type::analytic) {
        init_full_domain_from_analytic();
    }

    if (lev==0) {
        if (hires_init_level < 0) {
            if (solverChoice.ic_type == IC_Type::analytic) {
                init_analytic(lev);
            } else if (solverChoice.ic_type == IC_Type::netcdf) {
#ifdef REMORA_USE_NETCDF
                amrex::Print() << "Calling init_data_from_netcdf " << std::endl;
                init_data_from_netcdf(lev);
                bool apply_eminusp = false;
                set_zeta_to_Ztavg(lev, apply_eminusp);
                amrex::Print() << "Initial data loaded from netcdf file \n " << std::endl;
#endif
            } else {
                amrex::Abort("Unknown IC_Type");
            }
            // Biology last, and outside the ic_type branches: its source is
            // chosen independently by remora.biology_ic_type, and the analytic
            // profiles need the physical fields already in place.
            init_biology_ic(lev);
        } else {
            set_init_data_averaged_down(lev); // also sets biology data
            bool apply_eminusp = false;
            set_zeta_to_Ztavg(lev, apply_eminusp);
            // Since set_grid_scale is usually called from init_analytic for analytic problems
            if (solverChoice.ic_type == IC_Type::analytic) {
                set_grid_scale(lev);
            }
        }
    } else {
        if (lev > hires_init_level) {
            FillCoarsePatch(lev, time, cons_new[lev], cons_new[lev-1],BCVars::Temp_bc_comp,BdyVars::t);
            FillCoarsePatch(lev, time, xvel_new[lev], xvel_new[lev-1], xvel_bc(), BdyVars::u);
            FillCoarsePatch(lev, time, yvel_new[lev], yvel_new[lev-1], yvel_bc(), BdyVars::v);
            FillCoarsePatch(lev, time, zvel_new[lev], zvel_new[lev-1], zvel_bc(), BdyVars::null);
        } else {
            set_init_data_averaged_down(lev); // also sets biology data
            bool apply_eminusp = false;
            set_zeta_to_Ztavg(lev, apply_eminusp);
            if (solverChoice.ic_type == IC_Type::analytic) {
                // Since set_grid_scale is usually called from init_analytic for analytic problems
                set_grid_scale(lev);
            }
        }
    }

    // Ensure that the face-based data are the same on both sides of a periodic domain.
    // The data associated with the lower grid ID is considered the correct value.
    xvel_new[lev]->OverrideSync(geom[lev].periodicity());
    yvel_new[lev]->OverrideSync(geom[lev].periodicity());
    zvel_new[lev]->OverrideSync(geom[lev].periodicity());

    set_2darrays(lev);

    init_set_vmix(lev);
    set_hmixcoef(lev);
    set_coriolis(lev);

    // Previously set smflux here with OverrideSync:
//    set_smflux(lev);
//    prob->init_analytic_smflux(lev, geom[lev], solverChoice, *this, *vec_sustr[lev], *vec_svstr[lev]);
//    vec_sustr[lev]->OverrideSync(geom[lev].periodicity());
//    vec_svstr[lev]->OverrideSync(geom[lev].periodicity());

}

void
REMORA::ReadParameters ()
{
    BL_PROFILE("REMORA::ReadParameters()");
    {
        ParmParse pp;  // Traditionally, max_step and stop_time do not have prefix, so allow it for now.
        bool noprefix_max_step = pp.queryAdd("max_step", max_step);
        bool noprefix_stop_time = pp.queryAdd("stop_time", stop_time);
        bool remora_max_step = pp.queryAdd("remora.max_step", max_step);
        bool remora_stop_time = pp.queryAdd("remora.stop_time", stop_time);
        if (remora_max_step and noprefix_max_step) {
            Abort("remora.max_step and max_step are both specified. Please use only one!");
        }
        if (remora_stop_time and noprefix_stop_time) {
            Abort("remora.stop_time and stop_time are both specified. Please use only one!");
        }
    }

    ParmParse pp(pp_prefix);

    // Common physics and simulation parameters
    std::string biology_model_string = REMORABiology::biology_model_name(biology_model);
    pp.queryAdd("biology_model", biology_model_string);
    biology_model = REMORABiology::parse_biology_model(biology_model_string);

    if (REMORABiology::has_biology(biology_model)) {
        fennel_params.init_params(pp_prefix);

        // Source of the biology initial condition, independent of ic_type.
        // Default "follow" reproduces the previous behaviour exactly.
        std::string biology_ic_string = REMORABiology::biology_ic_type_name(biology_ic_type);
        pp.queryAdd("biology_ic_type", biology_ic_string);
        biology_ic_type = REMORABiology::parse_biology_ic_type(biology_ic_string);

        // Bridge-vs-native selection and diagnostic verbosity are runtime
        // controls so a parity comparison never requires a rebuild. Both
        // parse unconditionally; without USE_FENNEL_FORT there is no bridge
        // to select, so asking for it is an error rather than a silent
        // fallback to the path being validated.
        pp.queryAdd("use_biology_cpp_answer", use_biology_cpp_answer);
        pp.queryAdd("biology_debug", biology_debug);
        pp.queryAdd("biology_debug_i", biology_debug_i);
        pp.queryAdd("biology_debug_j", biology_debug_j);
#ifndef REMORA_USE_FENNEL_FORT
        if (use_biology_cpp_answer == 0) {
            amrex::Abort("remora.use_biology_cpp_answer = 0 selects the ROMS "
                         "Fennel Fortran bridge, which is not compiled in. "
                         "Rebuild with USE_FENNEL_FORT=TRUE (GNUmake) or "
                         "-DREMORA_ENABLE_FENNEL_FORT=ON (CMake).");
        }
#endif
#ifndef REMORA_USE_BIOLOGY_DIAG
        if (biology_debug > 0) {
            amrex::Abort("remora.biology_debug > 0 requests the Fennel parity "
                         "diagnostics, which are not compiled in. Rebuild with "
                         "USE_BIOLOGY_DIAG=TRUE (GNUmake) or "
                         "-DREMORA_ENABLE_BIOLOGY_DIAG=ON (CMake).");
        }
#endif

        // Biology tracers are counted separately from the passive scalars, so a run can
        // carry dye and biology at once.
        nbio = static_cast<int>(REMORABiology::tracer_names(biology_model, fennel_params).size());
    } else {
        nbio = 0;
    }
    // Dye is opt-in, biology or not: a component nothing asked for is one more thing to
    // advect, diffuse, and explain in every plotfile and boundary file.
    pp.queryAdd("nscalar", nscalar);
    if (nscalar < 0) {
        amrex::Abort("remora.nscalar must be non-negative");
    }
    Bio_comp = Tracer_comp + nscalar;
    ncons = Tracer_comp + nscalar + nbio;
    init_scalar_metadata();

    // remora.nscalar used to be required to equal the biology tracer count; it now counts
    // dye only and adds to it. Print the layout so a run carrying both is unmistakable,
    // and an input written against the old meaning is caught by eye rather than by a
    // surprising component count much later.
    if (nbio > 0 && nscalar > 0) {
        amrex::Print() << "Carrying " << nscalar << " passive scalar(s) and " << nbio
                       << " biology tracer(s), for " << ncons << " cell-centered components: ";
        for (int icomp = 0; icomp < ncons; ++icomp) {
            amrex::Print() << cons_names[icomp] << (icomp + 1 < ncons ? " " : "\n");
        }
    }

    pp.queryAdd("check_file", check_file);
    pp.queryAdd("check_int", check_int);
    pp.queryAdd("check_int_time", check_int_time);
    pp.queryAdd("expand_plotvars_to_unif_rr", expand_plotvars_to_unif_rr);
    pp.query("plotfile_fill_value", plotfile_fill_value);
    pp.query("netcdf_fill_value", netcdf_fill_value);
    pp.queryAdd("restart", restart_chkfile);
    pp.queryAdd("start_time", start_time);

    num_boxes_at_level.resize(max_level + 1, 0);
    boxes_at_level.resize(max_level + 1);
    num_boxes_at_level[0] = 1;
    boxes_at_level[0].resize(1);
    boxes_at_level[0][0] = geom[0].Domain();

    if (pp.contains("data_log")) {
        int num_datalogs = pp.countval("data_log");
        datalog.resize(num_datalogs);
        datalogname.resize(num_datalogs);
        pp.queryarr("data_log", datalogname, 0, num_datalogs);
        for (int i = 0; i < num_datalogs; i++)
            setRecordDataInfo(i, datalogname[i]);
    }

    pp.queryAdd("v", verbose);
    pp.queryAdd("sum_interval", sum_interval);
    pp.queryAdd("sum_period", sum_per);
    pp.queryAdd("file_min_digits", file_min_digits);

    if (file_min_digits < 0) {
        amrex::Abort("remora.file_min_digits must be non-negative");
    }

    pp.queryAdd("cfl", cfl);
    pp.queryAdd("change_max", change_max);
    pp.queryAdd("fixed_dt", fixed_dt);
    pp.queryAdd("fixed_fast_dt", fixed_fast_dt);
    pp.queryAdd("fixed_ndtfast_ratio", fixed_ndtfast_ratio);

    if (fixed_dt > zero && fixed_fast_dt > zero && fixed_ndtfast_ratio > 0) {
        if (fixed_dt / fixed_fast_dt != fixed_ndtfast_ratio) {
            amrex::Abort("Dt is over-specfied");
        }
    } else if (fixed_dt > zero && fixed_fast_dt > zero && fixed_ndtfast_ratio <= 0) {
        fixed_ndtfast_ratio = static_cast<int>(fixed_dt / fixed_fast_dt);
    }
    AMREX_ASSERT(cfl > zero || fixed_dt > zero);

    num_files_at_level.resize(max_level + 1, 0);
    num_boxes_at_level.resize(max_level + 1, 0);
    boxes_at_level.resize(max_level + 1);
    num_boxes_at_level[0] = 1;
    boxes_at_level[0].resize(1);
    boxes_at_level[0][0] = geom[0].Domain();

    pp.queryAdd("plot_file", plot_file_name);
    pp.queryAdd("plot_int", plot_int);
    pp.queryAdd("plot_int_time", plot_int_time);
    pp.query("plot_staggered_vels", plot_staggered_vels);
    pp.query("plot_nodal_data", plot_nodal_data);

    std::string plotfile_type_str = "amrex";
    pp.queryAdd("plotfile_type", plotfile_type_str);
    if (plotfile_type_str == "amrex") {
        plotfile_type = PlotfileType::amrex;
    } else if (plotfile_type_str == "netcdf" || plotfile_type_str == "NetCDF") {
        plotfile_type = PlotfileType::netcdf;
#ifdef REMORA_USE_NETCDF
        pp.queryAdd("write_history_file",write_history_file);
        pp.queryAdd("chunk_history_file",chunk_history_file);
        pp.queryAdd("steps_per_history_file",steps_per_history_file);
        // Estimate size of domain for one timestep of netcdf
        auto dom = geom[0].Domain();
        int nx = dom.length(0) + 2;
        int ny = dom.length(1) + 2;
        int nz = dom.length(2);
        Real two_gb = Real(1.6e10);
        Real double_bits = Real(64.0);
        if (write_history_file and chunk_history_file and (steps_per_history_file <= 0)) {
            // Estimate number of steps that will fit into a 2GB file.
            steps_per_history_file = int((two_gb - NCH2D * nx * ny * double_bits)
                    / (nx * ny * double_bits * (NC3D*nz + NC2D)));
            // If we calculate that a single step will exceed 2GB and the user has
            // requested automatic history file sizing, warn about a possible impending
            // error, and set steps_per_history_file = 1 to attempt output anyway.
            if (steps_per_history_file == 0) {
                amrex::Warning("NetCDF output for a single timestep appears to exceed 2GB. NetCDF output may not work. See Documentation for information about tested MPICH versions.");
                steps_per_history_file = 1;
            }
        } else if (write_history_file and !chunk_history_file) {
            // Estimate number of output steps we'll need
            int nt_out = int((max_step) / plot_int) + 1;
            Real est_hist_file_size = NCH2D * nx * ny * double_bits + nt_out * nx * ny * double_bits * (NC3D*nz + NC2D);
            if (est_hist_file_size > two_gb) {
                amrex::Warning("WARNING: NetCDF history file may be larger than 2GB limit. Consider setting remora.chunk_history_file=true");
            }
        }
        if (write_history_file and chunk_history_file) {
            Print() << "NetCDF history files will have " << steps_per_history_file << " steps per file." << std::endl;
        }
#endif
    } else {
        amrex::Print() << "User selected plotfile_type = " << plotfile_type_str << std::endl;
        amrex::Abort("Dont know this plotfile_type");
    }
#ifndef REMORA_USE_NETCDF
    if (plotfile_type == PlotfileType::netcdf)
    {
        amrex::Abort("Please compile with NetCDF in order to enable NetCDF plotfiles");
    }

#endif
#ifdef REMORA_USE_NETCDF
    nc_init_file.resize(max_level+1);
    nc_grid_file.resize(max_level+1);
    num_files_at_level.resize(max_level + 1, 0);

    boundary_series.resize(max_level+1);


    // NetCDF initialization files -- possibly multiple files at each of multiple levels
    //        but we always have exactly one file at level 0
    for (int lev = 0; lev <= max_level; lev++)
    {
        const std::string nc_file_names = amrex::Concatenate("nc_init_file_",lev,1);
        const std::string nc_bathy_file_names = amrex::Concatenate("nc_grid_file_",lev,1);

        if (pp.contains(nc_file_names.c_str()))
        {
            int num_files = pp.countval(nc_file_names.c_str());
            int num_bathy_files = pp.countval(nc_bathy_file_names.c_str());
            if (num_files != num_bathy_files) {
                amrex::Error("Must have same number of netcdf files for grid info as for solution");
            }

            num_files_at_level[lev] = num_files;
            nc_init_file[lev].resize(num_files);
            nc_grid_file[lev].resize(num_files);

            pp.queryarr(nc_file_names.c_str()      , nc_init_file[lev]     ,0,num_files);
            pp.queryarr(nc_bathy_file_names.c_str(), nc_grid_file[lev],0,num_files);
        }
    }

    pp.queryAdd("nc_grid_file_hires", nc_grid_file_hires);
    pp.queryAdd("nc_init_file_hires", nc_init_file_hires);

    // We only read boundary data at level 0
    pp.queryarr("nc_bdry_file", nc_bdry_file);

    // Also only read forcings at level 0 (for now)
    if (pp.contains("nc_frc_file")) {
        int num_files = pp.countval("nc_frc_file");
        nc_frc_file.resize(num_files);
        pp.queryarr("nc_frc_file", nc_frc_file, 0, num_files);
    }

    // Get river file
    if (pp.contains("nc_river_file")) {
        int num_files = pp.countval("nc_river_file");
        nc_riv_file.resize(num_files);
        pp.queryarr("nc_river_file", nc_riv_file, 0, num_files);
    }

    // Read in file names for climatology history and nudging weights
    if (pp.contains("nc_clim_his_file")) {
        int num_files = pp.countval("nc_clim_his_file");
        nc_clim_his_file.resize(num_files);
        pp.queryarr("nc_clim_his_file", nc_clim_his_file, 0, num_files);
    }
    pp.queryAdd("nc_clim_coeff_file", nc_clim_coeff_file);

    for (int i=0; i<BdyVars::NumTypes(ncons); i++) {
        bdry_time_name_byvar.push_back("");
    }
    pp.queryAdd("bdy_time_varname",bdry_time_varname);
    // Every tracer takes its time-axis name from its own variable name, so temp and salt
    // keep bdy_temp_time_varname / bdy_salt_time_varname and a biology tracer uses e.g.
    // bdy_NO3_time_varname
    for (int icomp = 0; icomp < ncons; ++icomp) {
        pp.queryAdd(("bdy_"+cons_names[icomp]+"_time_varname").c_str(),
                    bdry_time_name_byvar[BdyVars::cons(icomp)]);
    }
    pp.queryAdd("bdy_u_time_varname",bdry_time_name_byvar[BdyVars::u]);
    pp.queryAdd("bdy_v_time_varname",bdry_time_name_byvar[BdyVars::v]);
    pp.queryAdd("bdy_ubar_time_varname",bdry_time_name_byvar[BdyVars::ubar(ncons)]);
    pp.queryAdd("bdy_vbar_time_varname",bdry_time_name_byvar[BdyVars::vbar(ncons)]);
    pp.queryAdd("bdy_zeta_time_varname",bdry_time_name_byvar[BdyVars::zeta(ncons)]);

    // If not specified per variable, populate with the default
    for (int i=0; i<BdyVars::NumTypes(ncons); i++) {
        if (bdry_time_name_byvar[i] == "") {
            bdry_time_name_byvar[i] = bdry_time_varname;
        }
    }

    pp.queryAdd("frc_time_varname",frc_time_varname);

    pp.queryAdd("riv_time_varname",riv_time_varname);

    pp.queryAdd("clim_ubar_time_varname",clim_ubar_time_varname);
    pp.queryAdd("clim_vbar_time_varname",clim_vbar_time_varname);
    pp.queryAdd("clim_u_time_varname",clim_u_time_varname);
    pp.queryAdd("clim_v_time_varname",clim_v_time_varname);
    // As for the boundary data, each tracer takes its climatology time-axis name from its
    // own variable name, so temp and salt keep clim_temp_time_varname /
    // clim_salt_time_varname and a biology tracer uses e.g. clim_NO3_time_varname
    clim_cons_time_varname.assign(ncons, "ocean_time");
    for (int icomp = 0; icomp < ncons; ++icomp) {
        pp.queryAdd(("clim_"+cons_names[icomp]+"_time_varname").c_str(),
                    clim_cons_time_varname[icomp]);
    }

#endif
    // A hires level of 0 is not "level 0 is the hires level", it is a null pointer: the
    // full-domain arrays are only allocated for lev > 0 (see allocate_init_full_domain and
    // allocate_bathymetry_grid_vars_full_domain), while every consumer branch tests < 0 and
    // so would take the averaged-down path against an unallocated MultiFab. -1 means off.
    pp.queryAdd("hires_grid_level", hires_grid_level);
    if (hires_grid_level > max_level) {
        amrex::Abort("hires_grid_level must be less than or equal to amr.max_level");
    }
    if (hires_grid_level == 0) {
        amrex::Abort("hires_grid_level must be greater than 0; use -1 to specify grid data at level 0");
    }
    pp.queryAdd("hires_init_level", hires_init_level);
    if (hires_init_level > max_level) {
        amrex::Abort("hires_init_level must be less than or equal to amr.max_level");
    }
    if (hires_init_level == 0) {
        amrex::Abort("hires_init_level must be greater than 0; use -1 to specify initial data at level 0");
    }
#ifdef REMORA_USE_PARTICLES
    readTracersParams();
#endif

    {
        ParmParse pp_amr("amr");
        pp_amr.queryAdd("regrid_int", regrid_int);
        pp_amr.queryAdd("do_substep", do_substep);
        if (do_substep) {
            amrex::Abort("Time substepping is not yet implemented. amr.do_substep must be 0");
        }

    }
    solverChoice.init_params(ncons, nscalar, cons_names);

    // The biology IC source is chosen independently of ic_type, but only one of the two
    // mixed combinations works: NetCDF physics with analytic biology. The reverse has no
    // file to read from -- nc_init_file is only populated on the netcdf path -- so catch
    // it here instead of failing inside PnetCDF on an empty file name.
    if (REMORABiology::has_biology(biology_model) and
        biology_ic_type == REMORABiology::BiologyICType::netcdf and
        solverChoice.ic_type != IC_Type::netcdf) {
        amrex::Abort("remora.biology_ic_type = netcdf requires remora.ic_type = netcdf: the biology "
                     "initial data is read from the same files as the physical initial data, and no "
                     "such file is given for analytic initial conditions. Use "
                     "remora.biology_ic_type = analytic (or follow) instead.");
    }

    // NOTE: This feature is not yet implemented because it will require passing x,y,z to prob functions.
    // Currently these are accessed by passing a pointer to the REMORA class. However, this requires the
    // coordinates at hires_init_level to already exist (and specifically for the hires_init_level level
    // to already be initialized), which is generally not the case. A solution is to create a separate
    // coordinates object that is passed to the prob functions instead of the REMORA object. Then x,y,z
    // coordinates can be calculated at any level without the corresponding level having been created.
    if (hires_init_level >= 0 and solverChoice.ic_type == IC_Type::analytic) {
        amrex::Abort("Cannot do high-resolution initialization for analytic initial conditions. Not yet implemented");
    }

    // flat_bathymetry takes precedence at level 0 (see set_bathymetry), so a run asking for
    // both would silently ignore the high-resolution bathymetry it was given. The flat branch
    // of several analytic bathymetry hooks also writes h's second component, which the
    // one-component vec_h_full_domain does not have.
    if (hires_grid_level > 0 and solverChoice.flat_bathymetry) {
        amrex::Abort("remora.flat_bathymetry is incompatible with hires_grid_level > 0: flat bathymetry "
                     "would override the high-resolution bathymetry at level 0. Use one or the other");
    }

}


void
REMORA::AverageDown ()
{
    BL_PROFILE("REMORA::AverageDown()");
    for (int lev = finest_level-1; lev >= 0; --lev)
    {
        AverageDownTo(lev);
    }
}

/**
 * @param[in   ] crse_lev  level to average down to
 */
void
REMORA::AverageDownTo (int crse_lev)
{
    BL_PROFILE("REMORA::AverageDownTo()");
    average_down(*cons_new[crse_lev+1], *cons_new[crse_lev],
                 0, cons_new[crse_lev]->nComp(), refRatio(crse_lev));
    average_down(*vec_Zt_avg1[crse_lev+1].get(), *vec_Zt_avg1[crse_lev].get(),
                 0, vec_Zt_avg1[crse_lev]->nComp(), refRatio(crse_lev));

    Array<MultiFab*,AMREX_SPACEDIM>  faces_crse;
    Array<MultiFab*,AMREX_SPACEDIM>  faces_fine;
    faces_crse[0] = xvel_new[crse_lev];
    faces_crse[1] = yvel_new[crse_lev];
    faces_crse[2] = zvel_new[crse_lev];

    faces_fine[0] = xvel_new[crse_lev+1];
    faces_fine[1] = yvel_new[crse_lev+1];
    faces_fine[2] = zvel_new[crse_lev+1];

    average_down_faces(GetArrOfConstPtrs(faces_fine), faces_crse,
                       refRatio(crse_lev),geom[crse_lev]);
    stretch_transform(crse_lev);
}

/**
 * @param[in   ] crse_lev   level to average data down to
 * @param[inout] vec_mf     vector over levels of multifabs containing data to average
 */
void
REMORA::average_down_with_grow_cells (int crse_lev, Vector<std::unique_ptr<MultiFab>>& vec_mf)
{
    auto const& crsema = vec_mf[crse_lev]->arrays();
    auto const& finema = vec_mf[crse_lev+1]->const_arrays();
    auto ref_ratio_crse = refRatio(crse_lev);
    auto index_type = (vec_mf[crse_lev]->boxArray().ixType()).toIntVect();
    auto nghost_crse = cum_ref_ratios[crse_lev] - index_type;
    if (index_type[0]==0 and index_type[1]==0) {
        ParallelFor(*vec_mf[crse_lev], nghost_crse, vec_mf[crse_lev]->nComp(),
                [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
        {
            amrex_avgdown(i,j,k,n,crsema[box_no],finema[box_no],0,0,ref_ratio_crse);
        });
    } else if (index_type[0]==1 and index_type[1]==0) {
        ParallelFor(*vec_mf[crse_lev], nghost_crse, 1,
                [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
        {
            amrex_avgdown_faces(i,j,k,n,crsema[box_no],finema[box_no],0,0,ref_ratio_crse,0);
        });
    } else if (index_type[0]==0 and index_type[1]==1) {
        ParallelFor(*vec_mf[crse_lev], nghost_crse, 1,
                [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
        {
            amrex_avgdown_faces(i,j,k,n,crsema[box_no],finema[box_no],0,0,ref_ratio_crse,1);
        });
    } else {
        amrex::Abort("Unexpected nodality in average_down_with_grow_cells");
    }
    Gpu::streamSynchronize();
}

/**
 * @param[in   ] lev    level at which to get time
 */
amrex::Real REMORA::get_t_old(int lev) const
{
    return t_old[lev];
}
