/**
 * \file REMORA.cpp
 */

#include <REMORA_prob_common.H>
#include <REMORA.H>

#include <AMReX_buildInfo.H>

using namespace amrex;

amrex::Real REMORA::startCPUTime        = 0.0_rt;
amrex::Real REMORA::previousCPUTimeUsed = 0.0_rt;

Vector<AMRErrorTag> REMORA::ref_tags;

SolverChoice REMORA::solverChoice;

// Time step control
amrex::Real REMORA::cfl           =  0.8_rt;
amrex::Real REMORA::fixed_dt      = -1.0_rt;
amrex::Real REMORA::fixed_fast_dt = -1.0_rt;
amrex::Real REMORA::change_max    =  1.1_rt;

int   REMORA::fixed_ndtfast_ratio = 0;

// Dictate verbosity in screen output
int         REMORA::verbose       = 0;

// Frequency of diagnostic output
int         REMORA::sum_interval  = -1;
amrex::Real REMORA::sum_per       = -1.0_rt;

// Minimum number of digits in plotfile name
int         REMORA::file_min_digits = 5;

// Do we include staggered velocities in the plotfile?
int         REMORA::plot_staggered_vels = 0;

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

amrex::Vector<std::string> BCNames = {"xlo", "ylo", "zlo", "xhi", "yhi", "zhi"};

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
    const std::string& pv1 = "plot_vars"; setPlotVariables(pv1);

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

    t_new.resize(nlevs_max, 0.0_rt);
    t_old.resize(nlevs_max, -1.e100_rt);
    dt.resize(nlevs_max, 1.e100_rt);

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

REMORA::~REMORA ()
{
}

void
REMORA::Evolve ()
{
    BL_PROFILE("REMORA::Evolve()");
    Real cur_time = t_new[0];

    // Take one coarse timestep by calling timeStep -- which recursively calls timeStep
    //      for finer levels (with or without subcycling)
    for (int step = istep[0]; step < max_step && cur_time < stop_time; ++step)
    {
        amrex::Print() << "\nCoarse STEP " << step+1 << " starts ..." << std::endl;

        ComputeDt();

        int lev = 0;
        int iteration = 1;

        if (max_level == 0) {
            timeStep(lev, cur_time, iteration);
        }
        else {
            timeStepML(cur_time, iteration);
        }

        cur_time  += dt[0];

        amrex::Print() << "Coarse STEP " << step+1 << " ends." << " TIME = " << cur_time
                       << " DT = " << dt[0]  << std::endl;

        if ((plot_int > 0 && (step+1 - last_plot_file_step) == plot_int)
                || (plot_int_time > 0 && cur_time >= (last_plot_file_time + plot_int_time))) {
            last_plot_file_step = step+1;
            last_plot_file_time = cur_time;
            if (plotfile_type == PlotfileType::amrex) {

                WritePlotFile();
            }
#ifdef REMORA_USE_NETCDF
            else if (plotfile_type == PlotfileType::netcdf) {
                WriteNCPlotFile(step+1);
                history_count++;
            }
#endif
        }

        if ((check_int > 0 && (step+1 - last_check_file_step) == check_int)
                || (check_int_time > 0 && cur_time >= (last_check_file_time + check_int_time))) {
            last_check_file_step = step+1;
            last_check_file_time = cur_time;
            WriteCheckpointFile();
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

    if ((plot_int > 0 || plot_int_time > 0.0) && istep[0] > last_plot_file_step) {
        if (plotfile_type == PlotfileType::amrex) {
            WritePlotFile();
        }
#ifdef REMORA_USE_NETCDF
        if (plotfile_type == PlotfileType::netcdf) {
            WriteNCPlotFile(istep[0]);
            history_count++;
        }
#endif
    }

    if ((check_int > 0 || check_int_time > 0.0) && istep[0] > last_check_file_step) {
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
    // Initialize the start time for our CPU-time tracker
    startCPUTime = Real(ParallelDescriptor::second());

    // Map the words in the inputs file to BC types, then translate
    //     those types into what they mean for each variable
    init_bcs();

    // Init vertical stretching coeffs
    init_stretch_coeffs();

    last_plot_file_step = -1;
    last_check_file_step = -1;
    last_plot_file_time = -1.0_rt;
    last_check_file_time = -1.0_rt;

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
            advflux_reg[lev] = new YAFluxRegister(grids[lev], grids[lev-1],
                                                   dmap[lev],  dmap[lev-1],
                                                   geom[lev],  geom[lev-1],
                                              ref_ratio[lev-1], lev, NCONS);
        }
    }

    // Fill ghost cells/faces
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        if (lev > 0 && cf_width >= 0) {
            Construct_REMORAFillPatchers(lev);
        }

        if (restart_chkfile == "") {
            FillPatch(lev, t_new[lev], *cons_new[lev], cons_new, BCVars::cons_bc, BdyVars::t, 0, true, false,0,0,0.0,*cons_new[lev]);
            FillPatch(lev, t_new[lev], *xvel_new[lev], xvel_new, BCVars::xvel_bc, BdyVars::u, 0, true, false,0,0,0.0,*xvel_new[lev]);
            FillPatch(lev, t_new[lev], *yvel_new[lev], yvel_new, BCVars::yvel_bc, BdyVars::v, 0, true, false,0,0,0.0,*yvel_new[lev]);
            FillPatch(lev, t_new[lev], *zvel_new[lev], zvel_new, BCVars::zvel_bc, BdyVars::null, 0, true, false);

            // Copy from new into old just in case when initializing from scratch
            int ngs   = cons_new[lev]->nGrow();
            int ngvel = xvel_new[lev]->nGrow();
            MultiFab::Copy(*cons_old[lev],*cons_new[lev],0,0,NCONS,ngs);
            MultiFab::Copy(*xvel_old[lev],*xvel_new[lev],0,0,1,ngvel);
            MultiFab::Copy(*yvel_old[lev],*yvel_new[lev],0,0,1,ngvel);
            MultiFab::Copy(*zvel_old[lev],*zvel_new[lev],0,0,1,IntVect(ngvel,ngvel,0));
        }
    } // lev

    // Check for additional plotting variables that are available after
    // particle containers are setup.
    const std::string& pv1 = "plot_vars"; appendPlotVariables(pv1);

    if (restart_chkfile == "" && (check_int > 0 || check_int_time > 0.0_rt))
    {
        WriteCheckpointFile();
        last_check_file_step = 0;
    }

    if ( (restart_chkfile == "") ||
         (restart_chkfile != "" && plot_file_on_restart) )
    {
        if (plot_int > 0 || plot_int_time > 0.0)
        {
            if (plotfile_type == PlotfileType::amrex)
                WritePlotFile();
#ifdef REMORA_USE_NETCDF
            if (plotfile_type == PlotfileType::netcdf) {
                int step0 = 0;
                WriteNCPlotFile(step0);
                history_count++;
            }
#endif
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
    if (solverChoice.ic_type == IC_Type::analytic) {
        prob->init_analytic_zeta(lev, geom[lev], solverChoice, *this, *vec_zeta[lev]);

#ifdef REMORA_USE_NETCDF
    } else if (solverChoice.ic_type == IC_Type::netcdf) {
        amrex::Print() << "Calling init_zeta_from_netcdf on level " << lev << std::endl;
        init_zeta_from_netcdf(lev);
        amrex::Print() << "Sea surface height loaded from netcdf file \n " << std::endl;
#endif
    } else {
        Abort("Don't know this ic_type!");
    }
    vec_zeta[lev]->FillBoundary(geom[lev].periodicity());
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
    if (solverChoice.init_l0int_h) {
        if (lev==0) {
            if (solverChoice.ic_type == IC_Type::analytic) {
                if (!solverChoice.flat_bathymetry) {
                    prob->init_analytic_bathymetry(lev, geom[lev], solverChoice, *this, *vec_hOfTheConfusingName[lev]);
                } else {
                    init_flat_bathymetry(lev);
                }
#ifdef REMORA_USE_NETCDF
            } else if (solverChoice.ic_type == IC_Type::netcdf) {
                if (!solverChoice.flat_bathymetry) {
                    amrex::Print() << "Calling init_bathymetry_from_netcdf " << std::endl;
                    init_bathymetry_from_netcdf(lev);
                    amrex::Print() << "Bathymetry loaded from netcdf file \n " << std::endl;
                } else {
                    init_flat_bathymetry(lev);
                }
#endif
            } else {
                Abort("Don't know this ic_type!");
            }
            // Need FillBoundary to fill at grid-grid boundaries, and EnforcePeriodicity
            // to make sure ghost cells in the domain corners are consistent.
            vec_hOfTheConfusingName[lev]->FillBoundary(geom[lev].periodicity());
            vec_hOfTheConfusingName[lev]->EnforcePeriodicity(geom[lev].periodicity());
        } else {
            Real dummy_time = 0.0_rt;
            FillCoarsePatch(lev,dummy_time,vec_hOfTheConfusingName[lev].get(), vec_hOfTheConfusingName[lev-1].get(),BCVars::cons_bc);
        }
    } else if (solverChoice.init_ana_h) {
        if (solverChoice.ic_type == IC_Type::analytic) {
            if (!solverChoice.flat_bathymetry) {
                prob->init_analytic_bathymetry(lev, geom[lev], solverChoice, *this,*vec_hOfTheConfusingName[lev]);
            } else {
                init_flat_bathymetry(lev);
            }
#ifdef REMORA_USE_NETCDF
        } else if (solverChoice.ic_type == IC_Type::netcdf) {
            amrex::Print() << "Calling init_bathymetry_from_netcdf level " << lev << std::endl;
            init_bathymetry_from_netcdf(lev);
            amrex::Print() << "Bathymetry loaded from netcdf file \n " << std::endl;
#endif
        } else {
            Abort("Don't know this ic_type!");
        }
        // Need FillBoundary to fill at grid-grid boundaries, and EnforcePeriodicity
        // to make sure ghost cells in the domain corners are consistent.
        vec_hOfTheConfusingName[lev]->FillBoundary(geom[lev].periodicity());
        vec_hOfTheConfusingName[lev]->EnforcePeriodicity(geom[lev].periodicity());
    } else if (solverChoice.init_l1ad_h) {
        if (solverChoice.ic_type == IC_Type::analytic) {
            if (!solverChoice.flat_bathymetry) {
                prob->init_analytic_bathymetry(lev, geom[lev], solverChoice, *this, *vec_hOfTheConfusingName[lev]);
            } else {
                init_flat_bathymetry(lev);
            }
#ifdef REMORA_USE_NETCDF
        } else if (solverChoice.ic_type == IC_Type::netcdf) {
            if (!solverChoice.flat_bathymetry) {
                amrex::Print() << "Calling init_bathymetry_from_netcdf " << std::endl;
                init_bathymetry_from_netcdf(lev);
                amrex::Print() << "Bathymetry loaded from netcdf file \n " << std::endl;
            } else {
                init_flat_bathymetry(lev);
            }
#endif
        } else {
            Abort("Don't know this ic_type!");
        }
        // Need FillBoundary to fill at grid-grid boundaries, and EnforcePeriodicity
        // to make sure ghost cells in the domain corners are consistent.
        vec_hOfTheConfusingName[lev]->FillBoundary(geom[lev].periodicity());
        vec_hOfTheConfusingName[lev]->EnforcePeriodicity(geom[lev].periodicity());
    } else {
        amrex::Abort("Don't know this h init type");
    }
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
            amrex::Print() << "Calling init_coriolis_from_netcdf " << std::endl;
            init_coriolis_from_netcdf(lev);
            amrex::Print() << "Coriolis loaded from netcdf file \n" << std::endl;
#endif
        } else {
            Abort("Don't know this coriolis_type!");
        }

        Real time = 0.0_rt;
        FillPatch(lev, time, *vec_fcor[lev], GetVecOfPtrs(vec_fcor),BCVars::foextrap_bc);
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
    Real time = 0.0_rt;
    prob->init_analytic_vmix(lev, geom[lev], solverChoice, *this,*vec_Akv[lev], *vec_Akt[lev]);
    FillPatch(lev, time, *vec_Akv[lev], GetVecOfPtrs(vec_Akv),BCVars::zvel_bc,BdyVars::null,0,true,false);
    for (int n=0; n<NCONS;n++) {
        FillPatch(lev, time, *vec_Akt[lev], GetVecOfPtrs(vec_Akt),BCVars::zvel_bc,BdyVars::null,0,false,false);
    }
}

/**
 * @param[in   ] lev    level to operate on
 */
void
REMORA::set_hmixcoef(int lev)
{
    BL_PROFILE("REMORA::set_hmixcoef()");
    if (solverChoice.horiz_mixing_type == HorizMixingType::analytic) {
        prob->init_analytic_hmix(lev, geom[lev], solverChoice, *this, *vec_visc2_p[lev], *vec_visc2_r[lev], *vec_diff2[lev]);
    } else if (solverChoice.horiz_mixing_type == HorizMixingType::constant) {
        vec_visc2_p[lev]->setVal(solverChoice.visc2);
        vec_visc2_r[lev]->setVal(solverChoice.visc2);
        for (int n=0; n<NCONS; n++) {
            vec_diff2[lev]->setVal(solverChoice.tnu2[n],n,1);
        }
    } else {
        Abort("Don't know this horizontal mixing type");
    }
    Real time = 0.0_rt;
    FillPatch(lev, time, *vec_visc2_p[lev], GetVecOfPtrs(vec_visc2_p),BCVars::foextrap_periodic_bc);
    FillPatch(lev, time, *vec_visc2_r[lev], GetVecOfPtrs(vec_visc2_r),BCVars::foextrap_periodic_bc);
    for (int n=0; n<NCONS; n++) {
        FillPatch(lev, time, *vec_diff2[lev]  , GetVecOfPtrs(vec_diff2),BCVars::foextrap_periodic_bc,BdyVars::null,n,false);
    }
}

/**
 * @param[in   ] lev    level to operate on
 */
void
REMORA::init_flat_bathymetry(int lev)
{
    BL_PROFILE("REMORA::init_flat_bathymetry()");
    vec_hOfTheConfusingName[lev]->setVal(-geom[0].ProbLo()[2]);
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
        sustr_data_from_file->update_interpolated_to_time(t_old[lev]);
        svstr_data_from_file->update_interpolated_to_time(t_old[lev]);
        FillPatch(lev, t_old[lev], *vec_sustr[lev], GetVecOfPtrs(vec_sustr),BCVars::foextrap_periodic_bc,BdyVars::null,0,false);
        FillPatch(lev, t_old[lev], *vec_svstr[lev], GetVecOfPtrs(vec_svstr),BCVars::foextrap_periodic_bc,BdyVars::null,0,false);
#endif
    }
}

/**
 * @param[in   ] lev    level to operate on
 */
void
REMORA::set_wind(int lev)
{
    BL_PROFILE("REMORA::set_wind()");
    if (solverChoice.wind_type == WindType::analytic) {
        prob->init_analytic_wind(lev,geom[lev], solverChoice, *this, *vec_uwind[lev], *vec_vwind[lev]);
    } else if (solverChoice.wind_type == WindType::netcdf) {
#ifdef REMORA_USE_NETCDF
        Uwind_data_from_file->update_interpolated_to_time(t_old[lev]);
        Vwind_data_from_file->update_interpolated_to_time(t_old[lev]);
        FillPatch(lev, t_old[lev], *vec_uwind[lev], GetVecOfPtrs(vec_uwind),BCVars::foextrap_periodic_bc,BdyVars::null,0,false);
        FillPatch(lev, t_old[lev], *vec_vwind[lev], GetVecOfPtrs(vec_vwind),BCVars::foextrap_periodic_bc,BdyVars::null,0,false);
#endif
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
        amrex::Print() << "Calling init_masks_from_netcdf level " << lev << std::endl;
        init_masks_from_netcdf(lev);
        amrex::Print() << "Masks loaded from netcdf file \n " << std::endl;
#endif
    }
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
    t_old[lev] = time - 1.e200_rt;

    cons_new[lev]->setVal(0.0_rt);
    xvel_new[lev]->setVal(0.0_rt);
    yvel_new[lev]->setVal(0.0_rt);
    zvel_new[lev]->setVal(0.0_rt);

    xvel_old[lev]->setVal(0.0_rt);
    yvel_old[lev]->setVal(0.0_rt);
    zvel_old[lev]->setVal(0.0_rt);

    vec_ru[lev]->setVal(0.0_rt);
    vec_rv[lev]->setVal(0.0_rt);

    vec_ru2d[lev]->setVal(0.0_rt);
    vec_rv2d[lev]->setVal(0.0_rt);

    if (solverChoice.ic_type == IC_Type::analytic) {
        set_grid_scale(lev);
    }
    set_masks(lev);

#ifdef REMORA_USE_NETCDF
    if (solverChoice.ic_type == IC_Type::netcdf) {
        init_clim_nudg_coeff(lev);

        if (solverChoice.do_any_clim_nudg) {
            if (nc_clim_his_file.empty()) {
                amrex::Error("NetCDF climatology file name must be provided via input");
            }
            if (solverChoice.do_m2_clim_nudg) {
                ubar_clim_data_from_file = new NCTimeSeries(nc_clim_his_file, "ubar", clim_ubar_time_varname, geom[lev].Domain(),vec_ubar[lev].get(),true,true);
                vbar_clim_data_from_file = new NCTimeSeries(nc_clim_his_file, "vbar", clim_ubar_time_varname, geom[lev].Domain(),vec_vbar[lev].get(),true,true);
                ubar_clim_data_from_file->Initialize();
                vbar_clim_data_from_file->Initialize();
            }
            if (solverChoice.do_m3_clim_nudg) {
                u_clim_data_from_file = new NCTimeSeries(nc_clim_his_file, "u", clim_u_time_varname, geom[lev].Domain(),xvel_new[lev],false,true);
                v_clim_data_from_file = new NCTimeSeries(nc_clim_his_file, "v", clim_v_time_varname, geom[lev].Domain(),yvel_new[lev],false,true);
                u_clim_data_from_file->Initialize();
                v_clim_data_from_file->Initialize();
            }
            // Since the NCTimeSeries object isn't filling the cons_new MultiFab directly, we don't have to specify a component.
            // It just needs to know the shape of the MultiFab
            if (solverChoice.do_temp_clim_nudg) {
                temp_clim_data_from_file = new NCTimeSeries(nc_clim_his_file, "temp", clim_temp_time_varname,geom[lev].Domain(),cons_new[lev],false,true);
                temp_clim_data_from_file->Initialize();
            }
            if (solverChoice.do_salt_clim_nudg) {
                salt_clim_data_from_file = new NCTimeSeries(nc_clim_his_file, "salt", clim_salt_time_varname,geom[lev].Domain(),cons_new[lev],false,true);
                salt_clim_data_from_file->Initialize();
            }
        }
    }

    if (solverChoice.boundary_from_netcdf) {
        amrex::Print() << "Calling init_bdry_from_netcdf " << std::endl;
        init_bdry_from_netcdf();
        amrex::Print() << "Boundary data loaded from netcdf file \n " << std::endl;
    }

    // This will be a non-op if forcings specified analytically
    if (solverChoice.wind_type == WindType::netcdf) {
        if (nc_frc_file.empty()) {
            amrex::Error("NetCDF forcing file name must be provided via input for winds");
        }
        Uwind_data_from_file = new NCTimeSeries(nc_frc_file, "Uwind", frc_time_varname, geom[lev].Domain(),vec_uwind[lev].get(), true, false);
        Vwind_data_from_file = new NCTimeSeries(nc_frc_file, "Vwind", frc_time_varname, geom[lev].Domain(),vec_vwind[lev].get(), true, false);
        Uwind_data_from_file->Initialize();
        Vwind_data_from_file->Initialize();
    } else if (solverChoice.smflux_type == SMFluxType::netcdf) {
        if (nc_frc_file.empty()) {
            amrex::Error("NetCDF forcing file name must be provided via input for surface momentum fluxes");
        }
        sustr_data_from_file = new NCTimeSeries(nc_frc_file, "sustr", frc_time_varname, geom[lev].Domain(),vec_sustr[lev].get(), true, false);
        svstr_data_from_file = new NCTimeSeries(nc_frc_file, "svstr", frc_time_varname, geom[lev].Domain(),vec_svstr[lev].get(), true, false);
        sustr_data_from_file->Initialize();
        svstr_data_from_file->Initialize();
    }

    if (solverChoice.do_rivers) {
        auto dom = geom[0].Domain();
        int nz = dom.length(2);
        river_source_cons.resize(NCONS);
        Print() << solverChoice.do_rivers_cons[0] << std::endl;
        if ((bool) solverChoice.do_rivers_cons[Salt_comp]) {
            river_source_cons[Salt_comp] = new NCTimeSeriesRiver(nc_riv_file, "river_salt", riv_time_varname, nz);
            river_source_cons[Salt_comp]->Initialize();
        }
        if (solverChoice.do_rivers_cons[Temp_comp]) {
            river_source_cons[Temp_comp] = new NCTimeSeriesRiver(nc_riv_file, "river_temp", riv_time_varname, nz);
            river_source_cons[Temp_comp]->Initialize();
        }
        if (solverChoice.do_rivers_cons[Scalar_comp]) {
            river_source_cons[Scalar_comp] = new NCTimeSeriesRiver(nc_riv_file, "river_scalar", riv_time_varname, nz);
            river_source_cons[Scalar_comp]->Initialize();
        }
        river_source_transport = new NCTimeSeriesRiver(nc_riv_file, "river_transport", riv_time_varname, nz);
        river_source_transport->Initialize();
        river_source_transportbar = new NCTimeSeriesRiver(nc_riv_file, "river_transport", riv_time_varname, nz, 1);
        river_source_transportbar->Initialize();
        init_riv_pos_from_netcdf(lev);
    }
#else
    if (solverChoice.boundary_from_netcdf) {
        Abort("Not compiled with NetCDF, but selected boundary conditions require NetCDF");
    }
#endif

    set_bathymetry(lev);
    set_zeta(lev);
    stretch_transform(lev);

    if (solverChoice.init_l0int_T) {
        if (lev==0) {
            if (solverChoice.ic_type == IC_Type::analytic)
            {
                init_analytic(lev);
#ifdef REMORA_USE_NETCDF
            } else if (solverChoice.ic_type == IC_Type::netcdf) {
                amrex::Print() << "Calling init_data_from_netcdf " << std::endl;
                init_data_from_netcdf(lev);
                set_zeta_to_Ztavg(lev);
                amrex::Print() << "Initial data loaded from netcdf file \n " << std::endl;

#endif
            } else {
                Abort("Need to specify ic_type");
            }
        } else {
            FillCoarsePatch(lev, time, cons_new[lev], cons_new[lev-1],BCVars::Temp_bc_comp,BdyVars::t);
            FillCoarsePatch(lev, time, xvel_new[lev], xvel_new[lev-1],BCVars::xvel_bc,BdyVars::u);
            FillCoarsePatch(lev, time, yvel_new[lev], yvel_new[lev-1],BCVars::yvel_bc,BdyVars::v);
            FillCoarsePatch(lev, time, zvel_new[lev], zvel_new[lev-1],BCVars::zvel_bc,BdyVars::null);
        }
    } else if (solverChoice.init_ana_T || solverChoice.init_l1ad_T) {
        if (solverChoice.ic_type == IC_Type::analytic)
        {
            init_analytic(lev);
#ifdef REMORA_USE_NETCDF
        } else if (solverChoice.ic_type == IC_Type::netcdf) {
            amrex::Print() << "Calling init_data_from_netcdf on level " << lev << std::endl;
            init_data_from_netcdf(lev);
            set_zeta_to_Ztavg(lev);
            amrex::Print() << "Initial data loaded from netcdf file on level " << lev << std::endl;

#endif
        } else {
            Abort("Need to specify ic_type");
        }
    } else {
        amrex::Abort("Need to specify T init procedure");
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
        bool noprefix_max_step = pp.query("max_step", max_step);
        bool noprefix_stop_time = pp.query("stop_time", stop_time);
        bool remora_max_step = pp.query("remora.max_step", max_step);
        bool remora_stop_time = pp.query("remora.stop_time", stop_time);
        if (remora_max_step and noprefix_max_step) {
            Abort("remora.max_step and max_step are both specified. Please use only one!");
        }
        if (remora_stop_time and noprefix_stop_time) {
            Abort("remora.stop_time and stop_time are both specified. Please use only one!");
        }
    }

    ParmParse pp(pp_prefix);
    ParmParse pp_amr("amr");
    {
        pp_amr.query("regrid_int", regrid_int);
        pp.query("check_file", check_file);
        pp.query("check_int", check_int);
        pp_amr.query("check_int", check_int);
        pp.query("check_int_time", check_int_time);
        pp_amr.query("check_int_time", check_int_time);

        pp.query("restart", restart_chkfile);
        pp_amr.query("restart", restart_chkfile);
        pp.query("start_time",start_time);

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
        pp.query("file_min_digits", file_min_digits);

        if (file_min_digits < 0) {
            amrex::Abort("remora.file_min_digits must be non-negative");
        }

        // Time step controls
        pp.query("cfl", cfl);
        pp.query("change_max", change_max);

        pp.query("fixed_dt", fixed_dt);
        pp.query("fixed_fast_dt", fixed_fast_dt);

        pp.query("fixed_ndtfast_ratio", fixed_ndtfast_ratio);

        // If all three are specified, they must be consistent
        if (fixed_dt > 0. && fixed_fast_dt > 0. &&  fixed_ndtfast_ratio > 0)
        {
            if (fixed_dt / fixed_fast_dt != fixed_ndtfast_ratio)
            {
                amrex::Abort("Dt is over-specfied");
            }
        }
        // If two are specified, initialize fixed_ndtfast_ratio
        else if (fixed_dt > 0. && fixed_fast_dt > 0. &&  fixed_ndtfast_ratio <= 0)
        {
            fixed_ndtfast_ratio = static_cast<int>(fixed_dt / fixed_fast_dt);
        }

        AMREX_ASSERT(cfl > 0. || fixed_dt > 0.);

        pp_amr.query("do_substep", do_substep);
        if (do_substep) {
            amrex::Abort("Time substepping is not yet implemented. amr.do_substep must be 0");
        }

        // We use this to keep track of how many boxes we read in from WRF initialization
        num_files_at_level.resize(max_level+1,0);

        // We use this to keep track of how many boxes are specified thru the refinement indicators
        num_boxes_at_level.resize(max_level+1,0);
            boxes_at_level.resize(max_level+1);

        // We always have exactly one file at level 0
        num_boxes_at_level[0] = 1;
        boxes_at_level[0].resize(1);
        boxes_at_level[0][0] = geom[0].Domain();

        // Plotfile name and frequency
        pp.query("plot_file", plot_file_name);
        pp.query("plot_int", plot_int);
        pp.query("plot_int_time", plot_int_time);

        // Should we plot the staggered face velocities (without averaging to cell centers)
        pp.query("plot_staggered_vels", plot_staggered_vels);

        // Output format
        std::string plotfile_type_str = "amrex";
        pp.query("plotfile_type", plotfile_type_str);
        if (plotfile_type_str == "amrex") {
            plotfile_type = PlotfileType::amrex;
        } else if (plotfile_type_str == "netcdf" || plotfile_type_str == "NetCDF") {
            plotfile_type = PlotfileType::netcdf;
#ifdef REMORA_USE_NETCDF
            pp.query("write_history_file",write_history_file);
            pp.query("chunk_history_file",chunk_history_file);
            pp.query("steps_per_history_file",steps_per_history_file);
            // Estimate size of domain for one timestep of netcdf
            auto dom = geom[0].Domain();
            int nx = dom.length(0) + 2;
            int ny = dom.length(1) + 2;
            int nz = dom.length(2);
            if (write_history_file and chunk_history_file and (steps_per_history_file <= 0)) {
                // Estimate number of steps that will fit into a 2GB file.
                steps_per_history_file = int((1.6e10 - NCH2D * nx * ny * 64.0_rt)
                        / (nx * ny * 64.0_rt * (NC3D*nz + NC2D)));
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
                Real est_hist_file_size = NCH2D * nx * ny * 64.0_rt + nt_out * nx * ny * 64.0_rt * (NC3D*nz + NC2D);
                if (est_hist_file_size > 1.6e10) {
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
        // We only read boundary data at level 0
        pp.queryarr("nc_bdry_file", nc_bdry_file);

        // Also only read forcings at level 0 (for now)
        pp.query("nc_frc_file", nc_frc_file);

        // Get river file
        pp.query("nc_river_file", nc_riv_file);

        // Read in file names for climatology history and nudging weights
        pp.query("nc_clim_his_file", nc_clim_his_file);
        pp.query("nc_clim_coeff_file", nc_clim_coeff_file);

        pp.query("bdy_time_varname",bdry_time_varname);
        pp.query("frc_time_varname",frc_time_varname);
        pp.query("riv_time_varname",riv_time_varname);
        pp.query("clim_ubar_time_varname",clim_ubar_time_varname);
        pp.query("clim_vbar_time_varname",clim_vbar_time_varname);
        pp.query("clim_u_time_varname",clim_u_time_varname);
        pp.query("clim_v_time_varname",clim_v_time_varname);
        pp.query("clim_salt_time_varname",clim_salt_time_varname);
        pp.query("clim_temp_time_varname",clim_temp_time_varname);

#endif

#ifdef REMORA_USE_PARTICLES
        readTracersParams();
#endif
    }

    solverChoice.init_params();
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
}

/**
 * @param[in   ] lev    level at which to get time
 */
amrex::Real REMORA::get_t_old(int lev) const
{
    return t_old[lev];
}
