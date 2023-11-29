/**
 * \file ROMSX.cpp
 */

#include <prob_common.H>
#include <EOS.H>
#include <ROMSX.H>

#include <AMReX_buildInfo.H>

using namespace amrex;

amrex::Real ROMSX::startCPUTime        = 0.0;
amrex::Real ROMSX::previousCPUTimeUsed = 0.0;

Vector<AMRErrorTag> ROMSX::ref_tags;

SolverChoice ROMSX::solverChoice;

#ifdef ROMSX_USE_PARTICLES
ParticleData ROMSX::particleData;
#endif

// Time step control
amrex::Real ROMSX::cfl           =  0.8;
amrex::Real ROMSX::fixed_dt      = -1.0;
amrex::Real ROMSX::fixed_fast_dt = -1.0;
amrex::Real ROMSX::init_shrink   =  1.0;
amrex::Real ROMSX::change_max    =  1.1;
int   ROMSX::fixed_ndtfast_ratio = 0;

// Dictate verbosity in screen output
int         ROMSX::verbose       = 0;

// Frequency of diagnostic output
int         ROMSX::sum_interval  = -1;
amrex::Real ROMSX::sum_per       = -1.0;

// Native AMReX vs NetCDF
std::string ROMSX::plotfile_type    = "amrex";

// init_type:  "custom", "ideal", "real"
std::string ROMSX::init_type        = "custom";

#ifdef ROMSX_USE_NETCDF
// NetCDF wrfinput (initialization) file(s)
amrex::Vector<amrex::Vector<std::string>> ROMSX::nc_init_file = {{""}}; // Must provide via input
std::string ROMSX::nc_bdy_file;
#endif

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
        nsubsteps[lev] = do_substep ? MaxRefRatio(lev-1) : 1;
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

#ifdef ROMSX_USE_PARTICLES
    particleData.Redistribute();
#endif

    if (solverChoice.coupling_type == CouplingType::TwoWay)
    {
        for (int lev = finest_level-1; lev >= 0; lev--)
        {
            // This call refluxes from the lev/lev+1 interface onto lev
            getAdvFluxReg(lev+1)->Reflux(vars_new[lev][Vars::cons], 0, 0, NVAR);

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

        for (int lev = 0; lev <= finest_level; lev++)
            init_only(lev, time);

        if (solverChoice.coupling_type == CouplingType::TwoWay) {
            AverageDown();
        }

#ifdef ROMSX_USE_PARTICLES
        particleData.init_particles((amrex::ParGDBBase*)GetParGDB(),vec_z_phys_nd);
#endif

    } else { // Restart from a checkpoint

        restart();

    }

    // Initialize flux registers (whether we start from scratch or restart)
    if (solverChoice.coupling_type == CouplingType::TwoWay) {
        advflux_reg[0] = nullptr;
        for (int lev = 1; lev <= finest_level; lev++)
        {
            advflux_reg[lev] = new YAFluxRegister(grids[lev], grids[lev-1],
                                                   dmap[lev],  dmap[lev-1],
                                                   geom[lev],  geom[lev-1],
                                              ref_ratio[lev-1], lev, NVAR);
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

void
ROMSX::set_bathymetry(int lev)
{
    init_custom_bathymetry(geom[lev], *vec_hOfTheConfusingName[lev], *vec_Zt_avg1[lev], solverChoice);
}

void
ROMSX::set_vmix(int lev) {
    init_custom_vmix(geom[lev], *vec_Akv[lev], *vec_Akt[lev], *vec_z_w[lev], solverChoice);
}

void
ROMSX::set_hmixcoef(int lev) {
    init_custom_hmix(geom[lev], *vec_visc2_p[lev], *vec_visc2_r[lev], *vec_diff2_salt[lev],
            *vec_diff2_temp[lev], solverChoice);
}

void
ROMSX::set_smflux(int lev, Real time) {
    init_custom_smflux(geom[lev], time, *vec_sustr[lev], *vec_svstr[lev], solverChoice);
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

        pp.query("do_substep", do_substep);

        AMREX_ASSERT(cfl > 0. || fixed_dt > 0.);

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
            plotfile_type != "netcdf" && plotfile_type != "NetCDF")
        {
            amrex::Print() << "User selected plotfile_type = " << plotfile_type << std::endl;
            amrex::Abort("Dont know this plotfile_type");
        }
#ifndef ROMSX_USE_NETCDF
        if (plotfile_type == "netcdf" || plotfile_type == "NetCDF")
        {
            amrex::Abort("Please compile with NetCDF in order to enable NetCDF plotfiles");
        }
#endif
        pp.query("plot_file_1", plot_file_1);
        pp.query("plot_file_2", plot_file_2);
        pp.query("plot_int_1", plot_int_1);
        pp.query("plot_int_2", plot_int_2);

#ifdef ROMSX_USE_PARTICLES
        particleData.init_particle_params();
#endif
    }

    solverChoice.init_params();
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
