/**
 * \file REMORA.cpp
 */

#include <prob_common.H>
#include <EOS.H>
#include <REMORA.H>

#include <AMReX_buildInfo.H>

using namespace amrex;

amrex::Real REMORA::startCPUTime        = 0.0;
amrex::Real REMORA::previousCPUTimeUsed = 0.0;

Vector<AMRErrorTag> REMORA::ref_tags;

SolverChoice REMORA::solverChoice;

#ifdef REMORA_USE_PARTICLES
ParticleData REMORA::particleData;
#endif

// Time step control
amrex::Real REMORA::cfl           =  0.8;
amrex::Real REMORA::fixed_dt      = -1.0;
amrex::Real REMORA::fixed_fast_dt = -1.0;
amrex::Real REMORA::init_shrink   =  1.0;
amrex::Real REMORA::change_max    =  1.1;

int   REMORA::fixed_ndtfast_ratio = 0;

// Dictate verbosity in screen output
int         REMORA::verbose       = 0;

// Frequency of diagnostic output
int         REMORA::sum_interval  = -1;
amrex::Real REMORA::sum_per       = -1.0;

// Native AMReX vs NetCDF
PlotfileType REMORA::plotfile_type    = PlotfileType::amrex;

#ifdef REMORA_USE_NETCDF

int   REMORA::total_nc_plot_file_step = 1;

// Do we write one file per timestep (false) or one file for all timesteps (true)
bool  REMORA::write_history_file      = false;

// NetCDF initialization file
std::string REMORA::nc_bdry_file = ""; // Must provide via input
amrex::Vector<amrex::Vector<std::string>> REMORA::nc_init_file = {{""}}; // Must provide via input
amrex::Vector<amrex::Vector<std::string>> REMORA::nc_grid_file = {{""}}; // Must provide via input
#endif

amrex::Vector<std::string> BCNames = {"xlo", "ylo", "zlo", "xhi", "yhi", "zhi"};

// constructor - reads in parameters from inputs file
//             - sizes multilevel arrays and data structures
//             - initializes BCRec boundary condition object
REMORA::REMORA ()
{
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

    physbcs.resize(nlevs_max);

    t_new.resize(nlevs_max, 0.0);
    t_old.resize(nlevs_max, -1.e100);
    dt.resize(nlevs_max, 1.e100);

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

// advance solution to final time
void
REMORA::Evolve ()
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

        if (max_level == 0) {
            timeStep(lev, cur_time, iteration);
        }
        else {
            timeStepML(cur_time, iteration);
        }

        cur_time  += dt[0];

        amrex::Print() << "Coarse STEP " << step+1 << " ends." << " TIME = " << cur_time
                       << " DT = " << dt[0]  << std::endl;

        if (plot_int_1 > 0 && (step+1) % plot_int_1 == 0) {
            last_plot_file_step_1 = step+1;
            if (plotfile_type == PlotfileType::amrex)
                WritePlotFile(1,plot_var_names_1);
#ifdef REMORA_USE_NETCDF
            else if (plotfile_type == PlotfileType::netcdf)
                WriteNCPlotFile(step);
#endif
        }
        if (plot_int_2 > 0 && (step+1) % plot_int_2 == 0) {
            last_plot_file_step_2 = step+1;
            if (plotfile_type == PlotfileType::amrex)
                WritePlotFile(2,plot_var_names_2);
#ifdef REMORA_USE_NETCDF
            else if (plotfile_type == PlotfileType::netcdf)
                WriteNCPlotFile(step);
#endif
        }

        if (check_int > 0 && (step+1) % check_int == 0) {
            last_check_file_step = step+1;
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

    if (plot_int_1 > 0 && istep[0] > last_plot_file_step_1) {
        if (plotfile_type == PlotfileType::amrex)
            WritePlotFile(1,plot_var_names_1);
#ifdef REMORA_USE_NETCDF
        if (plotfile_type == PlotfileType::netcdf)
            WriteNCPlotFile(istep[0]);
#endif
    }
    if (plot_int_2 > 0 && istep[0] > last_plot_file_step_2) {
        if (plotfile_type == PlotfileType::amrex)
            WritePlotFile(2,plot_var_names_2);
#ifdef REMORA_USE_NETCDF
        else if (plotfile_type == PlotfileType::netcdf)
            WriteNCPlotFile(istep[0]);
#endif
    }

    if (check_int > 0 && istep[0] > last_check_file_step) {
        WriteCheckpointFile();
    }

}

// Called after every coarse timestep
void
REMORA::post_timestep (int nstep, Real time, Real dt_lev0)
{
    BL_PROFILE("REMORA::post_timestep()");

#ifdef REMORA_USE_PARTICLES
    particleData.Redistribute();
#endif

    if (solverChoice.coupling_type == CouplingType::TwoWay)
    {
        for (int lev = finest_level-1; lev >= 0; lev--)
        {
            // This call refluxes from the lev/lev+1 interface onto lev
            getAdvFluxReg(lev+1)->Reflux(*cons_new[lev], 0, 0, NCONS);

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
REMORA::InitData ()
{
    // Initialize the start time for our CPU-time tracker
    startCPUTime = ParallelDescriptor::second();

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

#ifdef REMORA_USE_PARTICLES
        particleData.init_particles((ParGDBBase*)GetParGDB(),vec_z_phys_nd);
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
                                              ref_ratio[lev-1], lev, NCONS);
        }
    }

    // Fill ghost cells/faces
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        FillPatch(lev, t_new[lev], *cons_new[lev], cons_new, BdyVars::t);
        FillPatch(lev, t_new[lev], *xvel_new[lev], xvel_new, BdyVars::u);
        FillPatch(lev, t_new[lev], *yvel_new[lev], yvel_new, BdyVars::v);
        FillPatch(lev, t_new[lev], *zvel_new[lev], zvel_new, BdyVars::null);

        //
        // Copy from new into old just in case
        //
        int ngs   = cons_new[lev]->nGrow();
        int ngvel = xvel_new[lev]->nGrow();
        MultiFab::Copy(*cons_old[lev],*cons_new[lev],0,0,NCONS,ngs);
        MultiFab::Copy(*xvel_old[lev],*xvel_new[lev],0,0,1,ngvel);
        MultiFab::Copy(*yvel_old[lev],*yvel_new[lev],0,0,1,ngvel);
        MultiFab::Copy(*zvel_old[lev],*zvel_new[lev],0,0,1,IntVect(ngvel,ngvel,0));
    } // lev
    if (restart_chkfile == "" && check_int > 0)
    {
        WriteCheckpointFile();
        last_check_file_step = 0;
    }

    if ( (restart_chkfile == "") ||
         (restart_chkfile != "" && plot_file_on_restart) )
    {
        if (plot_int_1 > 0)
        {
            if (plotfile_type == PlotfileType::amrex)
                WritePlotFile(1,plot_var_names_1);
#ifdef REMORA_USE_NETCDF
            if (plotfile_type == PlotfileType::netcdf) {
                int step0 = 0;
                WriteNCPlotFile(step0);
            }
#endif
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

}

void
REMORA::restart()
{
    ReadCheckpointFile();

    // We set this here so that we don't over-write the checkpoint file we just started from
    last_check_file_step = istep[0];
}

void
REMORA::set_bathymetry(int lev)
{
    if (solverChoice.ic_bc_type == IC_BC_Type::Custom) {
        init_custom_bathymetry(geom[lev], *vec_hOfTheConfusingName[lev], solverChoice);

#ifdef REMORA_USE_NETCDF
    } else if (solverChoice.ic_bc_type == IC_BC_Type::Real) {
        amrex::Print() << "Calling init_bathymetry_from_netcdf " << std::endl;
        init_bathymetry_from_netcdf(lev);
        amrex::Print() << "Bathymetry loaded from netcdf file \n " << std::endl;
#endif
    } else {
        Abort("Don't know this ic_bc_type!");
    }

    // HACK -- SHOULD WE ALWAYS DO THIS??
    vec_Zt_avg1[lev]->setVal(0.0);

    Real time = 0.0;
    FillPatch(lev, time, *vec_Zt_avg1[lev], GetVecOfPtrs(vec_Zt_avg1));
    FillPatch(lev, time, *vec_hOfTheConfusingName[lev], GetVecOfPtrs(vec_hOfTheConfusingName));
}

void
REMORA::set_coriolis(int lev) {
    if (solverChoice.use_coriolis) {
        if (solverChoice.coriolis_type == Cor_Type::Custom) {
            init_custom_coriolis(geom[lev], *vec_fcor[lev], solverChoice);
        } else if (solverChoice.coriolis_type == Cor_Type::Beta_Plane) {
            init_beta_plane_coriolis(lev);
#ifdef REMORA_USE_NETCDF
        } else if (solverChoice.coriolis_type == Cor_Type::Real) {
            amrex::Print() << "Calling init_coriolis_from_netcdf " << std::endl;
            init_coriolis_from_netcdf(lev);
            amrex::Print() << "Coriolis loaded from netcdf file \n" << std::endl;
#endif
        } else {
            Abort("Don't know this coriolis_type!");
        }

        Real time = 0.0;
        FillPatch(lev, time, *vec_fcor[lev], GetVecOfPtrs(vec_fcor));
    }
}

void
REMORA::set_vmix(int lev) {
    init_custom_vmix(geom[lev], *vec_Akv[lev], *vec_Akt[lev], *vec_z_w[lev], solverChoice);

    Real time = 0.0;
    FillPatch(lev, time, *vec_Akv[lev], GetVecOfPtrs(vec_Akv));
    FillPatch(lev, time, *vec_Akt[lev], GetVecOfPtrs(vec_Akt));
}

void
REMORA::set_hmixcoef(int lev)
{
    init_custom_hmix(geom[lev], *vec_visc2_p[lev], *vec_visc2_r[lev], *vec_diff2[lev], solverChoice);

    Real time = 0.0;
    FillPatch(lev, time, *vec_visc2_p[lev], GetVecOfPtrs(vec_visc2_p));
    FillPatch(lev, time, *vec_visc2_r[lev], GetVecOfPtrs(vec_visc2_r));
    FillPatch(lev, time, *vec_diff2[lev]  , GetVecOfPtrs(vec_diff2));
}

void
REMORA::set_smflux(int lev, Real time)
{
    init_custom_smflux(geom[lev], time, *vec_sustr[lev], *vec_svstr[lev], solverChoice);

    // FillPatch(lev, time, *vec_sustr[lev], GetVecOfPtrs(vec_sustr));
    // FillPatch(lev, time, *vec_svstr[lev], GetVecOfPtrs(vec_svstr));
}

void
REMORA::init_only(int lev, Real time)
{
    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    cons_new[lev]->setVal(0.0);
    xvel_new[lev]->setVal(0.0);
    yvel_new[lev]->setVal(0.0);
    zvel_new[lev]->setVal(0.0);

    if (solverChoice.ic_bc_type == IC_BC_Type::Custom)
    {
        init_custom(lev);
#ifdef REMORA_USE_NETCDF
    } else if (solverChoice.ic_bc_type == IC_BC_Type::Real) {

        amrex::Print() << "Calling init_data_from_netcdf " << std::endl;
        init_data_from_netcdf(lev);
        amrex::Print() << "Initial data loaded from netcdf file \n " << std::endl;

        amrex::Print() << "Calling init_bdry_from_netcdf " << std::endl;
        init_bdry_from_netcdf();
        amrex::Print() << "Boundary data loaded from netcdf file \n " << std::endl;
#endif
    } else {
        Abort("Need to specify ic_bc_type");
    }

    // Ensure that the face-based data are the same on both sides of a periodic domain.
    // The data associated with the lower grid ID is considered the correct value.
    xvel_new[lev]->OverrideSync(geom[lev].periodicity());
    yvel_new[lev]->OverrideSync(geom[lev].periodicity());
    zvel_new[lev]->OverrideSync(geom[lev].periodicity());
}

// read in some parameters from inputs file
void
REMORA::ReadParameters ()
{
    {
        ParmParse pp;  // Traditionally, max_step and stop_time do not have prefix.
        pp.query("max_step", max_step);
        pp.query("stop_time", stop_time);
    }

    ParmParse pp(pp_prefix);
    {
        pp.query("regrid_int", regrid_int);
        pp.query("check_file", check_file);
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

        // We use this to keep track of how many boxes we read in from WRF initialization
        num_files_at_level.resize(max_level+1,0);

        // We use this to keep track of how many boxes are specified thru the refinement indicators
        num_boxes_at_level.resize(max_level+1,0);
            boxes_at_level.resize(max_level+1);

        // We always have exactly one file at level 0
        num_boxes_at_level[0] = 1;
        boxes_at_level[0].resize(1);
        boxes_at_level[0][0] = geom[0].Domain();

        // Output format
        std::string plotfile_type_str = "amrex";
        pp.query("plotfile_type", plotfile_type_str);
        if (plotfile_type_str == "amrex") {
            plotfile_type = PlotfileType::amrex;
        } else if (plotfile_type_str == "netcdf" || plotfile_type_str == "NetCDF") {
            plotfile_type = PlotfileType::netcdf;
#ifdef REMORA_USE_NETCDF
            pp.query("write_history_file",write_history_file);
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
        pp.query("nc_bdry_file", nc_bdry_file);

        // Query the set and total widths for bdy interior ghost cells
        pp.query("bdy_width", bdy_width);
        pp.query("bdy_set_width", bdy_set_width);
        AMREX_ALWAYS_ASSERT(bdy_width >= 0);
        AMREX_ALWAYS_ASSERT(bdy_set_width >= 0);
        AMREX_ALWAYS_ASSERT(bdy_width >= bdy_set_width);
#endif
        pp.query("plot_file_1", plot_file_1);
        pp.query("plot_file_2", plot_file_2);
        pp.query("plot_int_1", plot_int_1);
        pp.query("plot_int_2", plot_int_2);

#ifdef REMORA_USE_PARTICLES
        particleData.init_particle_params();
#endif
    }

    solverChoice.init_params();
}

// Set covered coarse cells to be the average of overlying fine cells for all levels
void
REMORA::AverageDown ()
{
    for (int lev = finest_level-1; lev >= 0; --lev)
    {
        AverageDownTo(lev);
    }
}

// Set covered coarse cells to be the average of overlying fine cells at level crse_lev
void
REMORA::AverageDownTo (int crse_lev)
{
    average_down(*cons_new[crse_lev+1], *cons_new[crse_lev],
                 0, cons_new[crse_lev]->nComp(), refRatio(crse_lev));

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
