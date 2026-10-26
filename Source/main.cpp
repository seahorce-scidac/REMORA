#include <iostream>

#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_ParallelDescriptor.H>

#include <REMORA.H>

#if defined(AMREX_USE_FLOAT) && (defined(__SSE__) || defined(__x86_64__))
#include <xmmintrin.h>
#include <pmmintrin.h>
#define REMORA_CAN_FLUSH_SUBNORMALS
#endif

std::string inputs_name = "";

using namespace amrex;

// In single precision, small intermediate products (e.g. velocity times metric
// terms) can fall below FLT_MIN. Operations on subnormal floats take a slow
// microcode path on x86, which can make single precision slower than double.
// Flush them to zero (FTZ) and treat subnormal inputs as zero (DAZ) on every
// thread. Disable with remora.flush_subnormals = false.
void set_flush_subnormals () {
#ifdef REMORA_CAN_FLUSH_SUBNORMALS
    bool flush_subnormals = true;
    ParmParse pp("remora");
    pp.query("flush_subnormals",flush_subnormals);
    if (!flush_subnormals) { return; }

    auto set_mode = [] () {
        _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
        _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
    };
    set_mode();
#ifdef AMREX_USE_OMP
#pragma omp parallel
    set_mode();
#endif
#endif
}

// Set the refine_grid_layout flags to (NGROW-1,NGROW-1,0) by default
// since the REMORA default is different from the amrex default (NGROW-1,NGROW-1,NGROW-1)
// Also set max_grid_size to very large since the only reason for
// chopping grids is if Nprocs > Ngrids
void add_par () {
    ParmParse pp("amr");
    ParmParse pp_remora("remora");

    // Set the refine_grid_layout flags to (NGROW-1,NGROW-1,0) by default
    pp.add("refine_grid_layout_x",1);
    pp.add("refine_grid_layout_y",1);
    pp.add("refine_grid_layout_z",0);

    // n_proper is the minimum number of coarse cells between coarse-fine boundaries
    // between levels (ell and ell+1) and levels (ell-1 and ell).   We want this to be
    // greater than or equal to the stencil width (a function of spatial order) divided by
    // ref_ratio (which can be 2,3 or 4).  This ensures that fillpatch at level (ell)
    // does not need to reach beyond level (ell-1). Here to be conservative we set this to 2
    // (rather than the amrex default of 1).
    pp.add("n_proper",2);

    int max_grid_size = 2048;
    pp.queryAdd("max_grid_size",max_grid_size);

    // This will set the default value of blocking_factor to be 1, but will allow
    //     the user to override it in the inputs file or on command line
    int blocking_factor = 1;
    pp.queryAdd("blocking_factor",blocking_factor);

    // Default error buf of 0, user can override. Note that the whole column will be tagged
    // if any one cell is tagged, even if n_error_buf_z = 0
    int n_error_buf = 0;
    pp.queryAdd("n_error_buf",n_error_buf);

    // Inject variables into ParmParse database with the prefixes AMReX expects
    // This allows a REMORA user to use the remora. prefix for some AMReX-specific options
    amrex::Vector<int> n_cell;
    pp_remora.queryarr("n_cell",n_cell);
    pp.queryAdd("n_cell",n_cell);

    int max_level = 0;
    // Default max_level to 0 if not in file
    if (!pp.query("max_level",max_level)) {
        pp.queryAdd("max_level",max_level);
    }

    // If amr.v is not specified, set it equal to remora.v
    int amr_verbosity = 0;
    bool defined_amr_verbosity = pp.query("v",amr_verbosity);
    if (!defined_amr_verbosity) {
        int rem_verbosity = 0;
        pp_remora.query("v",rem_verbosity);
        pp.queryAdd("v",rem_verbosity);
    }

    ParmParse pp_geometry("geometry");
    amrex::Vector<amrex::Real> prob_lo;
    amrex::Vector<amrex::Real> prob_hi;
    amrex::Vector<int> is_periodic;
    pp_remora.queryarr("prob_lo",prob_lo);
    pp_remora.queryarr("prob_hi",prob_hi);
    pp_remora.queryarr("is_periodic",is_periodic);
    pp_geometry.queryAdd("prob_lo",prob_lo);
    pp_geometry.queryAdd("prob_hi",prob_hi);
    pp_geometry.queryAdd("is_periodic",is_periodic);

    ParmParse pp_fabarray("fabarray");
#if defined(REMORA_USE_GPU) || defined(AMREX_USE_OMP)
    // No tiling. On the GPU each box is one kernel launch. With OpenMP the threads share
    // the boxes of a level: the step kernels write shared scratch MultiFabs on tile boxes
    // grown by up to three cells, so tiles of one box would overlap and race (seen as
    // run-to-run differences in Upwelling_Fennel). Boxes do not overlap, so a rank needs at
    // least as many boxes per level as threads for them all to work.
    amrex::Vector<int> omp_tile_size = {1024000,1024000,1024000};
#else
    // Serial CPU: full rows in i, 16 cells in j, the whole column. The kernels recompute the
    // grown cells of every tile, so 8x8 tiles cost 1.7x the untiled time on a 128x128-column
    // box, while j-strips of 16 run 8% faster than no tiling.
    amrex::Vector<int> omp_tile_size = {1024000, 16, 1024000};
#endif
    pp_remora.queryarr("omp_tile_size",omp_tile_size,0,3);
#ifdef AMREX_USE_OMP
    // A tile smaller than the largest box in x or y splits boxes. Inputs commonly write
    // "1024 1024 1024" to mean no tiling, so compare with the box extent, not a fixed value.
    if (omp_get_max_threads() > 1) {
        amrex::Vector<int> box_extent = {max_grid_size, max_grid_size};
        pp.query("max_grid_size_x", box_extent[0]);
        pp.query("max_grid_size_y", box_extent[1]);
        if (n_cell.size() >= 2) {
            box_extent[0] = std::min(box_extent[0], n_cell[0]);
            box_extent[1] = std::min(box_extent[1], n_cell[1]);
        }
        if (omp_tile_size[0] < box_extent[0] || omp_tile_size[1] < box_extent[1]) {
            amrex::Abort("remora.omp_tile_size turns on tiling with more than one OpenMP "
                         "thread. The step kernels write shared scratch on grown tile boxes, "
                         "so tiles of one box race and the answers depend on thread timing. "
                         "Use the default (no tiling; threads share boxes) or run with "
                         "OMP_NUM_THREADS=1.");
        }
    }
#endif
    pp_fabarray.queryAdd("mfiter_tile_size",omp_tile_size);
}

int main(int argc, char* argv[])
{
{
#ifdef AMREX_USE_MPI
    MPI_Init(&argc, &argv);
#endif

    if (argc < 2) {
        // Print usage and exit with error code if no input file was provided.
        REMORA::print_usage(MPI_COMM_WORLD, std::cout);
        REMORA::print_error(
            MPI_COMM_WORLD, "No input file provided. Exiting!!");
        return 1;
    }

    // Look for "-h" or "--help" flag and print usage
    for (auto i = 1; i < argc; i++) {
        const std::string param(argv[i]);
        if ((param == "--help") || (param == "-h") || (param == "--usage")) {
            REMORA::print_banner(MPI_COMM_WORLD, std::cout);
            REMORA::print_usage(MPI_COMM_WORLD, std::cout);
            return 0;
        }
    }

    if (!amrex::FileSystem::Exists(std::string(argv[1]))) {
        // Print usage and exit with error code if we cannot find the input file
        REMORA::print_usage(MPI_COMM_WORLD, std::cout);
        REMORA::print_error(
            MPI_COMM_WORLD, "Input file does not exist = " +
                                std::string(argv[1]) + ". Exiting!!");
        return 1;
    }

  //  print_banner(MPI_COMM_WORLD, std::cout);
    // Check to see if the command line contains --describe
    if (argc >= 2) {
        for (auto i = 1; i < argc; i++) {
            if (std::string(argv[i]) == "--describe") {
                REMORA::writeBuildInfo(std::cout);
                return 0;
            }
        }
    } else if(argc < 2) {
        amrex::Print() << "REMORA is currently under development as a next-generation version of the Regional Ocean Modeling System (ROMS). See Docs for more information.\n"
        <<"  usage:\n"
        <<"    ./REMORA3d.xxx.yyy.ex inputs" <<std::endl;
//        std::cerr << "inputs should follow executable on command line" << std::endl;
        return -1;
    }

    amrex::Initialize(argc,argv,true,MPI_COMM_WORLD,add_par);

    set_flush_subnormals();

    // Save the inputs file name for later.
    if (!strchr(argv[1], '=')) {
      inputs_name = argv[1];
    }

    // timer for profiling
    BL_PROFILE_VAR("main()", pmain);

    // wallclock time
    const Real strt_total = Real(amrex::second());

    {
        BL_PROFILE_VAR("REMORA::InitData()", pinit);
        // constructor - reads in parameters from inputs file
        //             - sizes multilevel arrays and data structures
        REMORA remora;

        // initialize AMR data
        remora.InitData();

        BL_PROFILE_VAR_STOP(pinit);

        // advance solution to final time
        remora.Evolve();

        // wallclock time
        Real end_total = Real(amrex::second()) - strt_total;

        // print wallclock time
        ParallelDescriptor::ReduceRealMax(end_total ,ParallelDescriptor::IOProcessorNumber());
        if (remora.Verbose()) {
            amrex::Print() << "\nTotal Time: " << end_total << '\n';
        }
    }

    // destroy timer for profiling
    BL_PROFILE_VAR_STOP(pmain);
    amrex::Finalize();
#ifdef AMREX_USE_MPI
    MPI_Finalize();
#endif
}
}
