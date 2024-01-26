#include <iostream>

#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_ParallelDescriptor.H>

#include <REMORA.H>

std::string inputs_name = "";

using namespace amrex;

// Set the refine_grid_layout flags to (NGROW-1,NGROW-1,0) by default
// since the REMORA default is different from the amrex default (NGROW-1,NGROW-1,NGROW-1)
// Also set max_grid_size to very large since the only reason for
// chopping grids is if Nprocs > Ngrids
void add_par () {
   ParmParse pp("amr");

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

   pp.add("n_error_buf",0);
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

    // Save the inputs file name for later.
    if (!strchr(argv[1], '=')) {
      inputs_name = argv[1];
    }

    // timer for profiling
    BL_PROFILE_VAR("main()", pmain);

    // wallclock time
    const Real strt_total = Real(amrex::second());

    {
        // constructor - reads in parameters from inputs file
        //             - sizes multilevel arrays and data structures
        REMORA remora;

        // initialize AMR data
        remora.InitData();

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
    }
}
