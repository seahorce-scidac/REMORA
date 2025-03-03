#include <iomanip>
#include <iostream>
#include <string>
#include <ctime>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <AMReX_Utility.H>
#include <AMReX_buildInfo.H>
#include <AMReX_ParmParse.H>

#include "REMORA.H"
#include "REMORA_NCInterface.H"
#include "REMORA_NCPlotFile.H"
#include "REMORA_IndexDefines.H"

using namespace amrex;

void REMORA::WriteNCPlotFile(int which_step) {
    AMREX_ASSERT(max_level == 0);
    // For right now we assume single level -- we will generalize this later to multilevel
    int lev = 0;
    int which_subdomain = 0;
    int which_step_in_chunk = -1;

    // Create filename
    std::string plt_string;
    std::string plotfilename;
    if (REMORA::write_history_file) {
        plotfilename = plot_file_name + "_his";
    } else {
        plotfilename = Concatenate(plot_file_name, which_step, file_min_digits);
    }
    // If chunking, concatenate with which file we're in
    if (REMORA::write_history_file and REMORA::chunk_history_file) {
        int which_chunk = history_count / REMORA::steps_per_history_file;
        plotfilename = Concatenate(plotfilename, which_chunk, file_min_digits);
        which_step_in_chunk = history_count - which_chunk * REMORA::steps_per_history_file;
    }

    // Set the full IO path for NetCDF output
    std::string FullPath = plotfilename;
    if (lev == 0) {
        const std::string &extension = amrex::Concatenate("_d", lev + 1, 2);
        FullPath += extension + ".nc";
    } else {
        const std::string &extension = amrex::Concatenate("_d", lev + 1 + which_subdomain, 2);
        FullPath += extension + ".nc";
    }

    //
    // Check if this file/directory already exists and if so,
    //       have the IOProcessor move the existing
    //       file/directory to filename.old
    //
    if ((!REMORA::write_history_file) || (which_step == 0) || (which_step_in_chunk == 0)) {
        if (amrex::ParallelDescriptor::IOProcessor()) {
            if (amrex::FileExists(FullPath)) {
                std::string newoldname(FullPath + ".old." + amrex::UniqueString());
                amrex::Print() << "WriteNCPlotFile:  " << FullPath << " exists.  Renaming to:  " << newoldname << std::endl;
                if (std::rename(FullPath.c_str(), newoldname.c_str())) {
                    amrex::Abort("WriteNCPlotfile:: std::rename failed");
                }
            }
        }
        ParallelDescriptor::Barrier();
    }

    bool is_history;

    if (REMORA::write_history_file) {
        is_history = true;
        bool write_header = !(amrex::FileExists(FullPath));

        auto ncf =
                write_header ?
                        ncutils::NCFile::create(FullPath, NC_CLOBBER|NC_64BIT_DATA, amrex::ParallelContext::CommunicatorSub(), MPI_INFO_NULL) :
                        ncutils::NCFile::open(FullPath, NC_WRITE, amrex::ParallelContext::CommunicatorSub(), MPI_INFO_NULL);

        amrex::Print() << "Writing into level " << lev << " NetCDF history file " << FullPath << std::endl;

        WriteNCPlotFile_which(lev, which_subdomain, write_header, ncf, is_history);

    } else {

        is_history = false;
        bool write_header = true;

        // Open new netcdf file to write data
        auto ncf = ncutils::NCFile::create(FullPath, NC_CLOBBER|NC_64BIT_DATA, amrex::ParallelContext::CommunicatorSub(), MPI_INFO_NULL);
        amrex::Print() << "Writing level " << lev << " NetCDF plot file " << FullPath << std::endl;

        WriteNCPlotFile_which(lev, which_subdomain, write_header, ncf, is_history);
    }
}

void REMORA::WriteNCPlotFile_which(int lev, int which_subdomain, bool write_header, ncutils::NCFile &ncf, bool is_history) {
    // Number of cells in this "domain" at this level
    std::vector<int> n_cells;

    // We only do single-level writes when using NetCDF format
    int flev = 1; //max_level;

    Box subdomain;
    if (lev == 0) {
        subdomain = geom[lev].Domain();
    } else {
        subdomain = boxes_at_level[lev][which_subdomain];
    }

    int nx = subdomain.length(0);
    int ny = subdomain.length(1);
    int nz = subdomain.length(2);

    // unsigned long int nt= NC_UNLIMITED;
    if (is_history && max_step < 0)
        amrex::Abort("Need to know max_step if writing history file");
    long long int nt = is_history ? static_cast<long long int>(max_step / std::min(plot_int, max_step)) + 1 : 1;
    if (chunk_history_file) {
        // First index of the last history file
        int last_file_index = REMORA::steps_per_history_file * int(nt / REMORA::steps_per_history_file);
        if (history_count >= last_file_index) {
            nt = nt - last_file_index;
        } else {
            nt = REMORA::steps_per_history_file;
        }
    }

    n_cells.push_back(nx);
    n_cells.push_back(ny);
    n_cells.push_back(nz);

    int num_pts = nx * ny * nz;

    const std::string nt_name = "ocean_time";
    const std::string ndim_name = "num_geo_dimensions";

    const std::string flev_name = "FINEST_LEVEL";

    const std::string nx_name = "NX";
    const std::string ny_name = "NY";
    const std::string nz_name = "NZ";

    const std::string nx_r_name = "xi_rho";
    const std::string ny_r_name = "eta_rho";
    const std::string nz_r_name = "s_rho";

    const std::string nx_u_name = "xi_u";
    const std::string ny_u_name = "eta_u";

    const std::string nx_v_name = "xi_v";
    const std::string ny_v_name = "eta_v";

    const std::string nx_p_name = "xi_psi";
    const std::string ny_p_name = "eta_psi";
    const std::string nz_w_name = "s_w";

    const Real fill_value = 1.0e37_rt;

    if (write_header) {
        ncf.enter_def_mode();
        ncf.put_attr("title", "REMORA data ");
        ncf.def_dim(nt_name, nt);
        ncf.def_dim(ndim_name, AMREX_SPACEDIM);

        ncf.def_dim(nx_r_name, nx + 2);
        ncf.def_dim(ny_r_name, ny + 2);
        ncf.def_dim(nz_r_name, nz);

        ncf.def_dim(nx_u_name, nx + 1);
        ncf.def_dim(ny_u_name, ny + 2);

        ncf.def_dim(nx_v_name, nx + 2);
        ncf.def_dim(ny_v_name, ny + 1);

        ncf.def_dim(nx_p_name, nx + 1);
        ncf.def_dim(ny_p_name, ny + 1);

        ncf.def_dim(nz_w_name, nz + 1);

        ncf.def_dim(flev_name, flev);

        ncf.def_dim(nx_name, n_cells[0]);
        ncf.def_dim(ny_name, n_cells[1]);
        ncf.def_dim(nz_name, n_cells[2]);

        ncf.def_var("probLo", ncutils::NCDType::Real, { ndim_name });
        ncf.var("probLo").put_attr("long_name","Low side of problem domain in internal AMReX grid");
        ncf.var("probLo").put_attr("units","meter");
        ncf.def_var("probHi", ncutils::NCDType::Real, { ndim_name });
        ncf.var("probHi").put_attr("long_name","High side of problem domain in internal AMReX grid");
        ncf.var("probHi").put_attr("units","meter");

        ncf.def_var("Geom.smallend", NC_INT, { flev_name, ndim_name });
        ncf.var("Geom.smallend").put_attr("long_name","Low side of problem domain in index space");
        ncf.def_var("Geom.bigend", NC_INT, { flev_name, ndim_name });
        ncf.var("Geom.bigend").put_attr("long_name","High side of problem domain in index space");
        ncf.def_var("CellSize", ncutils::NCDType::Real, { flev_name, ndim_name });
        ncf.var("CellSize").put_attr("long_name","Cell size on internal AMReX grid");
        ncf.var("CellSize").put_attr("units","meter");

        ncf.def_var("theta_s",ncutils::NCDType::Real,{});
        ncf.var("theta_s").put_attr("long_name","S-coordinate surface control parameter");
        ncf.def_var("theta_b",ncutils::NCDType::Real,{});
        ncf.var("theta_b").put_attr("long_name","S-coordinate bottom control parameter");
        ncf.def_var("hc",ncutils::NCDType::Real,{});
        ncf.var("hc").put_attr("long_name","S-coordinate parameter, critical depth");
        ncf.var("hc").put_attr("units","meter");

        ncf.def_var("grid",NC_INT, {});
        ncf.var("grid").put_attr("cf_role","grid_topology");
        ncf.var("grid").put_attr("topology_dimension","2");
        ncf.var("grid").put_attr("node_dimensions", "xi_psi eta_psi");
        ncf.var("grid").put_attr("face_dimensions", "xi_rho: xi_psi (padding: both) eta_rho: eta_psi (padding: both)");
        ncf.var("grid").put_attr("edge1_dimensions", "xi_u: xi_psi eta_u: eta_psi (padding: both)");
        ncf.var("grid").put_attr("edge2_dimensions", "xi_v: xi_psi (padding: both) eta_v: eta_psi");
        ncf.var("grid").put_attr("node_coordinates", "x_psi y_psi");
        ncf.var("grid").put_attr("face_coordinates", "x_rho y_rho");
        ncf.var("grid").put_attr("edge1_coordinates", "x_u y_u");
        ncf.var("grid").put_attr("edge2_coordinates", "x_v y_v");
        ncf.var("grid").put_attr("vertical_dimensions", "s_rho: s_w (padding: none)");

        ncf.def_var("s_rho",ncutils::NCDType::Real, {nz_r_name});
        ncf.var("s_rho").put_attr("long_name","S-coordinate at RHO-points");
        ncf.var("s_rho").put_attr("field","s_rho, scalar");

        ncf.def_var("s_w",ncutils::NCDType::Real, {nz_w_name});
        ncf.var("s_w").put_attr("long_name","S-coordinate at W-points");
        ncf.var("s_w").put_attr("field","s_w, scalar");

        ncf.def_var("pm",ncutils::NCDType::Real, {ny_r_name, nx_r_name});
        ncf.var("pm").put_attr("long_name","curvilinear coordinate metric in XI");
        ncf.var("pm").put_attr("units","meter-1");
        ncf.var("pm").put_attr("grid","grid");
        ncf.var("pm").put_attr("location","face");
        ncf.var("pm").put_attr("coordinates","x_rho y_rho");
        ncf.var("pm").put_attr("field","pm, scalar");

        ncf.def_var("pn",ncutils::NCDType::Real, {ny_r_name, nx_r_name});
        ncf.var("pn").put_attr("long_name","curvilinear coordinate metric in ETA");
        ncf.var("pn").put_attr("units","meter-1");
        ncf.var("pn").put_attr("grid","grid");
        ncf.var("pn").put_attr("location","face");
        ncf.var("pn").put_attr("coordinates","x_rho y_rho");
        ncf.var("pn").put_attr("field","pn, scalar");

        ncf.def_var("f",ncutils::NCDType::Real, {ny_r_name, nx_r_name});
        ncf.var("f").put_attr("long_name","Coriolis parameter at RHO-points");
        ncf.var("f").put_attr("units","second-1");
        ncf.var("f").put_attr("grid","grid");
        ncf.var("f").put_attr("location","face");
        ncf.var("f").put_attr("coordinates","x_rho y_rho");
        ncf.var("f").put_attr("field","coriolis, scalar");

        ncf.def_var("x_rho",ncutils::NCDType::Real, {ny_r_name, nx_r_name});
        ncf.var("x_rho").put_attr("long_name","x-locations of RHO-points");
        ncf.var("x_rho").put_attr("units","meter");
        ncf.var("x_rho").put_attr("field","x_rho, scalar");

        ncf.def_var("y_rho",ncutils::NCDType::Real, {ny_r_name, nx_r_name});
        ncf.var("y_rho").put_attr("long_name","y-locations of RHO-points");
        ncf.var("y_rho").put_attr("units","meter");
        ncf.var("y_rho").put_attr("field","y_rho, scalar");

        ncf.def_var("x_u",ncutils::NCDType::Real, {ny_u_name, nx_u_name});
        ncf.var("x_u").put_attr("long_name","x-locations of U-points");
        ncf.var("x_u").put_attr("units","meter");
        ncf.var("x_u").put_attr("field","x_u, scalar");

        ncf.def_var("y_u",ncutils::NCDType::Real, {ny_u_name, nx_u_name});
        ncf.var("y_u").put_attr("long_name","y-locations of U-points");
        ncf.var("y_u").put_attr("units","meter");
        ncf.var("y_u").put_attr("field","y_u, scalar");

        ncf.def_var("x_v",ncutils::NCDType::Real, {ny_v_name, nx_v_name});
        ncf.var("x_v").put_attr("long_name","x-locations of V-points");
        ncf.var("x_v").put_attr("units","meter");
        ncf.var("x_v").put_attr("field","x_v, scalar");

        ncf.def_var("y_v",ncutils::NCDType::Real, {ny_v_name, nx_v_name});
        ncf.var("y_v").put_attr("long_name","y-locations of V-points");
        ncf.var("y_v").put_attr("units","meter");
        ncf.var("y_v").put_attr("field","y_v, scalar");

        ncf.def_var("x_psi",ncutils::NCDType::Real, {ny_p_name, nx_p_name});
        ncf.var("x_psi").put_attr("long_name","x-locations of PSI-points");
        ncf.var("x_psi").put_attr("units","meter");
        ncf.var("x_psi").put_attr("field","x_psi, scalar");

        ncf.def_var("y_psi",ncutils::NCDType::Real, {ny_p_name, nx_p_name});
        ncf.var("y_psi").put_attr("long_name","y-locations of PSI-points");
        ncf.var("y_psi").put_attr("units","meter");
        ncf.var("y_psi").put_attr("field","y_psi, scalar");

        ncf.def_var("ocean_time", ncutils::NCDType::Real, { nt_name });
        ncf.var("ocean_time").put_attr("long_name","time since initialization");
        ncf.var("ocean_time").put_attr("units","seconds since 0001-01-01 00:00:00");
        ncf.var("ocean_time").put_attr("field","time, scalar, series");

        ncf.def_var("h", ncutils::NCDType::Real, { ny_r_name, nx_r_name });
        ncf.var("h").put_attr("long_name","bathymetry at RHO-points");
        ncf.var("h").put_attr("units","meter");
        ncf.var("h").put_attr("grid","grid");
        ncf.var("h").put_attr("location","face");
        ncf.var("h").put_attr("coordinates","x_rho y_rho");
        ncf.var("h").put_attr("field","bath, scalar");

        ncf.def_var_fill("zeta", ncutils::NCDType::Real, { nt_name, ny_r_name, nx_r_name }, &fill_value);
        ncf.var("zeta").put_attr("long_name","free-surface");
        ncf.var("zeta").put_attr("units","meter");
        ncf.var("zeta").put_attr("time","ocean_time");
        ncf.var("zeta").put_attr("grid","grid");
        ncf.var("zeta").put_attr("location","face");
        ncf.var("zeta").put_attr("coordinates","x_rho y_rho ocean_time");
        ncf.var("zeta").put_attr("field","free-surface, scalar, series");

        ncf.def_var_fill("temp", ncutils::NCDType::Real, { nt_name, nz_r_name, ny_r_name, nx_r_name }, &fill_value);
        ncf.var("temp").put_attr("long_name","potential temperature");
        ncf.var("temp").put_attr("units","Celsius");
        ncf.var("temp").put_attr("time","ocean_time");
        ncf.var("temp").put_attr("grid","grid");
        ncf.var("temp").put_attr("location","face");
        ncf.var("temp").put_attr("coordinates","x_rho y_rho s_rho ocean_time");
        ncf.var("temp").put_attr("field","temperature, scalar, series");

        ncf.def_var_fill("salt", ncutils::NCDType::Real, { nt_name, nz_r_name, ny_r_name, nx_r_name }, &fill_value);
        ncf.var("salt").put_attr("long_name","salinity");
        ncf.var("salt").put_attr("time","ocean_time");
        ncf.var("salt").put_attr("grid","grid");
        ncf.var("salt").put_attr("location","face");
        ncf.var("salt").put_attr("coordinates","x_rho y_rho s_rho ocean_time");
        ncf.var("salt").put_attr("field","salinity, scalar, series");

        ncf.def_var_fill("tracer", ncutils::NCDType::Real, { nt_name, nz_r_name, ny_r_name, nx_r_name }, &fill_value);
        ncf.var("tracer").put_attr("long_name","passive tracer");
        ncf.var("tracer").put_attr("time","ocean_time");
        ncf.var("tracer").put_attr("grid","grid");
        ncf.var("tracer").put_attr("location","face");
        ncf.var("tracer").put_attr("coordinates","x_rho y_rho s_rho ocean_time");
        ncf.var("tracer").put_attr("field","tracer, scalar, series");

        ncf.def_var_fill("u", ncutils::NCDType::Real, { nt_name, nz_r_name, ny_u_name, nx_u_name }, &fill_value);
        ncf.var("u").put_attr("long_name","u-momentum component");
        ncf.var("u").put_attr("units","meter second-1");
        ncf.var("u").put_attr("time","ocean_time");
        ncf.var("u").put_attr("grid","grid");
        ncf.var("u").put_attr("location","edge1");
        ncf.var("u").put_attr("coordinates","x_u y_u s_rho ocean_time");
        ncf.var("u").put_attr("field","u-velocity, scalar, series");

        ncf.def_var_fill("v", ncutils::NCDType::Real, { nt_name, nz_r_name, ny_v_name, nx_v_name }, &fill_value);
        ncf.var("v").put_attr("long_name","v-momentum component");
        ncf.var("v").put_attr("units","meter second-1");
        ncf.var("v").put_attr("time","ocean_time");
        ncf.var("v").put_attr("grid","grid");
        ncf.var("v").put_attr("location","edge2");
        ncf.var("v").put_attr("coordinates","x_v y_v s_rho ocean_time");
        ncf.var("v").put_attr("field","v-velocity, scalar, series");

        ncf.def_var_fill("ubar", ncutils::NCDType::Real, { nt_name, ny_u_name, nx_u_name }, &fill_value);
        ncf.var("ubar").put_attr("long_name","vertically integrated u-momentum component");
        ncf.var("ubar").put_attr("units","meter second-1");
        ncf.var("ubar").put_attr("time","ocean_time");
        ncf.var("ubar").put_attr("grid","grid");
        ncf.var("ubar").put_attr("location","edge1");
        ncf.var("ubar").put_attr("coordinates","x_u y_u ocean_time");
        ncf.var("ubar").put_attr("field","ubar-velocity, scalar, series");

        ncf.def_var_fill("vbar", ncutils::NCDType::Real, { nt_name, ny_v_name, nx_v_name }, &fill_value);
        ncf.var("vbar").put_attr("long_name","vertically integrated v-momentum component");
        ncf.var("vbar").put_attr("units","meter second-1");
        ncf.var("vbar").put_attr("time","ocean_time");
        ncf.var("vbar").put_attr("grid","grid");
        ncf.var("vbar").put_attr("location","edge2");
        ncf.var("vbar").put_attr("coordinates","x_v y_v ocean_time");
        ncf.var("vbar").put_attr("field","vbar-velocity, scalar, series");

        ncf.def_var("sustr", ncutils::NCDType::Real, { nt_name, ny_u_name, nx_u_name });
        ncf.var("sustr").put_attr("long_name","surface u-momentum stress");
        ncf.var("sustr").put_attr("units","newton meter-2");
        ncf.var("sustr").put_attr("time","ocean_time");
        ncf.var("sustr").put_attr("grid","grid");
        ncf.var("sustr").put_attr("location","edge1");
        ncf.var("sustr").put_attr("coordinates","x_u y_u ocean_time");
        ncf.var("sustr").put_attr("field","surface u-momentum stress, scalar, series");

        ncf.def_var("svstr", ncutils::NCDType::Real, { nt_name, ny_v_name, nx_v_name });
        ncf.var("svstr").put_attr("long_name","surface v-momentum stress");
        ncf.var("svstr").put_attr("units","newton meter-2");
        ncf.var("svstr").put_attr("time","ocean_time");
        ncf.var("svstr").put_attr("grid","grid");
        ncf.var("svstr").put_attr("location","edge2");
        ncf.var("svstr").put_attr("coordinates","x_v y_v ocean_time");
        ncf.var("svstr").put_attr("field","surface v-momentum stress, scalar, series");

        Real time = 0.;

        // Right now this is hard-wired to {temp, salt, tracer, u, v}
        int n_data_items = 5;
//        ncf.put_attr("number_variables", std::vector<int> { n_data_items });
        ncf.put_attr("space_dimension", std::vector<int> { AMREX_SPACEDIM });
//        ncf.put_attr("current_time", std::vector<double> { time });
        ncf.put_attr("start_time", std::vector<double> { start_bdy_time });
        ncf.put_attr("CurrentLevel", std::vector<int> { flev });
        ncf.put_attr("DefaultGeometry", std::vector<int> { amrex::DefaultGeometry().Coord() });

        ncf.exit_def_mode();

        // We are doing single-level writes but it doesn't have to be level 0
        //
        // Write out the header information.
        //

        Real dx[AMREX_SPACEDIM];
        for (int i = 0; i < AMREX_SPACEDIM; i++) {
            dx[i] = geom[lev].CellSize()[i];
        }
        const auto *base = geom[lev].ProbLo();
        RealBox rb(subdomain, dx, base);

        amrex::Vector<Real> probLo;
        amrex::Vector<Real> probHi;
        for (int i = 0; i < AMREX_SPACEDIM; i++) {
            probLo.push_back(rb.lo(i));
            probHi.push_back(rb.hi(i));
        }

        //nc_probLo.par_access(NC_COLLECTIVE);
        // small variable data written by just the master proc
        ncmpi_begin_indep_data(ncf.ncid);
        if (amrex::ParallelDescriptor::IOProcessor()) // only master proc
        {
            auto nc_probLo = ncf.var("probLo");

            nc_probLo.put(probLo.data(), { 0 }, { AMREX_SPACEDIM });

            auto nc_probHi = ncf.var("probHi");
            //nc_probHi.par_access(NC_COLLECTIVE);
            nc_probHi.put(probHi.data(), { 0 }, { AMREX_SPACEDIM });

            amrex::Vector<int> smallend;
            amrex::Vector<int> bigend;
            for (int i = lev; i < flev; i++) {
                smallend.clear();
                bigend.clear();
                for (int j = 0; j < AMREX_SPACEDIM; j++) {
                    smallend.push_back(subdomain.smallEnd(j));
                    bigend.push_back(subdomain.bigEnd(j));
                }
                auto nc_Geom_smallend = ncf.var("Geom.smallend");
                //nc_Geom_smallend.par_access(NC_COLLECTIVE);
                nc_Geom_smallend.put(smallend.data(), { static_cast<long long int>(i - lev), 0 }, { 1,
                AMREX_SPACEDIM });

                auto nc_Geom_bigend = ncf.var("Geom.bigend");
                //nc_Geom_bigend.par_access(NC_COLLECTIVE);
                nc_Geom_bigend.put(bigend.data(), { static_cast<long long int>(i - lev), 0 }, { 1,
                AMREX_SPACEDIM });
            }

            amrex::Vector<Real> CellSize;
            for (int i = lev; i < flev; i++) {
                CellSize.clear();
                for (Real &j : dx) {
                    CellSize.push_back(amrex::Real(j));
                }
                auto nc_CellSize = ncf.var("CellSize");
                //nc_CellSize.par_access(NC_COLLECTIVE);
                nc_CellSize.put(CellSize.data(), { static_cast<long long int>(i - lev), 0 }, { 1,
                AMREX_SPACEDIM });
            }
            Real hc = solverChoice.tcline;
            Real theta_s = solverChoice.theta_s;
            Real theta_b = solverChoice.theta_b;
            ncf.var("hc").put(&hc);
            ncf.var("theta_s").put(&theta_s);
            ncf.var("theta_b").put(&theta_b);

        }
        ncmpi_end_indep_data(ncf.ncid);

    } // end if write_header

    ncmpi_begin_indep_data(ncf.ncid);
    //
    // We compute the offsets based on location of the box within the domain
    //
    long long adjusted_history_count = chunk_history_file ? history_count % steps_per_history_file : history_count;
    long long local_start_nt = (is_history ? static_cast<long long>(adjusted_history_count) : static_cast<long long>(0));
    long long local_nt = 1; // We write data for only one time

    {
        auto nc_plot_var = ncf.var("ocean_time");
        //nc_plot_var.par_access(NC_COLLECTIVE);
        nc_plot_var.put(&t_new[lev], { local_start_nt }, { local_nt });
    }
    // do all independent writes
    //ncmpi_end_indep_data(ncf.ncid);

    cons_new[lev]->FillBoundary(geom[lev].periodicity());

    mask_arrays_for_write(lev, (Real) fill_value);

    for (MFIter mfi(*cons_new[lev], false); mfi.isValid(); ++mfi) {
        auto bx = mfi.validbox();
        if (subdomain.contains(bx)) {
            //
            // We only include one grow cell at subdomain boundaries, not internal grid boundaries
            //
            Box tmp_bx(bx);
            if (tmp_bx.smallEnd()[0] == subdomain.smallEnd()[0])
                tmp_bx.growLo(0, 1);
            if (tmp_bx.smallEnd()[1] == subdomain.smallEnd()[1])
                tmp_bx.growLo(1, 1);
            if (tmp_bx.bigEnd()[0] == subdomain.bigEnd()[0])
                tmp_bx.growHi(0, 1);
            if (tmp_bx.bigEnd()[1] == subdomain.bigEnd()[1])
                tmp_bx.growHi(1, 1);
            // amrex::Print() << "    BX " << bx << std::endl;
            // amrex::Print() << "TMP_BX " << tmp_bx << std::endl;

            Box tmp_bx_2d(tmp_bx);
            tmp_bx_2d.makeSlab(2, 0);

            Box tmp_bx_1d(tmp_bx);
            tmp_bx_1d.makeSlab(0, 0);
            tmp_bx_1d.makeSlab(1, 0);

            //
            // These are the dimensions of the data we write for only this box
            //
            long long local_nx = tmp_bx.length()[0];
            long long local_ny = tmp_bx.length()[1];
            long long local_nz = tmp_bx.length()[2];

            // We do the "+1" because the offset needs to start at 0
            long long local_start_x = static_cast<long long>(tmp_bx.smallEnd()[0] + 1);
            long long local_start_y = static_cast<long long>(tmp_bx.smallEnd()[1] + 1);
            long long local_start_z = static_cast<long long>(tmp_bx.smallEnd()[2]);

            if (write_header) {
                // Only write out s_rho and s_w at x=0,y=0 to avoid NaNs
                if (bx.contains(IntVect(0,0,0)))
                {
                    {
                        FArrayBox tmp_srho;
                        tmp_srho.resize(tmp_bx_1d, 1, amrex::The_Pinned_Arena());

                        tmp_srho.template copy<RunOn::Device>((*vec_s_r[lev])[mfi.index()], 0, 0, 1);
                        Gpu::streamSynchronize();

                        auto nc_plot_var = ncf.var("s_rho");
                        //nc_plot_var.par_access(NC_INDEPENDENT);
                        nc_plot_var.put(tmp_srho.dataPtr(), { local_start_z }, { local_nz });
                    }
                    {
                        FArrayBox tmp_sw;
                        tmp_sw.resize(convert(tmp_bx_1d,IntVect(0,0,1)), 1, amrex::The_Pinned_Arena());

                        tmp_sw.template copy<RunOn::Device>((*vec_s_w[lev])[mfi.index()], 0, 0, 1);
                        Gpu::streamSynchronize();

                        auto nc_plot_var = ncf.var("s_w");
                        //nc_plot_var.par_access(NC_INDEPENDENT);
                        nc_plot_var.put(tmp_sw.dataPtr(), { local_start_z }, { local_nz + 1});
                    }
                }

                {
                    FArrayBox tmp_bathy;
                    tmp_bathy.resize(tmp_bx_2d, 1, amrex::The_Pinned_Arena());

                    tmp_bathy.template copy<RunOn::Device>((*vec_hOfTheConfusingName[lev])[mfi.index()], 0, 0, 1);
                    Gpu::streamSynchronize();

                    auto nc_plot_var = ncf.var("h");
                    //nc_plot_var.par_access(NC_INDEPENDENT);
                    nc_plot_var.put(tmp_bathy.dataPtr(), { local_start_y, local_start_x }, { local_ny, local_nx });
                }

                {
                    FArrayBox tmp_pm;
                    tmp_pm.resize(tmp_bx_2d, 1, amrex::The_Pinned_Arena());

                    tmp_pm.template copy<RunOn::Device>((*vec_pm[lev])[mfi.index()], 0, 0, 1);
                    Gpu::streamSynchronize();

                    auto nc_plot_var = ncf.var("pm");
                    //nc_plot_var.par_access(NC_INDEPENDENT);

                    nc_plot_var.put(tmp_pm.dataPtr(), { local_start_y, local_start_x }, { local_ny, local_nx });
                }

                {
                    FArrayBox tmp_pn;
                    tmp_pn.resize(tmp_bx_2d, 1, amrex::The_Pinned_Arena());

                    tmp_pn.template copy<RunOn::Device>((*vec_pn[lev])[mfi.index()], 0, 0, 1);
                    Gpu::streamSynchronize();

                    auto nc_plot_var = ncf.var("pn");
                    //nc_plot_var.par_access(NC_INDEPENDENT);
                    nc_plot_var.put(tmp_pn.dataPtr(), { local_start_y, local_start_x }, { local_ny, local_nx });
                }

                {
                    FArrayBox tmp_f;
                    tmp_f.resize(tmp_bx_2d, 1, amrex::The_Pinned_Arena());

                    tmp_f.template copy<RunOn::Device>((*vec_fcor[lev])[mfi.index()], 0, 0, 1);
                    Gpu::streamSynchronize();

                    auto nc_plot_var = ncf.var("f");
                    //nc_plot_var.par_access(NC_INDEPENDENT);
                    nc_plot_var.put(tmp_f.dataPtr(), { local_start_y, local_start_x }, { local_ny, local_nx });
                }

                {
                    FArrayBox tmp_xr;
                    tmp_xr.resize(tmp_bx_2d, 1, amrex::The_Pinned_Arena());

                    tmp_xr.template copy<RunOn::Device>((*vec_xr[lev])[mfi.index()], 0, 0, 1);
                    Gpu::streamSynchronize();

                    auto nc_plot_var = ncf.var("x_rho");
                    //nc_plot_var.par_access(NC_INDEPENDENT);

                    nc_plot_var.put(tmp_xr.dataPtr(), { local_start_y, local_start_x }, { local_ny, local_nx });
                }

                {
                    FArrayBox tmp_yr;
                    tmp_yr.resize(tmp_bx_2d, 1, amrex::The_Pinned_Arena());

                    tmp_yr.template copy<RunOn::Device>((*vec_yr[lev])[mfi.index()], 0, 0, 1);
                    Gpu::streamSynchronize();

                    auto nc_plot_var = ncf.var("y_rho");
                    //nc_plot_var.par_access(NC_INDEPENDENT);
                    nc_plot_var.put(tmp_yr.dataPtr(), { local_start_y, local_start_x }, { local_ny, local_nx });
                }
            }

            {
                FArrayBox tmp_zeta;
                tmp_zeta.resize(tmp_bx_2d, 1, amrex::The_Pinned_Arena());
                tmp_zeta.template copy<RunOn::Device>((*vec_Zt_avg1[lev])[mfi.index()], 0, 0, 1);
                Gpu::streamSynchronize();

                auto nc_plot_var = ncf.var("zeta");
                nc_plot_var.put(tmp_zeta.dataPtr(), { local_start_nt, local_start_y, local_start_x }, { local_nt, local_ny,
                        local_nx });
            }

            {
                FArrayBox tmp_temp;
                tmp_temp.resize(tmp_bx, 1, amrex::The_Pinned_Arena());
                tmp_temp.template copy<RunOn::Device>((*cons_new[lev])[mfi.index()], Temp_comp, 0, 1);
                Gpu::streamSynchronize();

                auto nc_plot_var = ncf.var("temp");
                nc_plot_var.put(tmp_temp.dataPtr(), { local_start_nt, local_start_z, local_start_y, local_start_x }, { local_nt,
                        local_nz, local_ny, local_nx });
            }

            {
                FArrayBox tmp_salt;
                tmp_salt.resize(tmp_bx, 1, amrex::The_Pinned_Arena());
                tmp_salt.template copy<RunOn::Device>((*cons_new[lev])[mfi.index()], Salt_comp, 0, 1);
                Gpu::streamSynchronize();

                auto nc_plot_var = ncf.var("salt");
                nc_plot_var.put(tmp_salt.dataPtr(), { local_start_nt, local_start_z, local_start_y, local_start_x }, { local_nt,
                        local_nz, local_ny, local_nx });
            }

            {
                FArrayBox tmp_tracer;
                tmp_tracer.resize(tmp_bx, 1, amrex::The_Pinned_Arena());
                tmp_tracer.template copy<RunOn::Device>((*cons_new[lev])[mfi.index()], Scalar_comp, 0, 1);
                Gpu::streamSynchronize();

                auto nc_plot_var = ncf.var("tracer");
                nc_plot_var.put(tmp_tracer.dataPtr(), { local_start_nt, local_start_z, local_start_y, local_start_x }, { local_nt,
                        local_nz, local_ny, local_nx });
            }
        } // subdomain
    } // mfi

    //ncf.wait_all(irq, &requests[0]);
    //requests.resize(0);
    //irq = 0;
    // Writing u (we loop over cons to get cell-centered box)
    for (MFIter mfi(*cons_new[lev], false); mfi.isValid(); ++mfi) {
        Box bx = mfi.validbox();

        if (subdomain.contains(bx)) {
            //
            // We only include one grow cell at subdomain boundaries, not internal grid boundaries
            //
            Box tmp_bx(bx);
            tmp_bx.surroundingNodes(0);
            if (tmp_bx.smallEnd()[1] == subdomain.smallEnd()[1])
                tmp_bx.growLo(1, 1);
            if (tmp_bx.bigEnd()[1] == subdomain.bigEnd()[1])
                tmp_bx.growHi(1, 1);
            Box tmp_bx_2d(tmp_bx);
            tmp_bx_2d.makeSlab(2, 0);

            //
            // These are the dimensions of the data we write for only this box
            //
            long long local_nx = tmp_bx.length()[0];
            long long local_ny = tmp_bx.length()[1];
            long long local_nz = tmp_bx.length()[2];

            // We do the "+1" because the offset needs to start at 0
            long long local_start_x = static_cast<long long>(tmp_bx.smallEnd()[0]);
            long long local_start_y = static_cast<long long>(tmp_bx.smallEnd()[1] + 1);
            long long local_start_z = static_cast<long long>(tmp_bx.smallEnd()[2]);

            if (write_header) {
                {
                    FArrayBox tmp;
                    tmp.resize(tmp_bx_2d, 1, amrex::The_Pinned_Arena());
                    tmp.template copy<RunOn::Device>((*vec_xu[lev])[mfi.index()], 0, 0, 1);
                    Gpu::streamSynchronize();

                    auto nc_plot_var = ncf.var("x_u");
                    //nc_plot_var.par_access(NC_INDEPENDENT);
                    nc_plot_var.put(tmp.dataPtr(), { local_start_y, local_start_x }, { local_ny, local_nx });
                }
                {
                    FArrayBox tmp;
                    tmp.resize(tmp_bx_2d, 1, amrex::The_Pinned_Arena());
                    tmp.template copy<RunOn::Device>((*vec_yu[lev])[mfi.index()], 0, 0, 1);
                    Gpu::streamSynchronize();

                    auto nc_plot_var = ncf.var("y_u");
                    //nc_plot_var.par_access(NC_INDEPENDENT);
                    nc_plot_var.put(tmp.dataPtr(), { local_start_y, local_start_x }, { local_ny, local_nx });
                }
            }

            {
                FArrayBox tmp;
                tmp.resize(tmp_bx, 1, amrex::The_Pinned_Arena());
                tmp.template copy<RunOn::Device>((*xvel_new[lev])[mfi.index()], 0, 0, 1);
                Gpu::streamSynchronize();

                auto nc_plot_var = ncf.var("u");
                nc_plot_var.put(tmp.dataPtr(), { local_start_nt, local_start_z, local_start_y, local_start_x }, { local_nt,
                        local_nz, local_ny, local_nx });
            }

            {
                FArrayBox tmp;
                tmp.resize(tmp_bx_2d, 1, amrex::The_Pinned_Arena());
                tmp.template copy<RunOn::Device>((*vec_ubar[lev])[mfi.index()], 0, 0, 1);
                Gpu::streamSynchronize();

                auto nc_plot_var = ncf.var("ubar");
                nc_plot_var.put(tmp.dataPtr(), { local_start_nt, local_start_y, local_start_x }, { local_nt, local_ny, local_nx });
            }
            {
                FArrayBox tmp;
                tmp.resize(tmp_bx_2d, 1, amrex::The_Pinned_Arena());
                tmp.template copy<RunOn::Device>((*vec_sustr[lev])[mfi.index()], 0, 0, 1);
                Gpu::streamSynchronize();

                auto nc_plot_var = ncf.var("sustr");
                nc_plot_var.put(tmp.dataPtr(), { local_start_nt, local_start_y, local_start_x }, { local_nt, local_ny, local_nx });
            }
        } // in subdomain
    } // mfi

    // Writing v (we loop over cons to get cell-centered box)
    for (MFIter mfi(*cons_new[lev], false); mfi.isValid(); ++mfi) {
        Box bx = mfi.validbox();

        if (subdomain.contains(bx)) {
            //
            // We only include one grow cell at subdomain boundaries, not internal grid boundaries
            //
            Box tmp_bx(bx);
            tmp_bx.surroundingNodes(1);
            if (tmp_bx.smallEnd()[0] == subdomain.smallEnd()[0])
                tmp_bx.growLo(0, 1);
            if (tmp_bx.bigEnd()[0] == subdomain.bigEnd()[0])
                tmp_bx.growHi(0, 1);
            // amrex::Print() << "    BX " << bx << std::endl;
            // amrex::Print() << "TMP_BX " << tmp_bx << std::endl;

            Box tmp_bx_2d(tmp_bx);
            tmp_bx_2d.makeSlab(2, 0);

            //
            // These are the dimensions of the data we write for only this box
            //
            long long local_nx = tmp_bx.length()[0];
            long long local_ny = tmp_bx.length()[1];
            long long local_nz = tmp_bx.length()[2];

            // We do the "+1" because the offset needs to start at 0
            long long local_start_x = static_cast<long long>(tmp_bx.smallEnd()[0] + 1);
            long long local_start_y = static_cast<long long>(tmp_bx.smallEnd()[1]);
            long long local_start_z = static_cast<long long>(tmp_bx.smallEnd()[2]);

            if (write_header) {
            {
                FArrayBox tmp;
                tmp.resize(tmp_bx_2d, 1, amrex::The_Pinned_Arena());
                tmp.template copy<RunOn::Device>((*vec_xv[lev])[mfi.index()], 0, 0, 1);
                Gpu::streamSynchronize();

                auto nc_plot_var = ncf.var("x_v");
                //nc_plot_var.par_access(NC_INDEPENDENT);
                nc_plot_var.put(tmp.dataPtr(), { local_start_y, local_start_x }, { local_ny, local_nx });
            }
            {
                FArrayBox tmp;
                tmp.resize(tmp_bx_2d, 1, amrex::The_Pinned_Arena());
                tmp.template copy<RunOn::Device>((*vec_yv[lev])[mfi.index()], 0, 0, 1);
                Gpu::streamSynchronize();

                auto nc_plot_var = ncf.var("y_v");
                //nc_plot_var.par_access(NC_INDEPENDENT);
                nc_plot_var.put(tmp.dataPtr(), { local_start_y, local_start_x }, { local_ny, local_nx });
            }

            }
            {
                FArrayBox tmp;
                tmp.resize(tmp_bx, 1, amrex::The_Pinned_Arena());
                tmp.template copy<RunOn::Device>((*yvel_new[lev])[mfi.index()], 0, 0, 1);
                Gpu::streamSynchronize();

                auto nc_plot_var = ncf.var("v");
                nc_plot_var.put(tmp.dataPtr(), { local_start_nt, local_start_z, local_start_y, local_start_x }, { local_nt,
                        local_nz, local_ny, local_nx });
            }

            {
                FArrayBox tmp;
                tmp.resize(tmp_bx_2d, 1, amrex::The_Pinned_Arena());
                tmp.template copy<RunOn::Device>((*vec_vbar[lev])[mfi.index()], 0, 0, 1);
                Gpu::streamSynchronize();

                auto nc_plot_var = ncf.var("vbar");
                nc_plot_var.put(tmp.dataPtr(), { local_start_nt, local_start_y, local_start_x }, { local_nt, local_ny, local_nx });
            }
            {
                FArrayBox tmp;
                tmp.resize(tmp_bx_2d, 1, amrex::The_Pinned_Arena());
                tmp.template copy<RunOn::Device>((*vec_svstr[lev])[mfi.index()], 0, 0, 1);
                Gpu::streamSynchronize();

                auto nc_plot_var = ncf.var("svstr");
                nc_plot_var.put(tmp.dataPtr(), { local_start_nt, local_start_y, local_start_x }, { local_nt, local_ny, local_nx });
            }

        } // in subdomain
    } // mfi

    for (MFIter mfi(*cons_new[lev], false); mfi.isValid(); ++mfi) {
        Box bx = mfi.validbox();

        if (subdomain.contains(bx)) {
            //
            // We only include one grow cell at subdomain boundaries, not internal grid boundaries
            //
            Box tmp_bx(bx);
            tmp_bx.surroundingNodes(0);
            tmp_bx.surroundingNodes(1);

            Box tmp_bx_2d(tmp_bx);
            tmp_bx_2d.makeSlab(2, 0);

            //
            // These are the dimensions of the data we write for only this box
            //
            long long local_nx = tmp_bx.length()[0];
            long long local_ny = tmp_bx.length()[1];
            long long local_nz = tmp_bx.length()[2];

            // We do the "+1" because the offset needs to start at 0
            long long local_start_x = static_cast<long long>(tmp_bx.smallEnd()[0]);
            long long local_start_y = static_cast<long long>(tmp_bx.smallEnd()[1]);
            long long local_start_z = static_cast<long long>(tmp_bx.smallEnd()[2]);

            if (write_header) {
            {
                FArrayBox tmp;
                tmp.resize(tmp_bx_2d, 1, amrex::The_Pinned_Arena());
                tmp.template copy<RunOn::Device>((*vec_xp[lev])[mfi.index()], 0, 0, 1);
                Gpu::streamSynchronize();

                auto nc_plot_var = ncf.var("x_psi");
                //nc_plot_var.par_access(NC_INDEPENDENT);
                nc_plot_var.put(tmp.dataPtr(), { local_start_y, local_start_x }, { local_ny, local_nx });
            }
            {
                FArrayBox tmp;
                tmp.resize(tmp_bx_2d, 1, amrex::The_Pinned_Arena());
                tmp.template copy<RunOn::Device>((*vec_yp[lev])[mfi.index()], 0, 0, 1);
                Gpu::streamSynchronize();

                auto nc_plot_var = ncf.var("y_psi");
                //nc_plot_var.par_access(NC_INDEPENDENT);
                nc_plot_var.put(tmp.dataPtr(), { local_start_y, local_start_x }, { local_ny, local_nx });
            }

            } // header
        } // in subdomain
    } // mfi

    mask_arrays_for_write(lev, 0.0_rt);

    ncf.close();

    REMORA::total_nc_plot_file_step += 1;
}

