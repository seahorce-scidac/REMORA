#include "REMORA_NCTimeSeriesBoundary.H"
#include "REMORA_NCFile.H"

#include "AMReX_ParallelDescriptor.H"

#include <string>

#ifdef REMORA_USE_NETCDF
/**
 * @param[in   ] a_file_name          file name to read from
 * @param[in   ] a_field_name         name of field to read in
 * @param[in   ] a_time_name          name of time variable in NetCDF file
 * @param[in   ] a_domain             simulation domain
 * @param[in   ] a_is2d               Whether the variable we're working with is 2D
 */
NCTimeSeriesBoundary::NCTimeSeriesBoundary (const amrex::Vector<std::string>& a_file_names, const std::string a_field_name,
                                            const std::string a_time_name,
                                            const amrex::Box& a_domain,
                                            const amrex::IntVect a_index_type,
                                            const amrex::GpuArray<bool, AMREX_SPACEDIM*2>* a_var_need_data,
                                            bool a_is2d) {
    file_names.assign(a_file_names.begin(), a_file_names.end());
    time_name = a_time_name;
    field_name = a_field_name;
    domain = a_domain;
    index_type = a_index_type;
    var_need_data = *a_var_need_data;
    is2d = a_is2d;
}

void NCTimeSeriesBoundary::Initialize() {
    // open file
    amrex::Print() << "Loading boundary data for " << field_name << " from NetCDF file " << std::endl;

    // The time field can have any number of names, depending on the field.
    // If not specified in input file (time_name.empty()) then set it by default
    if (time_name.empty())
    {
        time_name = "ocean_time";
    }

    for (int ifile = 0; ifile < file_names.size(); ifile++) {
        std::string file_name = file_names[ifile];
        // Check units of time stamps; should be days
        std::string unit_str = ReadNetCDFVarAttrStr(file_name, time_name, "units"); // works on proc 0
        if (amrex::ParallelDescriptor::IOProcessor())
        {
            if (unit_str.find("days") == std::string::npos) {
                amrex::Print() << "Units of ocean_time given as: " << unit_str << std::endl;
                amrex::Abort("Units must be in days.");
            }
        }
        // get times and put in array
        using RARRAY = NDArray<amrex::Real>;
        amrex::Vector<RARRAY> array_ts(1);
        ReadNetCDFFile(file_name, {time_name}, array_ts); // filled only on proc 0
        if (amrex::ParallelDescriptor::IOProcessor())
        {
            int ntimes_io = array_ts[0].get_vshape()[0];
            for (int nt(0); nt < ntimes_io; nt++)
            {
                // Convert ocean time from days to seconds
                bry_times.push_back((*(array_ts[0].get_data() + nt)) * amrex::Real(60.0) * amrex::Real(60.0) * amrex::Real(24.0));
                file_for_time.push_back(ifile);
                file_itime_offset.push_back(nt);
                // amrex::Print() << "TIMES " << bry_times[nt] << std::endl;
            }
        }
    }

    AMREX_ASSERT(std::is_sorted(bry_times.begin(), bry_times.end()));

    int ntimes = bry_times.size();
    int ioproc = amrex::ParallelDescriptor::IOProcessorNumber();
    amrex::ParallelDescriptor::Bcast(&ntimes,1,ioproc);
    if (!(amrex::ParallelDescriptor::IOProcessor())) {
        bry_times.resize(ntimes);
    }
    amrex::ParallelDescriptor::Bcast(bry_times.data(), bry_times.size(), ioproc);

    // Initialize Fabs
    amrex::Arena* Arena_Used = amrex::The_Arena();
#ifdef AMREX_USE_GPU
    Arena_Used = amrex::The_Pinned_Arena();
#endif
    const auto& lo = domain.loVect();
    const auto& hi = domain.hiVect();

    amrex::Box xlo_bx(amrex::IntVect(lo[0]+index_type[0]-1, lo[1]+index_type[1]-1, lo[2]),
                      amrex::IntVect(lo[0]+index_type[0]-1, hi[1]+1              , hi[2]), index_type);
    amrex::Box xhi_bx(amrex::IntVect(hi[0]+1              , lo[1]+index_type[1]-1, lo[2]),
                      amrex::IntVect(hi[0]+1              , hi[1]+1              , hi[2]), index_type);
    amrex::Box ylo_bx(amrex::IntVect(lo[0]+index_type[0]-1, lo[1]+index_type[1]-1, lo[2]),
                      amrex::IntVect(hi[0]+1              , lo[1]+index_type[1]-1, hi[2]), index_type);
    amrex::Box yhi_bx(amrex::IntVect(lo[0]+index_type[0]-1, hi[1]+1              , lo[2]),
                      amrex::IntVect(hi[0]+1              , hi[1]+1              , hi[2]), index_type);
    if (is2d) {
        xlo_bx.makeSlab(2,0);
        xhi_bx.makeSlab(2,0);
        ylo_bx.makeSlab(2,0);
        yhi_bx.makeSlab(2,0);
    }

    amrex::Print() << xlo_bx << " " << xhi_bx << " " << ylo_bx << " " << yhi_bx << std::endl;

    xlo_dat_before = amrex::FArrayBox(xlo_bx, 1, Arena_Used);
    xhi_dat_before = amrex::FArrayBox(xhi_bx, 1, Arena_Used);
    ylo_dat_before = amrex::FArrayBox(ylo_bx, 1, Arena_Used);
    yhi_dat_before = amrex::FArrayBox(yhi_bx, 1, Arena_Used);

    xlo_dat_after = amrex::FArrayBox(xlo_bx, 1, Arena_Used);
    xhi_dat_after = amrex::FArrayBox(xhi_bx, 1, Arena_Used);
    ylo_dat_after = amrex::FArrayBox(ylo_bx, 1, Arena_Used);
    yhi_dat_after = amrex::FArrayBox(yhi_bx, 1, Arena_Used);

    xlo_dat_interp = amrex::FArrayBox(xlo_bx, 1, Arena_Used);
    xhi_dat_interp = amrex::FArrayBox(xhi_bx, 1, Arena_Used);
    ylo_dat_interp = amrex::FArrayBox(ylo_bx, 1, Arena_Used);
    yhi_dat_interp = amrex::FArrayBox(yhi_bx, 1, Arena_Used);

    // dummy initialization
    i_time_before = -100;

    if (var_need_data[amrex::Orientation(amrex::Direction::x,amrex::Orientation::low)] == true) {
        nc_var_names.push_back(field_name + "_west");
    }
    if (var_need_data[amrex::Orientation(amrex::Direction::x,amrex::Orientation::high)] == true) {
        nc_var_names.push_back(field_name + "_east");
    }
    if (var_need_data[amrex::Orientation(amrex::Direction::y,amrex::Orientation::low)] == true) {
        nc_var_names.push_back(field_name + "_south");
    }
    if (var_need_data[amrex::Orientation(amrex::Direction::y,amrex::Orientation::high)] == true) {
        nc_var_names.push_back(field_name + "_north");
    }
}

/**
 * @param time   time to interpolate to
 */
void NCTimeSeriesBoundary::update_interpolated_to_time (amrex::Real time) {
    // Figure out time index:
    AMREX_ASSERT(time >= bry_times[0]);
    AMREX_ASSERT(time <= bry_times[bry_times.size()-1]);
    int i_time_before_old = i_time_before;
    for (int nt=0; nt < bry_times.size()-1; nt++) {
        if ((bry_times[nt] <= time) and (bry_times[nt+1] >= time)) {
            i_time_before = nt;
            time_before = bry_times[nt];
            time_after = bry_times[nt+1];
            break;
        }
    }

    int i_time_after = i_time_before + 1;
    if (i_time_before_old + 1 == i_time_before) {
        // swap multifabs so we only have to read in one MultiFab
        std::swap(xlo_dat_before, xlo_dat_after);
        std::swap(xhi_dat_before, xhi_dat_after);
        std::swap(ylo_dat_before, ylo_dat_after);
        std::swap(yhi_dat_before, yhi_dat_after);
        read_in_at_time(xlo_dat_after,xhi_dat_after, ylo_dat_after, yhi_dat_after, i_time_after);
    } else if (i_time_before_old != i_time_before) {
        read_in_at_time(xlo_dat_after,xhi_dat_after, ylo_dat_after, yhi_dat_after, i_time_after);
        read_in_at_time(xlo_dat_before,xhi_dat_before, ylo_dat_before,yhi_dat_before, i_time_before);
    }

    amrex::Real dt = time_after - time_before;

    amrex::Real time_before_copy = time_before;

    amrex::Array4<amrex::Real> xlo_interp_arr = xlo_dat_interp.array();
    amrex::Array4<amrex::Real> xhi_interp_arr = xhi_dat_interp.array();
    amrex::Array4<amrex::Real> ylo_interp_arr = ylo_dat_interp.array();
    amrex::Array4<amrex::Real> yhi_interp_arr = yhi_dat_interp.array();

    amrex::Array4<amrex::Real> xlo_before_arr = xlo_dat_before.array();
    amrex::Array4<amrex::Real> xhi_before_arr = xhi_dat_before.array();
    amrex::Array4<amrex::Real> ylo_before_arr = ylo_dat_before.array();
    amrex::Array4<amrex::Real> yhi_before_arr = yhi_dat_before.array();

    amrex::Array4<amrex::Real> xlo_after_arr = xlo_dat_after.array();
    amrex::Array4<amrex::Real> xhi_after_arr = xhi_dat_after.array();
    amrex::Array4<amrex::Real> ylo_after_arr = ylo_dat_after.array();
    amrex::Array4<amrex::Real> yhi_after_arr = yhi_dat_after.array();

    amrex::ParallelFor(xlo_dat_interp.box(), [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        xlo_interp_arr(i,j,k) = xlo_before_arr(i,j,k) + (time - time_before_copy) * (xlo_after_arr(i,j,k) - xlo_before_arr(i,j,k)) / dt;
    });

    amrex::ParallelFor(xhi_dat_interp.box(), [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        xhi_interp_arr(i,j,k) = xhi_before_arr(i,j,k) + (time - time_before_copy) * (xhi_after_arr(i,j,k) - xhi_before_arr(i,j,k)) / dt;
    });

    amrex::ParallelFor(ylo_dat_interp.box(), [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        ylo_interp_arr(i,j,k) = ylo_before_arr(i,j,k) + (time - time_before_copy) * (ylo_after_arr(i,j,k) - ylo_before_arr(i,j,k)) / dt;
    });

    amrex::ParallelFor(yhi_dat_interp.box(), [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        yhi_interp_arr(i,j,k) = yhi_before_arr(i,j,k) + (time - time_before_copy) * (yhi_after_arr(i,j,k) - yhi_before_arr(i,j,k)) / dt;
    });
}

/**
 * @param[inout] fab_xlo   fab to store xlo boundary data into
 * @param[inout] fab_xhi   fab to store xhi boundary data into
 * @param[inout] fab_ylo   fab to store ylo boundary data into
 * @param[inout] fab_yhi   fab to store yhi boundary data into
 * @param[in   ] itime     index of time step to read from file
 */
void NCTimeSeriesBoundary::read_in_at_time (amrex::FArrayBox& fab_xlo,
                                            amrex::FArrayBox& fab_xhi,
                                            amrex::FArrayBox& fab_ylo,
                                            amrex::FArrayBox& fab_yhi,
                                            int itime) {
    using RARRAY = NDArray<amrex::Real>;
    amrex::Vector<RARRAY> arrays(nc_var_names.size());

    // The width of the boundary region we need to read is 1
    int width = 1;

    amrex::Print() << "Reading in " << field_name << " at time " << bry_times[itime] << std::endl;
    std::string nc_bdry_file = file_names[file_for_time[itime]];
    int itime_offset = file_itime_offset[itime];
    ReadNetCDFFile(nc_bdry_file, nc_var_names, arrays, true, itime_offset); // does work on proc 0 only

    for (int iv=0; iv < nc_var_names.size(); iv++) {
        std::string  last4 = nc_var_names[iv].substr(nc_var_names[iv].size()-4, 4);
        std::string  last5 = nc_var_names[iv].substr(nc_var_names[iv].size()-5, 5);
        int nx, ny, nz, n_plane;
        int i, j, k, ioff, joff;
        if (last4 == "west") {
            amrex::Box my_box = fab_xlo.box();
            if (is2d) {
                nz = 1;
                ny = arrays[iv].get_vshape()[1];
            } else {
                nz = arrays[iv].get_vshape()[1];
                ny = arrays[iv].get_vshape()[2];
            }
            n_plane = ny * nz;
            AMREX_ALWAYS_ASSERT(my_box.numPts() == n_plane);

            i    = my_box.smallEnd()[0];
            joff = my_box.smallEnd()[1];

            amrex::Array4<amrex::Real> fab_arr = fab_xlo.array();
            for (int n(0); n < n_plane; n++) {
                k = n / ny;
                j = n - (k * ny);
                fab_arr(i, j+joff, k, 0) = static_cast<amrex::Real>(*(arrays[iv].get_data() + n));
            }
        } else if (last4 == "east") {
            amrex::Box my_box = fab_xhi.box();
            if (is2d) {
                nz = 1;
                ny = arrays[iv].get_vshape()[1];
            } else {
                nz = arrays[iv].get_vshape()[1];
                ny = arrays[iv].get_vshape()[2];
            }
            n_plane = ny * nz;
            AMREX_ALWAYS_ASSERT(my_box.numPts() == n_plane);

            i    = my_box.smallEnd()[0];
            joff = my_box.smallEnd()[1];

            amrex::Array4<amrex::Real> fab_arr = fab_xhi.array();
            for (int n(0); n < n_plane; n++) {
                k = n / ny;
                j = n - (k * ny);
                fab_arr(i, j+joff, k, 0) = static_cast<amrex::Real>(*(arrays[iv].get_data() + n));
            }
        } else if (last5 == "south") {
            amrex::Box my_box = fab_ylo.box();
            if (is2d) {
                nz = 1;
                nx = arrays[iv].get_vshape()[1];
            } else {
                nz = arrays[iv].get_vshape()[1];
                nx = arrays[iv].get_vshape()[2];
            }
            n_plane = nx * nz;
            AMREX_ALWAYS_ASSERT(my_box.numPts() == n_plane);

            j    = my_box.smallEnd()[1];
            ioff = my_box.smallEnd()[0];

            amrex::Array4<amrex::Real> fab_arr = fab_ylo.array();
            for (int n(0); n < n_plane; n++) {
                k = n / nx;
                i = n - (k * nx);
                fab_arr(i+ioff, j, k, 0) = static_cast<amrex::Real>(*(arrays[iv].get_data() + n));
            }
        } else if (last5 == "north") {
            amrex::Box my_box = fab_yhi.box();
            if (is2d) {
                nz = 1;
                nx = arrays[iv].get_vshape()[1];
            } else {
                nz = arrays[iv].get_vshape()[1];
                nx = arrays[iv].get_vshape()[2];
            }
            n_plane = nx * nz;
            AMREX_ALWAYS_ASSERT(my_box.numPts() == n_plane);

            j    = my_box.smallEnd()[1];
            ioff = my_box.smallEnd()[0];

            amrex::Array4<amrex::Real> fab_arr = fab_yhi.array();
            for (int n(0); n < n_plane; n++) {
                k = n / nx;
                i = n - (k * nx);
                fab_arr(i+ioff, j, k, 0) = static_cast<amrex::Real>(*(arrays[iv].get_data() + n));
            }
        }
    }

    amrex::ParallelDescriptor::Barrier();

    // When an FArrayBox is built, space is allocated on every rank.  However, we only
    //    filled the data in these FABs on the IOProcessor.  So here we broadcast
    //    the data to every rank.
    int ioproc = amrex::ParallelDescriptor::IOProcessorNumber();  // I/O rank
    amrex::ParallelDescriptor::Bcast(fab_xlo.dataPtr(),fab_xlo.box().numPts(),ioproc);
    amrex::ParallelDescriptor::Bcast(fab_xhi.dataPtr(),fab_xhi.box().numPts(),ioproc);
    amrex::ParallelDescriptor::Bcast(fab_ylo.dataPtr(),fab_ylo.box().numPts(),ioproc);
    amrex::ParallelDescriptor::Bcast(fab_yhi.dataPtr(),fab_yhi.box().numPts(),ioproc);

}
#endif // REMORA_USE_NETCDF
