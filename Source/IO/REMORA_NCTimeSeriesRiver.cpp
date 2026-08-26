#include "REMORA_NCTimeSeriesRiver.H"
#include "REMORA_NCFile.H"

#include "AMReX_ParallelDescriptor.H"

#include <string>

#ifdef REMORA_USE_NETCDF

/**
 * @param[in   ] a_file_names         vector of file name(s) to read from
 * @param[in   ] a_field_name         name of field to read in
 * @param[in   ] a_time_name          name of time variable in NetCDF file
 * @param[in   ] a_nz                 number of vertical levels in domain
 * @param[in   ] a_use_vert_integ     whether the data in the file is vertically integrated
 * @param[in   ] a_is_transport       whether the field is the river transport rather than a tracer
 */
NCTimeSeriesRiver::NCTimeSeriesRiver (const amrex::Vector<std::string>& a_file_names, const std::string a_field_name,
                                      const std::string a_time_name,
                                      const int a_nz, const int a_use_vert_integ,
                                      const int a_is_transport) {
    file_names.assign(a_file_names.begin(), a_file_names.end());
    time_name = a_time_name;
    field_name = a_field_name;
    nz   = a_nz;
    use_vert_integ = a_use_vert_integ;
    is_transport = a_is_transport;
}

void NCTimeSeriesRiver::Initialize() {
    // open file
    amrex::Print() << "Loading " << field_name << " from rivers NetCDF file(s)" << std::endl;

    // The time field can have any number of names, depending on the field.
    // If not specified in input file (time_name.empty()) then set it by default
    if (time_name.empty())
    {
        time_name = "river_time";
    }

    for (int ifile = 0; ifile < file_names.size(); ++ifile) {
        const std::string& file_name = file_names[ifile];

        // Check units of time stamps; should be days
        std::string unit_str = ReadNetCDFVarAttrStr(file_name, time_name, "units"); // works on proc 0
        if (amrex::ParallelDescriptor::IOProcessor())
        {
            if (unit_str.find("days") == std::string::npos) {
                amrex::Print() << "Units of river_time given as: " << unit_str << std::endl;
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
                // Convert river time from days to seconds
                river_times.push_back((*(array_ts[0].get_data() + nt)) * amrex::Real(60.0) * amrex::Real(60.0) * amrex::Real(24.0));
                file_for_time.push_back(ifile);
                file_itime_offset.push_back(nt);
            }
        }
    }
    int ntimes = river_times.size();
    // Only do checks on IO processors since river_times isn't populated on other ranks yet
    if (amrex::ParallelDescriptor::IOProcessor()) {
        AMREX_ASSERT(std::is_sorted(river_times.begin(), river_times.end()));
        if (ntimes <= 1) {
            amrex::Error("River data must be given at at least two times");
        }
    }
    int ioproc = amrex::ParallelDescriptor::IOProcessorNumber();
    amrex::ParallelDescriptor::Bcast(&ntimes,1,ioproc);
    if (!(amrex::ParallelDescriptor::IOProcessor())) {
        river_times.resize(ntimes);
        file_for_time.resize(ntimes);
        file_itime_offset.resize(ntimes);
    }
    amrex::ParallelDescriptor::Bcast(river_times.data(), river_times.size(), ioproc);
    amrex::ParallelDescriptor::Bcast(file_for_time.data(), file_for_time.size(), ioproc);
    amrex::ParallelDescriptor::Bcast(file_itime_offset.data(), file_itime_offset.size(), ioproc);

    auto ncf = ncutils::NCFile::open(file_names[0], NC_NOCLOBBER );
    if (amrex::ParallelDescriptor::IOProcessor())
    {
        std::vector<MPI_Offset> shape = ncf.var(field_name).shape();
        if (shape.size() == 2) {
            has_z = 0;
            nriv = shape[1];
        } else if (shape.size() == 3) {
            has_z = 1;
            nriv = shape[2];
        } else {
            amrex::Abort("River field shape not 2 or 3");
        }
    }
    amrex::ParallelDescriptor::Bcast(&has_z, 1, ioproc);
    amrex::ParallelDescriptor::Bcast(&nriv, 1, ioproc);

    // river_Vshape distributes the total transport in the vertical (ROMS Qsrc = Qbar*Qshape),
    // so it applies only to the transport. A tracer given as (river_time, river) is a
    // concentration, and is used unscaled at every level.
    if (!has_z && !is_transport) {
        amrex::Print() << "Warning: " << field_name << " has no s_rho dimension in "
                       << file_names[0] << "; the same value will be used at every "
                       << "vertical level. ROMS expects river tracers to be given as "
                       << "(river_time, s_rho, river)." << std::endl;
    }

    amrex::Box vshape_box(amrex::IntVect(0,0,0),amrex::IntVect(nriv,0,nz));
    if (!has_z && !use_vert_integ && is_transport) {
        amrex::Vector<amrex::FArrayBox*> NC_fabs;
        amrex::Vector<std::string> NC_names;
        amrex::Vector<enum NC_Data_Dims_Type> NC_dim_types;
        amrex::Print() << "Reading in river_Vshape from " << file_names[0] << std::endl;

        fab_vshape = new amrex::FArrayBox();
        NC_fabs.push_back(fab_vshape); NC_names.push_back("river_Vshape");
        NC_dim_types.push_back(NC_Data_Dims_Type::BT_Riv);
        BuildFABsFromNetCDFFile<amrex::FArrayBox,amrex::Real>(vshape_box, file_names[0], NC_names, NC_dim_types, NC_fabs);
    }

    nzbox = (use_vert_integ) ? 1 : nz;
    amrex::Box riv_box(amrex::IntVect(0,0,0),amrex::IntVect(nriv,0,nzbox));
#ifdef AMREX_USE_GPU
    // It's possible there should be a different arena
    fab_before = new amrex::FArrayBox(riv_box,1,amrex::The_Pinned_Arena());
    fab_after = new amrex::FArrayBox(riv_box,1,amrex::The_Pinned_Arena());
    fab_interp = new amrex::FArrayBox(riv_box,1,amrex::The_Pinned_Arena());
#else
    fab_before = new amrex::FArrayBox(riv_box,1);
    fab_after = new amrex::FArrayBox(riv_box,1);
    fab_interp = new amrex::FArrayBox(riv_box,1);
#endif

    // dummy initialization
    i_time_before = -100;
}

void NCTimeSeriesRiver::update_interpolated_to_time (amrex::Real time) {
    // Figure out time index:
    AMREX_ASSERT(time >= river_times[0]);
    AMREX_ASSERT(time <= river_times[river_times.size()-1]);
    int i_time_before_old = i_time_before;
    for (int nt=0; nt < river_times.size()-1; nt++) {
        if ((river_times[nt] <= time) and (river_times[nt+1] >= time)) {
            i_time_before = nt;
            time_before = river_times[nt];
            time_after = river_times[nt+1];
            break;
        }
    }

    int i_time_after = i_time_before + 1;
    if (i_time_before_old + 1 == i_time_before) {
        // swap data vectors so we only have to read in one MultiFab
        std::swap(fab_before, fab_after);
        read_in_at_time(fab_after, i_time_after);
    } else if (i_time_before_old != i_time_before) {
        read_in_at_time(fab_after,  i_time_after);
        read_in_at_time(fab_before, i_time_before);
    }

    amrex::Real dt = time_after - time_before;
    amrex::Real time_before_copy = time_before;

    amrex::Box fab_domain(amrex::IntVect(0,0,0), amrex::IntVect(nriv-1,0,nzbox-1));
    auto interp_array = fab_interp->array();
    auto before_array = fab_before->array();
    auto after_array = fab_after->array();
    amrex::ParallelFor(fab_domain, [=] AMREX_GPU_DEVICE (int r, int , int k) {
        interp_array(r,0,k) = before_array(r,0,k) + (time - time_before_copy) * (after_array(r,0,k) - before_array(r,0,k)) / dt;
    });
}

void NCTimeSeriesRiver::read_in_at_time (amrex::FArrayBox* fab_dat, int itime) {
    amrex::FArrayBox NC_fab;
    amrex::Vector<amrex::FArrayBox*> NC_fabs;
    amrex::Vector<std::string> NC_names;
    amrex::Vector<enum NC_Data_Dims_Type> NC_dim_types;
    // actual dims don't really matter here; only lower is used in call

    const std::string& file_name = file_names[file_for_time[itime]];
    const int itime_offset = file_itime_offset[itime];

    amrex::Print() << "Reading in " << field_name << " at time index " << itime
                   << " from " << file_name << std::endl;

    NC_fabs.push_back(&NC_fab) ; NC_names.push_back(field_name);

    if (has_z) {
        NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_Riv);
    } else {
        NC_dim_types.push_back(NC_Data_Dims_Type::Time_Riv);
    }

    amrex::Box riv_domain(amrex::IntVect(0,0,0), amrex::IntVect(nriv-1,0,nz-1));
    amrex::Box fab_domain(amrex::IntVect(0,0,0), amrex::IntVect(nriv-1,0,nzbox-1));

    BuildFABsFromNetCDFFile<amrex::FArrayBox,amrex::Real>(riv_domain, file_name, NC_names, NC_dim_types, NC_fabs, true, itime_offset);

    auto dat_array = fab_dat->array();
    auto tmp_array = NC_fabs[0]->const_array();
    if (has_z || use_vert_integ) {
        amrex::ParallelFor(fab_domain, [=] AMREX_GPU_DEVICE (int r, int , int k) {
            dat_array(r,0,k) = tmp_array(r,0,k);
        });
    } else if (is_transport) {
        // Distribute the vertically integrated transport over the levels
        auto array_vshape = fab_vshape->const_array();
        amrex::ParallelFor(fab_domain, [=] AMREX_GPU_DEVICE (int r, int , int k) {
            dat_array(r,0,k) = tmp_array(r,0,0) * array_vshape(r,0,k);
        });
    } else {
        // Concentration given only as (river_time, river); use it at every level
        amrex::ParallelFor(fab_domain, [=] AMREX_GPU_DEVICE (int r, int , int k) {
            dat_array(r,0,k) = tmp_array(r,0,0);
        });
    }
}
#endif // REMORA_USE_NETCDF
