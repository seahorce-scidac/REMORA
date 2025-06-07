#include "REMORA_NCTimeSeriesRiver.H"
#include "REMORA_NCFile.H"

#include "AMReX_ParallelDescriptor.H"

#include <string>

#ifdef REMORA_USE_NETCDF
NCTimeSeriesRiver::NCTimeSeriesRiver (const std::string a_file_name, const std::string a_field_name,
                                      const std::string a_time_name,
                                      const int a_nz) {
    file_name = a_file_name;
    time_name = a_time_name;
    field_name = a_field_name;
    nz   = a_nz;
}

void NCTimeSeriesRiver::Initialize() {
    // open file
    amrex::Print() << "Loading " << field_name << " from rivers NetCDF file " << file_name << std::endl;

    // The time field can have any number of names, depending on the field.
    // If not specified in input file (time_name.empty()) then set it by default
    if (time_name.empty())
    {
        time_name = "river_time";
    }

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
            // amrex::Print() << "TIMES " << river_times[nt] << std::endl;
        }
    }
    int ntimes = river_times.size();
    int ioproc = amrex::ParallelDescriptor::IOProcessorNumber();
    amrex::ParallelDescriptor::Bcast(&ntimes,1,ioproc);
    if (!(amrex::ParallelDescriptor::IOProcessor())) {
        river_times.resize(ntimes);
    }
    amrex::ParallelDescriptor::Bcast(river_times.data(), river_times.size(), ioproc);

    auto ncf = ncutils::NCFile::open(file_name, NC_NOCLOBBER );
    if (amrex::ParallelDescriptor::IOProcessor())
    {
        std::vector<MPI_Offset> shape = ncf.var(field_name).shape();
        if (shape.size() == 2) {
            has_z = 0;
        } else if (shape.size() == 3) {
            has_z = 1;
        } else {
            amrex::Abort("River field shape not 2 or 3");
        }
        nriv = shape[0];
    }
    amrex::ParallelDescriptor::Bcast(&has_z, 1, ioproc);
    amrex::ParallelDescriptor::Bcast(&nriv, 1, ioproc);

    amrex::Box vshape_box(amrex::IntVect(0,0,0),amrex::IntVect(nriv,0,nz));
    if (!has_z) {
        amrex::Vector<amrex::FArrayBox*> NC_fabs;
        amrex::Vector<std::string> NC_names;
        amrex::Vector<enum NC_Data_Dims_Type> NC_dim_types;
        amrex::Print() << "Reading in river_Vshape from " << file_name << std::endl;

        NC_fabs.push_back(fab_vshape); NC_names.push_back("river_Vshape");
        NC_dim_types.push_back(NC_Data_Dims_Type::BT_Riv);
        BuildFABsFromNetCDFFile<amrex::FArrayBox,amrex::Real>(vshape_box, file_name, NC_names, NC_dim_types, NC_fabs);
    }

    amrex::Box riv_box(amrex::IntVect(0,0,0),amrex::IntVect(nriv,0,nz));
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

    auto interp_array = fab_interp->array();
    auto before_array = fab_before->array();
    auto after_array = fab_after->array();
    for (int r=0; r < nriv; r++) {
        for (int k=0; k < nz; k++) {
            interp_array(r,0,k) = before_array(r,0,k) + (time - time_before_copy) * (after_array(r,0,k) - before_array(r,0,k)) / dt;
        }
    }
}

void NCTimeSeriesRiver::read_in_at_time (amrex::FArrayBox* fab_dat, int itime) {
    amrex::FArrayBox NC_fab;
    amrex::Vector<amrex::FArrayBox*> NC_fabs;
    amrex::Vector<std::string> NC_names;
    amrex::Vector<enum NC_Data_Dims_Type> NC_dim_types;
    // actual dims don't really matter here; only lower is used in call
    amrex::Box riv_domain(amrex::IntVect(0,0,0), amrex::IntVect(nriv,1,nz));

    amrex::Print() << "Reading in " << field_name << " at  time index " << itime << " from " << file_name << std::endl;

    NC_fabs.push_back(&NC_fab) ; NC_names.push_back(field_name);

    if (has_z) {
        NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_Riv);
    } else {
        NC_dim_types.push_back(NC_Data_Dims_Type::Time_Riv);
    }

    BuildFABsFromNetCDFFile<amrex::FArrayBox,amrex::Real>(riv_domain, file_name, NC_names, NC_dim_types, NC_fabs, true, itime);

    auto dat_array = fab_dat->array();
    auto tmp_array = NC_fabs[0]->array();
    for (int r=0; r < nriv; r++) {
        if (has_z) {
            for (int k=0; k < nz; k++) {
                dat_array(r,0,k) = tmp_array(r,0,k);
            }
        } else {
            auto array_vshape = fab_vshape->array();
            for (int k=0; k < nz; k++) {
                dat_array(r,0,k) = tmp_array(r,0,k) * array_vshape(r,0,k);
            }
        }
    }

//    using RARRAY = amrex::NDArray<Real>;
//    amrex::Vector<RARRAY> array_dat(1);
//
//    int nz_read; // Number of z points we will have read in
//    auto ncf = ncutils::NCFile::open(file_name, NC_NOCLOBBER );
//    ncmpi_begin_indep_data(ncf.ncid);
//    if (amrex::ParallelDescriptor::IOProcessor())
//    {
//        std::vector<MPI_Offset> shape = ncf.var(field_name).shape();
//        array_dat = NDArray<DType>(field_name, shape);
//        DType* dataPtr = array_dat.get_data();
//
//        std::vector<MPI_Offset> start(shape.size(), 0);
//        start[0] = itime;
//        shape[0] = 1;
//        // Don't have z-data
//        if (shape.size() == 2) {
//            nz_read = 1;
//        } else if (shape.size() == 3) {
//            nz_read = shape[1];
//        } else {
//            amrex::Abort("Unexpected river data dimension count");
//        }
//
//        ncf.var(field_name).get(dataPtr, start, shape);
//
//        if (nz_read == 1) {
//            RARRAY vshape_dat = NDArray<DType>(field_name, shape);
//        }
//
//        nriv = array_dat[0].get_vshape()[1];
//        nz = array_dat[0].get_vshape()[1];
//
//        if (nz_read == 1) {
//            for (int k=0; k < nz; k++) {
//                for (int r=0; r < nriv; r++) {
//                    vec[r][k] = *(array_dat.get_data() + r) * vshape[r][k];
//                }
//            }
//        } else {
//            for (int k=0; k < nz; k++) {
//                for (int r=0; r < nriv; r++) {
//                    vec[r][k] = *(array_dat.get_data() + r + k * nriv);
//                }
//            }
//        }
//    }
//    amrex::ParallelDescriptor::Bcast(vec.data(), vec.size(), ioproc);

}
#endif // REMORA_USE_NETCDF
