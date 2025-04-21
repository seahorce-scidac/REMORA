#include "REMORA_NCTimeSeries.H"
#include "REMORA_NCFile.H"

#include "AMReX_ParallelDescriptor.H"

#include <string>

#ifdef REMORA_USE_NETCDF
NCTimeSeries::NCTimeSeries (const std::string a_file_name, const std::string a_field_name,
                            const std::string a_time_name,
                            const amrex::Box& a_domain,
                            amrex::MultiFab* a_mf_var, bool a_is2d, bool a_save_interpolated) {
    file_name = a_file_name;
    time_name = a_time_name;
    field_name = a_field_name;
    domain = a_domain;
    mf_var = a_mf_var;
    is2d = a_is2d;
    save_interpolated = a_save_interpolated;
}

void NCTimeSeries::Initialize() {
    // open file
    amrex::Print() << "Loading " << field_name << " from NetCDF file " << file_name << std::endl;

    // The time field can have any number of names, depending on the field.
    // If not specified in input file (time_name.empty()) then set it by default
    if (time_name.empty())
    {
        if (field_name.find("wind") != std::string::npos) {
            time_name = "wind_time";
        } else if ((field_name.find("str") != std::string::npos) and (field_name[0] == 's')) {
            time_name = "sms_time";
        } else if ((field_name.find("str") != std::string::npos) and (field_name[0] == 'b')) {
            time_name = "bms_time";
        } else {
            time_name = "ocean_time";
        }
    }

    amrex::Print() << time_name << std::endl;

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
            ocean_times.push_back((*(array_ts[0].get_data() + nt)) * amrex::Real(60.0) * amrex::Real(60.0) * amrex::Real(24.0));
            // amrex::Print() << "TIMES " << ocean_times[nt] << std::endl;
        }
    }
    int ntimes = ocean_times.size();
    int ioproc = amrex::ParallelDescriptor::IOProcessorNumber();
    amrex::ParallelDescriptor::Bcast(&ntimes,1,ioproc);
    if (!(amrex::ParallelDescriptor::IOProcessor())) {
        ocean_times.resize(ntimes);
    }
    amrex::ParallelDescriptor::Bcast(ocean_times.data(), ocean_times.size(), ioproc);

    // Initialize MultiFabs
    // This assumes that the level 0 distribution map and box array won't
    mf_before = new amrex::MultiFab(mf_var->boxArray(), mf_var->DistributionMap(), 1, mf_var->nGrowVect());
    mf_after = new amrex::MultiFab(mf_var->boxArray(), mf_var->DistributionMap(), 1, mf_var->nGrowVect());
    if (save_interpolated) {
        mf_interpolated = new amrex::MultiFab(mf_var->boxArray(), mf_var->DistributionMap(), 1, mf_var->nGrowVect());
    }

    // dummy initialization
    i_time_before = -100;
}

void NCTimeSeries::update_interpolated_to_time (amrex::Real time) {
    // Figure out time index:
    AMREX_ASSERT(time >= ocean_times[0]);
    AMREX_ASSERT(time <= ocean_times[ocean_times.size()-1]);
    int i_time_before_old = i_time_before;
    for (int nt=0; nt < ocean_times.size()-1; nt++) {
        if ((ocean_times[nt] <= time) and (ocean_times[nt+1] >= time)) {
            i_time_before = nt;
            time_before = ocean_times[nt];
            time_after = ocean_times[nt+1];
            break;
        }
    }

    int i_time_after = i_time_before + 1;
    if (i_time_before_old + 1 == i_time_before) {
        // swap multifabs so we only have to read in one MultiFab
        std::swap(mf_before, mf_after);
        read_in_at_time(mf_after, i_time_after);
    } else if (i_time_before_old != i_time_before) {
        read_in_at_time(mf_after,  i_time_after);
        read_in_at_time(mf_before, i_time_before);
    }

    amrex::Real dt = time_after - time_before;

    auto nodality = mf_var->ixType();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(*mf_var,true); mfi.isValid(); ++mfi) {
        // Adjust box to match ROMS grid
        amrex::Box bx = mfi.growntilebox(amrex::IntVect(1-nodality[0],1-nodality[1],0));

        amrex::Real time_before_copy = time_before;

        // If we're saving the interpolated values, we fill mf_interpolated; otherwise, directly fill
        // the multifab from the REMORA class
        amrex::MultiFab* mf_to_fill = save_interpolated ? mf_interpolated : mf_var;
        amrex::Array4<amrex::Real> to_fill = mf_to_fill->array(mfi);
        amrex::Array4<const amrex::Real> before = mf_before->const_array(mfi);
        amrex::Array4<const amrex::Real> after  = mf_after->const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            to_fill(i,j,k) = before(i,j,k) + (time - time_before_copy) * (after(i,j,k) - before(i,j,k)) / dt;
        });
    }
}

void NCTimeSeries::read_in_at_time (amrex::MultiFab* mf, int itime) {
    // This all assumes that we're on level 0 with only one boxes_at_level
    amrex::FArrayBox NC_fab;
    amrex::Vector<amrex::FArrayBox*> NC_fabs;
    amrex::Vector<std::string> NC_names;
    amrex::Vector<enum NC_Data_Dims_Type> NC_dim_types;

    amrex::Print() << "Reading in " << field_name << " at  time index " << itime << " from " << file_name << std::endl;

    NC_fabs.push_back(&NC_fab) ; NC_names.push_back(field_name);

    if (is2d) {
        NC_dim_types.push_back(NC_Data_Dims_Type::Time_SN_WE);
    } else {
        NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE);
    }

    BuildFABsFromNetCDFFile<amrex::FArrayBox,amrex::Real>(domain, file_name, NC_names, NC_dim_types, NC_fabs, true, itime);

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    {
    // Don't tile this since we are operating on full FABs in this routine
    for ( amrex::MFIter mfi(*mf, false); mfi.isValid(); ++mfi )
    {
        amrex::FArrayBox &fab  = (*mf)[mfi];

        //
        // FArrayBox to FArrayBox copy does "copy on intersection"
        // This only works here because we have broadcast the FArrayBox of data from the netcdf file to all ranks

        fab.template    copy<amrex::RunOn::Device>(NC_fab);
    } // mf
    } // omp
}
#endif // REMORA_USE_NETCDF
