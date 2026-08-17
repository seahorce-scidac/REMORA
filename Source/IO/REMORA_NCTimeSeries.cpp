#include "REMORA_Constants.H"
#include "REMORA_NCTimeSeries.H"
#include "REMORA_NCFile.H"

#include "AMReX_FillPatchUtil.H"
#include "AMReX_Interpolater.H"
#include "AMReX_ParallelDescriptor.H"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#ifdef REMORA_USE_NETCDF
/**
 * @param[in   ] a_file_names         vector of file name(s) to read from
 * @param[in   ] a_field_name         name of field to read in
 * @param[in   ] a_time_name          name of time variable in NetCDF file
 * @param[in   ] a_domain             simulation domain
 * @param[inout] a_mf_var             MultiFab of data to either store into or reference for dimensions
 * @param[in   ] a_is2d               Whether the variable we're working with is 2D
 * @param[in   ] a_save_interpolated  Whether the interpolated value should be saved internally
 */
NCTimeSeries::NCTimeSeries (const amrex::Vector<std::string>& a_file_names, const std::string a_field_name,
                            const std::string a_time_name,
                            const amrex::Box& a_domain,
                            amrex::MultiFab* a_mf_var, bool a_is2d, bool a_save_interpolated) {
    file_names.assign(a_file_names.begin(), a_file_names.end());
    time_name = a_time_name;
    field_name = a_field_name;
    domain = a_domain;
    mf_var = a_mf_var;
    is2d = a_is2d;
    save_interpolated = a_save_interpolated;
}

void NCTimeSeries::Initialize() {
    // open file
    amrex::Print() << "Loading " << field_name << " from NetCDF file(s)" << std::endl;

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

    amrex::Vector<int> file_is_cycle;
    amrex::Vector<amrex::Real> file_cycle_length;
    file_is_spatially_uniform.clear();
    for (int ifile = 0; ifile < file_names.size(); ++ifile) {
        const std::string& file_name = file_names[ifile];

        // Check units of time stamps; should be days
        std::string unit_str = ReadNetCDFVarAttrStr(file_name, time_name, "units"); // works on proc 0
        if (amrex::ParallelDescriptor::IOProcessor())
        {
            if (unit_str.find("days") == std::string::npos) {
                amrex::Print() << "Units of ocean_time given as: " << unit_str << std::endl;
                amrex::Abort("Units must be in days.");
            }
        }

        // `cycle_length` is normally numeric, so use the NCVar utility directly.
        // By ROMS convention it is in the same units as the time variable.
        // REMORA currently requires time units in days, so convert to seconds below.
        bool l_is_cycle = false;
        amrex::Real l_cycle_length = 0.0;
        int l_is_spatially_uniform = 0;

        auto ncf = ncutils::NCFile::open(file_name, NC_NOCLOBBER);
        ncmpi_begin_indep_data(ncf.ncid);

        if (amrex::ParallelDescriptor::IOProcessor())
        {
            auto time_var = ncf.var(time_name);

            l_is_cycle = time_var.has_attr("cycle_length");

            if (l_is_cycle) {
                std::vector<double> cycle_attr;
                time_var.get_attr("cycle_length", cycle_attr);

                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                    cycle_attr.size() == 1,
                    "NetCDF time variable cycle_length attribute must be scalar");

                l_cycle_length = static_cast<amrex::Real>(cycle_attr[0])
                             * amrex::Real(60.0)
                             * amrex::Real(60.0)
                             * amrex::Real(24.0);
            }

            if (!ncf.has_var(field_name)) {
                amrex::Abort("NetCDF time series variable " + field_name +
                             " not found in " + file_name);
            }

            const int field_rank = ncf.var(field_name).ndim();
            const int expected_spatial_rank = is2d ? 3 : 4;
            l_is_spatially_uniform = (field_rank == 1) ? 1 : 0;

            if (field_rank != 1 && field_rank != expected_spatial_rank) {
                amrex::Abort("Unsupported NetCDF time series rank for variable " +
                             field_name + " in " + file_name + ": rank " +
                             std::to_string(field_rank) + ". Expected rank 1 for " +
                             "uniform time-only data or rank " +
                             std::to_string(expected_spatial_rank) +
                             " for spatial data.");
            }
        }

        ncf.close();

        const int ioproc = amrex::ParallelDescriptor::IOProcessorNumber();
        int is_cycle_int = l_is_cycle ? 1 : 0;
        amrex::ParallelDescriptor::Bcast(&is_cycle_int, 1, ioproc);
        amrex::ParallelDescriptor::Bcast(&l_cycle_length, 1, ioproc);
        amrex::ParallelDescriptor::Bcast(&l_is_spatially_uniform, 1, ioproc);

        file_is_cycle.push_back(is_cycle_int);
        file_cycle_length.push_back(l_cycle_length);
        file_is_spatially_uniform.push_back(l_is_spatially_uniform);

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
                file_for_time.push_back(ifile);
                file_itime_offset.push_back(nt);
            }
        }
    }

    bool all_files_cycle_equal = std::all_of(file_is_cycle.begin(), file_is_cycle.end(),
                               [&](const auto& x) { return x == file_is_cycle.front(); });
    bool all_files_cycle_length_equal = std::all_of(file_cycle_length.begin(), file_cycle_length.end(),
                               [&](const auto& x) { return x == file_cycle_length.front(); });

    if (!all_files_cycle_equal || !all_files_cycle_length_equal) {
        amrex::Abort("If one time series file in a set has a cycle, they all must, and cycle lengths must be equal");
    }

    // Store values to class members
    is_cycle = file_is_cycle[0];
    cycle_length = file_cycle_length[0];

    // Arrays will be padded if time series file gives a cycle
    int ntimes = (is_cycle) ? ocean_times.size() + 2 : ocean_times.size();
    // Only do checks on IO processor since ocean_times isn't populated on other ranks yet
    if (amrex::ParallelDescriptor::IOProcessor()) {
        AMREX_ASSERT(std::is_sorted(ocean_times.begin(), ocean_times.end()));
        if (ntimes <= 1) {
            amrex::Error("Time series data must be given at at least two times");
        }
    }
    int ioproc = amrex::ParallelDescriptor::IOProcessorNumber();
    amrex::ParallelDescriptor::Bcast(&ntimes,1,ioproc);
    if (!(amrex::ParallelDescriptor::IOProcessor())) {
        ocean_times.resize(ntimes);
        file_for_time.resize(ntimes);
        file_itime_offset.resize(ntimes);
    } else {
        // If we're in a cycle, hack the lists to close the loop on the IOProc rank
        if (is_cycle) {
            ocean_times.insert(ocean_times.begin(),ocean_times[ocean_times.size()-1]-cycle_length);
            // "first time" is now at index 1 because we already added to the front.
            ocean_times.insert(ocean_times.end(), ocean_times[1]+cycle_length);
            file_for_time.insert(file_for_time.begin(),file_for_time[file_for_time.size()-1]);
            file_for_time.insert(file_for_time.end(), file_for_time[1]);
            file_itime_offset.insert(file_itime_offset.begin(),file_itime_offset[file_itime_offset.size()-1]);
            file_itime_offset.insert(file_itime_offset.end(), file_itime_offset[1]);
        }
    }
    amrex::ParallelDescriptor::Bcast(ocean_times.data(), ocean_times.size(), ioproc);
    amrex::ParallelDescriptor::Bcast(file_for_time.data(), file_for_time.size(), ioproc);
    amrex::ParallelDescriptor::Bcast(file_itime_offset.data(), file_itime_offset.size(), ioproc);

    // Initialize MultiFabs
    // NetCDF data is always read and temporally interpolated on level 0.
    mf_before = new amrex::MultiFab(mf_var->boxArray(), mf_var->DistributionMap(), 1, mf_var->nGrowVect());
    mf_after = new amrex::MultiFab(mf_var->boxArray(), mf_var->DistributionMap(), 1, mf_var->nGrowVect());
    mf_interp_lev0 = new amrex::MultiFab(mf_var->boxArray(), mf_var->DistributionMap(), 1, mf_var->nGrowVect());

    // dummy initialization
    i_time_before = -100;
}

/**
 * @param time   time to interpolate to
 */
void NCTimeSeries::update_interpolated_to_time (amrex::Real time, int lev,
                                                amrex::MultiFab* mf_lev,
                                                const amrex::Vector<amrex::Geometry>& geom,
                                                const amrex::Vector<amrex::IntVect>& ref_ratio) {

    amrex::Real l_time = is_cycle ? std::fmod(time,cycle_length) : time;
    AMREX_ASSERT(l_time >= ocean_times[0]);
    AMREX_ASSERT(l_time <= ocean_times[ocean_times.size()-1]);
    // Figure out time index:
    int i_time_before_old = i_time_before;
    for (int nt=0; nt < ocean_times.size()-1; nt++) {
        if ((ocean_times[nt] <= l_time) and (ocean_times[nt+1] >= l_time)) {
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

    auto nodality = mf_interp_lev0->ixType();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(*mf_interp_lev0,true); mfi.isValid(); ++mfi) {
        // Adjust box to match ROMS grid
        amrex::Box bx = mfi.growntilebox(amrex::IntVect(1-nodality[0],1-nodality[1],0));

        amrex::Real time_before_copy = time_before;

        // Temporal interpolation is done once on level 0.
        amrex::MultiFab* mf_to_fill = mf_interp_lev0;
        amrex::Array4<amrex::Real> to_fill = mf_to_fill->array(mfi);
        amrex::Array4<const amrex::Real> before = mf_before->const_array(mfi);
        amrex::Array4<const amrex::Real> after  = mf_after->const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            to_fill(i,j,k) = before(i,j,k) + (l_time - time_before_copy) * (after(i,j,k) - before(i,j,k)) / dt;
        });
    }

    amrex::MultiFab* mf_to_fill_lev = mf_lev;
    if (save_interpolated) {
        if (mf_interpolated_lev.size() <= static_cast<amrex::Long>(lev)) {
            mf_interpolated_lev.resize(lev+1);
        }
        if (!mf_interpolated_lev[lev] ||
            mf_interpolated_lev[lev]->boxArray() != mf_lev->boxArray() ||
            mf_interpolated_lev[lev]->nGrowVect() != mf_lev->nGrowVect()) {
            mf_interpolated_lev[lev] = std::make_unique<amrex::MultiFab>(
                mf_lev->boxArray(), mf_lev->DistributionMap(), 1, mf_lev->nGrowVect());
        }
        mf_to_fill_lev = mf_interpolated_lev[lev].get();
    }

    if (lev == 0) {
        amrex::MultiFab::Copy(*mf_to_fill_lev, *mf_interp_lev0, 0, 0, 1, mf_to_fill_lev->nGrowVect());
        return;
    }

    amrex::PhysBCFunctNoOp null_bc_for_fill;
    amrex::Vector<amrex::BCRec> null_dom_bcs(1);
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        null_dom_bcs[0].setLo(dir, amrex::BCType::int_dir);
        null_dom_bcs[0].setHi(dir, amrex::BCType::int_dir);
    }

    amrex::Interpolater* mapper = nullptr;
    const auto& idx_type = mf_to_fill_lev->ixType();
    if (idx_type == amrex::IndexType(amrex::IntVect(0,0,0))) {
        mapper = &amrex::cell_cons_interp;
    } else {
        mapper = &amrex::face_cons_linear_interp;
    }

    amrex::InterpFromCoarseLevel(*mf_to_fill_lev, zero, *mf_interp_lev0,
                                 0, 0, 1,
                                 geom[0], geom[lev],
                                 null_bc_for_fill, 0, null_bc_for_fill, 0,
                                 cumulative_ref_ratio(lev, ref_ratio),
                                 mapper, null_dom_bcs, 0);

    mf_to_fill_lev->FillBoundary(geom[lev].periodicity());
}

amrex::IntVect
NCTimeSeries::cumulative_ref_ratio (int lev,
                                    const amrex::Vector<amrex::IntVect>& ref_ratio) const {
    amrex::IntVect rr(1,1,1);
    for (int l = 0; l < lev; ++l) {
        rr[0] *= ref_ratio[l][0];
        rr[1] *= ref_ratio[l][1];
        rr[2] *= ref_ratio[l][2];
    }
    return rr;
}

const amrex::MultiFab*
NCTimeSeries::get_interpolated_mf (int lev) const {
    AMREX_ALWAYS_ASSERT(save_interpolated);
    AMREX_ALWAYS_ASSERT(lev < static_cast<int>(mf_interpolated_lev.size()));
    AMREX_ALWAYS_ASSERT(mf_interpolated_lev[lev]);
    return mf_interpolated_lev[lev].get();
}

/**
 * @param[inout] mf        multifab to store time step data into
 * @param[in   ] itime     index of time step to read from file
 */
void NCTimeSeries::read_in_at_time (amrex::MultiFab* mf, int itime) {
    const int ifile = file_for_time[itime];
    const std::string& file_name = file_names[ifile];
    const int itime_offset = file_itime_offset[itime];

    amrex::Print() << "Reading in " << field_name << " at time index " << itime
                   << " from " << file_name << std::endl;

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        ifile < static_cast<int>(file_is_spatially_uniform.size()),
        "NetCDF time series file layout was not initialized");

    if (file_is_spatially_uniform[ifile]) {
        using RARRAY = NDArray<amrex::Real>;
        amrex::Vector<RARRAY> array_dat(1);
        ReadNetCDFFile(file_name, {field_name}, array_dat, true, itime_offset);

        amrex::Real uniform_value = 0.0;
        if (amrex::ParallelDescriptor::IOProcessor()) {
            const std::vector<MPI_Offset> shape = array_dat[0].get_vshape();
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                shape.size() == 1 && shape[0] == 1,
                "Uniform NetCDF time series read did not produce one value");
            uniform_value = *(array_dat[0].get_data());
        }

        const int ioproc = amrex::ParallelDescriptor::IOProcessorNumber();
        amrex::ParallelDescriptor::Bcast(&uniform_value, 1, ioproc);
        mf->setVal(uniform_value);
        return;
    }

    // This all assumes that we're on level 0 with only one boxes_at_level
    amrex::FArrayBox NC_fab;
    amrex::Vector<amrex::FArrayBox*> NC_fabs;
    amrex::Vector<std::string> NC_names;
    amrex::Vector<enum NC_Data_Dims_Type> NC_dim_types;

    NC_fabs.push_back(&NC_fab);
    NC_names.push_back(field_name);

    if (is2d) {
        NC_dim_types.push_back(NC_Data_Dims_Type::Time_SN_WE);
    } else {
        NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE);
    }

    BuildFABsFromNetCDFFile<amrex::FArrayBox,amrex::Real>(domain, file_name, NC_names, NC_dim_types,
                                                          NC_fabs, true, itime_offset);

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

        fab.template copy<amrex::RunOn::Device>(NC_fab);
    } // mf
    } // omp
}
#endif // REMORA_USE_NETCDF
