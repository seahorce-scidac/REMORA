#include "REMORA_NCTimeSeriesBoundary.H"
#include "REMORA_NCFile.H"
#include "REMORA_BoundaryInterpMath.H"

#include "AMReX_ParallelDescriptor.H"

#include <string>

#ifdef REMORA_USE_NETCDF
/**
 * @param[in   ] a_lev                level at which we will store the data
 * @param[in   ] a_geom               Vector of Geometry objects at all levels
 * @param[in   ] a_file_name          vector of file name(s) to read from
 * @param[in   ] a_field_name         name of field to read in
 * @param[in   ] a_time_name          name of time variable in NetCDF file
 * @param[in   ] a_index_type         nodality of data field
 * @param[in   ] a_var_need_data      array over boundaries of flags that indicate whether we need data for the variable
 * @param[in   ] a_is2d               Whether the variable we're working with is 2D
 * @param[in   ] a_rx                 refinement ratio in x relative to level 0
 * @param[in   ] a_ry                 refinement ratio in y relative to level 0
 * @param[in   ] a_time_ref           remora.time_ref, the reference date of the model clock
 * @param[in   ] a_verbose            REMORA's verbosity; 1 or more prints each record read
 */
NCTimeSeriesBoundary::NCTimeSeriesBoundary (int a_lev, const amrex::Vector<amrex::Geometry> a_geom,
                                            const amrex::Vector<std::string>& a_file_names, const std::string a_field_name,
                                            const std::string a_time_name,
                                            const amrex::IntVect a_index_type,
                                            const amrex::GpuArray<bool, AMREX_SPACEDIM*2>& a_var_need_data,
                                            bool a_is2d, int a_rx, int a_ry, double a_time_ref,
                                            int a_verbose)
{
    time_ref = a_time_ref;
    m_lev = a_lev;
    m_geom = a_geom;
    file_names.assign(a_file_names.begin(), a_file_names.end());
    time_name = a_time_name;
    field_name = a_field_name;
    domain = a_geom[a_lev].Domain();
    index_type = a_index_type;
    // Copied, not aliased: the caller's array (REMORA::phys_bc_need_data) is an
    // amrex::Vector element, whose address is not stable across a resize.
    var_need_data = a_var_need_data;
    is2d = a_is2d;
    m_rx = a_rx;
    m_ry = a_ry;
    m_verbose = a_verbose;
}

void NCTimeSeriesBoundary::Initialize()
{
    // Initialize Fabs
    amrex::Arena* Arena_Used = amrex::The_Arena();
#ifdef AMREX_USE_GPU
    Arena_Used = amrex::The_Pinned_Arena();
#endif

    // open file
    // Only announce the variables that will really be read. A series is constructed for
    // every boundary variable, but most of them typically sit on a local boundary
    // condition and never touch the file.
    bool any_side_needs_data = false;
    for (int ori = 0; ori < AMREX_SPACEDIM*2; ++ori) {
        any_side_needs_data = any_side_needs_data || var_need_data[ori];
    }
    if (any_side_needs_data && m_verbose >= 1) {
        amrex::Print() << "Setting up NetCDF boundary data for " << field_name
                       << std::endl;
    }

    // The time field can have any number of names, depending on the field.
    // If not specified in input file (time_name.empty()) then set it by default
    if (time_name.empty())
    {
        time_name = "ocean_time";
    }

    int ioproc = amrex::ParallelDescriptor::IOProcessorNumber();
    for (int ifile = 0; ifile < file_names.size(); ifile++) {
        std::string file_name = file_names[ifile];
        // The time stamps must count days, if they carry units at all; without
        // units they are assumed to be days on the model clock. A reference date
        // in the units is honored: this is the shift that puts them on the model
        // clock. Only announced for a series that will be read.
        const double offset_days = ReadNetCDFTimeOffsetDays(file_name, time_name, time_ref,
                                                            false, any_side_needs_data);
        // get times and put in array
        amrex::Vector<NDArray<double>> array_ts(1);
        ReadNetCDFFile(file_name, {time_name}, array_ts); // filled only on proc 0
        if (amrex::ParallelDescriptor::IOProcessor())
        {
            int ntimes_io = array_ts[0].get_vshape()[0];
            for (int nt(0); nt < ntimes_io; nt++)
            {
                // Onto the model clock, then from days to seconds
                bry_times.push_back((*(array_ts[0].get_data() + nt) + offset_days) * 60.0 * 60.0 * 24.0);
                file_for_time.push_back(ifile);
                file_itime_offset.push_back(nt);
            }
        }
    }

    int ntimes = bry_times.size();
    // Only do checks on IO processor since bry_times isn't populated on other ranks yet
    if (amrex::ParallelDescriptor::IOProcessor()) {
        AMREX_ASSERT(std::is_sorted(bry_times.begin(), bry_times.end()));
        if (ntimes <= 1) {
            amrex::Error("Time series of boundary data must be given at at least two times");
        }
    }
    amrex::ParallelDescriptor::Bcast(&ntimes,1,ioproc);
    if (!(amrex::ParallelDescriptor::IOProcessor())) {
        bry_times.resize(ntimes);
        file_for_time.resize(ntimes);
        file_itime_offset.resize(ntimes);
    }
    amrex::ParallelDescriptor::Bcast(bry_times.data(), bry_times.size(), ioproc);
    amrex::ParallelDescriptor::Bcast(file_for_time.data(), file_for_time.size(), ioproc);
    amrex::ParallelDescriptor::Bcast(file_itime_offset.data(), file_itime_offset.size(), ioproc);

    using amrex::Orientation;
    using amrex::Direction;
    for (const Orientation ori : {Orientation(Direction::x, Orientation::low),
                                  Orientation(Direction::x, Orientation::high),
                                  Orientation(Direction::y, Orientation::low),
                                  Orientation(Direction::y, Orientation::high)}) {
        // fill_from_bdyfiles takes the interpolated strip of every side, so each gets one; on a
        // side that takes no boundary data it holds 0. Only the sides read hold the two records.
        // The arena is host-accessible (pinned on GPU builds).
        m_dat_interp[ori] = amrex::FArrayBox(strip_box(domain, ori), 1, Arena_Used);
        m_dat_interp[ori].setVal<amrex::RunOn::Host>(amrex::Real(0.0));
        if (var_need_data[ori]) {
            m_sides.push_back(ori);
            m_dat_before[ori] = amrex::FArrayBox(strip_box(domain, ori), 1, Arena_Used);
            m_dat_after[ori]  = amrex::FArrayBox(strip_box(domain, ori), 1, Arena_Used);
        }
    }

    // dummy initialization
    i_time_before = -100;

    // Nothing to check for a variable that reads no boundary data at all, and skipping
    // keeps a run with many tracers from opening the file once per inert tracer
    if (m_sides.empty()) {
        return;
    }

    // Fail early and by name. read_in_at_time dereferences each variable unconditionally,
    // so a file that is missing e.g. NO3_west would otherwise fail deep inside the reader
    // with nothing to say which variable, which file, or which input asked for it.
    const amrex::Vector<std::string> names = var_names();
    for (const auto& file_name : file_names) {
        if (!QueryNetCDFHasVars(file_name, names)) {
            std::string msg = "Boundary file " + file_name + " does not contain all of:";
            for (const auto& var_name : names) {
                msg += " " + var_name;
            }
            msg += ". The boundary condition requested for " + field_name +
                   " reads from file (clamped, chapman, flather, or orlanski_rad_nudg)."
                   " Either add the missing variable(s) to the file, or give " + field_name +
                   " a boundary condition that needs no file data, such as outflow.";
            amrex::Abort(msg);
        }
    }
}

/**
 * @param time   time to interpolate to
 */
void NCTimeSeriesBoundary::update_interpolated_to_time (double time)
{
    BL_PROFILE("NCTimeSeriesBoundary::update_interpolated_to_time()");
    // Nothing in the file for this variable, so there is nothing to read or
    // interpolate. Returning here also keeps the reader from announcing that it is
    // reading a variable the file does not contain: the caller updates every tracer
    // in one sweep, and most of them are typically on a local boundary condition.
    if (m_sides.empty()) {
        return;
    }
    // The strips are already at this time; callers only read them
    if (have_interp && interp_time == time) {
        return;
    }

    // Initialize Fabs
    amrex::Arena* Arena_Used = amrex::The_Arena();
#ifdef AMREX_USE_GPU
    Arena_Used = amrex::The_Pinned_Arena();
#endif

    // Figure out time index:
    int i_time_before_old = i_time_before;
    int i_time_new = -1;
    for (int nt=0; nt < bry_times.size()-1; nt++) {
        if ((bry_times[nt] <= time) and (bry_times[nt+1] >= time)) {
            i_time_new = nt;
            break;
        }
    }
    // Bracketing must succeed: falling through would read the series at a stale or
    // uninitialized index. A runtime check rather than an assert, because the failure
    // is a mismatch between the file and the run.
    if (i_time_new < 0) {
        amrex::Abort("Time " + std::to_string(time) + " is not spanned by the boundary time"
                     + " series for '" + field_name + "', which covers ["
                     + std::to_string(bry_times[0]) + ", "
                     + std::to_string(bry_times[bry_times.size()-1])
                     + "]. Extend the boundary file to cover the run. The file's times were"
                     + " placed on the model clock from the reference date in their units"
                     + " attribute (the mapping is printed when the series is loaded), so"
                     + " also check remora.time_ref and remora.start_time.");
    }
    i_time_before = i_time_new;
    time_before = bry_times[i_time_before];
    time_after = bry_times[i_time_before+1];

    int i_time_after = i_time_before + 1;

    // Level 0 reads the file onto its own strips. A refined level reads it into level-0 strips,
    // made on the first read of this call and only for the sides that take data, and
    // interpolates them onto its own.
    std::array<amrex::FArrayBox, AMREX_SPACEDIM*2> crse;
    auto read_record = [&] (std::array<amrex::FArrayBox, AMREX_SPACEDIM*2>& dat, int itime)
    {
        if (m_lev == 0) {
            read_in_at_time(dat, itime);
            return;
        }
        for (const auto& ori : m_sides) {
            if (!crse[ori].isAllocated()) {
                crse[ori] = amrex::FArrayBox(strip_box(m_geom[0].Domain(), ori), 1, Arena_Used);
            }
        }
        read_in_at_time(crse, itime);
        for (const auto& ori : m_sides) {
            interp_fab(crse[ori], dat[ori], ori);
        }
    };

    bool read_any = true;
    if (i_time_before_old + 1 == i_time_before) {
        // The old after record becomes the before record, so only one record is read
        for (const auto& ori : m_sides) {
            std::swap(m_dat_before[ori], m_dat_after[ori]);
        }
        read_record(m_dat_after, i_time_after);
    } else if (i_time_before_old != i_time_before) {
        read_record(m_dat_before, i_time_before);
        read_record(m_dat_after, i_time_after);
    } else {
        read_any = false;
    }
    // The_Pinned_Arena is not stream-ordered, so wait for interp_fab's kernels before crse is
    // freed at return
    if (read_any && m_lev > 0) {
        amrex::Gpu::streamSynchronize();
    }

    // Both differences are small, so Real holds them even when the times do not.
    const amrex::Real dt = static_cast<amrex::Real>(time_after - time_before);
    const amrex::Real time_since_before = static_cast<amrex::Real>(time - time_before);

    for (const auto& ori : m_sides) {
        amrex::Array4<amrex::Real> interp_arr = m_dat_interp[ori].array();
        amrex::Array4<amrex::Real const> before_arr = m_dat_before[ori].const_array();
        amrex::Array4<amrex::Real const> after_arr = m_dat_after[ori].const_array();
        amrex::ParallelFor(m_dat_interp[ori].box(), [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            interp_arr(i,j,k) = before_arr(i,j,k)
                              + time_since_before * (after_arr(i,j,k) - before_arr(i,j,k)) / dt;
        });
    }
    have_interp = true;
    interp_time = time;
}

/**
 * @param[in   ] dom   domain whose side it is
 * @param[in   ] ori   horizontal side
 */
amrex::Box
NCTimeSeriesBoundary::strip_box (const amrex::Box& dom, amrex::Orientation ori) const
{
    // The domain and its boundary ring at this index type, then one point thick across ori
    const int d = ori.coordDir();
    amrex::IntVect lo(dom.smallEnd(0) + index_type[0] - 1, dom.smallEnd(1) + index_type[1] - 1,
                      dom.smallEnd(2));
    amrex::IntVect hi(dom.bigEnd(0) + 1, dom.bigEnd(1) + 1, dom.bigEnd(2));
    if (ori.isLow()) {
        hi[d] = lo[d];
    } else {
        lo[d] = hi[d];
    }
    amrex::Box bx(lo, hi, index_type);
    if (is2d) {
        bx.makeSlab(2,0);
    }
    return bx;
}

amrex::Vector<std::string>
NCTimeSeriesBoundary::var_names () const
{
    // ROMS names a side's boundary variable by its compass direction
    amrex::Vector<std::string> names;
    for (const auto& ori : m_sides) {
        const char* side = (ori.coordDir() == 0) ? (ori.isLow() ? "_west" : "_east")
                                                 : (ori.isLow() ? "_south" : "_north");
        names.push_back(field_name + side);
    }
    return names;
}

/**
 * Read the record at time index itime onto the strips of the sides that take boundary data.
 *
 * @param[inout] fabs    strips, indexed by Orientation; only those of m_sides are filled
 * @param[in   ] itime   index of time step to read from file
 */
void
NCTimeSeriesBoundary::read_in_at_time (std::array<amrex::FArrayBox, AMREX_SPACEDIM*2>& fabs,
                                       int itime)
{
    BL_PROFILE("NCTimeSeriesBoundary::read_in_at_time()");
    // The host writes below land in pinned fabs that queued kernels may still be
    // reading in place, so drain the stream first
    amrex::Gpu::streamSynchronize();

    const amrex::Vector<std::string> names = var_names();
    using RARRAY = NDArray<amrex::Real>;
    amrex::Vector<RARRAY> arrays(names.size());

    if (m_verbose >= 1) {
        amrex::Print() << "Reading boundary record " << itime << " of " << field_name << " (time "
                       << bry_times[itime] << " s)" << std::endl;
    }
    std::string nc_bdry_file = file_names[file_for_time[itime]];
    int itime_offset = file_itime_offset[itime];
    ReadNetCDFFile(nc_bdry_file, names, arrays, true, itime_offset); // does work on proc 0 only

    if (amrex::ParallelDescriptor::IOProcessor()) {
        for (int iv = 0; iv < names.size(); iv++) {
            const amrex::Orientation ori = m_sides[iv];
            amrex::FArrayBox& fab = fabs[ori];
            // The file holds ([s,] along) for this record, with the along-boundary index fastest
            const int along = 1 - ori.coordDir();
            const auto& shape = arrays[iv].get_vshape();
            const int n_along = static_cast<int>(is2d ? shape[1] : shape[2]);
            const int nz = is2d ? 1 : static_cast<int>(shape[1]);
            const int n_plane = n_along * nz;
            AMREX_ALWAYS_ASSERT(fab.box().numPts() == n_plane);

            const amrex::Dim3 lo = amrex::lbound(fab.box());
            amrex::Array4<amrex::Real> fab_arr = fab.array();
            for (int n = 0; n < n_plane; n++) {
                const int k = n / n_along;
                const int a = n - k * n_along;
                const int i = (along == 0) ? lo.x + a : lo.x;
                const int j = (along == 0) ? lo.y : lo.y + a;
                fab_arr(i, j, k, 0) = static_cast<amrex::Real>(*(arrays[iv].get_data() + n));
            }
        }
    }

    amrex::ParallelDescriptor::Barrier();

    // The strips were filled on the I/O rank only; send them to every rank
    int ioproc = amrex::ParallelDescriptor::IOProcessorNumber();
    for (const auto& ori : m_sides) {
        amrex::ParallelDescriptor::Bcast(fabs[ori].dataPtr(), fabs[ori].box().numPts(), ioproc);
    }

}

/**
 * \brief Interpolate a level-0 boundary strip onto this level's strip.
 *
 * Values are interpolated linearly along the strip between the coarse points on either side
 * of each fine point: cell centers where the strip is cell-centered along its length, nodes
 * where it is nodal. The point across the boundary is not interpolated.
 *
 * @param[in   ] dat_crse  fab of coarse data to interpolate from
 * @param[  out] dat_fine  fab of fine data to interpolate to
 * @param[in   ] ori       the side the strip lies on
 */
void
NCTimeSeriesBoundary::interp_fab (const amrex::FArrayBox& dat_crse, amrex::FArrayBox& dat_fine,
                                  amrex::Orientation ori)
{
    amrex::Array4<amrex::Real const> crse_arr = dat_crse.const_array();
    amrex::Array4<amrex::Real> fine_arr = dat_fine.array();

    // The strip is one point thick across the boundary and runs along the other direction
    AMREX_ASSERT(dat_crse.box().length(ori.coordDir()) == 1 &&
                 dat_fine.box().length(ori.coordDir()) == 1);
    const int along = 1 - ori.coordDir();
    const int across_crse = dat_crse.box().smallEnd(1 - along);
    const int ratio = (along == 0) ? m_rx : m_ry;
    const bool nodal = dat_crse.box().ixType()[along] == 1;

    amrex::ParallelFor(dat_fine.box(), [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        int c;
        amrex::Real w;
        const int f = (along == 0) ? i : j;
        if (nodal) {
            boundary_interp::nodal_weight(f, ratio, c, w);
        } else {
            boundary_interp::cell_centered_weight(f, ratio, c, w);
        }
        const int ic = (along == 0) ? c : across_crse;
        const int jc = (along == 0) ? across_crse : c;
        const amrex::Real lo = crse_arr(ic,jc,k);
        if (w == amrex::Real(0.0)) {
            fine_arr(i,j,k) = lo;
        } else {
            const amrex::Real hi = (along == 0) ? crse_arr(ic+1,jc,k) : crse_arr(ic,jc+1,k);
            fine_arr(i,j,k) = lo + w * (hi - lo);
        }
    });
}

#endif // REMORA_USE_NETCDF
