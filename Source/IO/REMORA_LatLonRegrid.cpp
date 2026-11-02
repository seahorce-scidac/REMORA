#include "REMORA_LatLonRegrid.H"
#include "REMORA_NCFile.H"

#include <AMReX_GpuContainers.H>
#include <AMReX_Utility.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParallelReduce.H>

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>

#ifdef REMORA_USE_NETCDF

using namespace amrex;

namespace {

/** Broadcast a string from the IO rank */
void bcast_string (std::string& s)
{
    amrex::BroadcastString(s, ParallelDescriptor::MyProc(), ParallelDescriptor::IOProcessorNumber(),
                           ParallelDescriptor::Communicator());
}

/** Broadcast a vector of doubles from the IO rank */
void bcast_vector (Vector<double>& v)
{
    amrex::BroadcastArray(v, ParallelDescriptor::MyProc(), ParallelDescriptor::IOProcessorNumber(),
                          ParallelDescriptor::Communicator());
}

/** Name of dimension idim of a variable */
std::string dim_name (const ncutils::NCFile& ncf, const ncutils::NCVar& var, int idim)
{
    std::vector<int> dimids(var.ndim());
    ncmpi_inq_vardimid(ncf.ncid, var.varid, dimids.data());
    char name[NC_MAX_NAME + 1];
    ncmpi_inq_dimname(ncf.ncid, dimids[idim], name);
    return std::string(name);
}

/** Whether a units attribute names degrees east (or north) */
bool units_say (const ncutils::NCVar& var, const std::string& direction)
{
    if (!var.has_attr("units")) { return false; }
    std::string u = var.get_attr("units");
    std::transform(u.begin(), u.end(), u.begin(), [] (unsigned char c) { return std::tolower(c); });
    return u.find("degree") != std::string::npos &&
           (u.find(direction) != std::string::npos ||
            u.find(direction == "east" ? "_e" : "_n") != std::string::npos);
}

/** The source grid of a variable: the names and values of its coordinate axes, its
 *  missing-data markers and its packing */
struct SourceGrid
{
    std::string lon_name, lat_name;
    Vector<double> lon, lat;
    int has_fill = 0, has_missing = 0;
    double fill_value = 0.0, missing_value = 0.0;
    double scale_factor = 1.0, add_offset = 0.0;
};

/** Whether dimension idim of a variable has a CF coordinate variable: a 1D variable of the
 *  dimension's own name */
bool has_dimension_coordinate (const ncutils::NCFile& ncf, const ncutils::NCVar& var, int idim)
{
    const std::string d = dim_name(ncf, var, idim);
    return ncf.has_var(d) && ncf.var(d).ndim() == 1 && dim_name(ncf, ncf.var(d), 0) == d;
}

/**
 * Find and read the coordinate axes of a (time, lat, lon) variable. The names come from
 * the inputs when given; otherwise from the variable's CF coordinates attribute, each name
 * identified by its units (degrees east or north) or, failing that, by its leading
 * "lon"/"lat" or, failing that, by order as ROMS's get_varcoords takes them; otherwise from
 * the CF coordinate variables of its last two dimensions (longitude(longitude), as ERA5 and
 * most reanalyses write them). Each axis must be 1D and dimensioned like the variable's
 * matching dimension. Collective; aborts on every rank with the reason, followed by
 * `context`, when the axes cannot be found.
 */
SourceGrid
find_source_grid (const std::string& file_name, const std::string& field_name,
                  const std::vector<long long>& field_shape,
                  const LatLonRegridOptions& options, const std::string& context)
{
    SourceGrid g;
    std::string error;

    auto ncf = ncutils::NCFile::open(file_name, NC_NOCLOBBER);
    ncmpi_begin_indep_data(ncf.ncid);
    if (ParallelDescriptor::IOProcessor())
    {
        const auto var = ncf.var(field_name);
        const int rank = var.ndim();

        g.lon_name = options.lon_name;
        g.lat_name = options.lat_name;
        if (g.lon_name.empty() || g.lat_name.empty()) {
            if (!var.has_attr("coordinates")) {
                if (has_dimension_coordinate(ncf, var, rank-1) &&
                    has_dimension_coordinate(ncf, var, rank-2)) {
                    if (g.lon_name.empty()) { g.lon_name = dim_name(ncf, var, rank-1); }
                    if (g.lat_name.empty()) { g.lat_name = dim_name(ncf, var, rank-2); }
                } else {
                    error = "it has no coordinates attribute naming its longitude and latitude, "
                            "its last two dimensions have no coordinate variables, and "
                            "remora.frc_lon_varname and remora.frc_lat_varname are not both set";
                }
            } else {
                std::istringstream tokens(var.get_attr("coordinates"));
                std::vector<std::string> names;
                for (std::string t; tokens >> t; ) {
                    if (ncf.has_var(t) && ncf.var(t).ndim() >= 1 &&
                        dim_name(ncf, ncf.var(t), 0) != dim_name(ncf, var, 0)) {
                        names.push_back(t);    // skip the time coordinate
                    }
                }
                std::string by_units_lon, by_units_lat, by_name_lon, by_name_lat;
                for (const auto& t : names) {
                    if (units_say(ncf.var(t), "east"))  { by_units_lon = t; }
                    if (units_say(ncf.var(t), "north")) { by_units_lat = t; }
                    if (t.rfind("lon", 0) == 0) { by_name_lon = t; }
                    if (t.rfind("lat", 0) == 0) { by_name_lat = t; }
                }
                std::string lon = !by_units_lon.empty() ? by_units_lon : by_name_lon;
                std::string lat = !by_units_lat.empty() ? by_units_lat : by_name_lat;
                if ((lon.empty() || lat.empty()) && names.size() == 2) {
                    lon = names[0];
                    lat = names[1];
                }
                if (g.lon_name.empty()) { g.lon_name = lon; }
                if (g.lat_name.empty()) { g.lat_name = lat; }
                if (g.lon_name.empty() || g.lat_name.empty()) {
                    error = "its coordinates attribute '" + var.get_attr("coordinates") +
                            "' does not name a longitude and a latitude variable; set "
                            "remora.frc_lon_varname and remora.frc_lat_varname";
                }
            }
        }

        // Each axis 1D, along the variable's last (lon) or next-to-last (lat) dimension
        const std::pair<std::string, int> axes[2] = {{g.lon_name, rank-1}, {g.lat_name, rank-2}};
        for (const auto& [name, idim] : axes) {
            if (!error.empty()) { break; }
            if (!ncf.has_var(name)) {
                error = "its coordinate variable '" + name + "' is not in the file";
            } else if (ncf.var(name).ndim() != 1) {
                error = "its coordinate variable '" + name + "' is " +
                        std::to_string(ncf.var(name).ndim()) + "D; only 1D (rectilinear) "
                        "longitude and latitude axes are supported";
            } else if (dim_name(ncf, ncf.var(name), 0) != dim_name(ncf, var, idim) ||
                       ncf.var(name).shape()[0] != field_shape[idim]) {
                error = "its coordinate variable '" + name + "' does not run along its " +
                        std::string(idim == rank-1 ? "last" : "next-to-last") + " dimension '" +
                        dim_name(ncf, var, idim) + "'; the variable must be (time, lat, lon)";
            }
        }

        int has_scale = 0, has_offset = 0;
        for (const auto& [attr, has, val] : { std::make_tuple("_FillValue", &g.has_fill, &g.fill_value),
                                              std::make_tuple("missing_value", &g.has_missing, &g.missing_value),
                                              std::make_tuple("scale_factor", &has_scale, &g.scale_factor),
                                              std::make_tuple("add_offset", &has_offset, &g.add_offset) })
        {
            if (var.has_attr(attr)) {
                std::vector<double> v;
                var.get_attr(attr, v);
                if (!v.empty()) { *has = 1; *val = v[0]; }
            }
        }
    }
    ncf.close();

    bcast_string(error);
    if (!error.empty()) {
        Abort("Cannot interpolate '" + field_name + "' in " + file_name +
              " onto the model grid: " + error + "." + context);
    }
    bcast_string(g.lon_name);
    bcast_string(g.lat_name);
    const int ioproc = ParallelDescriptor::IOProcessorNumber();
    ParallelDescriptor::Bcast(&g.has_fill, 1, ioproc);
    ParallelDescriptor::Bcast(&g.has_missing, 1, ioproc);
    ParallelDescriptor::Bcast(&g.fill_value, 1, ioproc);
    ParallelDescriptor::Bcast(&g.missing_value, 1, ioproc);
    ParallelDescriptor::Bcast(&g.scale_factor, 1, ioproc);
    ParallelDescriptor::Bcast(&g.add_offset, 1, ioproc);

    Vector<NDArray<double>> arrays(2);
    ReadNetCDFFile(file_name, {g.lon_name, g.lat_name}, arrays);
    if (ParallelDescriptor::IOProcessor()) {
        g.lon.assign(arrays[0].get_data(), arrays[0].get_data() + arrays[0].get_vshape()[0]);
        g.lat.assign(arrays[1].get_data(), arrays[1].get_data() + arrays[1].get_vshape()[0]);
    }
    bcast_vector(g.lon);
    bcast_vector(g.lat);
    return g;
}

/** Largest m with a[m] <= x in an increasing array, or -1 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int bracket (const Real* a, int n, Real x)
{
    int lo = -1, hi = n;    // a[lo] <= x < a[hi], with a[-1] = -inf and a[n] = +inf
    while (hi - lo > 1) {
        const int mid = (lo + hi) / 2;
        if (a[mid] <= x) { lo = mid; } else { hi = mid; }
    }
    return lo;
}

/** Whether a source value is data: apply() has already turned missing-data markers to NaN */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool is_valid (Real v)
{
    return (v == v) && std::abs(v) <= std::numeric_limits<Real>::max();
}

} // namespace

void
LatLonRegrid::define (const std::string& file_name, const std::string& field_name,
                      const std::vector<long long>& field_shape,
                      const MultiFab& target_lon, const MultiFab& target_lat,
                      const MultiFab& like, const Box& a_domain,
                      const LatLonRegridOptions& options, const std::string& context)
{
    const SourceGrid g = find_source_grid(file_name, field_name, field_shape, options, context);
    lon_name = g.lon_name;
    lat_name = g.lat_name;
    lon_file.assign(g.lon.begin(), g.lon.end());
    lat_file.assign(g.lat.begin(), g.lat.end());
    nlon = static_cast<int>(lon_file.size());
    nlat = static_cast<int>(lat_file.size());
    has_fill = g.has_fill;
    has_missing = g.has_missing;
    fill_value = g.fill_value;
    missing_value = g.missing_value;
    scale_factor = g.scale_factor;
    add_offset = g.add_offset;
    fill_max_rings = options.fill_max_rings;
    domain = makeSlab(a_domain, 2, 0);

    const std::string where = "'" + field_name + "' in " + file_name;
    if (nlon < 2 || nlat < 2) {
        Abort("Cannot interpolate " + where + ": its source grid is " + std::to_string(nlon) +
              " x " + std::to_string(nlat) + "; it needs at least two points along each axis.");
    }

    // Longitude made continuous, so an axis crossing the branch cut (170 .. 180, -180 .. -170)
    // is monotone; then both axes put in increasing order.
    lon_inc = lon_file;
    for (int m = 1; m < nlon; ++m) {
        lon_inc[m] -= 360.0 * std::round((lon_inc[m] - lon_inc[m-1]) / 360.0);
    }
    lat_inc = lat_file;
    for (auto [axis, flip, name] : { std::make_tuple(&lon_inc, &flip_lon, lon_name),
                                     std::make_tuple(&lat_inc, &flip_lat, lat_name) })
    {
        *flip = axis->back() < axis->front();
        if (*flip) { std::reverse(axis->begin(), axis->end()); }
        for (std::size_t m = 1; m < axis->size(); ++m) {
            if (!((*axis)[m] > (*axis)[m-1])) {
                Abort("Cannot interpolate " + where + ": its axis '" + name +
                      "' is not strictly monotone (at index " + std::to_string(m) + ").");
            }
        }
    }

    // A longitude axis whose spacing closes the circle covers the globe
    const double span = lon_inc.back() - lon_inc.front();
    const double dlon = span / (nlon - 1);
    periodic_lon = std::abs(span + dlon - 360.0) < 1.e-3 * dlon;
    if (span > 360.0 + 1.e-9) {
        Abort("Cannot interpolate " + where + ": its longitude axis spans " + std::to_string(span) +
              " degrees, more than once around the globe.");
    }

    weights.define(like.boxArray(), like.DistributionMap(), stencil_size, IntVect(1,1,0));
    indices.define(like.boxArray(), like.DistributionMap(), 3*stencil_size, IntVect(1,1,0));
    weights.setVal(0.0);
    indices.setVal(-1);

    // Target coordinates on the forcing layout, with the ring
    MultiFab lon(like.boxArray(), like.DistributionMap(), 1, IntVect(1,1,0));
    MultiFab lat(like.boxArray(), like.DistributionMap(), 1, IntVect(1,1,0));
    lon.ParallelCopy(target_lon, 0, 0, 1, IntVect(1,1,0), IntVect(1,1,0));
    lat.ParallelCopy(target_lat, 0, 0, 1, IntVect(1,1,0), IntVect(1,1,0));

    build_weights(lon, lat);

    amrex::Print() << "Interpolating " << field_name << " from " << summary() << std::endl;
}

void
LatLonRegrid::build_weights (const MultiFab& lon, const MultiFab& lat)
{
    Gpu::DeviceVector<Real> d_lon(nlon), d_lat(nlat);
    {
        std::vector<Real> h_lon(lon_inc.begin(), lon_inc.end());
        std::vector<Real> h_lat(lat_inc.begin(), lat_inc.end());
        Gpu::copy(Gpu::hostToDevice, h_lon.begin(), h_lon.end(), d_lon.begin());
        Gpu::copy(Gpu::hostToDevice, h_lat.begin(), h_lat.end(), d_lat.begin());
    }
    const Real* plon = d_lon.data();
    const Real* plat = d_lat.data();
    const int Nx = nlon, Ny = nlat;
    const bool periodic = periodic_lon, fx = flip_lon, fy = flip_lat;
    // A target longitude is taken into [lon_ref, lon_ref + 360): from the axis's first point
    // when it covers the globe, and otherwise centered on the axis, so that a point outside
    // it is reported on the side it misses.
    const Real lon_ref = periodic ? static_cast<Real>(lon_inc.front())
                                  : static_cast<Real>(0.5 * (lon_inc.front() + lon_inc.back()) - 180.0);
    const Real tol = Real(1.e-9);

    ReduceOps<ReduceOpSum, ReduceOpMin, ReduceOpMax, ReduceOpMin, ReduceOpMax> reduce_op;
    ReduceData<int, Real, Real, Real, Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

    for (MFIter mfi(weights); mfi.isValid(); ++mfi) {
        const Box bx = target_box(mfi);
        auto const& lon_a = lon.const_array(mfi);
        auto const& lat_a = lat.const_array(mfi);
        auto const& w   = weights.array(mfi);
        auto const& idx = indices.array(mfi);

        reduce_op.eval(bx, reduce_data, [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
        {
            // Into the source's longitude window
            const Real lt = lon_a(i,j,k) - Real(360.0) * std::floor((lon_a(i,j,k) - lon_ref) / Real(360.0));
            const Real pt = lat_a(i,j,k);

            int m, m2;
            Real p;
            bool inside = true;
            if (periodic) {
                m = amrex::max(0, bracket(plon, Nx, lt));
                m2 = (m + 1) % Nx;
                const Real right = (m == Nx-1) ? plon[0] + Real(360.0) : plon[m+1];
                p = (lt - plon[m]) / (right - plon[m]);
            } else {
                inside = (lt >= plon[0] - tol) && (lt <= plon[Nx-1] + tol);
                m = amrex::min(amrex::max(bracket(plon, Nx, lt), 0), Nx-2);
                m2 = m + 1;
                p = (lt - plon[m]) / (plon[m2] - plon[m]);
            }
            inside = inside && (pt >= plat[0] - tol) && (pt <= plat[Ny-1] + tol);
            const int n  = amrex::min(amrex::max(bracket(plat, Ny, pt), 0), Ny-2);
            const int n2 = n + 1;
            Real q = (pt - plat[n]) / (plat[n2] - plat[n]);
            p = amrex::min(amrex::max(p, Real(0.0)), Real(1.0));
            q = amrex::min(amrex::max(q, Real(0.0)), Real(1.0));

            // Stencil in the file's own index order
            auto fi = [=] (int a) { return fx ? Nx-1-a : a; };
            auto fj = [=] (int b) { return fy ? Ny-1-b : b; };
            const int    si[4] = { fi(m), fi(m2), fi(m2), fi(m)  };
            const int    sj[4] = { fj(n), fj(n),  fj(n2), fj(n2) };
            const Real   sw[4] = { (Real(1.0)-p)*(Real(1.0)-q), p*(Real(1.0)-q), p*q, (Real(1.0)-p)*q };
            for (int s = 0; s < 4; ++s) {
                w(i,j,k,s)       = inside ? sw[s] : Real(0.0);
                idx(i,j,k,3*s)   = inside ? si[s] : -1;
                idx(i,j,k,3*s+1) = inside ? sj[s] : -1;
                idx(i,j,k,3*s+2) = 0;
            }

            return { inside ? 0 : 1, lt, lt, pt, pt };
        });
    }

    ReduceTuple r = reduce_data.value(reduce_op);
    int  n_out   = amrex::get<0>(r);
    Real lon_min = amrex::get<1>(r), lon_max = amrex::get<2>(r);
    Real lat_min = amrex::get<3>(r), lat_max = amrex::get<4>(r);
    ParallelDescriptor::ReduceIntSum(n_out);
    ParallelDescriptor::ReduceRealMin(lon_min);
    ParallelDescriptor::ReduceRealMax(lon_max);
    ParallelDescriptor::ReduceRealMin(lat_min);
    ParallelDescriptor::ReduceRealMax(lat_max);

    if (n_out > 0) {
        std::ostringstream os;
        os << "The forcing grid does not contain the model grid: " << n_out << " rho points of"
           << " level 0 and its boundary ring lie outside it. Forcing grid: " << summary()
           << ". Model grid (in the forcing grid's longitude convention): lon " << lon_min
           << " to " << lon_max << ", lat " << lat_min << " to " << lat_max
           << ". Supply forcing that covers the domain with its boundary ring.";
        Abort(os.str());
    }
}

void
LatLonRegrid::apply (const std::string& file_name, const std::string& field_name, int itime,
                     MultiFab& dst) const
{
    Vector<NDArray<Real>> arrays(1);
    ReadNetCDFFile(file_name, {field_name}, arrays, true, itime);

    std::vector<Real> buf(static_cast<std::size_t>(nlon) * nlat);
    if (ParallelDescriptor::IOProcessor()) {
        const auto shape = arrays[0].get_vshape();
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(shape.size() == 3 && shape[1] == nlat && shape[2] == nlon,
            "Lat-lon forcing record of '" + field_name + "' has an unexpected shape");
        // Missing-data markers are compared exactly and in the stored (packed) values, as CF
        // defines them, then the data unpacked. A marker becomes NaN, which the interpolation
        // sets aside.
        const double nan = std::numeric_limits<double>::quiet_NaN();
        const Real* raw = arrays[0].get_data();
        for (std::size_t c = 0; c < buf.size(); ++c) {
            const double v = static_cast<double>(raw[c]);
            const bool marker = (has_fill && v == fill_value) || (has_missing && v == missing_value);
            buf[c] = static_cast<Real>(marker ? nan : v * scale_factor + add_offset);
        }
    }
    ParallelDescriptor::Bcast(buf.data(), buf.size(), ParallelDescriptor::IOProcessorNumber());

    // (lat, lon) row-major in the file is (i = lon, j = lat) in a FAB
    FArrayBox src(Box(IntVect(0,0,0), IntVect(nlon-1, nlat-1, 0)), 1);
    Gpu::copy(Gpu::hostToDevice, buf.begin(), buf.end(), src.dataPtr());

    interpolate(src, dst);

    // A target point with no valid source value within reach is left NaN
    ReduceOps<ReduceOpSum> op;
    ReduceData<int> data(op);
    for (MFIter mfi(dst); mfi.isValid(); ++mfi) {
        const Box bx = target_box(mfi);
        auto const& d = dst.const_array(mfi);
        op.eval(bx, data, [=] AMREX_GPU_DEVICE (int i, int j, int k) -> GpuTuple<int>
        {
            return { (d(i,j,k) == d(i,j,k)) ? 0 : 1 };
        });
    }
    int n_missing = amrex::get<0>(data.value(op));
    ParallelDescriptor::ReduceIntSum(n_missing);
    if (n_missing > 0) {
        Abort("Interpolating '" + field_name + "' record " + std::to_string(itime) + " of " +
              file_name + " onto the model grid: " + std::to_string(n_missing) +
              " rho points have no valid source value within " + std::to_string(fill_max_rings) +
              " cells (the source is _FillValue or missing_value there). Raise "
              "remora.frc_fill_max_rings or fill the forcing over the ocean.");
    }
}

void
LatLonRegrid::interpolate (const FArrayBox& src_fab, MultiFab& dst) const
{
    auto const& src = src_fab.const_array();
    const int Nx = nlon, Ny = nlat, R = fill_max_rings;
    const bool periodic = periodic_lon, fx = flip_lon, fy = flip_lat;
    const Real nan = std::numeric_limits<Real>::quiet_NaN();

    for (MFIter mfi(dst); mfi.isValid(); ++mfi) {
        const Box bx = target_box(mfi);
        auto const& d   = dst.array(mfi);
        auto const& w   = weights.const_array(mfi);
        auto const& idx = indices.const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            if (idx(i,j,k,0) < 0) { d(i,j,k) = nan; return; }

            Real sum = Real(0.0), wsum = Real(0.0);
            for (int s = 0; s < stencil_size; ++s) {
                const Real v = src(idx(i,j,k,3*s), idx(i,j,k,3*s+1), 0);
                if (w(i,j,k,s) > Real(0.0) && is_valid(v)) {
                    sum  += w(i,j,k,s) * v;
                    wsum += w(i,j,k,s);
                }
            }
            if (wsum > Real(0.0)) { d(i,j,k) = sum / wsum; return; }

            // Every corner with weight is missing: take the nearest valid source value. The
            // target's fractional source index follows from the stencil: corner 0 is its
            // lower-left cell, and the fractions are p = w1 + w2 and q = w2 + w3, counted
            // backward along an axis the file stores decreasing.
            const Real p = w(i,j,k,1) + w(i,j,k,2);
            const Real q = w(i,j,k,2) + w(i,j,k,3);
            const Real xf = Real(idx(i,j,k,0)) + (fx ? -p : p);
            const Real yf = Real(idx(i,j,k,1)) + (fy ? -q : q);
            const int ci = static_cast<int>(std::round(xf));
            const int cj = static_cast<int>(std::round(yf));
            Real best = nan, best_d2 = std::numeric_limits<Real>::max();
            for (int dj = -R; dj <= R; ++dj) {
                const int n = cj + dj;
                if (n < 0 || n >= Ny) { continue; }
                for (int di = -R; di <= R; ++di) {
                    int m = ci + di;
                    if (periodic) { m = ((m % Nx) + Nx) % Nx; }
                    else if (m < 0 || m >= Nx) { continue; }
                    const Real v = src(m, n, 0);
                    if (!is_valid(v)) { continue; }
                    Real dx = Real(m) - xf;
                    if (periodic) { dx -= Real(Nx) * std::round(dx / Real(Nx)); }
                    const Real dy = Real(n) - yf;
                    const Real d2 = dx*dx + dy*dy;
                    if (d2 < best_d2) { best_d2 = d2; best = v; }
                }
            }
            d(i,j,k) = best;
        });
    }
}

bool
LatLonRegrid::same_source_grid (const std::string& file_name, const std::string& field_name,
                                const std::vector<long long>& field_shape) const
{
    LatLonRegridOptions opts;
    opts.lon_name = lon_name;
    opts.lat_name = lat_name;
    const SourceGrid g = find_source_grid(file_name, field_name, field_shape, opts, "");
    auto same = [] (const std::vector<double>& a, const std::vector<double>& b) {
        if (a.size() != b.size()) { return false; }
        for (std::size_t m = 0; m < a.size(); ++m) {
            if (std::abs(a[m] - b[m]) > 1.e-9 * std::max(1.0, std::abs(a[m]))) { return false; }
        }
        return true;
    };
    return same(g.lon, lon_file) && same(g.lat, lat_file) &&
           g.has_fill == static_cast<int>(has_fill) && g.has_missing == static_cast<int>(has_missing) &&
           (!has_fill || g.fill_value == fill_value) && (!has_missing || g.missing_value == missing_value) &&
           g.scale_factor == scale_factor && g.add_offset == add_offset;
}

bool
LatLonRegrid::has_latlon_axes (const std::string& file_name, const std::string& field_name)
{
    int found = 0;
    auto ncf = ncutils::NCFile::open(file_name, NC_NOCLOBBER);
    ncmpi_begin_indep_data(ncf.ncid);
    if (ParallelDescriptor::IOProcessor() && ncf.has_var(field_name)) {
        const auto var = ncf.var(field_name);
        const int rank = var.ndim();
        if (rank >= 2 && has_dimension_coordinate(ncf, var, rank-1) &&
                         has_dimension_coordinate(ncf, var, rank-2)) {
            found = units_say(ncf.var(dim_name(ncf, var, rank-1)), "east") &&
                    units_say(ncf.var(dim_name(ncf, var, rank-2)), "north");
        }
    }
    ncf.close();
    ParallelDescriptor::Bcast(&found, 1, ParallelDescriptor::IOProcessorNumber());
    return found != 0;
}

std::string
LatLonRegrid::summary () const
{
    std::ostringstream os;
    os << "a " << nlon << " x " << nlat << " lat-lon grid (" << lon_name << " "
       << lon_file.front() << " to " << lon_file.back() << ", " << lat_name << " "
       << lat_file.front() << " to " << lat_file.back()
       << (periodic_lon ? ", periodic in longitude" : "") << ")";
    return os.str();
}

Box
LatLonRegrid::target_box (const MFIter& mfi) const
{
    Box bx = mfi.validbox();
    bx.grow(IntVect(1,1,0));
    bx &= amrex::grow(domain, IntVect(1,1,0));
    return bx;
}

#endif // REMORA_USE_NETCDF
