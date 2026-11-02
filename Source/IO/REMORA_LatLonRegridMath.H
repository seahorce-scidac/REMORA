#ifndef REMORA_LATLONREGRIDMATH_H_
#define REMORA_LATLONREGRIDMATH_H_

/**
 * \file REMORA_LatLonRegridMath.H
 *
 * The arithmetic of interpolating forcing from a rectilinear lat-lon grid onto REMORA's rho
 * points and rotating vectors onto the grid axes (LatLonRegrid, REMORA::rotate_to_grid_axes),
 * kept free of MultiFabs and NetCDF so that Tests/Unit can check it directly. Nothing here
 * depends on REMORA_USE_NETCDF.
 */

#include <AMReX_Algorithm.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace latlon_regrid {

/** Largest m with a[m] <= x in an increasing array of n values, or -1 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int bracket (const amrex::Real* a, int n, amrex::Real x)
{
    int lo = -1, hi = n;    // a[lo] <= x < a[hi], with a[-1] = -inf and a[n] = +inf
    while (hi - lo > 1) {
        const int mid = (lo + hi) / 2;
        if (a[mid] <= x) { lo = mid; } else { hi = mid; }
    }
    return lo;
}

/** Whether a source value is data rather than a missing-data marker turned to NaN */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool is_valid (amrex::Real v)
{
    return (v == v) && std::abs(v) <= std::numeric_limits<amrex::Real>::max();
}

/** A longitude taken into [ref, ref + 360) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real wrap_longitude (amrex::Real lon, amrex::Real ref)
{
    return lon - amrex::Real(360.0) * std::floor((lon - ref) / amrex::Real(360.0));
}

/** Four bilinear weights and the source cell (i, j) of each, in the file's index order:
 *  corner 0 is the lower-left of the target's cell in increasing coordinates, then
 *  counterclockwise. `inside` is false when the target lies outside the source grid. */
struct Stencil
{
    int i[4];
    int j[4];
    amrex::Real w[4];
    bool inside;
};

/**
 * The bilinear stencil of one target point.
 *
 * @param lon       longitude axis, increasing and continuous (see unwrap_longitudes)
 * @param nlon      its length
 * @param lat       latitude axis, increasing
 * @param nlat      its length
 * @param periodic  whether the longitude axis closes the circle, its last interval running
 *                  from lon[nlon-1] to lon[0] + 360
 * @param flip_lon  whether the file stores longitude decreasing, so that index a in the
 *                  increasing axis is nlon-1-a in the file
 * @param flip_lat  likewise for latitude
 * @param lt        target longitude, already in the axis's window (see wrap_longitude)
 * @param pt        target latitude
 * @param tol       how far outside the axes a target may lie and still count as inside
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Stencil bilinear_stencil (const amrex::Real* lon, int nlon, const amrex::Real* lat, int nlat,
                          bool periodic, bool flip_lon, bool flip_lat,
                          amrex::Real lt, amrex::Real pt, amrex::Real tol)
{
    using amrex::Real;
    int m, m2;
    Real p;
    bool inside = true;
    if (periodic) {
        m = amrex::max(0, bracket(lon, nlon, lt));
        m2 = (m + 1) % nlon;
        const Real right = (m == nlon-1) ? lon[0] + Real(360.0) : lon[m+1];
        p = (lt - lon[m]) / (right - lon[m]);
    } else {
        inside = (lt >= lon[0] - tol) && (lt <= lon[nlon-1] + tol);
        m = amrex::min(amrex::max(bracket(lon, nlon, lt), 0), nlon-2);
        m2 = m + 1;
        p = (lt - lon[m]) / (lon[m2] - lon[m]);
    }
    inside = inside && (pt >= lat[0] - tol) && (pt <= lat[nlat-1] + tol);
    const int n  = amrex::min(amrex::max(bracket(lat, nlat, pt), 0), nlat-2);
    const int n2 = n + 1;
    Real q = (pt - lat[n]) / (lat[n2] - lat[n]);
    p = amrex::min(amrex::max(p, Real(0.0)), Real(1.0));
    q = amrex::min(amrex::max(q, Real(0.0)), Real(1.0));

    const int fm  = flip_lon ? nlon-1-m  : m;
    const int fm2 = flip_lon ? nlon-1-m2 : m2;
    const int fn  = flip_lat ? nlat-1-n  : n;
    const int fn2 = flip_lat ? nlat-1-n2 : n2;

    Stencil s;
    s.i[0] = fm;  s.i[1] = fm2; s.i[2] = fm2; s.i[3] = fm;
    s.j[0] = fn;  s.j[1] = fn;  s.j[2] = fn2; s.j[3] = fn2;
    s.w[0] = (Real(1.0)-p)*(Real(1.0)-q);
    s.w[1] = p*(Real(1.0)-q);
    s.w[2] = p*q;
    s.w[3] = (Real(1.0)-p)*q;
    s.inside = inside;
    return s;
}

/**
 * Rotate a vector given as east and north components onto grid axes at angle a from east,
 * counterclockwise (ROMS's `angle`), as ROMS's set_data does:
 *     xi = E cos a + N sin a,   eta = -E sin a + N cos a.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void rotate_to_grid (amrex::Real east, amrex::Real north, amrex::Real angle,
                     amrex::Real& xi, amrex::Real& eta)
{
    const amrex::Real c = std::cos(angle);
    const amrex::Real s = std::sin(angle);
    xi  =  east * c + north * s;
    eta = -east * s + north * c;
}

/** Make a longitude axis continuous across the branch cut (170 .. 180, -180 .. -170 becomes
 *  170 .. 190), each step taken modulo 360 */
inline void unwrap_longitudes (std::vector<double>& lon)
{
    for (std::size_t m = 1; m < lon.size(); ++m) {
        lon[m] -= 360.0 * std::round((lon[m] - lon[m-1]) / 360.0);
    }
}

/** The span in degrees of an increasing, continuous longitude axis */
inline double longitude_span (const std::vector<double>& lon_inc)
{
    return lon_inc.empty() ? 0.0 : lon_inc.back() - lon_inc.front();
}

/** Whether an increasing, continuous longitude axis closes the circle: its span plus one
 *  spacing is 360 degrees */
inline bool is_periodic (const std::vector<double>& lon_inc)
{
    if (lon_inc.size() < 2) { return false; }
    const double span = longitude_span(lon_inc);
    const double dlon = span / static_cast<double>(lon_inc.size() - 1);
    return std::abs(span + dlon - 360.0) < 1.e-3 * dlon;
}

/** The start of the longitude window targets are wrapped into: the axis's first point when
 *  it covers the globe, and otherwise centered on the axis, so that a point outside it is
 *  reported on the side it misses */
inline double longitude_reference (const std::vector<double>& lon_inc, bool periodic)
{
    return periodic ? lon_inc.front() : 0.5 * (lon_inc.front() + lon_inc.back()) - 180.0;
}

} // namespace latlon_regrid

#endif // REMORA_LATLONREGRIDMATH_H_
