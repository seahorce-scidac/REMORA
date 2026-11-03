/**
 * \file REMORA_LatLonRegrid_test.cpp
 *
 * Unit tests for Source/IO/REMORA_LatLonRegridMath.H, the arithmetic of interpolating forcing
 * from a lat-lon grid onto REMORA's rho points and rotating vectors onto the grid axes.
 *
 * The expectations come from properties rather than from the implementation:
 *
 *   - Bilinear weights are a partition of unity, so a constant is reproduced exactly, and
 *     they reproduce any a + b lon + c lat + d lon lat exactly on a rectilinear grid, even a
 *     non-uniform one.
 *   - Storing an axis decreasing moves the stencil's indices to the mirrored ones and leaves
 *     the interpolated value unchanged.
 *   - A longitude is the same point modulo 360, so wrapping, unwrapping and the periodic
 *     interval of a global grid must give the value an unwrapped grid gives.
 *   - Rotation by the grid angle is the ROMS formula: rotating east by pi/2 gives -eta, and
 *     a rotation keeps the magnitude.
 *
 * Every case runs in either build precision; tolerances scale with machine epsilon.
 */

#include "REMORA_LatLonRegridMath.H"

#include <cmath>
#include <cstdio>
#include <limits>
#include <vector>

using amrex::Real;
using namespace latlon_regrid;

namespace {

int failures = 0;
int checks = 0;

const Real eps = std::numeric_limits<Real>::epsilon();

void check (const char* what, bool ok)
{
    ++checks;
    if (!ok) { ++failures; std::printf("FAIL %s\n", what); }
}

void check_int (const char* what, int got, int expect)
{
    ++checks;
    if (got != expect) { ++failures; std::printf("FAIL %s: got %d, expected %d\n", what, got, expect); }
}

void check_near (const char* what, double got, double expect, double tol)
{
    ++checks;
    if (!(std::abs(got - expect) <= tol)) {
        ++failures;
        std::printf("FAIL %s: got %.17g, expected %.17g (tol %.3g)\n", what, got, expect, tol);
    }
}

//! Interpolate a field sampled on (lon, lat), stored in the file's index order, with a stencil
Real apply (const Stencil& s, const std::vector<Real>& f, int nlon)
{
    Real v = 0;
    for (int k = 0; k < 4; ++k) { v += s.w[k] * f[s.j[k] * nlon + s.i[k]]; }
    return v;
}

Real bilinear_field (Real x, Real y) { return Real(3.0) - Real(0.2) * x + Real(0.5) * y + Real(0.01) * x * y; }

void test_bracket ()
{
    const Real a[5] = {0, 1, 2.5, 4, 7};
    check_int("bracket below the axis", bracket(a, 5, Real(-1)), -1);
    check_int("bracket on the first node", bracket(a, 5, Real(0)), 0);
    check_int("bracket inside", bracket(a, 5, Real(3)), 2);
    check_int("bracket on an interior node", bracket(a, 5, Real(2.5)), 2);
    check_int("bracket on the last node", bracket(a, 5, Real(7)), 4);
    check_int("bracket above the axis", bracket(a, 5, Real(9)), 4);
}

void test_wrap_and_unwrap ()
{
    check_near("wrap -190 into [-180, 180)", wrap_longitude(Real(-190), Real(-180)), 170, 1e3*eps);
    check_near("wrap 190 into [-180, 180)",  wrap_longitude(Real(190),  Real(-180)), -170, 1e3*eps);
    check_near("wrap -0.5 into [0, 360)",    wrap_longitude(Real(-0.5), Real(0)), 359.5, 1e3*eps);
    check_near("wrap 720.25 into [0, 360)",  wrap_longitude(Real(720.25), Real(0)), 0.25, 1e3*eps);
    check_near("wrap leaves a point in the window", wrap_longitude(Real(10), Real(-180)), 10, 0);

    std::vector<double> lon = {170, 175, 180, -175, -170};
    unwrap_longitudes(lon);
    const double expect[5] = {170, 175, 180, 185, 190};
    bool ok = true;
    for (int m = 0; m < 5; ++m) { ok = ok && lon[m] == expect[m]; }
    check("unwrap continues an axis across the dateline", ok);

    std::vector<double> lon360 = {350, 355, 0, 5};
    unwrap_longitudes(lon360);
    check("unwrap continues an axis across Greenwich in [0, 360)", lon360[2] == 360 && lon360[3] == 365);

    std::vector<double> global, regional, closed;
    for (int m = 0; m < 144; ++m) { global.push_back(2.5 * m); }
    for (int m = 0; m <= 20; ++m) { regional.push_back(-130 + m); }
    for (int m = 0; m <= 144; ++m) { closed.push_back(2.5 * m); }
    check("a 0..357.5 axis at 2.5 degrees is periodic", is_periodic(global));
    check("a regional axis is not periodic", !is_periodic(regional));
    check("an axis repeating its first point at 360 is not periodic", !is_periodic(closed));
    check_near("reference of a periodic axis is its first point", longitude_reference(global, true), 0, 0);
    check_near("reference of a regional axis centers it", longitude_reference(regional, false), -300, 0);
}

void test_stencil_properties ()
{
    // Non-uniform axes
    const std::vector<Real> lon = {-3.0, -1.0, 0.5, 3.0, 4.0};
    const std::vector<Real> lat = {10.0, 10.5, 12.0, 15.0};
    const int nx = 5, ny = 4;
    std::vector<Real> f(nx*ny);
    for (int j = 0; j < ny; ++j) { for (int i = 0; i < nx; ++i) { f[j*nx+i] = bilinear_field(lon[i], lat[j]); } }

    const Real pts[6][2] = {{-2.2, 10.1}, {0.5, 12.0}, {3.9, 14.9}, {-3.0, 10.0}, {4.0, 15.0}, {1.7, 11.3}};
    double worst_sum = 0, worst_exact = 0;
    bool all_inside = true, all_nonneg = true;
    for (const auto& pt : pts) {
        const Stencil s = bilinear_stencil(lon.data(), nx, lat.data(), ny, false, false, false,
                                           pt[0], pt[1], Real(1e-9));
        all_inside = all_inside && s.inside;
        Real sum = 0;
        for (int k = 0; k < 4; ++k) { sum += s.w[k]; all_nonneg = all_nonneg && s.w[k] >= 0; }
        worst_sum = std::max(worst_sum, std::abs(double(sum) - 1.0));
        worst_exact = std::max(worst_exact, std::abs(double(apply(s, f, nx)) - double(bilinear_field(pt[0], pt[1]))));
    }
    check("points on and inside a non-uniform grid are inside", all_inside);
    check("bilinear weights are not negative", all_nonneg);
    check_near("bilinear weights sum to one", worst_sum, 0, 10*eps);
    check_near("a bilinear field is reproduced exactly", worst_exact, 0, 200*eps);

    // Outside, on each side
    const Real out[4][2] = {{-3.5, 11}, {4.5, 11}, {0, 9.5}, {0, 15.5}};
    bool none_inside = true;
    for (const auto& pt : out) {
        none_inside = none_inside && !bilinear_stencil(lon.data(), nx, lat.data(), ny, false, false, false,
                                                       pt[0], pt[1], Real(1e-9)).inside;
    }
    check("points beyond each side of the grid are outside", none_inside);
}

void test_flipped_axes ()
{
    // The same grid stored with both axes decreasing, as ERA5 stores latitude
    const std::vector<Real> lon = {0, 1, 2, 3}, lat = {0, 2, 4};
    const int nx = 4, ny = 3;
    std::vector<Real> f(nx*ny), g(nx*ny);
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            f[j*nx+i] = bilinear_field(lon[i], lat[j]);
            g[(ny-1-j)*nx + (nx-1-i)] = f[j*nx+i];        // the file's order when both are flipped
        }
    }
    const Real x = 1.3, y = 2.7;
    const Stencil a = bilinear_stencil(lon.data(), nx, lat.data(), ny, false, false, false, x, y, Real(1e-9));
    const Stencil b = bilinear_stencil(lon.data(), nx, lat.data(), ny, false, true, true, x, y, Real(1e-9));
    bool mirrored = true;
    for (int k = 0; k < 4; ++k) { mirrored = mirrored && b.i[k] == nx-1-a.i[k] && b.j[k] == ny-1-a.j[k] && b.w[k] == a.w[k]; }
    check("flipping both axes mirrors the stencil's indices and keeps its weights", mirrored);
    check_near("flipped storage interpolates to the same value", apply(b, g, nx), apply(a, f, nx), 0);
}

void test_periodic ()
{
    // Global 30-degree grid, field linear in a longitude continuous across Greenwich
    std::vector<Real> lon;
    for (int m = 0; m < 12; ++m) { lon.push_back(Real(30.0) * m); }
    const std::vector<Real> lat = {-10, 10};
    const int nx = 12, ny = 2;
    auto L = [] (Real x) { return x >= 180 ? x - 360 : x; };          // continuous near 0
    std::vector<Real> f(nx*ny);
    for (int j = 0; j < ny; ++j) { for (int i = 0; i < nx; ++i) { f[j*nx+i] = bilinear_field(L(lon[i]), lat[j]); } }

    // -10 degrees east lies in the interval from 330 to 360 = 0
    const Real lt = wrap_longitude(Real(-10), lon[0]);
    const Stencil s = bilinear_stencil(lon.data(), nx, lat.data(), ny, true, false, false, lt, Real(3), Real(1e-9));
    check("a target west of Greenwich on a global grid is inside", s.inside);
    check("the periodic interval joins the last and first columns",
          s.i[0] == nx-1 && s.i[1] == 0 && s.i[2] == 0 && s.i[3] == nx-1);
    check_near("the periodic interval interpolates as an unwrapped grid would",
               apply(s, f, nx), bilinear_field(Real(-10), Real(3)), 200*eps);

    const Stencil e = bilinear_stencil(lon.data(), nx, lat.data(), ny, true, false, false,
                                       wrap_longitude(Real(0), lon[0]), Real(0), Real(1e-9));
    check_near("a target on the first column takes its value", apply(e, f, nx), bilinear_field(0, 0), 200*eps);
}

void test_rotation ()
{
    const Real pi = std::acos(Real(-1));
    Real xi, eta;
    rotate_to_grid(Real(1), Real(0), Real(0), xi, eta);
    check("rotation by zero is the identity", xi == 1 && eta == 0);

    rotate_to_grid(Real(1), Real(0), pi/2, xi, eta);
    check_near("east on a grid turned by pi/2: xi", xi, 0, 10*eps);
    check_near("east on a grid turned by pi/2: eta", eta, -1, 10*eps);

    rotate_to_grid(Real(0), Real(1), pi/2, xi, eta);
    check_near("north on a grid turned by pi/2: xi", xi, 1, 10*eps);
    check_near("north on a grid turned by pi/2: eta", eta, 0, 10*eps);

    const Real e = 3.7, n = -1.2, a = 0.83;
    rotate_to_grid(e, n, a, xi, eta);
    check_near("rotation keeps the magnitude", std::hypot(xi, eta), std::hypot(e, n), 20*eps);
    Real back_e, back_n;
    rotate_to_grid(xi, eta, -a, back_e, back_n);
    check_near("rotating back recovers east", back_e, e, 50*eps);
    check_near("rotating back recovers north", back_n, n, 50*eps);
}

} // namespace

int main ()
{
    test_bracket();
    test_wrap_and_unwrap();
    test_stencil_properties();
    test_flipped_axes();
    test_periodic();
    test_rotation();

    std::printf("LatLonRegrid: %d of %d checks passed\n", checks - failures, checks);
    return failures == 0 ? 0 : 1;
}
