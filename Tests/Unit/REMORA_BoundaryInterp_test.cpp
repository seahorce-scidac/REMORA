/**
 * Checks the interpolation of boundary data along a boundary strip onto a refined level.
 * Linear interpolation between coarse points reproduces a linear profile exactly, so every
 * fine point must get the profile's value there. Positions are in coarse cell widths. A
 * cell-centered strip runs over coarse cells -1..n and fine cells -1..r*n, with centers at
 * c + 0.5 and (f + 0.5)/r. A nodal strip runs over coarse nodes 0..n and fine nodes 0..r*n,
 * at c and f/r. These are the extents of the boundary strips.
 */
#include "REMORA_BoundaryInterpMath.H"

#include <cmath>
#include <cstdio>
#include <limits>

int main ()
{
    int failures = 0;
    const int n = 7;
    const double tol = 64 * std::numeric_limits<amrex::Real>::epsilon();
    auto profile = [] (double x) { return 3.0 - 0.75 * x; };

    // Check one fine point: the coarse points it uses lie in [first, n] and give the profile
    auto check = [&] (const char* kind, int r, int f, int c, amrex::Real w, int first,
                      double crse_offset, double expected) {
        const bool in_strip = c >= first && (c + (w > 0 ? 1 : 0)) <= n;
        const double crse_lo = profile(c + crse_offset);
        const double crse_hi = profile(c + 1 + crse_offset);
        const double value = crse_lo + w * (crse_hi - crse_lo);
        if (!in_strip || !(w >= 0 && w < 1) || std::abs(value - expected) > tol) {
            std::printf("%s, r = %d, fine point %d: coarse point %d, weight %g, value %.17g,"
                        " expected %.17g\n", kind, r, f, c, static_cast<double>(w), value,
                        expected);
            ++failures;
        }
    };

    for (int r = 1; r <= 4; ++r) {
        for (int f = -1; f <= r*n; ++f) {
            int c;
            amrex::Real w;
            boundary_interp::cell_centered_weight(f, r, c, w);
            check("cell-centered", r, f, c, w, -1, 0.5, profile((f + 0.5) / r));
        }
        for (int f = 0; f <= r*n; ++f) {
            int c;
            amrex::Real w;
            boundary_interp::nodal_weight(f, r, c, w);
            check("nodal", r, f, c, w, 0, 0.0, profile(static_cast<double>(f) / r));
        }
    }
    std::printf("Boundary interpolation: %d failures\n", failures);
    return failures != 0;
}
