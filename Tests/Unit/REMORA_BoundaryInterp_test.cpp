/**
 * Checks the interpolation of cell-centered boundary data along a boundary strip onto a
 * refined level. Linear interpolation between coarse cell centers reproduces a linear
 * profile exactly, so every fine cell center must get the profile's value there. The strip
 * runs over coarse cells -1..n and fine cells -1..r*n, as the boundary strips do.
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
    for (int r = 1; r <= 4; ++r) {
        // A linear profile in position x, in coarse cell widths: coarse cell c is centered at
        // c + 0.5 and fine cell f at (f + 0.5)/r.
        auto profile = [] (double x) { return 3.0 - 0.75 * x; };
        for (int f = -1; f <= r*n; ++f) {
            int c;
            amrex::Real w;
            boundary_interp::cell_centered_weight(f, r, c, w);
            const bool in_strip = c >= -1 && (c + (w > 0 ? 1 : 0)) <= n;
            const double crse_lo = profile(c + 0.5);
            const double crse_hi = profile(c + 1.5);
            const double value = (1.0 - w) * crse_lo + w * crse_hi;
            const double expected = profile((f + 0.5) / r);
            if (!in_strip || !(w >= 0 && w < 1) || std::abs(value - expected) > tol) {
                std::printf("r = %d, fine cell %d: coarse cell %d, weight %g, value %.17g,"
                            " expected %.17g\n", r, f, c, static_cast<double>(w), value,
                            expected);
                ++failures;
            }
        }
    }
    std::printf("Boundary interpolation: %d failures\n", failures);
    return failures != 0;
}
