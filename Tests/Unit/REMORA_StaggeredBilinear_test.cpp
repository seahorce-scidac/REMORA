/**
 * Checks the bilinear value StaggeredBilinear computes from four coarse values. Without a period
 * a bilinear field is reproduced exactly, extrapolation included. With a period, a longitude or
 * angle that crosses a branch cut must give the value of the continuous field, taken to the
 * branch of the coarse value with the largest weight, so the cut stays between the same coarse
 * points.
 */
#include "REMORA_StaggeredBilinear.H"

#include <cmath>
#include <cstdio>
#include <limits>

int main ()
{
    using amrex::Real;
    int failures = 0;
    const double eps = std::numeric_limits<Real>::epsilon();
    const double pi = 3.14159265358979323846;

    auto check = [&] (const char* what, double got, double expected, double tol) {
        if (!(std::abs(got - expected) <= tol)) {
            std::printf("%s: got %.17g, expected %.17g\n", what, got, expected);
            ++failures;
        }
    };
    auto value = [] (double c00, double c10, double c01, double c11, double wx, double wy,
                     double period) {
        return static_cast<double>(staggered_bilinear_value(Real(c00), Real(c10), Real(c01),
                                                            Real(c11), Real(wx), Real(wy),
                                                            Real(period)));
    };

    // A bilinear field, inside the cell and extrapolated past it
    auto field = [] (double x, double y) { return 2.0 + 0.5 * x - 1.25 * y + 0.75 * x * y; };
    for (const double wx : {-0.75, -0.25, 0.0, 0.25, 0.5, 0.75, 1.25}) {
        for (const double wy : {-0.25, 0.25, 0.75, 1.75}) {
            check("bilinear field, no period",
                  value(field(0, 0), field(1, 0), field(0, 1), field(1, 1), wx, wy, 0),
                  field(wx, wy), 64 * eps);
        }
    }

    // Longitude 179.5 then -179.5 along x, in [-180, 180): each fine value keeps the branch of
    // its nearer coarse value, inside the cell and extrapolated past it
    const double ltol = 64 * eps * 360;
    check("date line, near the west point", value(179.5, -179.5, 179.5, -179.5, 0.25, 0.5, 360),
          179.75, ltol);
    check("date line, near the east point", value(179.5, -179.5, 179.5, -179.5, 0.75, 0.5, 360),
          -179.75, ltol);
    check("date line, extrapolated west", value(179.5, -179.5, 179.5, -179.5, -0.25, 0.5, 360),
          179.25, ltol);
    check("date line, extrapolated east", value(179.5, -179.5, 179.5, -179.5, 1.25, 0.5, 360),
          -179.25, ltol);

    // The same crossing at 0 in [0, 360)
    check("prime meridian, near the west point", value(359.5, 0.5, 359.5, 0.5, 0.25, 0.5, 360),
          359.75, ltol);
    check("prime meridian, near the east point", value(359.5, 0.5, 359.5, 0.5, 0.75, 0.5, 360),
          0.25, ltol);

    // A cut crossing the cell diagonally, at weights nearest c11: the continuous field is
    // 179 + 2 wx + 0.5 wy + 0.5 wx wy (c10 = 181, c01 = 179.5, c11 = 182 unwrapped), taken to
    // c11's branch
    {
        const double wx = 0.75, wy = 0.75;
        const double continuous = 179.0 + 2.0 * wx + 0.5 * wy + 0.5 * wx * wy;
        check("diagonal cut", value(179.0, -179.0, 179.5, -178.0, wx, wy, 360),
              continuous - 360.0, ltol);
    }

    // An angle crossing pi
    check("angle across pi",
          value(pi - 0.05, -pi + 0.05, pi - 0.05, -pi + 0.05, 0.25, 0.25, 2 * pi),
          pi - 0.025, 64 * eps * 2 * pi);
    check("angle across pi, near the far side",
          value(pi - 0.05, -pi + 0.05, pi - 0.05, -pi + 0.05, 0.75, 0.25, 2 * pi),
          -pi + 0.025, 64 * eps * 2 * pi);

    std::printf("Staggered bilinear interpolation: %d failures\n", failures);
    return failures != 0;
}
