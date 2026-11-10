#include "REMORA_TideMath.H"

#include <cmath>
#include <cstdio>

int main ()
{
    int failures = 0;
    const double pi = std::acos(-1.0);
    auto check = [&] (const char* name, double actual, double expected) {
        if (std::abs(actual - expected) > 1.0e-12) {
            std::printf("FAIL %s: %.17g != %.17g\n", name, actual, expected);
            ++failures;
        }
    };
    double cosine, sine;
    remora_tide_coefficients(2.0, 90.0, cosine, sine);
    check("90 degree lag at the epoch", cosine, 0.0);
    check("90 degree lag a quarter period later", sine, 2.0);
    remora_tide_coefficients(3.0, -180.0, cosine, sine);
    check("negative phase", cosine, -3.0);
    check("negative phase sine", sine, 0.0);
    remora_tide_coefficients(3.0, 540.0, cosine, sine);
    check("phase beyond one cycle", cosine, -3.0);
    check("quarter period", remora_tide_argument(63713433600.25, 63713433600.0, 1.0), pi/2);
    check("before the epoch", remora_tide_argument(99.0, 100.0, 4.0), -pi/2);
    check("many cycles", std::cos(remora_tide_argument(1.0e9+30.0, 30.0, 4.0)), 1.0);
    double c0, s0, c1, s1;
    remora_tide_coefficients(1.0, 359.0, c0, s0);
    remora_tide_coefficients(1.0, 1.0, c1, s1);
    check("phase wrap interpolation cosine", (c0+c1)/2, std::cos(pi/180));
    check("phase wrap interpolation sine", (s0+s1)/2, 0.0);
    return failures != 0;
}
