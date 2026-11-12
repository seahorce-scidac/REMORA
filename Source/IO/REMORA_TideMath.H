#ifndef REMORA_TideMath_H
#define REMORA_TideMath_H

#include "REMORA_Constants.H"

#include <cmath>

/**
 * \brief Convert a phase lag in degrees to the cosine/sine coefficients of
 * A*cos(omega*t-phase), as in ROMS Nonlinear/set_tides.F.
 *
 * @param[in   ] amplitude  elevation or velocity amplitude
 * @param[in   ] phase      phase lag in degrees
 * @param[  out] cosine     coefficient multiplying cos(omega*t)
 * @param[  out] sine       coefficient multiplying sin(omega*t)
 */
inline void
remora_tide_coefficients (double amplitude, double phase, double& cosine, double& sine)
{
    const double angle = std::remainder(phase, 360.0) * static_cast<double>(PI) / 180.0;
    cosine = amplitude * std::cos(angle);
    sine = amplitude * std::sin(angle);
}

/**
 * \brief Harmonic time argument on the model clock, as in ROMS set_tides.F.
 *
 * @param[in   ] time    time in seconds on the model clock
 * @param[in   ] epoch   phase-reference time in seconds on the model clock
 * @param[in   ] period  positive period in seconds
 */
inline double
remora_tide_argument (double time, double epoch, double period)
{
    return 2.0 * static_cast<double>(PI) * (std::remainder(time - epoch, period) / period);
}

#endif
