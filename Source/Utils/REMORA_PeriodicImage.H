#ifndef REMORA_PERIODIC_IMAGE_H_
#define REMORA_PERIODIC_IMAGE_H_

#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

#include <cmath>

/**
 * The image of a periodic value (a longitude, period 360, or an angle, period 2 pi) nearest
 * a reference value: x shifted by a whole number of periods to within half a period of ref.
 *
 * @param[in   ] x       the value
 * @param[in   ] ref     the reference value
 * @param[in   ] period  the period, positive
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
nearest_periodic_image (amrex::Real x, amrex::Real ref, amrex::Real period)
{
    return x - period * std::round((x - ref) / period);
}

#endif
