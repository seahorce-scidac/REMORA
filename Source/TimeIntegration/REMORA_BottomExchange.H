#ifndef REMORA_BOTTOM_EXCHANGE_H_
#define REMORA_BOTTOM_EXCHANGE_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>

namespace BottomExchange {

/** Bottom fluxes are positive into the ocean. Incoming transport is supplied
 * separately so refinement preserves Q*C, including where C is zero (fresh water).
 * The donor concentration is used only on outflow. */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real tracer_flux (amrex::Real volume, amrex::Real incoming,
                         amrex::Real ocean) noexcept
{
    return incoming + ((volume < 0) ? volume * ocean : amrex::Real(0));
}

/** Mechanical seepage measures; these are diagnostics, not closure coefficients.
 * d is a resolved near-bed variation length; eps is the GLS dissipation rate.
 * A zero physical stress has no meaningful wall ratio (reported separately). */
struct Scales {
    amrex::Real transport;
    amrex::Real adjustment;
    amrex::Real wall;
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Scales scales (amrex::Real speed, amrex::Real d, amrex::Real km,
               amrex::Real kk, amrex::Real kp, amrex::Real k,
               amrex::Real dissipation, amrex::Real stress,
               amrex::Real horizontal_speed) noexcept
{
    const amrex::Real minimum = (km < kk) ? ((km < kp) ? km : kp)
                                        : ((kk < kp) ? kk : kp);
    return {speed*d/minimum, speed*k/(d*dissipation),
            (stress > 0) ? speed*horizontal_speed/stress : amrex::Real(-1)};
}

} // namespace BottomExchange
#endif
