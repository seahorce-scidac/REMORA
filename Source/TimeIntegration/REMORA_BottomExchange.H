#ifndef REMORA_BOTTOM_EXCHANGE_H_
#define REMORA_BOTTOM_EXCHANGE_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>

namespace BottomExchange {

/** Advective tracer flux, with incoming water positive into the ocean.
 * @param[in   ] volume integrated water transport (m^3/s)
 * @param[in   ] incoming prescribed incoming tracer transport (concentration m^3/s)
 * @param[in   ] ocean ocean-side donor concentration
 * @return signed tracer transport (concentration m^3/s)
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
tracer_flux (amrex::Real volume, amrex::Real incoming, amrex::Real ocean) noexcept
{
    return incoming + (volume < amrex::Real(0) ? volume * ocean : amrex::Real(0));
}

/** Dimensionless mechanical seepage measures, independent of buoyancy strength. */
struct Scales {
    amrex::Real transport;   ///< Maximum advective-to-diffusive transport ratio.
    amrex::Real adjustment;  ///< Advective rate relative to GLS dissipation.
    amrex::Real wall;        ///< Momentum flux relative to wall stress; -1 at zero stress.
};

/** Compare seepage with GLS mixing, dissipation and wall stress.
 * @param[in   ] speed magnitude of water flux per horizontal area (m/s)
 * @param[in   ] d near-bed variation length (m)
 * @param[in   ] km momentum diffusivity (m^2/s)
 * @param[in   ] kk TKE diffusivity (m^2/s)
 * @param[in   ] kp GLS diffusivity (m^2/s)
 * @param[in   ] k turbulent kinetic energy (m^2/s^2)
 * @param[in   ] dissipation TKE dissipation rate (m^2/s^3)
 * @param[in   ] stress kinematic wall-stress magnitude (m^2/s^2)
 * @param[in   ] horizontal_speed near-bed horizontal speed (m/s)
 * @return ratios to assess the documented weak-seepage approximation
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Scales
scales (amrex::Real speed, amrex::Real d, amrex::Real km, amrex::Real kk,
        amrex::Real kp, amrex::Real k, amrex::Real dissipation, amrex::Real stress,
        amrex::Real horizontal_speed) noexcept
{
    const amrex::Real minimum = km < kk ? (km < kp ? km : kp) : (kk < kp ? kk : kp);
    return {speed * d / minimum, speed * k / (d * dissipation),
            stress > amrex::Real(0) ? speed * horizontal_speed / stress : amrex::Real(-1)};
}

} // namespace BottomExchange
#endif
