#ifndef REMORA_BIOLOGY_GAS_EXCHANGE_H_
#define REMORA_BIOLOGY_GAS_EXCHANGE_H_

#include <AMReX_Array4.H>
#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>
#include <cmath>

namespace REMORABiology {

/**
 * Air-sea O2 and CO2 exchange shared by the Fennel and ECB kernels, whose ROMS
 * sources (fennel.h, ecb.h) carry identical copies of these formulas.
 *
 * Each coefficient set belongs to one published relation, so the Schmidt
 * polynomial and the leading rate coefficient of the transfer velocity are
 * selected together; ROMS makes the same pairings with RW14_OXYGEN_SC,
 * OCMIP_OXYGEN_SC and RW14_CO2_SC. Pick a set on the host and capture it by
 * value: the kernels never read these namespace-scope objects on the device.
 */
struct SchmidtCoefficients {
    amrex::Real A;
    amrex::Real B;
    amrex::Real C;
    amrex::Real D;
    amrex::Real E;
    amrex::Real rate;   //!< leading coefficient of the transfer velocity [cm/hr]
};

namespace GasExchange {
//! O2, Wanninkhof (1992)
inline constexpr SchmidtCoefficients O2Wanninkhof1992 {
    amrex::Real(1953.4), amrex::Real(128.0), amrex::Real(3.9918),
    amrex::Real(0.050091), amrex::Real(0.0), amrex::Real(0.31)};
//! O2, Wanninkhof (2014)
inline constexpr SchmidtCoefficients O2Wanninkhof2014 {
    amrex::Real(1920.4), amrex::Real(135.6), amrex::Real(5.2122),
    amrex::Real(0.10939), amrex::Real(0.00093777), amrex::Real(0.251)};
//! O2, Keeling et al. (1998), as OCMIP. ROMS pairs it with the 1992 rate.
inline constexpr SchmidtCoefficients O2OCMIP {
    amrex::Real(1638.0), amrex::Real(81.83), amrex::Real(1.483),
    amrex::Real(0.008004), amrex::Real(0.0), amrex::Real(0.31)};
//! CO2, Wanninkhof (1992)
inline constexpr SchmidtCoefficients CO2Wanninkhof1992 {
    amrex::Real(2073.1), amrex::Real(125.62), amrex::Real(3.6276),
    amrex::Real(0.043219), amrex::Real(0.0), amrex::Real(0.31)};
//! CO2, Wanninkhof (2014)
inline constexpr SchmidtCoefficients CO2Wanninkhof2014 {
    amrex::Real(2116.8), amrex::Real(136.25), amrex::Real(4.7353),
    amrex::Real(0.092307), amrex::Real(0.0007555), amrex::Real(0.251)};
}

/** \brief Schmidt number at temperature \p T [degC]. */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real schmidt_number (SchmidtCoefficients const& sc, amrex::Real T) noexcept
{
    return sc.A - T * (sc.B - T * (sc.C - T * (sc.D - T * sc.E)));
}

/** \brief Transfer-velocity rate factor for one step of \p dtdays days. */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real transfer_rate (SchmidtCoefficients const& sc, amrex::Real dtdays) noexcept
{
    return dtdays * sc.rate * amrex::Real(24.0) / amrex::Real(100.0);
}

/** \brief Squared 10 m wind speed [m2/s2] at rho point (i,j).
 *
 * From the winds under bulk fluxes, otherwise from the surface stress
 * averaged to the rho point, as ROMS does without BULK_FLUXES.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real u10_squared (bool bulk_fluxes,
                         amrex::Array4<amrex::Real const> const& uwind,
                         amrex::Array4<amrex::Real const> const& vwind,
                         amrex::Array4<amrex::Real const> const& sustr,
                         amrex::Array4<amrex::Real const> const& svstr,
                         int i, int j, amrex::Real rho0) noexcept
{
    using amrex::Real;
    if (bulk_fluxes) {
        return uwind(i,j,0) * uwind(i,j,0) + vwind(i,j,0) * vwind(i,j,0);
    }
    const Real cff1 = rho0 * Real(550.0);
    return cff1 *
        std::sqrt((Real(0.5) * (sustr(i,j,0) + sustr(i+1,j,0))) *
                  (Real(0.5) * (sustr(i,j,0) + sustr(i+1,j,0))) +
                  (Real(0.5) * (svstr(i,j,0) + svstr(i,j+1,0))) *
                  (Real(0.5) * (svstr(i,j,0) + svstr(i,j+1,0))));
}

/** \brief Gas transfer velocity times dt, for rate factor \p cff2 from
 *         transfer_rate, squared wind \p u10squ and Schmidt number \p Sc. */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real transfer_velocity (amrex::Real cff2, amrex::Real u10squ, amrex::Real Sc) noexcept
{
    return cff2 * u10squ * std::sqrt(amrex::Real(660.0) / Sc);
}

/** \brief O2 saturation [mmol/m3] from temperature [degC] and salinity, by the
 *         Garcia and Gordon (1992) formula converted from ml/l. */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real o2_saturation (amrex::Real T, amrex::Real S) noexcept
{
    using amrex::Real;
    constexpr Real OA0 = Real(2.00907);
    constexpr Real OA1 = Real(3.22014);
    constexpr Real OA2 = Real(4.05010);
    constexpr Real OA3 = Real(4.94457);
    constexpr Real OA4 = Real(-0.256847);
    constexpr Real OA5 = Real(3.88767);
    constexpr Real OB0 = Real(-0.00624523);
    constexpr Real OB1 = Real(-0.00737614);
    constexpr Real OB2 = Real(-0.0103410);
    constexpr Real OB3 = Real(-0.00817083);
    constexpr Real OC0 = Real(-0.000000488682);
    constexpr Real l2mol = Real(1000.0) / Real(22.3916); // liter to mol

    const Real TS = std::log((Real(298.15) - T) / (Real(273.15) + T));
    const Real AA = OA0 + TS * (OA1 + TS * (OA2 + TS * (OA3 + TS * (OA4 + TS * OA5)))) +
                    S * (OB0 + TS * (OB1 + TS * (OB2 + TS * OB3))) +
                    OC0 * S * S;
    return l2mol * std::exp(AA);
}

/** \brief CO2 solubility [mol/(kg atm)] from temperature [degC] and salinity,
 *         by the Weiss (1974) formula. */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real co2_solubility (amrex::Real T, amrex::Real S) noexcept
{
    using amrex::Real;
    constexpr Real A1 = Real(-60.2409);
    constexpr Real A2 = Real(93.4517);
    constexpr Real A3 = Real(23.3585);
    constexpr Real B1 = Real(0.023517);
    constexpr Real B2 = Real(-0.023656);
    constexpr Real B3 = Real(0.0047036);

    const Real TempK = Real(0.01) * (T + Real(273.15));
    return std::exp(A1 +
                    A2 / TempK +
                    A3 * std::log(TempK) +
                    S * (B1 + TempK * (B2 + B3 * TempK)));
}


/** \brief The air-sea exchange factors for one surface cell. */
struct SurfaceGasExchange {
    amrex::Real O2_transfer = amrex::Real(0.0);   //!< O2 transfer velocity times dt
    amrex::Real O2satu = amrex::Real(0.0);        //!< O2 saturation [mmol/m3]
    amrex::Real CO2_transfer = amrex::Real(0.0);  //!< CO2 transfer velocity times dt
    amrex::Real CO2_sol = amrex::Real(0.0);       //!< CO2 solubility [mol/(kg atm)]
};

/** \brief Air-sea exchange factors at rho point (i,j) from the surface
 *         temperature \p T [degC] and salinity \p S and the wind or stress.
 *
 * They depend only on those, so a kernel computes them once per column, not
 * once per BioIter iteration. The O2 factors are computed only with \p oxygen
 * and the CO2 ones only with \p carbon; the others stay zero, and with
 * neither the wind and stress views are not read.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
SurfaceGasExchange surface_gas_exchange (bool oxygen, bool carbon,
                                         SchmidtCoefficients const& o2_schmidt,
                                         SchmidtCoefficients const& co2_schmidt,
                                         amrex::Real dtdays, bool bulk_fluxes,
                                         amrex::Array4<amrex::Real const> const& uwind,
                                         amrex::Array4<amrex::Real const> const& vwind,
                                         amrex::Array4<amrex::Real const> const& sustr,
                                         amrex::Array4<amrex::Real const> const& svstr,
                                         int i, int j, amrex::Real rho0,
                                         amrex::Real T, amrex::Real S) noexcept
{
    SurfaceGasExchange gas;
    if (!oxygen && !carbon) {
        return gas;
    }
    // Transfer velocity from u10squared (u10 in m/s).
    const amrex::Real u10squ = u10_squared(bulk_fluxes, uwind, vwind, sustr, svstr, i, j, rho0);
    if (oxygen) {
        gas.O2_transfer = transfer_velocity(transfer_rate(o2_schmidt, dtdays), u10squ,
                                            schmidt_number(o2_schmidt, T));
        gas.O2satu = o2_saturation(T, S);
    }
    if (carbon) {
        gas.CO2_transfer = transfer_velocity(transfer_rate(co2_schmidt, dtdays), u10squ,
                                             schmidt_number(co2_schmidt, T));
        gas.CO2_sol = co2_solubility(T, S);
    }
    return gas;
}

}
#endif
