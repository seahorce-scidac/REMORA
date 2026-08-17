#ifndef REMORA_INDEX_DEFINES_H_
#define REMORA_INDEX_DEFINES_H_

#include <AMReX_REAL.H>
#include <AMReX_Arena.H>
#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>

// Cell-centered state variables
#define Temp_comp      0
#define Salt_comp      1
#define Tracer_comp    2

// Number of active tracers, i.e. the ones that feed back on the density field:
// temperature and salinity. ROMS carries vertical mixing coefficients for these two
// alone, so Tracer_comp doubles as the count of them.
#define NAT            Tracer_comp

/**
 * The vertical diffusivity component that cons component icomp is mixed with.
 *
 * Akt is stored for the active tracers only, and every passive tracer -- dye and biology
 * alike -- mixes with the salinity coefficient. This is ROMS's ltrc=MIN(itrc,NAT) in
 * pre_step3d and step3d_t, and it is why vec_Akt has NAT components rather than ncons.
 * Horizontal diffusivity is a separate story: tnu2 really is per tracer.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int akt_comp (int icomp) noexcept { return (icomp < NAT) ? icomp : Salt_comp; }

#define NGROW          2

// Number of 3D/2D variables output in netcdf files
#define NC3D           5
#define NC2D           5
// Number of 2D variables in a netcdf file's header
#define NCH2D          12

/**
 * \brief Boundary fill types
 */
namespace BCVars {
    static constexpr int cons_bc         = 0;
    static constexpr int Temp_bc_comp    = 0;
    static constexpr int Salt_bc_comp    = 1;
    static constexpr int Scalar_bc_comp  = 2;

    inline int xvel_bc (int ncons) noexcept { return ncons; }
    inline int yvel_bc (int ncons) noexcept { return ncons + 1; }
    inline int zvel_bc (int ncons) noexcept { return ncons + 2; }
    inline int ubar_bc (int ncons) noexcept { return ncons + 3; }
    inline int vbar_bc (int ncons) noexcept { return ncons + 4; }
    inline int zeta_bc (int ncons) noexcept { return ncons + 5; }
    inline int tke_bc  (int ncons) noexcept { return ncons + 6; }
    inline int foextrap_periodic_bc (int ncons) noexcept { return ncons + 7; }
    inline int foextrap_bc          (int ncons) noexcept { return ncons + 8; }
    inline int u2d_simple_bc        (int ncons) noexcept { return ncons + 9; }
    inline int v2d_simple_bc        (int ncons) noexcept { return ncons + 10; }
    inline int NumTypes             (int ncons) noexcept { return ncons + AMREX_SPACEDIM + 8; }
}

/**
 * \brief These are only for the variables read in from netcdf boundary condition files
 *
 * The cell-centered tracers occupy a contiguous block of slots starting at t, laid
 * out in the same order as the components of the cons MultiFab. So slot t is temp,
 * slot t+1 is salt, and slot t+n is cons component n for every additional scalar
 * (passive tracers, biology tracers). Use cons(icomp) rather than hardcoding this,
 * and note that the logic in fill_from_bdyfiles relies on the block being contiguous
 * and in component order.
 *
 * Like BCVars, everything above the tracer block is a function of ncons rather than a
 * compile-time constant. With ncons == 2 these reduce to the historical fixed values
 * (ubar 4, vbar 5, zeta 6, NumTypes 7).
 */
namespace BdyVars {
    static constexpr int null =  -1;
    static constexpr int u    =   0;
    static constexpr int v    =   1;
    static constexpr int t    =   2; ///< cons component Temp_comp
    static constexpr int s    =   3; ///< cons component Salt_comp

    inline int cons (int icomp) noexcept { return t + icomp; }
    inline int ubar (int ncons) noexcept { return t + ncons;     }
    inline int vbar (int ncons) noexcept { return t + ncons + 1; }
    inline int zeta (int ncons) noexcept { return t + ncons + 2; }
    inline int NumTypes (int ncons) noexcept { return t + ncons + 3; }
}

/**
 * \brief Index constants for driver_atmos_state_from_driver and the
 *        matching atmospheric forcing MultiFab vector passed to
 *        REMORA::ApplyAtmosphericStates.
 *
 * Order mirrors AtmosToOceanIndex (state lanes 0-8) from
 * ERFRemoraCouplingContract.H; Pair uses REMORA's internal naming
 * convention rather than the contract's Patm.
 */
namespace AtmosState {
    enum {
        Uwind = 0, ///< 10-m zonal wind [m/s]
        Vwind,     ///< 10-m meridional wind [m/s]
        Pair,      ///< atmospheric pressure [Pa from driver, mb in REMORA]
        Qair,      ///< specific humidity [kg/kg]
        Tair,      ///< air temperature [K from driver, degC in REMORA]
        Cloud,     ///< cloud fraction [0-1]
        Rain,      ///< precipitation rate [kg/m^2/s]
        SWrad,     ///< downward shortwave radiation [W/m^2]
        LWrad,     ///< downward longwave radiation [W/m^2]
        NumTypes
    };
}

/**
 * \brief Index constants for the driver ATM->OCN flux-passing active view.
 *
 * Order mirrors AtmosToOceanActiveIndices(FluxPassing) from
 * ERFRemoraCouplingContract.H:
 * tau_x, tau_y, SHflux, LHflux, SWrad, LWrad, rain, evap.
 */
namespace AtmosFluxes {
    enum {
        TauX = 0,  ///< surface zonal stress lane
        TauY,      ///< surface meridional stress lane
        SHflux,    ///< sensible heat flux lane
        LHflux,    ///< latent heat flux lane
        SWrad,     ///< shortwave radiation lane
        LWrad,     ///< longwave flux lane
        Rain,      ///< precipitation rate lane
        Evap,      ///< evaporation rate lane
        NumTypes
    };
}

/**
 * \brief Source type for bulk-flux atmospheric forcing variables.
 *
 * computed is used for diagnostics/parameterizations that are evaluated by
 * REMORA rather than prescribed as an input field, e.g. Berliand longwave or
 * bulk evaporation minus precipitation.
 */
enum class BulkForcingType {
    constant,
    analytic,
    netcdf,
    computed
};

/**
 * \brief Index constants for SolverChoice::bulk_flux_type.
 */
namespace BulkFlux {
    enum {
        Uwind = 0,  ///< 10-m zonal wind [m/s]
        Vwind,     ///< 10-m meridional wind [m/s]
        Pair,      ///< atmospheric pressure [mb]
        Qair,      ///< specific humidity or relative humidity [kg/kg or fraction]
        Tair,      ///< air temperature [degC]
        Cloud,     ///< cloud fraction [0-1]
        Rain,      ///< precipitation rate [kg/m^2/s]
        SWrad,     ///< downward shortwave radiation [W/m^2]
        LWrad,     ///< longwave radiation [W/m^2]
        EminusP,   ///< evaporation minus precipitation [m/s]
        NumTypes
    };
}

enum struct REMORA_BC {
    symmetry, inflow, outflow, no_slip_wall, slip_wall, periodic,
    clamped, chapman, flather, orlanski_rad, orlanski_rad_nudge, undefined
};

/**
 * NOTE: the first of these must match up with the BCType enum
 * in amrex/Src/Base/AMReX_BC_TYPES.H.  We had extras at
 * the end to use locally
 */
namespace REMORABCType {
enum mathematicalBndryTypes : int { bogus        = -666,
    reflect_odd        = -1,
    int_dir            =  0,
    reflect_even       =  1,
    foextrap           =  2,
    ext_dir            =  3,
    clamped            =  4,
    chapman            =  5,
    flather            =  6,
    orlanski_rad       =  7,
    orlanski_rad_nudge =  8
};
}
#endif
