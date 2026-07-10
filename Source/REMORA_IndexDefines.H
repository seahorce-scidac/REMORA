#ifndef REMORA_INDEX_DEFINES_H_
#define REMORA_INDEX_DEFINES_H_

#include <AMReX_REAL.H>
#include <AMReX_Arena.H>

// Cell-centered state variables
#define Temp_comp      0
#define Salt_comp      1
#define Tracer_comp    2

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
 * Note that the salt component must be one greater than the temp component for the
 * logic in fill_from_bdyfiles to work
 */
namespace BdyVars {
    enum {
       null  = -1,
       u     =  0,
       v     =  1,
       t     =  2,
       s     =  3,
       ubar      ,
       vbar      ,
       zeta      ,
       NumTypes
    };
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
