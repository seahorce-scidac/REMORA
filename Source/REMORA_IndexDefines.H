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
