#ifndef REMORA_INDEX_DEFINES_H_
#define REMORA_INDEX_DEFINES_H_

#include <AMReX_REAL.H>
#include <AMReX_Arena.H>

// Cell-centered state variables
#define Temp_comp      0
#define Salt_comp      1
#define Scalar_comp    2
#define NCONS          3

#define NGROW          2

namespace BCVars {
    enum {
        cons_bc = 0,
        Temp_bc_comp = 0,
        Salt_bc_comp,
        Scalar_bc_comp,
        xvel_bc = NCONS,
        yvel_bc,
        zvel_bc,
        foextrap_bc,
        NumTypes
    };
}

// These are only for the variables read in from netcdf boundary condition files
// Note that the salt component must be one greater than the temp component for the
//      logic in fill_from_bdyfiles to work
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
    symmetry, inflow, outflow, no_slip_wall, slip_wall, periodic, undefined
};

// NOTE: the first of these must match up with the BCType enum
//     in amrex/Src/Base/AMReX_BC_TYPES.H.  We had extras at
//     the end to use locally
namespace REMORABCType {
enum mathematicalBndryTypes : int { bogus        = -666,
    reflect_odd  = -1,
    int_dir      =  0,
    reflect_even =  1,
    foextrap     =  2,
    ext_dir      =  3
};
}
#endif
