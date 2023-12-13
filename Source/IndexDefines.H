#ifndef ROMSX_INDEX_DEFINES_H_
#define ROMSX_INDEX_DEFINES_H_

#include <AMReX_REAL.H>
#include <AMReX_Arena.H>

// Cell-centered state variables
#define Temp_comp      0
#define Scalar_comp    1
#define Salt_comp      2
#define NCONS          3

#define NGROW          2

namespace BCVars {
    enum {
        cons_bc = 0,
        Temp_bc_comp = 0,
        Scalar_bc_comp,
        xvel_bc = NCONS,
        yvel_bc,
        zvel_bc,
        NumTypes
    };
}

enum struct ROMSX_BC {
    symmetry, inflow, outflow, no_slip_wall, slip_wall, periodic, MOST, undefined
};

// NOTE: the first of these must match up with the BCType enum
//     in amrex/Src/Base/AMReX_BC_TYPES.H.  We had extras at
//     the end to use locally
namespace ROMSXBCType {
enum mathematicalBndryTypes : int { bogus        = -666,
    reflect_odd  = -1,
    int_dir      =  0,
    reflect_even =  1,
    foextrap     =  2,
    ext_dir      =  3,
    MOST         =  101,
    ext_dir_ingested =  102
};
}
#endif
