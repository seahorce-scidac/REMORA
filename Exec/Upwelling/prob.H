#ifndef _PROB_H_
#define _PROB_H_

#include "AMReX_REAL.H"

struct ProbParm {
  // Linear equation of State parameters
  amrex::Real R0 = 1028; // background density value (Kg/m3) used in Linear Equation of State
  amrex::Real S0 = 35.0; // background salinity (nondimensional) constant
  amrex::Real T0 = 5.0; // background potential temperature (Celsius) constant

  /*
      // Background salinity (nondimensional) constant.
    amrex::Real S0;

    // Background Temperature (Celsius) constant.
    amrex::Real T0;
  */

    // mod_scalars values

    // Mean density (Kg/m3) used when the Boussinesq approximation is inferred.
    amrex::Real rho0 = 1025.0;
  
    // physical constant for gravity
#ifdef SOLITON
    amrex::Real g = 1.0;                // non-dimensional
#elif defined WBC_1 || defined WBC_2 || defined WBC_3
    amrex::Real g = 9.8;                // m/s2
#elif defined CIRCLE
    amrex::Real g = 3.92e-2;            // m/s2
#else
    amrex::Real g = 9.81;               // m/s2
#endif
  

  /*
  amrex::Real U_0 = 0.0;
  amrex::Real V_0 = 0.0;
  amrex::Real x_c = 0.0; // center of thermal perturbation
  amrex::Real z_c = 3000.0;
  amrex::Real x_r = 4000.0;
  amrex::Real z_r = 2000.0;
  amrex::Real T_pert = -15.0; // perturbation temperature
  // overridden physical constants
  amrex::Real C_p = 1004.0;
  */
}; // namespace ProbParm

extern ProbParm parms;

#endif

