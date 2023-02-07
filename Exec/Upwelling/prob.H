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
  

  amrex::Real rho_0   = 0.0;
  amrex::Real Theta_0 = 1.0;
  amrex::Real A_0     = 1.0;
  amrex::Real B_0     = 0.0;
  amrex::Real u_0     = 0.0;
  amrex::Real v_0     = 0.0;
  amrex::Real rad_0   = 0.0;
  amrex::Real z0      = 0.1;  // Surface Roughness
  amrex::Real zRef    = 80.0; // Reference Height
  amrex::Real uRef    = 0.0;  // Reference Wind Speed

  amrex::Real xc_frac = 0.5; // Location of "center" of scalar (multiplies domain length)
  amrex::Real yc_frac = 0.5; // Location of "center" of scalar (multiplies domain length)
  amrex::Real zc_frac = 0.5; // Location of "center" of scalar (multiplies domain length)

  int prob_type = -1;
}; // namespace ProbParm

extern ProbParm parms;

#endif

