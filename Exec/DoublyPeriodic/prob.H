#ifndef _PROB_H_
#define _PROB_H_

#include "AMReX_REAL.H"

struct ProbParm {

  amrex::Real u_0     = 0.0;
  amrex::Real v_0     = 0.0;
  amrex::Real z0      = 0.1;  // Surface Roughness
  amrex::Real zRef    = 80.0; // Reference Height
  amrex::Real uRef    = 0.0;  // Reference Wind Speed

}; // namespace ProbParm

extern ProbParm parms;

#endif

