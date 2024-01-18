#ifndef REMORA_CONSTANTS_H_
#define REMORA_CONSTANTS_H_

#include <AMReX_REAL.H>

constexpr amrex::Real PI = amrex::Real(3.14159265358979323846264338327950288);
//Note that amrex source code uses fewer digits of pi:
//     constexpr Real pi = Real(3.1415926535897932);
//TODO: Change these types of macros to 'const'
constexpr amrex::Real R_d      = amrex::Real(287.0);   // gas constant for dry air [J/(kg-K)]
constexpr amrex::Real R_v      = amrex::Real(461.6);   // gas constant for dry air [J/(kg-K)]
constexpr amrex::Real c_p      = amrex::Real(1004.5);   // specific heat at constant pressure for dry air [J/(kg-K)]
constexpr amrex::Real p_0        = amrex::Real(1.0e5); // reference surface pressure [Pa]
constexpr amrex::Real Gamma      = amrex::Real(1.4);   // C_p/C_v [-]
constexpr amrex::Real KAPPA      = amrex::Real(0.41);  // von Karman constant

// Derived Constants
constexpr amrex::Real ip_0    = amrex::Real(1./p_0);
constexpr amrex::Real iR_d    = amrex::Real(1./R_d);
constexpr amrex::Real iGamma  = amrex::Real(1./Gamma);
constexpr amrex::Real rdOcp   = amrex::Real(R_d/c_p);

#endif
