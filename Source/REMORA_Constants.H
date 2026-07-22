#ifndef REMORA_CONSTANTS_H_
#define REMORA_CONSTANTS_H_

#include <AMReX_REAL.H>

constexpr amrex::Real zero         = amrex::Real(0.0);
constexpr amrex::Real one          = amrex::Real(1.0);
constexpr amrex::Real two          = amrex::Real(2.0);
constexpr amrex::Real three        = amrex::Real(3.0);
constexpr amrex::Real four         = amrex::Real(4.0);
constexpr amrex::Real half         = amrex::Real(0.5);
constexpr amrex::Real fourth       = amrex::Real(0.25);
constexpr amrex::Real third        = amrex::Real(1.0/3.0);
constexpr amrex::Real eightthirds  = amrex::Real(8.0/3.0);

//Note that amrex source code uses fewer digits of pi:
//     constexpr Real pi = Real(3.1415926535897932);
//
constexpr amrex::Real PI = amrex::Real(3.14159265358979323846264338327950288);

// Safe overly large value to fill data so that it is not uninitialized
//      but we can easily identify if the value is used in a computation
// Would like to use the following but std::sqrt will not be allowable
//       in constexpr until C++26
// const amrex::Real bogus_large_value = std::sqrt(std::numeric_limits<amrex::Real>::max());
#ifdef AMREX_USE_FLOAT
constexpr amrex::Real bogus_large_value = amrex::Real(1.e18);
#else
constexpr amrex::Real bogus_large_value = amrex::Real(1.e150);
#endif

constexpr amrex::Real Cp           = amrex::Real(3985.0);              // Specific heat of seawater [Joules/kg/degC]
constexpr amrex::Real Csolar       = amrex::Real(1353.0);              // Solar irradiation constant [1360-1380 W/m2]
constexpr amrex::Real Eradius      = amrex::Real(6371315.0);           // Earth equatorial radius [m]
constexpr amrex::Real StefBo       = amrex::Real(5.670374419e-8);      // Stefan-Boltzmann constant [Watts/m2/K4]
constexpr amrex::Real emmiss       = amrex::Real(0.97);                // Infrared emissivity [non_dimensional]
constexpr amrex::Real rhow         = amrex::Real(1000.0);              // fresh water density [kg/m3]
constexpr amrex::Real g            = amrex::Real(9.80665);             // acceleration due to gravity [m/s2]
constexpr amrex::Real vonKar       = amrex::Real(0.41);                // von Karman constant

//-----------------------------------------------------------------------
// Constants used in surface fluxes bulk parameterization.
//-----------------------------------------------------------------------
//
//    blk_Cpa       Specific heat capacity for dry air (J/kg/K).
//    blk_Cpw       Specific heat capacity for seawater (J/kg/K).
//    blk_Rgas      Gas constant for dry air (J/kg/K).
//    blk_Zabl      Height (m) of atmospheric boundary layer.
//    blk_ZQ        Height (m) of surface air humidity measurement.
//    blk_ZT        Height (m) of surface air temperature measurement.
//    blk_ZW        Height (m) of surface winds measurement.
//    blk_beta      Beta parameter evaluated from Fairall low windspeed
//                    turbulence data.
//    blk_dter      Temperature change.
//    blk_tcw       Thermal conductivity of water (W/m/K).
//    blk_visw      Kinematic viscosity water (m2/s).

constexpr amrex::Real blk_Cpa  = amrex::Real(1004.67);    // (J/kg/K), Businger 1982
constexpr amrex::Real blk_Cpw  = amrex::Real(4000.0);     // (J/kg/K)
constexpr amrex::Real blk_Rgas = amrex::Real(287.1);      // (J/kg/K)
constexpr amrex::Real blk_Zabl = amrex::Real(600.0);      // (m)
constexpr amrex::Real blk_beta = amrex::Real(1.2);        // non-dimensional
constexpr amrex::Real blk_dter = amrex::Real(0.3);        // (K)
constexpr amrex::Real blk_tcw  = amrex::Real(0.6);        // (W/m/K)
constexpr amrex::Real blk_visw = amrex::Real(0.000001);   // (m2/s)


//  Polynomial  expansion  coefficients for the computation of
//  "in situ" density  and other associated quantities via the
//  nonlinear equation of state for seawater  as a function of
//  potential temperature, salinity, and pressure (Jackett and
//  McDougall, 1992).

constexpr amrex::Real A00 = amrex::Real(1.909256e+04);
constexpr amrex::Real A01 = amrex::Real(2.098925e+02);
constexpr amrex::Real A02 = amrex::Real(-3.041638e+00);
constexpr amrex::Real A03 = amrex::Real(-1.852732e-03);
constexpr amrex::Real A04 = amrex::Real(-1.361629e-05);
constexpr amrex::Real B00 = amrex::Real(1.044077e+02);
constexpr amrex::Real B01 = amrex::Real(-6.500517e+00);
constexpr amrex::Real B02 = amrex::Real(1.553190e-01);
constexpr amrex::Real B03 = amrex::Real(2.326469e-04);
constexpr amrex::Real D00 = amrex::Real(-5.587545e+00);
constexpr amrex::Real D01 = amrex::Real(7.390729e-01);
constexpr amrex::Real D02 = amrex::Real(-1.909078e-02);
constexpr amrex::Real E00 = amrex::Real(4.721788e-01);
constexpr amrex::Real E01 = amrex::Real(1.028859e-02);
constexpr amrex::Real E02 = amrex::Real(-2.512549e-04);
constexpr amrex::Real E03 = amrex::Real(-5.939910e-07);
constexpr amrex::Real F00 = amrex::Real(-1.571896e-02);
constexpr amrex::Real F01 = amrex::Real(-2.598241e-04);
constexpr amrex::Real F02 = amrex::Real(7.267926e-06);
constexpr amrex::Real G00 = amrex::Real(2.042967e-03);
constexpr amrex::Real G01 = amrex::Real(1.045941e-05);
constexpr amrex::Real G02 = amrex::Real(-5.782165e-10);
constexpr amrex::Real G03 = amrex::Real(1.296821e-07);
constexpr amrex::Real H00 = amrex::Real(-2.595994e-07);
constexpr amrex::Real H01 = amrex::Real(-1.248266e-09);
constexpr amrex::Real H02 = amrex::Real(-3.508914e-09);
constexpr amrex::Real Q00 = amrex::Real(9.99842594e+02);
constexpr amrex::Real Q01 = amrex::Real(6.793952e-02);
constexpr amrex::Real Q02 = amrex::Real(-9.095290e-03);
constexpr amrex::Real Q03 = amrex::Real(1.001685e-04);
constexpr amrex::Real Q04 = amrex::Real(-1.120083e-06);
constexpr amrex::Real Q05 = amrex::Real(6.536332e-09);
constexpr amrex::Real U00 = amrex::Real(8.24493e-01);
constexpr amrex::Real U01 = amrex::Real(-4.08990e-03);
constexpr amrex::Real U02 = amrex::Real(7.64380e-05);
constexpr amrex::Real U03 = amrex::Real(-8.24670e-07);
constexpr amrex::Real U04 = amrex::Real(5.38750e-09);
constexpr amrex::Real V00 = amrex::Real(-5.72466e-03);
constexpr amrex::Real V01 = amrex::Real(1.02270e-04);
constexpr amrex::Real V02 = amrex::Real(-1.65460e-06);
constexpr amrex::Real W00 = amrex::Real(4.8314e-04);

#endif
