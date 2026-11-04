#ifndef REMORA_BIOLOGY_CARBONATE_H_
#define REMORA_BIOLOGY_CARBONATE_H_

#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>
#include <cmath>

namespace REMORABiology {

using amrex::Real;

/**
 * Surface-water pCO2 [ppmv] from temperature [degC], salinity, total
 * inorganic carbon and total alkalinity [mmol/m3], by bracketing and
 * bisection for [H+]. This is ROMS pCO2_water with phosphate and silicate
 * set to zero, shared by the Fennel and ECB kernels, whose ROMS sources carry
 * identical copies.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
pco2_water(
  amrex::Real T, amrex::Real S, amrex::Real TIC, amrex::Real TAlk) noexcept
{
    constexpr Real zero = Real(0.0);
    constexpr int IbrackMax = 30;

    // Determine coefficients for surface carbon chemistry.
    const Real Tk = T + Real(273.15);
    const Real centiTk = Real(0.01) * Tk;
    const Real invTk = Real(1.0) / Tk;
    const Real logTk = std::log(Tk);
    const Real sqrtS = std::sqrt(S);
    const Real SO4 = Real(19.924) * S / (Real(1000.0) - Real(1.005) * S);
    const Real sqrtSO4 = std::sqrt(SO4);
    const Real scl = S / Real(1.80655);

    const Real alk = TAlk * Real(0.000001);
    const Real dic = TIC * Real(0.000001);
    const Real phos = zero;
    const Real sili = zero;

    // Correction term for non-ideality, ff=k0*(1-pH2O). Equation 13 with
    // table 6 values from Weiss and Price (1980, Mar. Chem., 8, 347-359).
    const Real ff = std::exp(
      -Real(162.8301) + Real(218.2968) / centiTk +
      std::log(centiTk) * Real(90.9241) - centiTk * centiTk * Real(1.47696) +
      S * (Real(0.025695) -
           centiTk * (Real(0.025225) - centiTk * Real(0.0049867))));

    // Compute first (K1) and second (K2) dissociation constant of carboinic
    // acid:
    //          K1 = [H][HCO3]/[H2CO3]
    //          K2 = [H][CO3]/[HCO3]
    // From Millero (1995; page 664) using Mehrbach et al. (1973) data on
    // seawater scale.
    const Real K1 = std::pow(
      Real(10.0), Real(62.008) - invTk * Real(3670.7) - logTk * Real(9.7944) +
                    S * (Real(0.0118) - S * Real(0.000116)));
    const Real K2 = std::pow(
      Real(10.0), -Real(4.777) - invTk * Real(1394.7) +
                    S * (Real(0.0184) - S * Real(0.000118)));

    // Compute dissociation constant of boric acid, Kb=[H][BO2]/[HBO2].
    // From Millero (1995; page 669) using data from Dickson (1990).
    const Real Kb = std::exp(
      -invTk *
        (Real(8966.90) +
         sqrtS * (Real(2890.53) +
                  sqrtS * (Real(77.942) -
                           sqrtS * (Real(1.728) - sqrtS * Real(0.0996))))) -
      logTk * (Real(24.4344) + sqrtS * (Real(25.085) + sqrtS * Real(0.2474))) +
      Tk * (sqrtS * Real(0.053105)) + Real(148.0248) +
      sqrtS * (Real(137.1942) + sqrtS * Real(1.62142)));

    // Compute first (K1p), second (K2p), and third (K3p) dissociation
    // constant of phosphoric acid:
    //          K1p = [H][H2PO4]/[H3PO4]
    //          K2p = [H][HPO4]/[H2PO4]
    //          K3p = [H][PO4]/[HPO4]
    // From DOE (1994) equations 7.2.20, 7.2.23, and 7.2.26, respectively.
    // With footnote using data from Millero (1974).
    const Real K1p = std::exp(
      Real(115.525) - invTk * Real(4576.752) - logTk * Real(18.453) +
      sqrtS * (Real(0.69171) - invTk * Real(106.736)) -
      S * (Real(0.01844) + invTk * Real(0.65643)));
    const Real K2p = std::exp(
      Real(172.0883) - invTk * Real(8814.715) - logTk * Real(27.927) +
      sqrtS * (Real(1.3566) - invTk * Real(160.340)) -
      S * (Real(0.05778) - invTk * Real(0.37335)));
    const Real K3p = std::exp(
      -Real(18.141) - invTk * Real(3070.75) +
      sqrtS * (Real(2.81197) + invTk * Real(17.27039)) -
      S * (Real(0.09984) + invTk * Real(44.99486)));

    // Compute dissociation constant of silica, Ksi=[H][SiO(OH)3]/[Si(OH)4].
    // From Millero (1995; page 671) using data from Yao and Millero (1995).
    const Real Ksi = std::exp(
      Real(117.385) - invTk * Real(8904.2) - logTk * Real(19.334) +
      sqrtSO4 * (Real(3.5913) - invTk * Real(458.79)) -
      SO4 * (Real(1.5998) - invTk * Real(188.74) -
             SO4 * (Real(0.07871) - invTk * Real(12.1652))) +
      std::log(Real(1.0) - Real(0.001005) * S));

    // Compute ion product of whater, Kw = [H][OH].
    // From Millero (1995; page 670) using composite data.
    const Real Kw = std::exp(
      Real(148.9652) - invTk * Real(13847.26) - logTk * Real(23.6521) -
      sqrtS * (Real(5.977) - invTk * Real(118.67) - logTk * Real(1.0495)) -
      S * Real(0.01615));

    // Compute salinity constant of hydrogen sulfate, Ks = [H][SO4]/[HSO4].
    // From Dickson (1990, J. chem. Thermodynamics 22, 113).
    const Real Ks = std::exp(
      Real(141.328) - invTk * Real(4276.1) - logTk * Real(23.093) +
      sqrtSO4 * (Real(324.57) - invTk * Real(13856.0) - logTk * Real(47.986) -
                 SO4 * invTk * Real(2698.0)) -
      SO4 * (Real(771.54) - invTk * Real(35474.0) - logTk * Real(114.723) -
             SO4 * invTk * Real(1776.0)) +
      std::log(Real(1.0) - Real(0.001005) * S));

    // Compute stability constant of hydrogen fluorid, Kf = [H][F]/[HF].
    // From Dickson and Riley (1979) -- change pH scale to total.
    const Real Kf = std::exp(
      -Real(12.641) + invTk * Real(1590.2) + sqrtSO4 * Real(1.525) +
      std::log(Real(1.0) - Real(0.001005) * S) +
      std::log(Real(1.0) + Real(0.1400) * scl / (Real(96.062) * Ks)));

    // Calculate concentrations for borate (Uppstrom, 1974), sulfate (Morris
    // and Riley, 1966), and fluoride (Riley, 1965).
    const Real borate = Real(0.000232) * scl / Real(10.811);
    const Real sulfate = Real(0.14) * scl / Real(96.062);
    const Real fluoride = Real(0.000067) * scl / Real(18.9984);

    // Bracket and bisection method.
    Real X_lo = std::pow(Real(10.0), -Real(10.0));
    Real X_hi = std::pow(Real(10.0), -Real(5.0));
    Real X_mid = Real(0.5) * (X_lo + X_hi);
    Real X = X_mid;
    const Real K12 = K1 * K2;
    const Real K12p = K1p * K2p;
    const Real K123p = K12p * K3p;
    const Real invKb = Real(1.0) / Kb;
    const Real invKs = Real(1.0) / Ks;
    const Real invKsi = Real(1.0) / Ksi;
    Real fni1 = zero;
    Real fni3 = zero;

    for (int Ibrack = 0; Ibrack < IbrackMax; ++Ibrack) {
        for (int Hstep = 1; Hstep <= 3; ++Hstep) {
            if (Hstep == 1) {
                X = X_hi;
            }
            if (Hstep == 2) {
                X = X_lo;
            }
            if (Hstep == 3) {
                X = X_mid;
            }

            // Set some common combinations of parameters used in the iterative
            // [H+] solver.
            const Real X2 = X * X;
            const Real X3 = X2 * X;
            const Real invX = Real(1.0) / X;
            const Real A = X * (K12p + X * (K1p + X)) + K123p;
            const Real B = X * (K1 + X) + K12;
            const Real C = Real(1.0) / (Real(1.0) + sulfate * invKs);

            // Evaluate f([H+]) for bracketing and mid-value cases.
            const Real fni = dic * (K1 * X + Real(2.0) * K12) / B +
                             borate / (Real(1.0) + X * invKb) + Kw * invX +
                             phos * (K12p * X + Real(2.0) * K123p - X3) / A +
                             sili / (Real(1.0) + X * invKsi) - X * C -
                             sulfate / (Real(1.0) + Ks * invX * C) -
                             fluoride / (Real(1.0) + Kf * invX) - alk;
            if (Hstep == 1) {
                fni1 = fni;
            }
            if (Hstep == 3) {
                fni3 = fni;
            }
        }

        // Now, bracket solution within two of three.
        if (fni3 == zero) {
            break;
        }
        const Real ftest = fni1 / fni3;
        if (ftest > zero) {
            X_hi = X_mid;
        } else {
            X_lo = X_mid;
        }
        X_mid = Real(0.5) * (X_lo + X_hi);
    }

    // Last iteration gives value.
    X = X_mid;

    // Determine pCO2. Total Hydrogen ion concentration, Htotal = [H+].
    const Real Htotal = X;
    const Real Htotal2 = Htotal * Htotal;

    // Calculate [CO2*] (mole/m3) as defined in DOE Methods Handbook 1994
    // Version 2, ORNL/CDIAC-74, Dickson and Goyet, Eds. (Chapter 2,
    // page 10, Eq A.49).
    const Real CO2star = dic * Htotal2 / (Htotal2 + K1 * Htotal + K1 * K2);

    // Add output argument for storing pCO2surf.
    return CO2star * Real(1000000.0) / ff;
}

} // namespace REMORABiology
#endif
