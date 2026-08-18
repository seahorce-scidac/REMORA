#include <REMORA.H>

#include <algorithm>
#include <cmath>
#include <string>

#include <AMReX.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Utility.H>

#include <REMORA_DateClock.H>

using namespace amrex;


namespace {

#ifdef REMORA_USE_BIOLOGY_DIAG
/**
 * Emit one Path B diagnostic record.
 *
 * Format is contractual and must stay byte-identical to fennel_tag_line in
 * Source/Biology/Fortran/REMORA_fennel_roms.F, whose Fortran edit descriptor
 * is ('FENNEL-FORT ',a10,' iter=',i3,' i=',i5,' j=',i5,' k=',i4,' ',a4,' ',
 * ES26.17E2). Both paths print k in REMORA's 0-based convention and iter
 * 1-based, so the two logs diff without any index map.
 *
 * Compiled only under REMORA_USE_BIOLOGY_DIAG (GNUmake USE_BIOLOGY_DIAG=TRUE,
 * CMake REMORA_ENABLE_BIOLOGY_DIAG=ON). This is called from inside the device
 * lambda, so in a default build the printf and its device-side runtime are
 * absent rather than merely unreached.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void remora_fennel_tag (const char* tag, int iter, int i, int j, int k,
                        const char* var, Real val) noexcept
{
    printf("FENNEL-CPP  %-10s iter=%3d i=%5d j=%5d k=%4d %-4s %26.17E\n",
           tag, iter, i, j, k, var, double(val));
}
#endif

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real fennel_pco2_water (Real T, Real S, Real TIC, Real TAlk) noexcept
{
    constexpr Real zero = Real(0.0);
    constexpr Real half = Real(0.5);
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
    const Real ff = std::exp(-Real(162.8301) +
                             Real(218.2968) / centiTk +
                             std::log(centiTk) * Real(90.9241) -
                             centiTk * centiTk * Real(1.47696) +
                             S * (Real(0.025695) -
                                  centiTk * (Real(0.025225) -
                                             centiTk * Real(0.0049867))));

    // Compute first (K1) and second (K2) dissociation constant of carboinic
    // acid:
    //          K1 = [H][HCO3]/[H2CO3]
    //          K2 = [H][CO3]/[HCO3]
    // From Millero (1995; page 664) using Mehrbach et al. (1973) data on
    // seawater scale.
    const Real K1 = std::pow(Real(10.0), Real(62.008) -
                             invTk * Real(3670.7) -
                             logTk * Real(9.7944) +
                             S * (Real(0.0118) - S * Real(0.000116)));
    const Real K2 = std::pow(Real(10.0), -Real(4.777) -
                             invTk * Real(1394.7) +
                             S * (Real(0.0184) - S * Real(0.000118)));

    // Compute dissociation constant of boric acid, Kb=[H][BO2]/[HBO2].
    // From Millero (1995; page 669) using data from Dickson (1990).
    const Real Kb = std::exp(-invTk * (Real(8966.90) +
                                       sqrtS * (Real(2890.53) +
                                                sqrtS * (Real(77.942) -
                                                         sqrtS * (Real(1.728) -
                                                                  sqrtS * Real(0.0996))))) -
                             logTk * (Real(24.4344) +
                                      sqrtS * (Real(25.085) + sqrtS * Real(0.2474))) +
                             Tk * (sqrtS * Real(0.053105)) +
                             Real(148.0248) +
                             sqrtS * (Real(137.1942) + sqrtS * Real(1.62142)));

    // Compute first (K1p), second (K2p), and third (K3p) dissociation
    // constant of phosphoric acid:
    //          K1p = [H][H2PO4]/[H3PO4]
    //          K2p = [H][HPO4]/[H2PO4]
    //          K3p = [H][PO4]/[HPO4]
    // From DOE (1994) equations 7.2.20, 7.2.23, and 7.2.26, respectively.
    // With footnote using data from Millero (1974).
    const Real K1p = std::exp(Real(115.525) -
                              invTk * Real(4576.752) -
                              logTk * Real(18.453) +
                              sqrtS * (Real(0.69171) - invTk * Real(106.736)) -
                              S * (Real(0.01844) + invTk * Real(0.65643)));
    const Real K2p = std::exp(Real(172.0883) -
                              invTk * Real(8814.715) -
                              logTk * Real(27.927) +
                              sqrtS * (Real(1.3566) - invTk * Real(160.340)) -
                              S * (Real(0.05778) - invTk * Real(0.37335)));
    const Real K3p = std::exp(-Real(18.141) -
                              invTk * Real(3070.75) +
                              sqrtS * (Real(2.81197) + invTk * Real(17.27039)) -
                              S * (Real(0.09984) + invTk * Real(44.99486)));

    // Compute dissociation constant of silica, Ksi=[H][SiO(OH)3]/[Si(OH)4].
    // From Millero (1995; page 671) using data from Yao and Millero (1995).
    const Real Ksi = std::exp(Real(117.385) -
                              invTk * Real(8904.2) -
                              logTk * Real(19.334) +
                              sqrtSO4 * (Real(3.5913) - invTk * Real(458.79)) -
                              SO4 * (Real(1.5998) - invTk * Real(188.74) -
                                     SO4 * (Real(0.07871) - invTk * Real(12.1652))) +
                              std::log(Real(1.0) - Real(0.001005) * S));

    // Compute ion product of whater, Kw = [H][OH].
    // From Millero (1995; page 670) using composite data.
    const Real Kw = std::exp(Real(148.9652) -
                             invTk * Real(13847.26) -
                             logTk * Real(23.6521) -
                             sqrtS * (Real(5.977) - invTk * Real(118.67) -
                                      logTk * Real(1.0495)) -
                             S * Real(0.01615));

    // Compute salinity constant of hydrogen sulfate, Ks = [H][SO4]/[HSO4].
    // From Dickson (1990, J. chem. Thermodynamics 22, 113).
    const Real Ks = std::exp(Real(141.328) -
                             invTk * Real(4276.1) -
                             logTk * Real(23.093) +
                             sqrtSO4 * (Real(324.57) - invTk * Real(13856.0) -
                                        logTk * Real(47.986) - SO4 * invTk * Real(2698.0)) -
                             SO4 * (Real(771.54) - invTk * Real(35474.0) -
                                    logTk * Real(114.723) - SO4 * invTk * Real(1776.0)) +
                             std::log(Real(1.0) - Real(0.001005) * S));

    // Compute stability constant of hydrogen fluorid, Kf = [H][F]/[HF].
    // From Dickson and Riley (1979) -- change pH scale to total.
    const Real Kf = std::exp(-Real(12.641) +
                             invTk * Real(1590.2) +
                             sqrtSO4 * Real(1.525) +
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
    Real X_mid = half * (X_lo + X_hi);
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
            if (Hstep == 1) { X = X_hi; }
            if (Hstep == 2) { X = X_lo; }
            if (Hstep == 3) { X = X_mid; }

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
                             borate / (Real(1.0) + X * invKb) +
                             Kw * invX +
                             phos * (K12p * X + Real(2.0) * K123p - X3) / A +
                             sili / (Real(1.0) + X * invKsi) -
                             X * C -
                             sulfate / (Real(1.0) + Ks * invX * C) -
                             fluoride / (Real(1.0) + Kf * invX) - alk;
            if (Hstep == 1) { fni1 = fni; }
            if (Hstep == 3) { fni3 = fni; }
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
        X_mid = half * (X_lo + X_hi);
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

/**
 * Atmospheric pCO2 (ppmv) for the surface CO2 gas exchange.
 *
 * \p time_seconds is the current model time and \p time_ref the reference date
 * and calendar it is measured against, both handed to remora_caldate. The two
 * time-dependent forms are ROMS's PCO2AIR_DATA and PCO2AIR_SECULAR; neither
 * varies in space, so this is evaluated once per call on the host and the
 * result handed to the kernel.
 */
Real
fennel_pco2_air (REMORABiology::PCO2AirType type, Real time_ref, Real time_seconds,
                 Real pco2air_constant) noexcept
{
    if (type == REMORABiology::PCO2AirType::constant) {
        return pco2air_constant;
    }

    constexpr Real pi2 = Real(6.2831853071796);

    int year = 0;
    Real yday = Real(0.0);
    remora_caldate(time_ref, time_seconds / Real(86400.0), year, yday);

    if (type == REMORABiology::PCO2AirType::data) {
        // Annual climatology of Laurent et al. (2017).
        return Real(380.464) + Real(9.321) *
               std::sin(pi2 * yday / Real(365.25) + Real(1.068));
    }

    // Secular trend. ROMS names this pmonth but computes years since 1951 and
    // multiplies by 12 in the linear term; keep both, so the coefficients are
    // the published ones.
    constexpr Real D0 = Real(282.6);
    constexpr Real D1 = Real(0.125);
    constexpr Real D2 = Real(-7.18);
    constexpr Real D3 = Real(0.86);
    constexpr Real D4 = Real(-0.99);
    constexpr Real D5 = Real(0.28);
    constexpr Real D6 = Real(-0.80);
    constexpr Real D7 = Real(0.06);

    const Real pmonth = Real(year) - Real(1951.0) + yday / Real(365.0);
    return D0 + D1 * pmonth * Real(12.0) +
           D2 * std::sin(pi2 * pmonth + D3) +
           D4 * std::sin(pi2 * pmonth + D5) +
           D6 * std::sin(pi2 * pmonth + D7);
}

}

namespace REMORABiology {

void
FennelParameters::init_params (const std::string& remora_prefix)
{
    ParmParse pp(remora_prefix + ".fennel");

    pp.queryAdd("BioIter", BioIter);
    pp.queryAdd("AttSW", AttSW);
    pp.queryAdd("AttChl", AttChl);
    pp.queryAdd("PARfrac", PARfrac);
    pp.queryAdd("Vp0", Vp0);
    pp.queryAdd("I_thNH4", I_thNH4);
    pp.queryAdd("D_p5NH4", D_p5NH4);
    pp.queryAdd("NitriR", NitriR);
    pp.queryAdd("K_NO3", K_NO3);
    pp.queryAdd("K_NH4", K_NH4);
    pp.queryAdd("K_PO4", K_PO4);
    pp.queryAdd("K_Phy", K_Phy);
    pp.queryAdd("Chl2C_m", Chl2C_m);
    pp.queryAdd("ChlMin", ChlMin);
    pp.queryAdd("PhyCN", PhyCN);
    pp.queryAdd("R_P2N", R_P2N);
    pp.queryAdd("PhyIP", PhyIP);
    pp.queryAdd("PhyIS", PhyIS);
    pp.queryAdd("PhyMin", PhyMin);
    pp.queryAdd("PhyMR", PhyMR);
    pp.queryAdd("ZooAE_N", ZooAE_N);
    pp.queryAdd("ZooCN", ZooCN);
    pp.queryAdd("ZooBM", ZooBM);
    pp.queryAdd("ZooER", ZooER);
    pp.queryAdd("ZooGR", ZooGR);
    pp.queryAdd("ZooMin", ZooMin);
    pp.queryAdd("ZooMR", ZooMR);
    pp.queryAdd("LDeRRN", LDeRRN);
    pp.queryAdd("LDeRRC", LDeRRC);
    pp.queryAdd("CoagR", CoagR);
    pp.queryAdd("SDeRRN", SDeRRN);
    pp.queryAdd("SDeRRC", SDeRRC);
    pp.queryAdd("RDeRRN", RDeRRN);
    pp.queryAdd("RDeRRC", RDeRRC);
    pp.queryAdd("wPhy", wPhy);
    pp.queryAdd("wLDet", wLDet);
    pp.queryAdd("wSDet", wSDet);
    pp.queryAdd("pCO2air", pCO2air);
    pp.queryAdd("po4", po4);
    pp.queryAdd("carbon", carbon);
    pp.queryAdd("oxygen", oxygen);
    pp.queryAdd("odu", odu);
    pp.queryAdd("denitrification", denitrification);
    pp.queryAdd("bio_sediment", bio_sediment);
    pp.queryAdd("river_don", river_don);
    pp.queryAdd("talk_nonconserv", talk_nonconserv);

    // Alkalinity only exists as a tracer under carbon, so every term this
    // option adds would have nowhere to go. ROMS would not compile in this
    // combination; say so rather than silently ignoring the request.
    if (talk_nonconserv && !carbon) {
        amrex::Abort("remora.fennel.talk_nonconserv requires remora.fennel.carbon: "
                     "alkalinity is a carbon-block tracer");
    }

    static std::string pco2air_type_string = "constant";
    pp.queryAdd("pco2air_type", pco2air_type_string);
    const std::string pco2air_type_ci = amrex::toLower(pco2air_type_string);
    if (pco2air_type_ci == "constant") {
        pco2air_type = PCO2AirType::constant;
    } else if (pco2air_type_ci == "data") {
        pco2air_type = PCO2AirType::data;
    } else if (pco2air_type_ci == "secular") {
        pco2air_type = PCO2AirType::secular;
    } else {
        amrex::Abort("Unknown remora.fennel.pco2air_type: " + pco2air_type_string +
                     ". Expected constant, data, or secular.");
    }

    static std::string co2_schmidt_string = "wanninkhof1992";
    pp.queryAdd("co2_schmidt", co2_schmidt_string);
    const std::string co2_schmidt_ci = amrex::toLower(co2_schmidt_string);
    if (co2_schmidt_ci == "wanninkhof1992" || co2_schmidt_ci == "w92") {
        co2_schmidt = CO2SchmidtType::wanninkhof1992;
    } else if (co2_schmidt_ci == "wanninkhof2014" || co2_schmidt_ci == "rw14") {
        co2_schmidt = CO2SchmidtType::wanninkhof2014;
    } else {
        amrex::Abort("Unknown remora.fennel.co2_schmidt: " + co2_schmidt_string +
                     ". Expected wanninkhof1992 or wanninkhof2014.");
    }

    static std::string o2_schmidt_string = "wanninkhof1992";
    pp.queryAdd("oxygen_schmidt", o2_schmidt_string);
    const std::string o2_schmidt_ci = amrex::toLower(o2_schmidt_string);
    if (o2_schmidt_ci == "wanninkhof1992" || o2_schmidt_ci == "w92") {
        o2_schmidt = O2SchmidtType::wanninkhof1992;
    } else if (o2_schmidt_ci == "wanninkhof2014" || o2_schmidt_ci == "rw14") {
        o2_schmidt = O2SchmidtType::wanninkhof2014;
    } else if (o2_schmidt_ci == "ocmip") {
        o2_schmidt = O2SchmidtType::ocmip;
    } else {
        amrex::Abort("Unknown remora.fennel.oxygen_schmidt: " + o2_schmidt_string +
                     ". Expected wanninkhof1992, wanninkhof2014, or ocmip.");
    }
}

BiologyModel
parse_biology_model (const std::string& name)
{
    const std::string model = amrex::toLower(name);
    if (model == "none" || model == "off") {
        return BiologyModel::none;
    }
    if (model == "fennel") {
        return BiologyModel::fennel;
    }
    amrex::Abort("Unknown remora.biology_model: " + name);
    return BiologyModel::none;
}

std::string
biology_model_name (BiologyModel model)
{
    switch (model) {
    case BiologyModel::none:
        return "none";
    case BiologyModel::fennel:
        return "fennel";
    }
    amrex::Abort("Invalid biology model");
    return "none";
}

BiologyICType
parse_biology_ic_type (const std::string& name)
{
    std::string lower = amrex::toLower(name);
    if (lower == "follow" || lower == "follow_ic_type" || lower == "default") {
        return BiologyICType::follow_ic_type;
    } else if (lower == "analytic") {
        return BiologyICType::analytic;
    } else if (lower == "netcdf") {
        return BiologyICType::netcdf;
    }
    amrex::Abort("remora.biology_ic_type must be one of: follow, analytic, netcdf");
    return BiologyICType::follow_ic_type;
}

std::string
biology_ic_type_name (BiologyICType type)
{
    switch (type) {
    case BiologyICType::follow_ic_type:
        return "follow";
    case BiologyICType::analytic:
        return "analytic";
    case BiologyICType::netcdf:
        return "netcdf";
    }
    amrex::Abort("Invalid biology IC type");
    return "follow";
}

bool
has_biology (BiologyModel model) noexcept
{
    return model != BiologyModel::none;
}

Vector<std::string>
tracer_names (BiologyModel model, FennelParameters const& fennel_parameters)
{
    switch (model) {
    case BiologyModel::none:
        return {};
    case BiologyModel::fennel: {
        Vector<std::string> names = {"NO3", "NH4", "chlorophyll", "phytoplankton",
                                     "zooplankton", "LdetritusN", "SdetritusN"};
        if (fennel_parameters.river_don) {
            names.emplace_back("RdetritusN");
        }
        if (fennel_parameters.po4) {
            names.emplace_back("PO4");
        }
        if (fennel_parameters.carbon) {
            names.emplace_back("LdetritusC");
            names.emplace_back("SdetritusC");
            names.emplace_back("TIC");
            names.emplace_back("alkalinity");
            if (fennel_parameters.river_don) {
                names.emplace_back("RdetritusC");
            }
        }
        if (fennel_parameters.oxygen) {
            names.emplace_back("oxygen");
        }
        if (fennel_parameters.odu) {
            names.emplace_back("ODU");
        }
        return names;
    }
    }
    amrex::Abort("Invalid biology model");
    return {};
}

}

void
REMORA::advance_biology (int lev, MultiFab const& mf_cons_old, MultiFab& mf_cons_new,
                         int N, Real dt_lev)
{
    // Add biological Source/Sink terms. Avoid computing source/sink
    // terms if no biological iterations are requested.
    if (!REMORABiology::has_biology(biology_model) || fennel_params.BioIter <= 0) {
        return;
    }

    if (biology_model != REMORABiology::BiologyModel::fennel) {
        amrex::Abort("advance_biology only supports fennel");
    }

#ifdef REMORA_USE_FENNEL_FORT
    // Path A: the ROMS Fortran kernel is the oracle until native parity is
    // proven for the active scope. Selected at run time; see
    // Source/Biology/Fortran/tag_map.md for the comparison protocol.
    if (use_biology_cpp_answer == 0) {
        advance_biology_fortran(lev, mf_cons_old, mf_cons_new, N, dt_lev);
        return;
    }
#endif

    const auto parms = fennel_params;

    const bool need_wind = (parms.oxygen || parms.carbon) && solverChoice.bulk_fluxes;
    const bool need_stress = (parms.oxygen || parms.carbon) && !solverChoice.bulk_fluxes;

    if (vec_Hz[lev] == nullptr || vec_z_w[lev] == nullptr || vec_mskr[lev] == nullptr ||
        vec_srflx[lev] == nullptr ||
        (need_stress && (vec_sustr[lev] == nullptr || vec_svstr[lev] == nullptr)) ||
        (need_wind && (vec_uwind[lev] == nullptr || vec_vwind[lev] == nullptr))) {
        amrex::Abort("Fennel biology requires Hz, z_w, rmask, srflx, and, when carbon or "
                     "oxygen is active, either uwind/vwind (bulk_fluxes) or sustr/svstr");
    }

    // Set time-stepping according to the number of iterations.
    const Real dtdays = dt_lev / Real(86400.0) / static_cast<Real>(parms.BioIter);
    const Real rho0 = solverChoice.rho0;
    const bool use_po4 = parms.po4;
    const bool use_carbon = parms.carbon;
    const bool use_oxygen = parms.oxygen;
    const bool use_odu = parms.odu;
    const bool use_denitrification = parms.denitrification;
    const bool use_river_don = parms.river_don;
    const bool use_river_don_c = parms.river_don && parms.carbon;
    const bool use_talk_nonconserv = parms.talk_nonconserv;
    const bool use_salt = use_oxygen || use_carbon;
    const bool do_bulk_flux = solverChoice.bulk_fluxes;

    // Set biological tracer component identifiers. The biology block starts after the
    // passive scalars, so this is Tracer_comp only when the run carries no dye.
    const auto bio_comp = REMORABiology::Fennel::components(parms, Bio_comp);

    // Schmidt-number coefficients and the leading rate coefficient for the gas
    // transfer velocity. Each pair belongs to one published relation, so they
    // are selected together rather than independently; ROMS makes the same
    // pairings with RW14_OXYGEN_SC, OCMIP_OXYGEN_SC and RW14_CO2_SC.
    Real A_O2 = Real(1953.4);       // Schmidt number
    Real B_O2 = Real(128.0);        // coefficients from
    Real C_O2 = Real(3.9918);       // Wanninkhof (1992)
    Real D_O2 = Real(0.050091);
    Real E_O2 = Real(0.0);
    Real o2_rate = Real(0.31);
    if (parms.o2_schmidt == REMORABiology::O2SchmidtType::wanninkhof2014) {
        A_O2 = Real(1920.4);
        B_O2 = Real(135.6);
        C_O2 = Real(5.2122);
        D_O2 = Real(0.10939);
        E_O2 = Real(0.00093777);
        o2_rate = Real(0.251);
    } else if (parms.o2_schmidt == REMORABiology::O2SchmidtType::ocmip) {
        // Keeling et al. (1998); Sc is slightly smaller up to about 35C. ROMS
        // pairs this set with the 1992 rate coefficient, not the 2014 one.
        A_O2 = Real(1638.0);
        B_O2 = Real(81.83);
        C_O2 = Real(1.483);
        D_O2 = Real(0.008004);
        E_O2 = Real(0.0);
    }

    Real A_CO2 = Real(2073.1);      // Schmidt number
    Real B_CO2 = Real(125.62);      // coefficients from
    Real C_CO2 = Real(3.6276);      // Wanninkhof (1992)
    Real D_CO2 = Real(0.043219);
    Real E_CO2 = Real(0.0);
    Real co2_rate = Real(0.31);
    if (parms.co2_schmidt == REMORABiology::CO2SchmidtType::wanninkhof2014) {
        A_CO2 = Real(2116.8);
        B_CO2 = Real(136.25);
        C_CO2 = Real(4.7353);
        D_CO2 = Real(0.092307);
        E_CO2 = Real(0.0007555);
        co2_rate = Real(0.251);
    }

    // Atmospheric pCO2 varies in time but not in space, so evaluate it once
    // here rather than once per column. t_old is the time at the start of the
    // step, which is the state the kernel reads.
    const Real pco2air = fennel_pco2_air(parms.pco2air_type, solverChoice.time_ref,
                                         t_old[lev], parms.pCO2air);

    // The secular fit is anchored to 1951 and extrapolates to a negative partial
    // pressure well before it, so a run whose clock sits near the origin of its
    // calendar silently draws CO2 out of the ocean for its whole length. That is
    // the default with remora.time_ref = 0, whose epoch is 0001-01-01.
    if (pco2air <= zero) {
        amrex::Abort("remora.fennel.pco2air_type = secular gives a non-positive atmospheric"
                     " pCO2 (" + std::to_string(pco2air) + " ppmv) at this model time. The"
                     " secular trend is fitted around 1951, so put the run in a real year:"
                     " set remora.time_ref to the reference date and remora.start_time to"
                     " the offset from it.");
    }

#ifdef REMORA_USE_BIOLOGY_DIAG
    // Path B half of the frozen diagnostic contract. Tag names, field order
    // and numeric format must stay byte-identical to the Fortran emitters in
    // Source/Biology/Fortran/REMORA_fennel_roms.F; see tag_map.md.
    const int dbg_level = biology_debug;
    const int dbg_i     = biology_debug_i;
    const int dbg_j     = biology_debug_j;
#endif

    // Scratch arrays mirror the ROMS Bio, qc, bR, bL, WR, WL, and
    // inverse-thickness work arrays, but use per-tile FArrayBox storage.
    // Optional tracer scratch is only allocated for active runtime options.
    int ncell = 0;
    const int sc_no3 = ncell++;
    const int sc_nh4 = ncell++;
    const int sc_chlo = ncell++;
    const int sc_phyt = ncell++;
    const int sc_zoop = ncell++;
    const int sc_lden = ncell++;
    const int sc_sden = ncell++;
    const int sc_rden = use_river_don ? ncell++ : -1;
    const int sc_po4 = use_po4 ? ncell++ : -1;
    const int sc_ldec = use_carbon ? ncell++ : -1;
    const int sc_sdec = use_carbon ? ncell++ : -1;
    const int sc_tic = use_carbon ? ncell++ : -1;
    const int sc_talk = use_carbon ? ncell++ : -1;
    const int sc_rdec = use_river_don_c ? ncell++ : -1;
    const int sc_oxyg = use_oxygen ? ncell++ : -1;
    const int sc_odu = use_odu ? ncell++ : -1;
    const int sc_temp = ncell++;
    const int sc_salt = use_salt ? ncell++ : -1;
    const int sc_inv_hz = ncell++;
    const int sc_inv_hz2 = ncell++;
    const int sc_inv_hz3 = ncell++;
    const int sc_qc = ncell++;
    const int sc_bR = ncell++;
    const int sc_bL = ncell++;
    const int sc_WR = ncell++;
    const int sc_WL = ncell++;

    int nw = 0;
    const int sw_FC = nw++;

    for (MFIter mfi(mf_cons_new, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        Box bx = mfi.tilebox();
        Box bx2d = bx;
        bx2d.makeSlab(2, 0);
        Box wbx = surroundingNodes(bx, 2);

        FArrayBox fab_cell(bx, ncell, amrex::The_Async_Arena());
        fab_cell.template setVal<RunOn::Device>(zero);
        FArrayBox fab_w(wbx, nw, amrex::The_Async_Arena());
        fab_w.template setVal<RunOn::Device>(zero);

        auto no3 = fab_cell.array(sc_no3);
        auto nh4 = fab_cell.array(sc_nh4);
        auto chlo = fab_cell.array(sc_chlo);
        auto phyt = fab_cell.array(sc_phyt);
        auto zoop = fab_cell.array(sc_zoop);
        auto lden = fab_cell.array(sc_lden);
        auto sden = fab_cell.array(sc_sden);
        Array4<Real> rden;
        if (use_river_don) {
            rden = fab_cell.array(sc_rden);
        }
        Array4<Real> po4;
        if (use_po4) {
            po4 = fab_cell.array(sc_po4);
        }
        Array4<Real> ldec;
        Array4<Real> sdec;
        Array4<Real> tic;
        Array4<Real> talk;
        if (use_carbon) {
            ldec = fab_cell.array(sc_ldec);
            sdec = fab_cell.array(sc_sdec);
            tic = fab_cell.array(sc_tic);
            talk = fab_cell.array(sc_talk);
        }
        Array4<Real> rdec;
        if (use_river_don_c) {
            rdec = fab_cell.array(sc_rdec);
        }
        Array4<Real> oxyg;
        if (use_oxygen) {
            oxyg = fab_cell.array(sc_oxyg);
        }
        Array4<Real> odu;
        if (use_odu) {
            odu = fab_cell.array(sc_odu);
        }
        auto temp = fab_cell.array(sc_temp);
        Array4<Real> salt;
        if (use_salt) {
            salt = fab_cell.array(sc_salt);
        }
        auto inv_hz = fab_cell.array(sc_inv_hz);
        auto inv_hz2 = fab_cell.array(sc_inv_hz2);
        auto inv_hz3 = fab_cell.array(sc_inv_hz3);
        auto qc = fab_cell.array(sc_qc);
        auto bR = fab_cell.array(sc_bR);
        auto bL = fab_cell.array(sc_bL);
        auto WR = fab_cell.array(sc_WR);
        auto WL = fab_cell.array(sc_WL);
        auto FC = fab_w.array(sw_FC);

        Array4<Real const> const& state_old = mf_cons_old.const_array(mfi);
        Array4<Real> const& state_new = mf_cons_new.array(mfi);
        Array4<Real const> const& Hz = vec_Hz[lev]->const_array(mfi);
        Array4<Real const> const& z_w = vec_z_w[lev]->const_array(mfi);
        Array4<Real const> const& srflx = vec_srflx[lev]->const_array(mfi);
        Array4<Real const> const& mskr = vec_mskr[lev]->const_array(mfi);
        // ROMS takes the gas-exchange transfer velocity from Uwind/Vwind
        // under BULK_FLUXES and from the surface stress otherwise; only the
        // selected pair is read.
        Array4<Real const> sustr;
        Array4<Real const> svstr;
        Array4<Real const> uwind;
        Array4<Real const> vwind;
        if (use_oxygen || use_carbon) {
            if (do_bulk_flux) {
                uwind = vec_uwind[lev]->const_array(mfi);
                vwind = vec_vwind[lev]->const_array(mfi);
            } else {
                sustr = vec_sustr[lev]->const_array(mfi);
                svstr = vec_svstr[lev]->const_array(mfi);
            }
        }

        ParallelFor(bx2d, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
            constexpr Real zero = Real(0.0);
            // Land columns contribute nothing: the increment below is scaled by rmask.
            // Skip them outright rather than evaluating the chemistry and multiplying the
            // result away -- a masked column may hold a NetCDF fill value, and the logs
            // and square roots in the pCO2 and oxygen-saturation blocks would turn that
            // into a NaN that survives the rmask factor (NaN * 0 is NaN) and lands in
            // cons_new. ROMS guards pCO2_water with #ifdef MASKING for the same reason.
            if (mskr(i,j,0) == zero) return;

            constexpr Real one = Real(1.0);
            constexpr Real two = Real(2.0);
            constexpr Real eps = Real(1.0e-20);
            constexpr Real minval = Real(1.0e-6);
            constexpr Real cff_weno = Real(1.0e-14);

            constexpr Real OA0 = Real(2.00907);       // Oxygen
            constexpr Real OA1 = Real(3.22014);       // saturation
            constexpr Real OA2 = Real(4.05010);       // coefficients
            constexpr Real OA3 = Real(4.94457);
            constexpr Real OA4 = Real(-0.256847);
            constexpr Real OA5 = Real(3.88767);
            constexpr Real OB0 = Real(-0.00624523);
            constexpr Real OB1 = Real(-0.00737614);
            constexpr Real OB2 = Real(-0.0103410);
            constexpr Real OB3 = Real(-0.00817083);
            constexpr Real OC0 = Real(-0.000000488682);
            constexpr Real rOxNO3 = Real(8.625);       // 138/16
            constexpr Real rOxNH4 = Real(6.625);       // 106/16
            constexpr Real rOxNH4Denit = Real(115.0) / Real(16.0);
            constexpr Real denitrification_NH4_fraction = Real(4.0) / Real(16.0);
            constexpr Real l2mol = Real(1000.0) / Real(22.3916); // liter to mol

            constexpr Real A1 = Real(-60.2409);       // surface
            constexpr Real A2 = Real(93.4517);        // CO2
            constexpr Real A3 = Real(23.3585);        // solubility
            constexpr Real B1 = Real(0.023517);       // coefficients
            constexpr Real B2 = Real(-0.023656);
            constexpr Real B3 = Real(0.0047036);

            // Compute inverse thickness to avoid repeated divisions.
            //
            // Extract biological variables from tracer arrays, place them
            // into scratch arrays, and restrict their values to be positive
            // definite. At input, tracers are read from the nstp state
            // (cons_old); the final increment is applied to nnew
            // (cons_new) below.
            for (int k = 0; k <= N; ++k) {
                inv_hz(i,j,k) = one / Hz(i,j,k);
                no3(i,j,k) = amrex::max(zero, state_old(i,j,k,bio_comp.no3));
                nh4(i,j,k) = amrex::max(zero, state_old(i,j,k,bio_comp.nh4));
                chlo(i,j,k) = amrex::max(zero, state_old(i,j,k,bio_comp.chlo));
                phyt(i,j,k) = amrex::max(zero, state_old(i,j,k,bio_comp.phyt));
                zoop(i,j,k) = amrex::max(zero, state_old(i,j,k,bio_comp.zoop));
                lden(i,j,k) = amrex::max(zero, state_old(i,j,k,bio_comp.lden));
                sden(i,j,k) = amrex::max(zero, state_old(i,j,k,bio_comp.sden));
                if (use_river_don) {
                    rden(i,j,k) = amrex::max(zero, state_old(i,j,k,bio_comp.rden));
                }
                if (use_po4) {
                    po4(i,j,k) = amrex::max(zero, state_old(i,j,k,bio_comp.po4));
                }
                if (use_carbon) {
                    ldec(i,j,k) = amrex::max(zero, state_old(i,j,k,bio_comp.ldec));
                    sdec(i,j,k) = amrex::max(zero, state_old(i,j,k,bio_comp.sdec));
                    tic(i,j,k) = amrex::max(zero, state_old(i,j,k,bio_comp.tic));
                    tic(i,j,k) = amrex::min(tic(i,j,k), Real(3000.0));
                    tic(i,j,k) = amrex::max(tic(i,j,k), Real(400.0));
                    talk(i,j,k) = amrex::max(zero, state_old(i,j,k,bio_comp.talk));
                    if (use_river_don) {
                        rdec(i,j,k) = amrex::max(zero, state_old(i,j,k,bio_comp.rdec));
                    }
                }
                if (use_oxygen) {
                    oxyg(i,j,k) = amrex::max(zero, state_old(i,j,k,bio_comp.oxyg));
                }
                if (use_odu) {
                    odu(i,j,k) = amrex::max(zero, state_old(i,j,k,bio_comp.odu));
                }
                // Extract potential temperature and salinity.
                temp(i,j,k) = amrex::min(state_old(i,j,k,Temp_comp), Real(35.0));
                if (use_salt) {
                    salt(i,j,k) = amrex::max(state_old(i,j,k,Salt_comp), zero);
                }
            }
            for (int k = 0; k < N; ++k) {
                inv_hz2(i,j,k) = one / (Hz(i,j,k) + Hz(i,j,k+1));
            }
            for (int k = 1; k < N; ++k) {
                inv_hz3(i,j,k) = one / (Hz(i,j,k-1) + Hz(i,j,k) + Hz(i,j,k+1));
            }

            // Calculate surface Photosynthetically Available Radiation
            // (PAR). REMORA stores srflx in Watts/m2, so only PARfrac is
            // applied here.
            const Real PARsur = parms.PARfrac * srflx(i,j,0);

            // Emit one Tier-1 tag for this column: every level, every active
            // tracer, in the frozen field order. Inactive tracers are skipped
            // on both paths so the two logs stay line-aligned.
            //
            // Two definitions rather than an #ifdef'd body: the disabled form
            // captures nothing, so the closure and its ~14 Array4 captures
            // disappear from the device lambda entirely instead of surviving
            // as dead weight. Call sites are identical either way.
#ifndef REMORA_USE_BIOLOGY_DIAG
            auto tag_state = [] (const char*, int) noexcept {};
#else
            auto tag_state = [=] (const char* tag, int it) noexcept
            {
                if (dbg_level <= 0) return;
                if (dbg_level == 1 && (i != dbg_i || j != dbg_j)) return;
                for (int k = 0; k <= N; ++k) {
                    remora_fennel_tag(tag, it, i, j, k, "NO3 ", no3(i,j,k));
                    remora_fennel_tag(tag, it, i, j, k, "NH4 ", nh4(i,j,k));
                    remora_fennel_tag(tag, it, i, j, k, "CHLO", chlo(i,j,k));
                    remora_fennel_tag(tag, it, i, j, k, "PHYT", phyt(i,j,k));
                    remora_fennel_tag(tag, it, i, j, k, "ZOOP", zoop(i,j,k));
                    remora_fennel_tag(tag, it, i, j, k, "LDEN", lden(i,j,k));
                    remora_fennel_tag(tag, it, i, j, k, "SDEN", sden(i,j,k));
                    if (use_river_don) {
                        remora_fennel_tag(tag, it, i, j, k, "RDEN", rden(i,j,k));
                    }
                    if (use_po4) {
                        remora_fennel_tag(tag, it, i, j, k, "PO4 ", po4(i,j,k));
                    }
                    if (use_carbon) {
                        remora_fennel_tag(tag, it, i, j, k, "LDEC", ldec(i,j,k));
                        remora_fennel_tag(tag, it, i, j, k, "SDEC", sdec(i,j,k));
                        remora_fennel_tag(tag, it, i, j, k, "TIC ", tic(i,j,k));
                        remora_fennel_tag(tag, it, i, j, k, "TALK", talk(i,j,k));
                        if (use_river_don) {
                            remora_fennel_tag(tag, it, i, j, k, "RDEC", rdec(i,j,k));
                        }
                    }
                    if (use_oxygen) {
                        remora_fennel_tag(tag, it, i, j, k, "OXYG", oxyg(i,j,k));
                    }
                    if (use_odu) {
                        remora_fennel_tag(tag, it, i, j, k, "ODU ", odu(i,j,k));
                    }
                }
            };
#endif

            tag_state("G00_PRE   ", 0);

            // Start internal iterations to achieve convergence of the
            // nonlinear backward-implicit solution.
            //
            // During the iterative procedure a series of fractional time
            // steps are performed in a chained mode, splitting different
            // biological conversion processes in sequence of the food chain.
            // In all stages the concentration of the component being consumed
            // is treated implicitly, so the algorithm guarantees
            // non-negative values. The overall algorithm, as well as any
            // stage of it, is formulated in conservative form except for
            // explicit sinking.
            //
            // The iterative loop is to iterate toward a backward-Euler
            // treatment of all terms; it damps numerical oscillations but
            // does not improve formal accuracy.
            for (int iter = 0; iter < parms.BioIter; ++iter) {

                // Light-limited computations.
                //
                // Compute attenuation coefficient based on chlorophyll-a in
                // each grid box. Then attenuate surface PAR down into the
                // water column. Thus, PAR at a depth depends on the whole
                // distribution of chlorophyll-a above.
                if (PARsur > zero) {
                    Real PAR = PARsur;
                    for (int k = N; k >= 0; --k) {
                        // Compute average light attenuation for each grid
                        // cell. Other attenuation contributions like
                        // suspended sediment or CDOM can be added here.
                        const Real thickness = z_w(i,j,k+1) - z_w(i,j,k);
                        const Real Att = (parms.AttSW + parms.AttChl * chlo(i,j,k)) * thickness;
                        const Real ExpAtt = std::exp(-Att);
                        const Real Itop = PAR;
                        const Real PARavg = Itop * (one - ExpAtt) / Att;

                        // Compute Chlorophyll-a phytoplankton ratio,
                        // [mg Chla / mg C].
                        const Real cff = parms.PhyCN * Real(12.0);
                        const Real Chl2C = amrex::min(chlo(i,j,k) / (phyt(i,j,k) * cff + eps), parms.Chl2C_m);

                        // Temperature-limited and light-limited growth rate
                        // (Eppley, R.W., 1972, Fishery Bulletin,
                        // 70: 1063-1085; here 0.59=ln(2)*0.851).
                        const Real Vp = parms.Vp0 * Real(0.59) * std::pow(Real(1.066), temp(i,j,k));
                        const Real fac1_light = PARavg * parms.PhyIS;
                        const Real Epp = Vp / std::sqrt(Vp * Vp + fac1_light * fac1_light);
                        const Real t_PPmax = Epp * fac1_light;

                        // Nutrient-limitation terms (Parker 1993 Ecol
                        // Mod., 66, 113-120).
                        const Real cff1 = nh4(i,j,k) * parms.K_NH4;
                        const Real cff2 = no3(i,j,k) * parms.K_NO3;
                        const Real inhNH4 = one / (one + cff1);
                        const Real L_NH4 = cff1 / (one + cff1);
                        const Real L_NO3 = cff2 * inhNH4 / (one + cff2);
                        const Real LTOT = L_NO3 + L_NH4;
                        Real LMIN = LTOT;
                        if (use_po4) {
                            const Real cff3 = po4(i,j,k) * parms.K_PO4;
                            const Real L_PO4 = cff3 / (one + cff3);
                            LMIN = amrex::min(LTOT, L_PO4);
                        }

                        // Nitrate and ammonium uptake by Phytoplankton.
                        Real cff4 = zero;
                        Real cff5 = zero;
                        Real cff6 = zero;
                        if (use_po4) {
                            const Real mu = dtdays * t_PPmax * LMIN;
                            cff4 = mu * phyt(i,j,k) * L_NO3 /
                                   amrex::max(minval, LTOT) / amrex::max(minval, no3(i,j,k));
                            cff5 = mu * phyt(i,j,k) * L_NH4 /
                                   amrex::max(minval, LTOT) / amrex::max(minval, nh4(i,j,k));
                            cff6 = parms.R_P2N * mu * phyt(i,j,k) /
                                   amrex::max(minval, po4(i,j,k));
                        } else {
                            const Real fac1 = dtdays * t_PPmax;
                            cff4 = fac1 * parms.K_NO3 * inhNH4 / (one + cff2) * phyt(i,j,k);
                            cff5 = fac1 * parms.K_NH4 / (one + cff1) * phyt(i,j,k);
                        }
                        no3(i,j,k) = no3(i,j,k) / (one + cff4);
                        nh4(i,j,k) = nh4(i,j,k) / (one + cff5);
                        if (use_po4) {
                            po4(i,j,k) = po4(i,j,k) / (one + cff6);
                        }
                        const Real N_Flux_NewProd = no3(i,j,k) * cff4;
                        const Real N_Flux_RegProd = nh4(i,j,k) * cff5;
                        phyt(i,j,k) += N_Flux_NewProd + N_Flux_RegProd;
                        chlo(i,j,k) += (dtdays * t_PPmax * t_PPmax * LMIN * LMIN *
                                        parms.Chl2C_m * chlo(i,j,k)) /
                                       (parms.PhyIS * amrex::max(Chl2C, eps) * PARavg + eps);
                        if (use_oxygen) {
                            oxyg(i,j,k) += N_Flux_NewProd * rOxNO3 + N_Flux_RegProd * rOxNH4;
                        }
                        if (use_carbon) {
                            // Total inorganic carbon (CO2) uptake during
                            // phytoplankton growth.
                            const Real cff1_carbon = parms.PhyCN *
                                                     (N_Flux_NewProd + N_Flux_RegProd);
                            tic(i,j,k) -= cff1_carbon;
                            if (use_talk_nonconserv) {
                                // Account for the uptake of NO3 on total alkalinity.
                                talk(i,j,k) += N_Flux_NewProd - N_Flux_RegProd;
                            }
                        }

                        // The Nitrification of NH4 ==> NO3 is thought to
                        // occur only in dark and only in aerobic water
                        // (Olson, R. J., 1981, JMR: 39, 227-238).
                        //
                        // NH4+ + 3/2 O2 ==> NO2- + H2O, via Nitrosomonas
                        // bacteria; NO2- + 1/2 O2 ==> NO3-, via Nitrobacter
                        // bacteria.
                        //
                        // Note that the entire process has a total loss of
                        // two moles of O2 per mole of NH4. If we were to
                        // resolve NO2 profiles, this is where we would change
                        // the code to split out the differential effects of
                        // the two different bacteria types. If OXYGEN is
                        // defined, nitrification is inhibited at low oxygen
                        // concentrations using a Michaelis-Menten term.
                        Real nitri_fac = dtdays * parms.NitriR;
                        if (use_oxygen) {
                            const Real fac2 = amrex::max(oxyg(i,j,k), zero); // O2 max
                            const Real fac3 = amrex::max(fac2 / (Real(3.0) + fac2), zero); // MM for O2 dependence
                            nitri_fac *= fac3;
                        }
                        const Real light_ratio = (PARavg - parms.I_thNH4) /
                                                 (parms.D_p5NH4 + PARavg - two * parms.I_thNH4);
                        const Real light_inhibit = one - amrex::max(zero, light_ratio);
                        const Real nitri = nitri_fac * light_inhibit;
                        nh4(i,j,k) = nh4(i,j,k) / (one + nitri);
                        const Real N_Flux_Nitrifi = nh4(i,j,k) * nitri;
                        no3(i,j,k) += N_Flux_Nitrifi;
                        if (use_oxygen) {
                            oxyg(i,j,k) -= two * N_Flux_Nitrifi;
                        }
                        if (use_talk_nonconserv) {
                            talk(i,j,k) -= two * N_Flux_Nitrifi;
                        }

                        // Light attenuation at the bottom of the grid cell.
                        // It is the starting PAR value for the next deeper
                        // vertical grid cell.
                        PAR = Itop * ExpAtt;
                    }
                } else {
                    // If PARsur=0, nitrification occurs at the maximum rate
                    // (NitriR).
                    const Real nitri = dtdays * parms.NitriR;
                    for (int k = N; k >= 0; --k) {
                        nh4(i,j,k) = nh4(i,j,k) / (one + nitri);
                        const Real N_Flux_Nitrifi = nh4(i,j,k) * nitri;
                        no3(i,j,k) += N_Flux_Nitrifi;
                        if (use_oxygen) {
                            oxyg(i,j,k) -= two * N_Flux_Nitrifi;
                        }
                        if (use_talk_nonconserv) {
                            talk(i,j,k) -= two * N_Flux_Nitrifi;
                        }
                    }
                }

                tag_state("G01_LIGHT ", iter + 1);

                // Phytoplankton grazing by zooplankton (rate: ZooGR),
                // phytoplankton assimilated to zooplankton (fraction:
                // ZooAE_N) and egested to small detritus, and
                // phytoplankton mortality (rate: PhyMR) to small detritus
                // [Landry 1993 L and O 38:468-472].
                const Real grazing_fac = dtdays * parms.ZooGR;
                const Real phyt_mort_fac = dtdays * parms.PhyMR;
                for (int k = 0; k <= N; ++k) {
                    // Phytoplankton grazing by zooplankton.
                    const Real phy2 = phyt(i,j,k) * phyt(i,j,k);
                    const Real graze = grazing_fac * zoop(i,j,k) * phyt(i,j,k) /
                                       (parms.K_Phy + phy2);
                    const Real consume = one / (one + graze);
                    phyt(i,j,k) *= consume;
                    chlo(i,j,k) *= consume;
                    // Phytoplankton assimilated to zooplankton and
                    // egested to small detritus.
                    const Real N_Flux_Assim = graze * phyt(i,j,k) * parms.ZooAE_N;
                    const Real N_Flux_Egest = phyt(i,j,k) * graze * (one - parms.ZooAE_N);
                    zoop(i,j,k) += N_Flux_Assim;
                    sden(i,j,k) += N_Flux_Egest;

                    // Phytoplankton mortality, limited by a
                    // phytoplankton minimum.
                    const Real N_Flux_Pmortal = phyt_mort_fac * amrex::max(phyt(i,j,k) - parms.PhyMin, zero);
                    phyt(i,j,k) -= N_Flux_Pmortal;
                    chlo(i,j,k) -= phyt_mort_fac * amrex::max(chlo(i,j,k) - parms.ChlMin, zero);
                    sden(i,j,k) += N_Flux_Pmortal;
                    if (use_carbon) {
                        sdec(i,j,k) += parms.PhyCN * (N_Flux_Egest + N_Flux_Pmortal) +
                                       (parms.PhyCN - parms.ZooCN) * N_Flux_Assim;
                    }
                }

                tag_state("G02_GRAZE ", iter + 1);

                // Zooplankton basal metabolism to NH4 (rate: ZooBM),
                // zooplankton mortality to small detritus (rate: ZooMR),
                // and zooplankton ingestion-related excretion (rate: ZooER).
                const Real zoo_metab_fac = dtdays * parms.ZooBM;
                const Real zoo_mort_fac = dtdays * parms.ZooMR;
                const Real zoo_excrete_fac = dtdays * parms.ZooER;
                for (int k = 0; k <= N; ++k) {
                    const Real phy2 = phyt(i,j,k) * phyt(i,j,k);
                    const Real ingestion_excretion = zoo_excrete_fac * phy2 /
                                                     (parms.K_Phy + phy2);
                    const Real cff2 = zoo_mort_fac * zoop(i,j,k);
                    const Real cff3 = ingestion_excretion * parms.ZooAE_N;
                    zoop(i,j,k) = zoop(i,j,k) / (one + cff2 + cff3);
                    // Zooplankton mortality and excretion.
                    const Real N_Flux_Zmortal = cff2 * zoop(i,j,k);
                    const Real N_Flux_Zexcret = cff3 * zoop(i,j,k);
                    nh4(i,j,k) += N_Flux_Zexcret;
                    if (use_po4) {
                        po4(i,j,k) += parms.R_P2N * N_Flux_Zexcret;
                    }
                    sden(i,j,k) += N_Flux_Zmortal;

                    // Zooplankton basal metabolism, limited by a
                    // zooplankton minimum.
                    const Real N_Flux_Zmetabo = zoo_metab_fac * amrex::max(zoop(i,j,k) - parms.ZooMin, zero);
                    zoop(i,j,k) -= N_Flux_Zmetabo;
                    nh4(i,j,k) += N_Flux_Zmetabo;
                    if (use_po4) {
                        po4(i,j,k) += parms.R_P2N * N_Flux_Zmetabo;
                    }
                    if (use_oxygen) {
                        oxyg(i,j,k) -= rOxNH4 * (N_Flux_Zmetabo + N_Flux_Zexcret);
                    }
                    if (use_carbon) {
                        sdec(i,j,k) += parms.ZooCN * N_Flux_Zmortal;
                        tic(i,j,k) += parms.ZooCN * (N_Flux_Zmetabo + N_Flux_Zexcret);
                        if (use_talk_nonconserv) {
                            talk(i,j,k) += N_Flux_Zmetabo + N_Flux_Zexcret;
                        }
                    }
                }

                tag_state("G03_ZMETAB", iter + 1);

                // Coagulation of phytoplankton and small detritus to
                // large detritus.
                const Real coag_fac = dtdays * parms.CoagR;
                for (int k = 0; k <= N; ++k) {
                    const Real coag = coag_fac * (sden(i,j,k) + phyt(i,j,k));
                    const Real consume = one / (one + coag);
                    phyt(i,j,k) *= consume;
                    chlo(i,j,k) *= consume;
                    sden(i,j,k) *= consume;
                    const Real N_Flux_CoagP = phyt(i,j,k) * coag;
                    const Real N_Flux_CoagD = sden(i,j,k) * coag;
                    lden(i,j,k) += N_Flux_CoagP + N_Flux_CoagD;
                    if (use_carbon) {
                        sdec(i,j,k) -= parms.PhyCN * N_Flux_CoagD;
                        ldec(i,j,k) += parms.PhyCN * (N_Flux_CoagP + N_Flux_CoagD);
                    }
                }

                tag_state("G04_COAG  ", iter + 1);

                // Detritus recycling to NH4, remineralization.
                if (use_oxygen) {
                    for (int k = 0; k <= N; ++k) {
                        const Real fac1 = amrex::max(oxyg(i,j,k) - Real(6.0), zero); // O2 off max
                        const Real fac2 = amrex::max(fac1 / (Real(3.0) + fac1), zero); // MM for O2 dependence
                        const Real remin_s = dtdays * parms.SDeRRN * fac2;
                        const Real remin_s_consume = one / (one + remin_s);
                        const Real remin_l = dtdays * parms.LDeRRN * fac2;
                        const Real remin_l_consume = one / (one + remin_l);
                        sden(i,j,k) *= remin_s_consume;
                        lden(i,j,k) *= remin_l_consume;
                        const Real N_Flux_RemineS = sden(i,j,k) * remin_s;
                        const Real N_Flux_RemineL = lden(i,j,k) * remin_l;
                        nh4(i,j,k) += N_Flux_RemineS + N_Flux_RemineL;
                        if (use_po4) {
                            po4(i,j,k) += parms.R_P2N * (N_Flux_RemineS + N_Flux_RemineL);
                        }
                        oxyg(i,j,k) -= (N_Flux_RemineS + N_Flux_RemineL) * rOxNH4;
                        if (use_talk_nonconserv) {
                            talk(i,j,k) += N_Flux_RemineS + N_Flux_RemineL;
                        }
                        if (use_river_don) {
                            const Real remin_r = dtdays * parms.RDeRRN * fac2;
                            const Real remin_r_consume = one / (one + remin_r);
                            rden(i,j,k) *= remin_r_consume;
                            const Real N_Flux_RemineR = rden(i,j,k) * remin_r;
                            nh4(i,j,k) += N_Flux_RemineR;
                            if (use_po4) {
                                po4(i,j,k) += parms.R_P2N * N_Flux_RemineR;
                            }
                            oxyg(i,j,k) -= N_Flux_RemineR * rOxNH4;
                            if (use_talk_nonconserv) {
                                talk(i,j,k) += N_Flux_RemineR;
                            }
                        }
                    }
                } else {
                    const Real remin_s = dtdays * parms.SDeRRN;
                    const Real remin_l = dtdays * parms.LDeRRN;
                    const Real remin_s_consume = one / (one + remin_s);
                    const Real remin_l_consume = one / (one + remin_l);
                    const Real remin_r = dtdays * parms.RDeRRN;
                    const Real remin_r_consume = one / (one + remin_r);
                    for (int k = 0; k <= N; ++k) {
                        sden(i,j,k) *= remin_s_consume;
                        lden(i,j,k) *= remin_l_consume;
                        const Real N_Flux_RemineS = sden(i,j,k) * remin_s;
                        const Real N_Flux_RemineL = lden(i,j,k) * remin_l;
                        nh4(i,j,k) += N_Flux_RemineS + N_Flux_RemineL;
                        if (use_po4) {
                            po4(i,j,k) += parms.R_P2N * (N_Flux_RemineS + N_Flux_RemineL);
                        }
                        if (use_talk_nonconserv) {
                            talk(i,j,k) += N_Flux_RemineS + N_Flux_RemineL;
                        }
                        if (use_river_don) {
                            rden(i,j,k) *= remin_r_consume;
                            const Real N_Flux_RemineR = rden(i,j,k) * remin_r;
                            nh4(i,j,k) += N_Flux_RemineR;
                            if (use_po4) {
                                po4(i,j,k) += parms.R_P2N * N_Flux_RemineR;
                            }
                            if (use_talk_nonconserv) {
                                talk(i,j,k) += N_Flux_RemineR;
                            }
                        }
                    }
                }

                tag_state("G05_REMINN", iter + 1);

                if (use_oxygen) {
                    // Surface O2 gas exchange.
                    //
                    // Compute surface O2 gas exchange.
                    const Real cff1 = rho0 * Real(550.0);
                    const Real cff2 = dtdays * o2_rate * Real(24.0) / Real(100.0);
                    const int k = N;

                    Real u10squ;
                    // Compute O2 transfer velocity: u10squared (u10 in m/s).
                    if (do_bulk_flux) {
                        u10squ = uwind(i,j,0) * uwind(i,j,0) + vwind(i,j,0) * vwind(i,j,0);
                    } else {
                        u10squ = cff1 *
                        std::sqrt((Real(0.5) * (sustr(i,j,0) + sustr(i+1,j,0))) *
                                  (Real(0.5) * (sustr(i,j,0) + sustr(i+1,j,0))) +
                                  (Real(0.5) * (svstr(i,j,0) + svstr(i,j+1,0))) *
                                  (Real(0.5) * (svstr(i,j,0) + svstr(i,j+1,0))));
                    }
                    const Real SchmidtN_Ox = A_O2 - temp(i,j,k) *
                        (B_O2 - temp(i,j,k) *
                         (C_O2 - temp(i,j,k) *
                          (D_O2 - temp(i,j,k) * E_O2)));
                    const Real cff3 = cff2 * u10squ * std::sqrt(Real(660.0) / SchmidtN_Ox);

                    // Calculate O2 saturation concentration using Garcia and
                    // Gordon L and O (1992) formula, (EXP(AA) is in ml/l).
                    const Real TS = std::log((Real(298.15) - temp(i,j,k)) /
                                             (Real(273.15) + temp(i,j,k)));
                    const Real AA = OA0 + TS * (OA1 + TS * (OA2 + TS * (OA3 + TS * (OA4 + TS * OA5)))) +
                                    salt(i,j,k) * (OB0 + TS * (OB1 + TS * (OB2 + TS * OB3))) +
                                    OC0 * salt(i,j,k) * salt(i,j,k);

                    // Convert from ml/l to mmol/m3.
                    const Real O2satu = l2mol * std::exp(AA);

                    // Add in O2 gas exchange.
                    const Real O2_Flux = cff3 * (O2satu - oxyg(i,j,k));
                    oxyg(i,j,k) += O2_Flux * inv_hz(i,j,k);

                    tag_state("G06_O2FLX ", iter + 1);
                }

                if (use_carbon) {
                    // Allow different remineralization rates for detrital C
                    // and detrital N.
                    const Real remin_sc = dtdays * parms.SDeRRC;
                    const Real remin_sc_consume = one / (one + remin_sc);
                    const Real remin_lc = dtdays * parms.LDeRRC;
                    const Real remin_lc_consume = one / (one + remin_lc);
                    const Real remin_rc = dtdays * parms.RDeRRC;
                    const Real remin_rc_consume = one / (one + remin_rc);
                    for (int k = 0; k <= N; ++k) {
                        sdec(i,j,k) *= remin_sc_consume;
                        ldec(i,j,k) *= remin_lc_consume;
                        const Real C_Flux_RemineS = sdec(i,j,k) * remin_sc;
                        const Real C_Flux_RemineL = ldec(i,j,k) * remin_lc;
                        tic(i,j,k) += C_Flux_RemineS + C_Flux_RemineL;
                        if (use_river_don_c) {
                            rdec(i,j,k) *= remin_rc_consume;
                            tic(i,j,k) += rdec(i,j,k) * remin_rc;
                        }
                    }

                    if (!use_talk_nonconserv) {
                        // Alkalinity is treated as a diagnostic variable. TAlk =
                        // f(S[PSU]) following Brewer et al. (1986). Under
                        // talk_nonconserv it is prognostic instead, and the
                        // increments accumulated above are what carry it, so
                        // overwriting here would discard every one of them.
                        for (int k = 0; k <= N; ++k) {
                            talk(i,j,k) = Real(587.05) + Real(50.56) * salt(i,j,k);
                        }
                    }

                    tag_state("G07_REMINC", iter + 1);

                    // Surface CO2 gas exchange.
                    //
                    // Compute equilibrium partial pressure inorganic carbon
                    // (ppmv) at the surface.
                    const int k = N;
                    const Real pCO2 = fennel_pco2_water(temp(i,j,k), salt(i,j,k),
                                                        tic(i,j,k), talk(i,j,k));

                    // Compute surface CO2 gas exchange.
                    const Real cff1 = rho0 * Real(550.0);
                    const Real cff2 = dtdays * co2_rate * Real(24.0) / Real(100.0);

                    // Compute CO2 transfer velocity: u10squared (u10 in m/s).
                    Real u10squ;
                    if (do_bulk_flux) {
                        u10squ = uwind(i,j,0) * uwind(i,j,0) + vwind(i,j,0) * vwind(i,j,0);
                    } else {
                        u10squ = cff1 *
                        std::sqrt((Real(0.5) * (sustr(i,j,0) + sustr(i+1,j,0))) *
                                  (Real(0.5) * (sustr(i,j,0) + sustr(i+1,j,0))) +
                                  (Real(0.5) * (svstr(i,j,0) + svstr(i,j+1,0))) *
                                  (Real(0.5) * (svstr(i,j,0) + svstr(i,j+1,0))));
                    }
                    const Real SchmidtN = A_CO2 - temp(i,j,k) *
                        (B_CO2 - temp(i,j,k) *
                         (C_CO2 - temp(i,j,k) *
                          (D_CO2 - temp(i,j,k) * E_CO2)));
                    const Real cff3 = cff2 * u10squ * std::sqrt(Real(660.0) / SchmidtN);

                    // Calculate CO2 solubility [mol/(kg.atm)] using Weiss
                    // (1974) formula.
                    const Real TempK = Real(0.01) * (temp(i,j,k) + Real(273.15));
                    const Real CO2_sol = std::exp(A1 +
                                                  A2 / TempK +
                                                  A3 * std::log(TempK) +
                                                  salt(i,j,k) * (B1 + TempK * (B2 + B3 * TempK)));

                    // Add in CO2 gas exchange.
                    const Real CO2_Flux = cff3 * CO2_sol * (pco2air - pCO2);
                    tic(i,j,k) += CO2_Flux * inv_hz(i,j,k);

                    tag_state("G08_CO2FLX", iter + 1);
                }

                // Vertical sinking terms.
                //
                // Set vertical sinking identification and sinking velocity
                // vectors in the same order: phytoplankton, chlorophyll,
                // small nitrogen-detritus, large nitrogen-detritus, and,
                // when CARBON is active, small and large carbon-detritus.
                const int nsink = use_carbon ? 6 : 4;
                for (int isink = 0; isink < nsink; ++isink) {
                    Array4<Real> bio = phyt;
                    Real wbio = parms.wPhy;
                    int sediment_type = 1;
                    const bool phyt_sink = (isink == 0);
                    if (isink == 1) {
                        bio = chlo;
                        wbio = parms.wPhy;
                        sediment_type = 0;
                    } else if (isink == 2) {
                        bio = sden;
                        wbio = parms.wSDet;
                    } else if (isink == 3) {
                        bio = lden;
                        wbio = parms.wLDet;
                    } else if (isink == 4) {
                        bio = sdec;
                        wbio = parms.wSDet;
                        sediment_type = 2;
                    } else if (isink == 5) {
                        bio = ldec;
                        wbio = parms.wLDet;
                        sediment_type = 2;
                    }

                    // Copy concentration of biological particulates into
                    // scratch array qc, restricted to positive values.
                    for (int k = 0; k <= N; ++k) {
                        qc(i,j,k) = bio(i,j,k);
                        bR(i,j,k) = qc(i,j,k);
                        bL(i,j,k) = qc(i,j,k);
                        WR(i,j,k) = zero;
                        WL(i,j,k) = zero;
                    }
                    for (int k = 0; k <= N + 1; ++k) {
                        FC(i,j,k) = zero;
                    }

                    // Reconstruct vertical profiles as parabolic segments
                    // within each grid box.
                    for (int iface = 1; iface <= N; ++iface) {
                        FC(i,j,iface) = (qc(i,j,iface) - qc(i,j,iface-1)) * inv_hz2(i,j,iface-1);
                    }
                    for (int k = 1; k < N; ++k) {
                        Real dltR = Hz(i,j,k) * FC(i,j,k+1);
                        Real dltL = Hz(i,j,k) * FC(i,j,k);
                        const Real cff = Hz(i,j,k-1) + two * Hz(i,j,k) + Hz(i,j,k+1);
                        const Real cffR = cff * FC(i,j,k+1);
                        const Real cffL = cff * FC(i,j,k);

                        // Apply PPM monotonicity constraint to prevent
                        // oscillations within the grid box.
                        if ((dltR * dltL) <= zero) {
                            dltR = zero;
                            dltL = zero;
                        } else if (std::abs(dltR) > std::abs(cffL)) {
                            dltR = cffL;
                        } else if (std::abs(dltL) > std::abs(cffR)) {
                            dltL = cffR;
                        }

                        // Compute right and left side values of parabolic
                        // segments; WR and WL are measures of quadratic
                        // variations.
                        const Real curv = (dltR - dltL) * inv_hz3(i,j,k);
                        dltR -= curv * Hz(i,j,k+1);
                        dltL += curv * Hz(i,j,k-1);
                        bR(i,j,k) = qc(i,j,k) + dltR;
                        bL(i,j,k) = qc(i,j,k) - dltL;
                        WR(i,j,k) = (two * dltR - dltL) * (two * dltR - dltL);
                        WL(i,j,k) = (dltR - two * dltL) * (dltR - two * dltL);
                    }
                    // Reconcile interfacial values using a WENO procedure
                    // so the whole profile remains monotonic.
                    for (int k = 1; k < N - 1; ++k) {
                        const Real dltL = amrex::max(cff_weno, WL(i,j,k));
                        const Real dltR = amrex::max(cff_weno, WR(i,j,k+1));
                        bR(i,j,k) = (dltR * bR(i,j,k) + dltL * bL(i,j,k+1)) / (dltR + dltL);
                        bL(i,j,k+1) = bR(i,j,k);
                    }
                    FC(i,j,N + 1) = zero;
                    bR(i,j,N) = qc(i,j,N);
                    bL(i,j,N) = qc(i,j,N);
                    bR(i,j,N - 1) = qc(i,j,N);
                    bL(i,j,1) = qc(i,j,0);
                    bR(i,j,0) = qc(i,j,0);
                    bL(i,j,0) = qc(i,j,0);

                    // Apply monotonicity constraint again, since
                    // reconciled interfacial values may cause non-monotonic
                    // behavior inside the grid box.
                    for (int k = 0; k <= N; ++k) {
                        Real dltR = bR(i,j,k) - qc(i,j,k);
                        Real dltL = qc(i,j,k) - bL(i,j,k);
                        const Real cffR = two * dltR;
                        const Real cffL = two * dltL;
                        if ((dltR * dltL) < zero) {
                            dltR = zero;
                            dltL = zero;
                        } else if (std::abs(dltR) > std::abs(cffL)) {
                            dltR = cffL;
                        } else if (std::abs(dltL) > std::abs(cffR)) {
                            dltL = cffR;
                        }
                        bR(i,j,k) = qc(i,j,k) + dltR;
                        bL(i,j,k) = qc(i,j,k) - dltL;
                    }

                    // After reconstruction, compute vertical advective
                    // fluxes. The algorithm is free of a CFL criterion by
                    // allowing semi-Lagrangian integration bounds to use as
                    // many upstream grid boxes as necessary.
                    //
                    // WL is the z-coordinate of the departure point for grid
                    // box interface z_w. FC is the finite-volume flux.
                    const Real sink_distance = dtdays * std::abs(wbio);
                    for (int k = 0; k <= N; ++k) {
                        FC(i,j,k) = zero;
                        WL(i,j,k) = z_w(i,j,k) + sink_distance;
                        WR(i,j,k) = Hz(i,j,k) * qc(i,j,k);
                    }
                    FC(i,j,N + 1) = zero;
                    for (int k = 0; k <= N; ++k) {
                        int ksource = k;
                        Real flux = zero;
                        for (int ks = k; ks < N; ++ks) {
                            if (WL(i,j,k) > z_w(i,j,ks+1)) {
                                ksource = ks + 1;
                                flux += WR(i,j,ks);
                            }
                        }
                        // Finalize flux computation by adding the
                        // fractional part from the source grid box.
                        const Real cu = amrex::min(one, (WL(i,j,k) - z_w(i,j,ksource)) *
                                                   inv_hz(i,j,ksource));
                        flux += Hz(i,j,ksource) * cu *
                                (bL(i,j,ksource) + cu * (Real(0.5) * (bR(i,j,ksource) - bL(i,j,ksource)) -
                                 (Real(1.5) - cu) * (bR(i,j,ksource) + bL(i,j,ksource) -
                                  two * qc(i,j,ksource))));
                        FC(i,j,k) = flux;
                    }
                    for (int k = 0; k <= N; ++k) {
                        bio(i,j,k) = qc(i,j,k) +
                                     (FC(i,j,k+1) - FC(i,j,k)) * inv_hz(i,j,k);
                    }
                    // Particulate flux reaching the seafloor is
                    // remineralized and returned to the dissolved nitrate
                    // pool. Without this conversion, particulate material
                    // falls out of the system. This is a temporary fix to
                    // restore total nitrogen conservation. It will be replaced
                    // later by a parameterization that includes the time delay
                    // of remineralization and dissolved oxygen.
                    if (parms.bio_sediment && sediment_type != 0) {
                        const Real cff1 = FC(i,j,0) * inv_hz(i,j,0);
                        if (sediment_type == 1) {
                            if (use_denitrification) {
                                nh4(i,j,0) += cff1 * denitrification_NH4_fraction;
                                if (use_oxygen) {
                                    oxyg(i,j,0) -= cff1 * rOxNH4Denit;
                                }
                            } else {
                                nh4(i,j,0) += cff1;
                                if (use_oxygen) {
                                    oxyg(i,j,0) -= cff1 * rOxNH4;
                                }
                                if (use_talk_nonconserv) {
                                    talk(i,j,0) += cff1;
                                }
                            }
                            if (use_po4) {
                                po4(i,j,0) += cff1 * parms.R_P2N;
                            }
                            if (use_carbon && phyt_sink) {
                                tic(i,j,0) += cff1 * parms.PhyCN;
                            }
                        } else if (use_carbon) {
                            tic(i,j,0) += cff1;
                        }
                    }
                }

                tag_state("G09_SINK  ", iter + 1);
            }

            // Emitted at the same semantic point as the Fortran path: right
            // after ITER_LOOP and before anything else touches the scratch
            // state, so it brackets the post-loop TIC clamp below. ROMS applies
            // that same clamp in the same place (fennel.h, just before "Update
            // global tracer variables"), so a divergence that appears between
            // G10 and G11 is in the clamp or the increment, not in whether the
            // clamp belongs there.
            tag_state("G10_POST  ", parms.BioIter);

            if (use_carbon) {
                for (int k = 0; k <= N; ++k) {
                    tic(i,j,k) = amrex::min(tic(i,j,k), Real(3000.0));
                    tic(i,j,k) = amrex::max(tic(i,j,k), Real(400.0));
                }
            }

            // Update global tracer variables: add increment due to BGC
            // processes to tracer array in time index nnew. Index nnew is
            // the solution after advection and mixing and has transport
            // units (m Tunits), hence the increment is multiplied by Hz.
            // Subtract the original values at the top of the routine to
            // account only for concentrations affected by BGC processes.
            // If Bio were unchanged by BGC processes, the increment would
            // be exactly zero. Final tracer values are not bounded >=0 so
            // total inventory can be preserved even when advection causes
            // tracer concentration to go negative.
            const Real rmask = mskr(i,j,0);
            for (int k = 0; k <= N; ++k) {
                const Real old_no3 = amrex::max(zero, state_old(i,j,k,bio_comp.no3));
                const Real old_nh4 = amrex::max(zero, state_old(i,j,k,bio_comp.nh4));
                Real old_po4 = zero;
                if (use_po4) {
                    old_po4 = amrex::max(zero, state_old(i,j,k,bio_comp.po4));
                }
                const Real old_chlo = amrex::max(zero, state_old(i,j,k,bio_comp.chlo));
                const Real old_phyt = amrex::max(zero, state_old(i,j,k,bio_comp.phyt));
                const Real old_zoop = amrex::max(zero, state_old(i,j,k,bio_comp.zoop));
                const Real old_lden = amrex::max(zero, state_old(i,j,k,bio_comp.lden));
                const Real old_sden = amrex::max(zero, state_old(i,j,k,bio_comp.sden));
                Real old_rden = zero;
                if (use_river_don) {
                    old_rden = amrex::max(zero, state_old(i,j,k,bio_comp.rden));
                }
                Real old_rdec = zero;
                if (use_river_don_c) {
                    old_rdec = amrex::max(zero, state_old(i,j,k,bio_comp.rdec));
                }
                Real old_ldec = zero;
                Real old_sdec = zero;
                Real old_tic = zero;
                Real old_talk = zero;
                if (use_carbon) {
                    old_ldec = amrex::max(zero, state_old(i,j,k,bio_comp.ldec));
                    old_sdec = amrex::max(zero, state_old(i,j,k,bio_comp.sdec));
                    old_tic = amrex::max(zero, state_old(i,j,k,bio_comp.tic));
                    old_tic = amrex::min(old_tic, Real(3000.0));
                    old_tic = amrex::max(old_tic, Real(400.0));
                    old_talk = amrex::max(zero, state_old(i,j,k,bio_comp.talk));
                }
                Real old_oxyg = zero;
                if (use_oxygen) {
                    old_oxyg = amrex::max(zero, state_old(i,j,k,bio_comp.oxyg));
                }
                Real old_odu = zero;
                if (use_odu) {
                    old_odu = amrex::max(zero, state_old(i,j,k,bio_comp.odu));
                }

                state_new(i,j,k,bio_comp.no3) += (no3(i,j,k) - old_no3) * rmask * Hz(i,j,k);
                state_new(i,j,k,bio_comp.nh4) += (nh4(i,j,k) - old_nh4) * rmask * Hz(i,j,k);
                if (use_po4) {
                    state_new(i,j,k,bio_comp.po4) += (po4(i,j,k) - old_po4) * rmask * Hz(i,j,k);
                }
                state_new(i,j,k,bio_comp.chlo) += (chlo(i,j,k) - old_chlo) * rmask * Hz(i,j,k);
                state_new(i,j,k,bio_comp.phyt) += (phyt(i,j,k) - old_phyt) * rmask * Hz(i,j,k);
                state_new(i,j,k,bio_comp.zoop) += (zoop(i,j,k) - old_zoop) * rmask * Hz(i,j,k);
                state_new(i,j,k,bio_comp.lden) += (lden(i,j,k) - old_lden) * rmask * Hz(i,j,k);
                state_new(i,j,k,bio_comp.sden) += (sden(i,j,k) - old_sden) * rmask * Hz(i,j,k);
                if (use_river_don) {
                    state_new(i,j,k,bio_comp.rden) += (rden(i,j,k) - old_rden) * rmask * Hz(i,j,k);
                }
                if (use_carbon) {
                    state_new(i,j,k,bio_comp.ldec) += (ldec(i,j,k) - old_ldec) * rmask * Hz(i,j,k);
                    state_new(i,j,k,bio_comp.sdec) += (sdec(i,j,k) - old_sdec) * rmask * Hz(i,j,k);
                    state_new(i,j,k,bio_comp.tic) += (tic(i,j,k) - old_tic) * rmask * Hz(i,j,k);
                    state_new(i,j,k,bio_comp.talk) += (talk(i,j,k) - old_talk) * rmask * Hz(i,j,k);
                    if (use_river_don) {
                        state_new(i,j,k,bio_comp.rdec) += (rdec(i,j,k) - old_rdec) * rmask * Hz(i,j,k);
                    }
                }
                if (use_oxygen) {
                    state_new(i,j,k,bio_comp.oxyg) += (oxyg(i,j,k) - old_oxyg) * rmask * Hz(i,j,k);
                }
                if (use_odu) {
                    state_new(i,j,k,bio_comp.odu) += (odu(i,j,k) - old_odu) * rmask * Hz(i,j,k);
                }
            }

#ifdef REMORA_USE_BIOLOGY_DIAG
            // G11_UPDATE reads the updated global tracer array rather than
            // the Bio scratch, matching fennel_tag_t on the Fortran path.
            // Guarded separately from tag_state because it reads state_new
            // rather than the Bio scratch and so cannot go through it.
            if (dbg_level > 0 &&
                (dbg_level == 2 || (i == dbg_i && j == dbg_j))) {
                const int it = parms.BioIter;
                for (int k = 0; k <= N; ++k) {
                    remora_fennel_tag("G11_UPDATE", it, i, j, k, "NO3 ",
                                      state_new(i,j,k,bio_comp.no3));
                    remora_fennel_tag("G11_UPDATE", it, i, j, k, "NH4 ",
                                      state_new(i,j,k,bio_comp.nh4));
                    remora_fennel_tag("G11_UPDATE", it, i, j, k, "CHLO",
                                      state_new(i,j,k,bio_comp.chlo));
                    remora_fennel_tag("G11_UPDATE", it, i, j, k, "PHYT",
                                      state_new(i,j,k,bio_comp.phyt));
                    remora_fennel_tag("G11_UPDATE", it, i, j, k, "ZOOP",
                                      state_new(i,j,k,bio_comp.zoop));
                    remora_fennel_tag("G11_UPDATE", it, i, j, k, "LDEN",
                                      state_new(i,j,k,bio_comp.lden));
                    remora_fennel_tag("G11_UPDATE", it, i, j, k, "SDEN",
                                      state_new(i,j,k,bio_comp.sden));
                    if (use_river_don) {
                        remora_fennel_tag("G11_UPDATE", it, i, j, k, "RDEN",
                                          state_new(i,j,k,bio_comp.rden));
                    }
                    if (use_po4) {
                        remora_fennel_tag("G11_UPDATE", it, i, j, k, "PO4 ",
                                          state_new(i,j,k,bio_comp.po4));
                    }
                    if (use_carbon) {
                        remora_fennel_tag("G11_UPDATE", it, i, j, k, "LDEC",
                                          state_new(i,j,k,bio_comp.ldec));
                        remora_fennel_tag("G11_UPDATE", it, i, j, k, "SDEC",
                                          state_new(i,j,k,bio_comp.sdec));
                        remora_fennel_tag("G11_UPDATE", it, i, j, k, "TIC ",
                                          state_new(i,j,k,bio_comp.tic));
                        remora_fennel_tag("G11_UPDATE", it, i, j, k, "TALK",
                                          state_new(i,j,k,bio_comp.talk));
                        if (use_river_don) {
                            remora_fennel_tag("G11_UPDATE", it, i, j, k, "RDEC",
                                              state_new(i,j,k,bio_comp.rdec));
                        }
                    }
                    if (use_oxygen) {
                        remora_fennel_tag("G11_UPDATE", it, i, j, k, "OXYG",
                                          state_new(i,j,k,bio_comp.oxyg));
                    }
                    if (use_odu) {
                        remora_fennel_tag("G11_UPDATE", it, i, j, k, "ODU ",
                                          state_new(i,j,k,bio_comp.odu));
                    }
                }
            }
#endif
        });
    }
}
