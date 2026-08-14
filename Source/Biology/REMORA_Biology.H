#ifndef REMORA_BIOLOGY_H_
#define REMORA_BIOLOGY_H_

#include <string>

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <REMORA_IndexDefines.H>

namespace REMORABiology {


enum class BiologyModel {
    none,
    fennel
};

/**
 * Source of the biology tracer initial condition, selected independently of
 * remora.ic_type.
 *
 * The two are decoupled because a NetCDF initial file supplies whichever
 * tracers it happens to contain, and validating an option set the file does
 * not cover (PO4, ODU) would otherwise require regenerating binary input.
 * With `analytic`, physical fields still come from NetCDF while biology comes
 * from the problem's analytic profile, which is reproducible from source.
 */
enum class BiologyICType {
    follow_ic_type,  //!< default: NetCDF when ic_type is netcdf, else analytic
    analytic,
    netcdf
};

struct FennelParameters {
    // Maximum number of iterations to achieve convergence of the nonlinear solution.
    int BioIter = 1;

    // Light attenuation due to seawater [1/m].
    amrex::Real AttSW = amrex::Real(0.04);

    // Light attenuation by chlorophyll [1/(mg_Chl m2)].
    amrex::Real AttChl = amrex::Real(0.02486);

    // Fraction of shortwave radiation that is photosynthetically active.
    amrex::Real PARfrac = amrex::Real(0.43);

    // Eppley temperature-limited and light-limited growth tuning parameter.
    amrex::Real Vp0 = amrex::Real(1.0);

    // Radiation threshold for nitrification inhibition [Watts/m2].
    amrex::Real I_thNH4 = amrex::Real(0.0095);

    // Half-saturation radiation for nitrification inhibition [Watts/m2].
    amrex::Real D_p5NH4 = amrex::Real(0.1);

    // Nitrification rate: oxidation of NH4 to NO3 [day-1].
    amrex::Real NitriR = amrex::Real(0.05);

    // Inverse half-saturation for phytoplankton NO3 uptake [1/(millimole_N m-3)].
    amrex::Real K_NO3 = amrex::Real(2.0);

    // Inverse half-saturation for phytoplankton NH4 uptake [1/(millimole_N m-3)].
    amrex::Real K_NH4 = amrex::Real(2.0);

    // Inverse half-saturation for phytoplankton PO4 uptake [1/(millimole_P m-3)].
    amrex::Real K_PO4 = amrex::Real(32.0);

    // Zooplankton half-saturation squared constant for ingestion [millimole_N m-3]^2.
    amrex::Real K_Phy = amrex::Real(1.0);

    // Maximum chlorophyll to carbon ratio [mg_Chl/mg_C].
    amrex::Real Chl2C_m = amrex::Real(0.0535);

    // Chlorophyll minimum threshold value [mg_Chl/m3].
    amrex::Real ChlMin = amrex::Real(0.01);

    // Phytoplankton Carbon:Nitrogen ratio [mole_C/mole_N].
    amrex::Real PhyCN = amrex::Real(6.625);

    // Phytoplankton Phosphate:Nitrogen ratio [mole_P/mole_N].
    amrex::Real R_P2N = amrex::Real(0.0625);

    // Phytoplankton NH4 inhibition parameter [1/(millimole_N)].
    amrex::Real PhyIP = amrex::Real(1.5);

    // Phytoplankton initial slope of P-I curve [1/(Watts m-2 day)].
    amrex::Real PhyIS = amrex::Real(0.025);

    // Phytoplankton minimum threshold value [millimole_N/m3].
    amrex::Real PhyMin = amrex::Real(0.01);

    // Phytoplankton mortality rate [1/day].
    amrex::Real PhyMR = amrex::Real(0.15);

    // Zooplankton Nitrogen assimilation efficiency.
    amrex::Real ZooAE_N = amrex::Real(0.75);

    // Zooplankton Carbon:Nitrogen ratio [mole_C/mole_N].
    amrex::Real ZooCN = amrex::Real(6.625);

    // Basal metabolism of zooplankton [1/day].
    amrex::Real ZooBM = amrex::Real(0.1);

    // Zooplankton specific excretion rate [1/day].
    amrex::Real ZooER = amrex::Real(0.1);

    // Zooplankton maximum growth rate [1/day].
    amrex::Real ZooGR = amrex::Real(0.6);

    // Zooplankton minimum threshold value [millimole_N/m3].
    amrex::Real ZooMin = amrex::Real(0.01);

    // Zooplankton mortality rate [1/day].
    amrex::Real ZooMR = amrex::Real(0.025);

    // Large nitrogen detritus remineralization rate [day-1].
    amrex::Real LDeRRN = amrex::Real(0.01);

    // Large detritus remineralization rate C-fraction [day-1].
    amrex::Real LDeRRC = amrex::Real(0.01);

    // Coagulation rate: aggregation rate of SDeN + Phyt ==> LDeN [day-1].
    amrex::Real CoagR = amrex::Real(0.005);

    // Small nitrogen detritus remineralization rate [day-1].
    amrex::Real SDeRRN = amrex::Real(0.1);

    // Small detritus remineralization rate C-fraction [day-1].
    amrex::Real SDeRRC = amrex::Real(0.03);

    // Vertical sinking velocity for phytoplankton [m/day].
    amrex::Real wPhy = amrex::Real(0.1);

    // Vertical sinking velocity for large detritus [m/day].
    amrex::Real wLDet = amrex::Real(1.0);

    // Vertical sinking velocity for small detritus [m/day].
    amrex::Real wSDet = amrex::Real(0.1);

    // CO2 partial pressure in the air (parts per million by volume).
    amrex::Real pCO2air = amrex::Real(370.0);

    // Enable ROMS PO4 option, adding phosphate as an active biological tracer.
    bool po4 = false;

    // Enable ROMS CARBON option, adding large and small carbon detritus,
    // total inorganic carbon, and alkalinity tracers.
    bool carbon = false;

    // Enable one additional biological tracer for dissolved oxygen. OXYGEN can
    // be activated independently of the CARBON option in ROMS.
    bool oxygen = false;

    // Enable ROMS ODU option, adding dissolved oxygen demand units.
    bool odu = false;

    // Runtime analogue of ROMS DENITRIFICATION bottom remineralization option.
    bool denitrification = false;

    // Runtime analogue of ROMS BIO_SEDIMENT bottom remineralization option.
    bool bio_sediment = true;

    void init_params (const std::string& remora_prefix);
};

namespace Fennel {
static constexpr int NumBaseTracers = 7;
static constexpr int NumCarbonTracers = 4;

struct Components {
    int no3 = -1;
    int nh4 = -1;
    int chlo = -1;
    int phyt = -1;
    int zoop = -1;
    int lden = -1;
    int sden = -1;
    int po4 = -1;
    int ldec = -1;
    int sdec = -1;
    int tic = -1;
    int talk = -1;
    int oxyg = -1;
    int odu = -1;
};

inline Components components (FennelParameters const& parameters) noexcept
{
    int comp = Tracer_comp;
    Components c;
    c.no3 = comp++;
    c.nh4 = comp++;
    c.chlo = comp++;
    c.phyt = comp++;
    c.zoop = comp++;
    c.lden = comp++;
    c.sden = comp++;
    c.po4 = parameters.po4 ? comp++ : -1;
    c.ldec = parameters.carbon ? comp++ : -1;
    c.sdec = parameters.carbon ? comp++ : -1;
    c.tic = parameters.carbon ? comp++ : -1;
    c.talk = parameters.carbon ? comp++ : -1;
    c.oxyg = parameters.oxygen ? comp++ : -1;
    c.odu = parameters.odu ? comp++ : -1;
    return c;
}

inline int num_tracers (FennelParameters const& parameters) noexcept
{
    return NumBaseTracers + (parameters.po4 ? 1 : 0) +
           (parameters.carbon ? NumCarbonTracers : 0) +
           (parameters.oxygen ? 1 : 0) + (parameters.odu ? 1 : 0);
}
}

BiologyModel parse_biology_model (const std::string& name);
std::string biology_model_name (BiologyModel model);
BiologyICType parse_biology_ic_type (const std::string& name);
std::string biology_ic_type_name (BiologyICType type);
bool has_biology (BiologyModel model) noexcept;
amrex::Vector<std::string> tracer_names (BiologyModel model,
                                         FennelParameters const& fennel_parameters);

}

#endif
