#ifndef REMORA_ECB_BIOLOGY_H_
#define REMORA_ECB_BIOLOGY_H_

#include <REMORA_Biology.H>

namespace REMORABiology {

/**
 * Runtime parameters of the ECB (Estuarine Carbon Biogeochemistry) model, read
 * from remora.ecb.*. Units and defaults follow ROMS External/bio_ECB.in at the
 * commit given in Docs/sphinx_doc/ECB.rst; the defaults are the values that
 * file sets, not the suggestions in its comments.
 */
struct ECBParameters {
    // Maximum number of iterations to achieve convergence of the nonlinear solution.
    int BioIter = 1;

    // Light attenuation due to seawater [1/m].
    amrex::Real AttSW = amrex::Real(0.04);

    // Light attenuation by chlorophyll [1/(mg_Chl m2)].
    amrex::Real AttChl = amrex::Real(0.02486);

    // Fraction of shortwave radiation that is photosynthetically active.
    amrex::Real PARfrac = amrex::Real(0.43);

    // Eppley temperature-limited and light-limited growth tuning parameter.
    amrex::Real Vp0 = amrex::Real(0.55);

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
    amrex::Real K_Phy = amrex::Real(2.0);

    // Maximum chlorophyll to carbon ratio [mg_Chl/mg_C].
    amrex::Real Chl2C_m = amrex::Real(0.0535);

    // Chlorophyll minimum threshold value [mg_Chl/m3].
    amrex::Real ChlMin = amrex::Real(0.001);

    // Phytoplankton Carbon:Nitrogen ratio [mole_C/mole_N].
    amrex::Real PhyCN = amrex::Real(6.625);

    // Phytoplankton Phosphate:Nitrogen ratio [mole_P/mole_N].
    amrex::Real R_P2N = amrex::Real(0.0625);

    // Phytoplankton initial slope of P-I curve [1/(Watts m-2 day)].
    amrex::Real PhyIS = amrex::Real(0.04);

    // Phytoplankton minimum threshold value [millimole_N/m3].
    amrex::Real PhyMin = amrex::Real(0.001);

    // Phytoplankton mortality rate [1/day].
    amrex::Real PhyMR = amrex::Real(0.05);

    // Phytoplankton exudation rate of labile DON [nondimensional].
    amrex::Real ElDON = amrex::Real(0.0);

    // Phytoplankton exudation rate of semilabile DON [nondimensional].
    amrex::Real EsDON = amrex::Real(0.07);

    // Fraction of semilabile DON to total DON in phytoplankton cell.
    amrex::Real deltN = amrex::Real(0.15);

    // Fraction of semilabile DOC to total DOC in phytoplankton cell.
    amrex::Real deltC = amrex::Real(0.275);

    // Parameter of carbon excess-based DOC exudation [nondimensional].
    amrex::Real gammaC = amrex::Real(0.20);

    // Fraction of semilabile DOC produced by carbon excess-based exudation.
    amrex::Real slCexc = amrex::Real(0.45);

    // Zooplankton Nitrogen assimilation efficiency.
    amrex::Real ZooAE_N = amrex::Real(0.75);

    // Basal metabolism of zooplankton [1/day].
    amrex::Real ZooBM = amrex::Real(0.1);

    // Zooplankton Carbon:Nitrogen ratio [mole_C/mole_N].
    amrex::Real ZooCN = amrex::Real(6.625);

    // Zooplankton specific excretion rate [1/day].
    amrex::Real ZooER = amrex::Real(0.1);

    // Zooplankton maximum growth rate [1/day].
    amrex::Real ZooGR = amrex::Real(0.03);

    // Zooplankton minimum threshold value [millimole_N/m3].
    amrex::Real ZooMin = amrex::Real(0.001);

    // Zooplankton mortality rate [1/day].
    amrex::Real ZooMR = amrex::Real(0.025);

    // Large nitrogen detritus remineralization rate at 0C [day-1].
    amrex::Real LDeRRN = amrex::Real(0.09);

    // Large carbon detritus remineralization rate at 0C [day-1].
    amrex::Real LDeRRC = amrex::Real(0.09);

    // Coagulation rate: aggregation rate of SDeN + Phyt ==> LDeN [m3 (millimole_N day)-1].
    amrex::Real CoagR = amrex::Real(0.008);

    // Small nitrogen detritus remineralization rate at 0C [day-1].
    amrex::Real SDeRRN = amrex::Real(0.09);

    // Small carbon detritus remineralization rate at 0C [day-1].
    amrex::Real SDeRRC = amrex::Real(0.09);

    // Semilabile DON remineralization rate at 0C [day-1].
    amrex::Real RDeRRN = amrex::Real(0.01);

    // Semilabile DOC remineralization rate at 0C [day-1].
    amrex::Real RDeRRC = amrex::Real(0.01);

    // Vertical sinking velocity for phytoplankton [m/day].
    amrex::Real wPhy = amrex::Real(0.1);

    // Vertical sinking velocity for large detritus [m/day].
    amrex::Real wLDet = amrex::Real(5.0);

    // Vertical sinking velocity for small detritus [m/day].
    amrex::Real wSDet = amrex::Real(0.1);

    // CO2 partial pressure in the air (parts per million by volume).
    amrex::Real pCO2air = amrex::Real(370.0);

    void init_params (const std::string& remora_prefix);
};

/**
 * ECB ecosystem with the fixed ROMS option set CARBON, OXYGEN, BIO_SEDIMENT,
 * DENITRIFICATION, PO4, TALK_NONCONSERV, RW14_OXYGEN_SC and RW14_CO2_SC.
 * Polymorphism is host-only; kernels capture plain data.
 */
class ECBBiology final : public Biology {
public:
    std::string name () const override { return "ecb"; }
    void configure (const std::string& remora_prefix) override;
    const BiologyDeclaration& declaration () const noexcept override { return m_declaration; }
    void advance (const BiologyAdvanceContext& context) const override;

private:
    ECBParameters m_parameters;
    BiologyDeclaration m_declaration;
};

}
#endif
