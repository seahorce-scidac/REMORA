#ifndef REMORA_ECB_BIOLOGY_H_
#define REMORA_ECB_BIOLOGY_H_
#include <REMORA_Biology.H>
namespace REMORABiology {
// Independent ECB parameters, in ROMS input units. All must be supplied.
struct ECBParameters {
    int BioIter = 1;
    amrex::Real AttSW = 0;
    amrex::Real AttChl = 0;
    amrex::Real Chl2C_m = 0;
    amrex::Real ChlMin = 0;
    amrex::Real CoagR = 0;
    amrex::Real D_p5NH4 = 0;
    amrex::Real I_thNH4 = 0;
    amrex::Real K_NH4 = 0;
    amrex::Real K_NO3 = 0;
    amrex::Real K_PO4 = 0;
    amrex::Real K_Phy = 0;
    amrex::Real LDeRRN = 0;
    amrex::Real LDeRRC = 0;
    amrex::Real NitriR = 0;
    amrex::Real PARfrac = 0;
    amrex::Real PhyCN = 0;
    amrex::Real R_P2N = 0;
    amrex::Real PhyIS = 0;
    amrex::Real PhyMin = 0;
    amrex::Real PhyMR = 0;
    amrex::Real SDeRRN = 0;
    amrex::Real SDeRRC = 0;
    amrex::Real RDeRRN = 0;
    amrex::Real RDeRRC = 0;
    amrex::Real Vp0 = 0;
    amrex::Real wLDet = 0;
    amrex::Real wPhy = 0;
    amrex::Real wSDet = 0;
    amrex::Real ZooAE_N = 0;
    amrex::Real ZooBM = 0;
    amrex::Real ZooCN = 0;
    amrex::Real ZooER = 0;
    amrex::Real ZooGR = 0;
    amrex::Real ZooMin = 0;
    amrex::Real ZooMR = 0;
    amrex::Real ElDON = 0;
    amrex::Real EsDON = 0;
    amrex::Real deltN = 0;
    amrex::Real deltC = 0;
    amrex::Real gammaC = 0;
    amrex::Real slCexc = 0;
    amrex::Real pCO2air = 0;
};
class ECBBiology final : public Biology {
public:
    std::string name () const override { return "ecb"; }
    std::string checkpoint_identity () const override {
        return "ecb/roms-31cabd318c9379429c365c308b9e0f2ee166a6fc/po4-talk-rw14-v1";
    }
    void configure (const std::string& remora_prefix) override;
    const BiologyDeclaration& declaration () const noexcept override { return m_declaration; }
    void advance (const BiologyAdvanceContext& context) const override;
private:
    ECBParameters m_parameters;
    BiologyDeclaration m_declaration;
};
}
#endif
