#ifndef REMORA_BIO_FENNEL_H_
#define REMORA_BIO_FENNEL_H_

#include "REMORA_biology.H"
#include "REMORA_IndexDefines.H"

namespace REMORA_Biology {

class FennelBiology final : public BiologyModel {
public:
    explicit FennelBiology(int ncons, const TracerIndexMap& tracer_indices)
        : m_ncons(ncons), m_tracer_indices(tracer_indices)
    {}

    void ComputeTendencies(
        const int& tile,
        const int& ng,
        const amrex::Box& bx,
        const amrex::Array4<const amrex::Real>& state,
        const amrex::Array4<const amrex::Real>& Hz,
        const amrex::Array4<const amrex::Real>& pm,
        const amrex::Array4<const amrex::Real>& pn,
        const amrex::Array4<const amrex::Real>& rmask,
        const amrex::Array4<const amrex::Real>& rmask_wet,
        const amrex::Array4<const amrex::Real>& forcing_data,
        amrex::Array4<amrex::Real>& rhs_tracers,
        int ncons
    ) override;

    int GetNumTracers() const override {
        return BiologyTracers::NumTracerTypes;
    }

    void ComputeConservationBudget(
        const amrex::Array4<const amrex::Real>& state,
        std::map<std::string, amrex::Real>& budgets
    ) override;

private:
    struct FennelParams {
        // Phase A core defaults (placeholders; to be aligned with ROMS fennel inputs)
        amrex::Real mu0 = amrex::Real(1.0);             // max phytoplankton growth rate [1/day]
        amrex::Real k_no3 = amrex::Real(0.5);           // nitrate half-sat [mmol N m-3]
        amrex::Real k_nh4 = amrex::Real(0.1);           // ammonium half-sat [mmol N m-3]
        amrex::Real phy_mort = amrex::Real(0.03);       // phytoplankton mortality [1/day]
        amrex::Real phy_mort_to_ldet = amrex::Real(0.5);// mortality split to large detritus [0-1]
        amrex::Real zoo_graz = amrex::Real(0.2);        // zooplankton grazing coefficient [1/day]
        amrex::Real zoo_graz_halfsat = amrex::Real(0.25);// phyt^2 half-saturation [mmol N m-3]^2
        amrex::Real zoo_assim_n = amrex::Real(0.75);    // assimilation efficiency [0-1]
        amrex::Real zoo_mort = amrex::Real(0.02);       // zooplankton mortality coefficient [1/day]
        amrex::Real zoo_excr = amrex::Real(0.02);       // ingestion-related excretion [1/day]
        amrex::Real zoo_basal_metab = amrex::Real(0.03);// basal metabolism [1/day]
        amrex::Real zoo_min = amrex::Real(0.01);        // minimum zooplankton threshold [mmol N m-3]
        amrex::Real remin = amrex::Real(0.02);          // bulk remineralization [1/day]
        amrex::Real sdet_remin = amrex::Real(0.02);     // small detritus remin [1/day]
        amrex::Real ldet_remin = amrex::Real(0.01);     // large detritus remin [1/day]
        amrex::Real nitrif = amrex::Real(0.03);         // nitrification [1/day]
        amrex::Real nitrif_light_inhib = amrex::Real(0.5); // light inhibition factor [0-1]
        amrex::Real uptake_light_factor = amrex::Real(0.5); // simple growth light factor [0-1]
        amrex::Real par_surface = amrex::Real(100.0);   // surface PAR proxy [arbitrary]
        amrex::Real par_att_sw = amrex::Real(0.10);     // background attenuation coefficient
        amrex::Real par_att_chl = amrex::Real(0.01);    // chlorophyll-dependent attenuation coefficient
        amrex::Real uptake_par_halfsat = amrex::Real(30.0); // light half-sat for uptake
        amrex::Real nitrif_par_halfsat = amrex::Real(30.0); // light half-sat for nitrification inhibition
        bool use_oxygen = false;                        // enable oxygen coupling if oxygen tracer exists
        amrex::Real o2_per_nitrif = amrex::Real(2.0);  // O2 sink per NH4 nitrified
        amrex::Real o2_per_nh4_remin = amrex::Real(6.625); // O2 sink per NH4 regenerated
        bool use_carbon = false;                        // enable TIC coupling if TIC tracer exists
        bool use_talk_nonconserv = false;               // enable non-conservative alkalinity coupling
        amrex::Real c_per_n_remin = amrex::Real(6.625); // TIC source per regenerated NH4
        amrex::Real talk_per_n_remin = amrex::Real(1.0);// TAlk source per regenerated NH4
        amrex::Real talk_per_n_nitrif = amrex::Real(2.0);// TAlk sink per nitrified NH4
        amrex::Real chlo_to_phy = amrex::Real(1.0);     // chlorophyll produced per phyt growth [chl units / N units]
        amrex::Real chlo_relax = amrex::Real(0.0);      // relaxation to target Chl:Phy [1/day]
        amrex::Real chlo_min = amrex::Real(0.0);        // minimum chlorophyll floor [chl units]
        bool diagnostics_enable = false;                // enable lightweight process-rate diagnostics
        int diagnostics_stride = 1;                     // sample every N calls to ComputeTendencies
        bool phase_a_stub_active = true;                // keep tendencies zero until full port
    };

    void LoadParams();
    struct FennelCompMap {
        int no3 = -1;
        int nh4 = -1;
        int phyt = -1;
        int chlo = -1;
        int zoop_s = -1;
        int zoop_l = -1;
        int sdet_n = -1;
        int ldet_n = -1;
        int don = -1;
        int oxy = -1;
        int tic = -1;
        int talk = -1;
    };

    bool BuildCompMap(int ncons, FennelCompMap& comp) const;

    int m_ncons;
    TracerIndexMap m_tracer_indices;
    FennelParams m_params;
    bool m_params_loaded = false;
    bool m_missing_required_tracers = false;
    long long m_calls = 0;

    long long m_diag_calls = 0;
    amrex::Real m_diag_uptake_no3 = 0.0;
    amrex::Real m_diag_uptake_nh4 = 0.0;
    amrex::Real m_diag_nitrif = 0.0;
    amrex::Real m_diag_remin = 0.0;
    amrex::Real m_diag_grazing = 0.0;
};

} // namespace REMORA_Biology

#endif // REMORA_BIO_FENNEL_H_
