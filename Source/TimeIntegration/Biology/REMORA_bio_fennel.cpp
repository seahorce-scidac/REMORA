#include "REMORA_bio_fennel.H"

#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <algorithm>
#include <cmath>
#include <vector>

namespace REMORA_Biology {

void FennelBiology::LoadParams()
{
    if (m_params_loaded) {
        return;
    }

    amrex::ParmParse pp("remora");
    pp.query("fennel_mu0", m_params.mu0);
    pp.query("fennel_k_no3", m_params.k_no3);
    pp.query("fennel_k_nh4", m_params.k_nh4);
    pp.query("fennel_phy_mort", m_params.phy_mort);
    pp.query("fennel_phy_mort_to_ldet", m_params.phy_mort_to_ldet);
    pp.query("fennel_zoo_graz", m_params.zoo_graz);
    pp.query("fennel_zoo_graz_halfsat", m_params.zoo_graz_halfsat);
    pp.query("fennel_zoo_assim_n", m_params.zoo_assim_n);
    pp.query("fennel_zoo_mort", m_params.zoo_mort);
    pp.query("fennel_zoo_excr", m_params.zoo_excr);
    pp.query("fennel_zoo_basal_metab", m_params.zoo_basal_metab);
    pp.query("fennel_zoo_min", m_params.zoo_min);
    pp.query("fennel_remin", m_params.remin);
    pp.query("fennel_sdet_remin", m_params.sdet_remin);
    pp.query("fennel_ldet_remin", m_params.ldet_remin);
    pp.query("fennel_nitrif", m_params.nitrif);
    pp.query("fennel_nitrif_light_inhib", m_params.nitrif_light_inhib);
    pp.query("fennel_uptake_light_factor", m_params.uptake_light_factor);
    pp.query("fennel_par_surface", m_params.par_surface);
    pp.query("fennel_par_att_sw", m_params.par_att_sw);
    pp.query("fennel_par_att_chl", m_params.par_att_chl);
    pp.query("fennel_uptake_par_halfsat", m_params.uptake_par_halfsat);
    pp.query("fennel_nitrif_par_halfsat", m_params.nitrif_par_halfsat);
    pp.query("fennel_use_oxygen", m_params.use_oxygen);
    pp.query("fennel_o2_per_nitrif", m_params.o2_per_nitrif);
    pp.query("fennel_o2_per_nh4_remin", m_params.o2_per_nh4_remin);
        pp.query("fennel_use_carbon", m_params.use_carbon);
        pp.query("fennel_use_talk_nonconserv", m_params.use_talk_nonconserv);
        pp.query("fennel_c_per_n_remin", m_params.c_per_n_remin);
        pp.query("fennel_talk_per_n_remin", m_params.talk_per_n_remin);
        pp.query("fennel_talk_per_n_nitrif", m_params.talk_per_n_nitrif);
    pp.query("fennel_chlo_to_phy", m_params.chlo_to_phy);
    pp.query("fennel_chlo_relax", m_params.chlo_relax);
    pp.query("fennel_chlo_min", m_params.chlo_min);
    pp.query("fennel_diagnostics_enable", m_params.diagnostics_enable);
    pp.query("fennel_diagnostics_stride", m_params.diagnostics_stride);
    pp.query("fennel_phase_a_stub_active", m_params.phase_a_stub_active);

    // If explicit detrital rates are not provided, fall back to the bulk rate.
    if (m_params.sdet_remin <= amrex::Real(0.0)) {
        m_params.sdet_remin = m_params.remin;
    }
    if (m_params.ldet_remin <= amrex::Real(0.0)) {
        m_params.ldet_remin = m_params.remin;
    }

    m_params.nitrif_light_inhib = amrex::max(amrex::Real(0.0),
                                             amrex::min(amrex::Real(1.0), m_params.nitrif_light_inhib));
    m_params.uptake_light_factor = amrex::max(amrex::Real(0.0),
                                              amrex::min(amrex::Real(1.0), m_params.uptake_light_factor));
    m_params.phy_mort_to_ldet = amrex::max(amrex::Real(0.0),
                                           amrex::min(amrex::Real(1.0), m_params.phy_mort_to_ldet));
    m_params.chlo_to_phy = amrex::max(amrex::Real(0.0), m_params.chlo_to_phy);
    m_params.chlo_relax = amrex::max(amrex::Real(0.0), m_params.chlo_relax);
    m_params.chlo_min = amrex::max(amrex::Real(0.0), m_params.chlo_min);
    m_params.diagnostics_stride = std::max(1, m_params.diagnostics_stride);
    m_params.par_surface = amrex::max(amrex::Real(0.0), m_params.par_surface);
    m_params.par_att_sw = amrex::max(amrex::Real(0.0), m_params.par_att_sw);
    m_params.par_att_chl = amrex::max(amrex::Real(0.0), m_params.par_att_chl);
    m_params.uptake_par_halfsat = amrex::max(amrex::Real(1.0e-12), m_params.uptake_par_halfsat);
    m_params.nitrif_par_halfsat = amrex::max(amrex::Real(1.0e-12), m_params.nitrif_par_halfsat);
    m_params.o2_per_nitrif = amrex::max(amrex::Real(0.0), m_params.o2_per_nitrif);
    m_params.o2_per_nh4_remin = amrex::max(amrex::Real(0.0), m_params.o2_per_nh4_remin);
    m_params.c_per_n_remin = amrex::max(amrex::Real(0.0), m_params.c_per_n_remin);
    m_params.talk_per_n_remin = amrex::max(amrex::Real(0.0), m_params.talk_per_n_remin);
    m_params.talk_per_n_nitrif = amrex::max(amrex::Real(0.0), m_params.talk_per_n_nitrif);
    m_params.zoo_assim_n = amrex::max(amrex::Real(0.0),
                                      amrex::min(amrex::Real(1.0), m_params.zoo_assim_n));
    m_params.zoo_graz_halfsat = amrex::max(amrex::Real(1.0e-12), m_params.zoo_graz_halfsat);
    m_params.zoo_min = amrex::max(amrex::Real(0.0), m_params.zoo_min);

    m_params_loaded = true;
}

bool FennelBiology::BuildCompMap(int ncons, FennelCompMap& comp) const
{
    auto get_index = [&](const std::string& key) -> int {
        auto it = m_tracer_indices.find(key);
        if (it == m_tracer_indices.end()) {
            return -1;
        }
        int idx = it->second;
        return (idx >= 0 && idx < ncons) ? idx : -1;
    };

    comp.no3 = get_index("NO3");
    comp.nh4 = get_index("NH4");
    comp.phyt = get_index("Phyt");
    comp.chlo = get_index("Chlo");
    comp.zoop_s = get_index("ZoopS");
    comp.zoop_l = get_index("ZoopL");
    comp.sdet_n = get_index("SDet_N");
    comp.ldet_n = get_index("LDet_N");
    comp.don = get_index("DON");

    comp.oxy = get_index("Oxygen");
    comp.tic = get_index("TIC");
    comp.talk = get_index("TAlk");

    return (comp.no3 >= 0 && comp.nh4 >= 0 && comp.phyt >= 0 && comp.chlo >= 0 &&
            comp.zoop_s >= 0 && comp.zoop_l >= 0 && comp.sdet_n >= 0 &&
            comp.ldet_n >= 0 && comp.don >= 0);
}

void FennelBiology::ComputeTendencies(
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
)
{
    static bool warned_once = false;
    static bool warned_missing = false;
    static bool warned_missing_oxygen = false;
    static bool warned_missing_carbon = false;

    LoadParams();
    ++m_calls;

    (void)tile;
    (void)ng;
    (void)rmask_wet;
    (void)forcing_data;
    (void)m_ncons;

    FennelCompMap comp;
    if (!BuildCompMap(ncons, comp)) {
        m_missing_required_tracers = true;
        if (!warned_missing) {
            warned_missing = true;
            amrex::Warning("[REMORA] FennelBiology requires packed tracer components through DON, "
                           "but ncons is too small. Skipping fennel tendency update for this run");
        }
        return;
    }

    m_missing_required_tracers = false;

    if (!warned_once) {
        warned_once = true;
        amrex::Print() << "[REMORA] FennelBiology Phase A scaffold active (component mapping + tile kernel)."
                       << " Params: mu0=" << m_params.mu0
                       << ", k_no3=" << m_params.k_no3
                       << ", k_nh4=" << m_params.k_nh4
                       << ", phy_mort=" << m_params.phy_mort
                       << ", phy_mort_to_ldet=" << m_params.phy_mort_to_ldet
                       << ", zoo_graz=" << m_params.zoo_graz
                       << ", zoo_graz_halfsat=" << m_params.zoo_graz_halfsat
                       << ", zoo_assim_n=" << m_params.zoo_assim_n
                       << ", zoo_mort=" << m_params.zoo_mort
                       << ", zoo_excr=" << m_params.zoo_excr
                       << ", zoo_basal_metab=" << m_params.zoo_basal_metab
                       << ", zoo_min=" << m_params.zoo_min
                       << ", remin=" << m_params.remin
                       << ", sdet_remin=" << m_params.sdet_remin
                       << ", ldet_remin=" << m_params.ldet_remin
                       << ", nitrif=" << m_params.nitrif
                       << ", nitrif_light_inhib=" << m_params.nitrif_light_inhib
                       << ", uptake_light_factor=" << m_params.uptake_light_factor
                       << ", par_surface=" << m_params.par_surface
                       << ", par_att_sw=" << m_params.par_att_sw
                       << ", par_att_chl=" << m_params.par_att_chl
                       << ", uptake_par_halfsat=" << m_params.uptake_par_halfsat
                       << ", nitrif_par_halfsat=" << m_params.nitrif_par_halfsat
                       << ", use_oxygen=" << m_params.use_oxygen
                       << ", o2_per_nitrif=" << m_params.o2_per_nitrif
                       << ", o2_per_nh4_remin=" << m_params.o2_per_nh4_remin
                       << ", use_carbon=" << m_params.use_carbon
                       << ", use_talk_nonconserv=" << m_params.use_talk_nonconserv
                       << ", c_per_n_remin=" << m_params.c_per_n_remin
                       << ", talk_per_n_remin=" << m_params.talk_per_n_remin
                       << ", talk_per_n_nitrif=" << m_params.talk_per_n_nitrif
                       << ", chlo_to_phy=" << m_params.chlo_to_phy
                       << ", chlo_relax=" << m_params.chlo_relax
                       << ", chlo_min=" << m_params.chlo_min
                       << ", diagnostics_enable=" << m_params.diagnostics_enable
                       << ", diagnostics_stride=" << m_params.diagnostics_stride
                       << ". Insert ROMS fennel source/sink equations in REMORA_bio_fennel.cpp.\n";
    }

    // Phase A scaffold:
    // 1) Gather model-state variables at each cell
    // 2) Compute local source/sink tendencies (currently placeholders)
    // 3) Accumulate into rhs_tracers
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        using amrex::Real;

        if (rmask(i,j,0) <= Real(0.0)) {
            return;
        }

        const Real inv_hz = (Hz(i,j,k) > Real(0.0)) ? Real(1.0) / Hz(i,j,k) : Real(0.0);
        const Real cell_metric = pm(i,j,0) * pn(i,j,0);

        const Real no3 = amrex::max(state(i,j,k,comp.no3), Real(0.0));
        const Real nh4 = amrex::max(state(i,j,k,comp.nh4), Real(0.0));
        const Real phyt = amrex::max(state(i,j,k,comp.phyt), Real(0.0));
        const Real chlo = amrex::max(state(i,j,k,comp.chlo), Real(0.0));
        const Real zoop_s = amrex::max(state(i,j,k,comp.zoop_s), Real(0.0));
        const Real zoop_l = amrex::max(state(i,j,k,comp.zoop_l), Real(0.0));
        const Real sdet_n = amrex::max(state(i,j,k,comp.sdet_n), Real(0.0));
        const Real ldet_n = amrex::max(state(i,j,k,comp.ldet_n), Real(0.0));
        const Real don = amrex::max(state(i,j,k,comp.don), Real(0.0));
        const bool has_oxygen = (comp.oxy >= 0 && comp.oxy < ncons);
        const bool has_tic = (comp.tic >= 0 && comp.tic < ncons);
        const bool has_talk = (comp.talk >= 0 && comp.talk < ncons);

        // Phase A process port (incremental):
        //   (1) NO3/NH4 uptake partitioned by ROMS-like nutrient limitation
        //   (2) Zooplankton grazing + assimilation/egestion + mortality/metabolism/excretion
        //   (3) Phytoplankton mortality routed to SDet_N/LDet_N
        //   (4) Detritus remineralization: SDet_N/LDet_N -> NH4
        //   (5) Nitrification: NH4 -> NO3 (with simple light inhibition factor)
        // Rates are converted from [1/day] to [1/s].
        const Real day_to_sec = Real(1.0) / Real(86400.0);

        const Real k_uptake = m_params.mu0 * day_to_sec;
        const Real k_phy_mort = m_params.phy_mort * day_to_sec;
        const Real k_zoo_graz = m_params.zoo_graz * day_to_sec;
        const Real k_zoo_mort = m_params.zoo_mort * day_to_sec;
        const Real k_zoo_excr = m_params.zoo_excr * day_to_sec;
        const Real k_zoo_metab = m_params.zoo_basal_metab * day_to_sec;
        const Real k_remin_s = m_params.sdet_remin * day_to_sec;
        const Real k_remin_l = m_params.ldet_remin * day_to_sec;
        const Real k_nitrif = m_params.nitrif * day_to_sec;
        const Real k_chlo_relax = m_params.chlo_relax * day_to_sec;
        // Approximate depth-dependent light profile (PAR):
        // - use tile-normalized vertical coordinate as depth proxy
        // - attenuate by seawater background + chlorophyll self-shading
        const int klo = bx.smallEnd(2);
        const int khi = bx.bigEnd(2);
        const Real kspan = amrex::max(Real(1.0), Real(khi - klo));
        const Real depth_frac = amrex::max(Real(0.0), amrex::min(Real(1.0), Real(k - klo) / kspan));
        const Real att_total = m_params.par_att_sw + m_params.par_att_chl * amrex::max(chlo, Real(0.0));
        const Real par = m_params.par_surface * std::exp(-att_total * depth_frac);

        // Light modulation (ROMS-like behavior, simplified):
        // uptake increases with PAR; nitrification inhibition increases with PAR.
        const Real growth_light = m_params.uptake_light_factor *
            (par / (par + m_params.uptake_par_halfsat));
        const Real nitrif_light = par / (par + m_params.nitrif_par_halfsat);
        const Real light_factor = Real(1.0) - m_params.nitrif_light_inhib * nitrif_light;

        // ROMS-like partition logic:
        // cff1 = NH4*K_NH4, cff2 = NO3*K_NO3
        // L_NH4 = cff1/(1+cff1), L_NO3 = cff2*(1/(1+cff1))/(1+cff2)
        // LTOT = L_NH4 + L_NO3
        const Real cff1 = nh4 * m_params.k_nh4;
        const Real cff2 = no3 * m_params.k_no3;
        const Real inh_nh4 = Real(1.0) / (Real(1.0) + cff1);
        const Real l_nh4 = cff1 / (Real(1.0) + cff1);
        const Real l_no3 = cff2 * inh_nh4 / (Real(1.0) + cff2);
        const Real ltot = l_nh4 + l_no3;

        const Real uptake_total = k_uptake * growth_light * phyt * ltot;
        const Real uptake_no3 = (ltot > Real(0.0)) ? uptake_total * (l_no3 / ltot) : Real(0.0);
        const Real uptake_nh4 = (ltot > Real(0.0)) ? uptake_total * (l_nh4 / ltot) : Real(0.0);

        const Real phy2 = phyt * phyt;
        const Real graze_lim = phy2 / (m_params.zoo_graz_halfsat + phy2);
        const Real grazing_s = k_zoo_graz * graze_lim * zoop_s;
        const Real grazing_l = k_zoo_graz * graze_lim * zoop_l;
        const Real grazing_total = grazing_s + grazing_l;

        const Real assim_s = grazing_s * m_params.zoo_assim_n;
        const Real assim_l = grazing_l * m_params.zoo_assim_n;
        const Real egest_s = grazing_s * (Real(1.0) - m_params.zoo_assim_n);
        const Real egest_l = grazing_l * (Real(1.0) - m_params.zoo_assim_n);
        const Real egest_total = egest_s + egest_l;

        const Real zmort_s = k_zoo_mort * zoop_s * zoop_s;
        const Real zmort_l = k_zoo_mort * zoop_l * zoop_l;
        const Real zexcr_s = k_zoo_excr * graze_lim * zoop_s * m_params.zoo_assim_n;
        const Real zexcr_l = k_zoo_excr * graze_lim * zoop_l * m_params.zoo_assim_n;
        const Real zmet_s = k_zoo_metab * amrex::max(zoop_s - m_params.zoo_min, Real(0.0));
        const Real zmet_l = k_zoo_metab * amrex::max(zoop_l - m_params.zoo_min, Real(0.0));

        const Real phy_mort = k_phy_mort * phyt;
        const Real mort_to_ldet = phy_mort * m_params.phy_mort_to_ldet;
        const Real mort_to_sdet = phy_mort - mort_to_ldet;

        const Real remine_s = k_remin_s * sdet_n;
        const Real remine_l = k_remin_l * ldet_n;
        const Real nitrif = k_nitrif * nh4 * light_factor;

        // Chlorophyll coupling:
        // - gains tied to phyt growth with configurable Chl:Phy ratio
        // - losses tied to phyt mortality and grazing pressure
        // - optional relaxation toward target ratio m_params.chlo_to_phy * phyt
        const Real eps = Real(1.0e-12);
        const Real chlo_growth = m_params.chlo_to_phy * (uptake_no3 + uptake_nh4);
        const Real phy_loss_rate = (phy_mort + grazing_total) / amrex::max(phyt, eps);
        const Real chlo_above_min = amrex::max(chlo - m_params.chlo_min, Real(0.0));
        const Real chlo_loss = phy_loss_rate * chlo_above_min;
        const Real chlo_target = m_params.chlo_to_phy * phyt;
        const Real chlo_relax = k_chlo_relax * (chlo_target - chlo);

        const Real active = m_params.phase_a_stub_active ? Real(0.0) : Real(1.0);

        const Real d_no3 = active * (nitrif - uptake_no3) * inv_hz * cell_metric;
        const Real d_nh4 = active * (remine_s + remine_l - nitrif - uptake_nh4 + zexcr_s + zexcr_l + zmet_s + zmet_l)
                   * inv_hz * cell_metric;
        const Real d_phyt = active * (uptake_no3 + uptake_nh4 - phy_mort - grazing_total) * inv_hz * cell_metric;
        const Real d_chlo = active * (chlo_growth - chlo_loss + chlo_relax) * inv_hz * cell_metric;
        const Real d_zoop_s = active * (assim_s - zmort_s - zexcr_s - zmet_s) * inv_hz * cell_metric;
        const Real d_zoop_l = active * (assim_l - zmort_l - zexcr_l - zmet_l) * inv_hz * cell_metric;
        const Real d_sdet_n = active * (mort_to_sdet + egest_total + zmort_s + zmort_l - remine_s)
                     * inv_hz * cell_metric;
        const Real d_ldet_n = active * (mort_to_ldet - remine_l) * inv_hz * cell_metric;
        const Real d_don = Real(0.0);
        Real d_oxy = Real(0.0);
        if (m_params.use_oxygen && has_oxygen) {
            const Real nh4_regen = remine_s + remine_l + zexcr_s + zexcr_l + zmet_s + zmet_l;
            const Real oxy_sink = m_params.o2_per_nh4_remin * nh4_regen + m_params.o2_per_nitrif * nitrif;
            d_oxy = active * (-oxy_sink) * inv_hz * cell_metric;
        }

        Real d_tic = Real(0.0);
        if (m_params.use_carbon && has_tic) {
            const Real nh4_regen = remine_s + remine_l + zexcr_s + zexcr_l + zmet_s + zmet_l;
            d_tic = active * (m_params.c_per_n_remin * nh4_regen) * inv_hz * cell_metric;
        }

        Real d_talk = Real(0.0);
        if (m_params.use_talk_nonconserv && has_talk) {
            const Real nh4_regen = remine_s + remine_l + zexcr_s + zexcr_l + zmet_s + zmet_l;
            const Real talk_src = m_params.talk_per_n_remin * nh4_regen;
            const Real talk_sink = m_params.talk_per_n_nitrif * nitrif;
            d_talk = active * (talk_src - talk_sink) * inv_hz * cell_metric;
        }

        rhs_tracers(i,j,k,comp.no3) += d_no3;
        rhs_tracers(i,j,k,comp.nh4) += d_nh4;
        rhs_tracers(i,j,k,comp.phyt) += d_phyt;
        rhs_tracers(i,j,k,comp.chlo) += d_chlo;
        rhs_tracers(i,j,k,comp.zoop_s) += d_zoop_s;
        rhs_tracers(i,j,k,comp.zoop_l) += d_zoop_l;
        rhs_tracers(i,j,k,comp.sdet_n) += d_sdet_n;
        rhs_tracers(i,j,k,comp.ldet_n) += d_ldet_n;
        rhs_tracers(i,j,k,comp.don) += d_don;
        if (m_params.use_oxygen && has_oxygen) {
            rhs_tracers(i,j,k,comp.oxy) += d_oxy;
        }
        if (m_params.use_carbon && has_tic) {
            rhs_tracers(i,j,k,comp.tic) += d_tic;
        }
        if (m_params.use_talk_nonconserv && has_talk) {
            rhs_tracers(i,j,k,comp.talk) += d_talk;
        }
    });

    if (m_params.use_oxygen && (comp.oxy < 0 || comp.oxy >= ncons) && !warned_missing_oxygen) {
        warned_missing_oxygen = true;
        amrex::Warning("[REMORA] fennel_use_oxygen=true but oxygen tracer component is not available in state. "
                       "Oxygen coupling is skipped.");
    }

    if ((m_params.use_carbon && (comp.tic < 0 || comp.tic >= ncons) ||
         m_params.use_talk_nonconserv && (comp.talk < 0 || comp.talk >= ncons)) &&
        !warned_missing_carbon) {
        warned_missing_carbon = true;
        amrex::Warning("[REMORA] fennel_use_carbon/fennel_use_talk_nonconserv enabled but TIC/TAlk tracer "
                       "components are not available in state. Carbon/alkalinity coupling is skipped.");
    }

    if (m_params.diagnostics_enable) {
        static bool warned_reductions_unavailable = false;
        ++m_diag_calls;
        if ((m_diag_calls % m_params.diagnostics_stride) == 0) {
            // Keep this path lightweight and compiler-portable; detailed spatial
            // reductions can be reintroduced with FabArray-based reductions.
            if (!warned_reductions_unavailable) {
                warned_reductions_unavailable = true;
                amrex::Print() << "[REMORA] Fennel diagnostics reductions are in compatibility mode; "
                               << "only fennel_diag_calls is accumulated in this build.\n";
            }
        }
    }
}

void FennelBiology::ComputeConservationBudget(
    const amrex::Array4<const amrex::Real>& state,
    std::map<std::string, amrex::Real>& budgets)
{
    (void)state;

    const auto requirements = GetModelTracerRequirements("fennel");
    AddScaffoldBudgetFields("fennel",
                            m_ncons,
                            m_tracer_indices,
                            m_calls,
                            requirements.required_tracer_keys,
                            requirements.optional_tracer_keys,
                            m_missing_required_tracers,
                            budgets);

    budgets["fennel_missing_required_tracers"] =
        m_missing_required_tracers ? amrex::Real(1.0) : amrex::Real(0.0);
    budgets["fennel_diag_calls"] = static_cast<amrex::Real>(m_diag_calls);
    budgets["fennel_uptake_no3_sum"] = m_diag_uptake_no3;
    budgets["fennel_uptake_nh4_sum"] = m_diag_uptake_nh4;
    budgets["fennel_nitrif_sum"] = m_diag_nitrif;
    budgets["fennel_remin_sum"] = m_diag_remin;
    budgets["fennel_grazing_sum"] = m_diag_grazing;
}

} // namespace REMORA_Biology
