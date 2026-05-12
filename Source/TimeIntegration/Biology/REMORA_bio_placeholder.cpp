#include "REMORA_bio_placeholder.H"

#include <AMReX_Print.H>

namespace REMORA_Biology {

void PlaceholderBiology::ComputeTendencies(
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
    int ncons)
{
    amrex::ignore_unused(tile, ng, bx, state, Hz, pm, pn, rmask, rmask_wet, forcing_data, rhs_tracers, ncons);
    amrex::ignore_unused(m_ncons);

    ++m_calls;
    if (!m_warned_once) {
        m_warned_once = true;
        amrex::Warning("[REMORA] Biology model '" + m_model_name +
                       "' is a placeholder plugin; no biology tendencies are applied");
    }
}

void PlaceholderBiology::ComputeConservationBudget(
    const amrex::Array4<const amrex::Real>& state,
    std::map<std::string, amrex::Real>& budgets)
{
    amrex::ignore_unused(state);
    const auto requirements = GetModelTracerRequirements(m_model_name);
    const bool missing_required = MissingAnyRequiredTracers(m_tracer_indices,
                                                            m_ncons,
                                                            requirements.required_tracer_keys);

    AddScaffoldBudgetFields(m_model_name,
                            m_ncons,
                            m_tracer_indices,
                            m_calls,
                            requirements.required_tracer_keys,
                            requirements.optional_tracer_keys,
                            missing_required,
                            budgets);

    budgets["placeholder_diag_calls"] = static_cast<amrex::Real>(m_calls);
    budgets["placeholder_ncons"] = static_cast<amrex::Real>(m_ncons);
    budgets["placeholder_tracer_map_size"] = static_cast<amrex::Real>(m_tracer_indices.size());
}

} // namespace REMORA_Biology