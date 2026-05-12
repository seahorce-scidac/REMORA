#ifndef REMORA_BIO_PLACEHOLDER_H_
#define REMORA_BIO_PLACEHOLDER_H_

#include "REMORA_biology.H"

#include <string>

namespace REMORA_Biology {

class PlaceholderBiology final : public BiologyModel {
public:
    PlaceholderBiology(std::string model_name, int ncons, const TracerIndexMap& tracer_indices)
        : m_model_name(std::move(model_name)), m_ncons(ncons), m_tracer_indices(tracer_indices)
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

    void ComputeConservationBudget(
        const amrex::Array4<const amrex::Real>& state,
        std::map<std::string, amrex::Real>& budgets
    ) override;

    int GetNumTracers() const override {
        return 0;
    }

private:
    std::string m_model_name;
    int m_ncons;
    TracerIndexMap m_tracer_indices;
    mutable bool m_warned_once = false;
    long long m_calls = 0;
};

} // namespace REMORA_Biology

#endif // REMORA_BIO_PLACEHOLDER_H_