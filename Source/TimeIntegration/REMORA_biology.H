#ifndef REMORA_BIOLOGY_H_
#define REMORA_BIOLOGY_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Array.H>
#include <memory>
#include <map>
#include <string>
#include <vector>

/**
 * \file REMORA_biology.H
 * 
 * \brief Generic biology framework for REMORA
 *
 * This module provides the pluggable biology framework for REMORA.
 * Individual models (fennel, ecb, ecosim, etc.) inherit from BiologyModel
 * and implement model-specific tendencies.
 *
 * Design:
 * - One generic biology driver (REMORA_biology.cpp)
 * - One abstract base class (BiologyModel) that all plugins inherit
 * - Model registry and factory for runtime selection
 * - Each plugin (REMORA_bio_*.H/cpp) implements the interface
 */

namespace REMORA_Biology {

using TracerIndexMap = std::map<std::string, int>;

struct ModelTracerRequirements {
    std::string normalized_model_name;
    std::vector<std::string> required_tracer_keys;
    std::vector<std::string> optional_tracer_keys;
    bool uses_packed_don_layout = false;
    int packed_min_nscalar = 0;
};

std::string NormalizeBiologyModelName(const std::string& model_name);
ModelTracerRequirements GetModelTracerRequirements(const std::string& model_name);
TracerIndexMap BuildDefaultTracerMapForModel(const std::string& model_name, int base_comp);
bool MissingAnyRequiredTracers(const TracerIndexMap& tracer_indices,
                               int ncons,
                               const std::vector<std::string>& required_keys);

/**
 * \class BiologyModel
 * 
 * \brief Abstract base class for all biology model plugins
 *
 * Derived classes (Fennel, ECB, ECOSIM, etc.) implement the pure virtual
 * methods to compute their model-specific tendencies and diagnostics.
 */
class BiologyModel {
public:
    virtual ~BiologyModel() = default;

    /**
     * \brief Compute tracer source/sink tendencies for one tile
     *
     * This is the main interface called by the time integrator.
     * Called once per tile per time step, after advection/mixing,
     * before vertical viscosity.
     *
     * \param[in]     tile              Tile index
     * \param[in]     ng                Grid number
     * \param[in]     bx                Tile box
     * \param[in]     state             Temperature, salinity, tracers [ncons]
     * \param[in]     Hz                Vertical grid spacing
     * \param[in]     grid_metrics      pm, pn, pmon_u, pnom_v, etc.
     * \param[in]     masks             rmask, rmask_wet, etc.
     * \param[in]     forcing           Wind stress, PAR, etc. (if needed by model)
     * \param[out]    rhs_tracers       RHS for each tracer (accumulated)
     * \param[in]     solverChoice      Runtime options and flags
     */
    virtual void ComputeTendencies(
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
    ) = 0;

    /**
     * \brief Compute optional diagnostics
     *
     * Model-specific diagnostic fields (e.g., GPP, NPP, grazing rate).
     * By default, no diagnostics are computed.
     *
     * \param[in]  tile             Tile index
     * \param[in]  bx               Tile box
     * \param[in]  state            Tracer state
     * \param[out] diag_output      Diagnostic fields (model-specific)
     */
    virtual void ComputeDiagnostics(
        const int& tile,
        const amrex::Box& bx,
        const amrex::Array4<const amrex::Real>& state,
        amrex::Array4<amrex::Real>& diag_output
    ) {}

    /**
     * \brief Compute and report conservation budgets
     *
     * Optional: compute nitrogen, carbon, phosphorus inventories and
     * net changes per time step. Useful for validation and debugging.
     *
     * \param[in]  state             Tracer state
     * \param[out] budgets           Map of budget names -> values
     */
    virtual void ComputeConservationBudget(
        const amrex::Array4<const amrex::Real>& state,
        std::map<std::string, amrex::Real>& budgets
    ) {}

    /**
     * \brief Return the number of biology tracers used by this model
     */
    virtual int GetNumTracers() const = 0;

    /**
     * \brief Return the number of optional diagnostics computed by this model
     */
    virtual int GetNumDiagnostics() const { return 0; }

protected:
    void AddScaffoldBudgetFields(
        const std::string& model_name,
        int ncons,
        const TracerIndexMap& tracer_indices,
        long long calls,
        const std::vector<std::string>& required_tracer_keys,
        const std::vector<std::string>& optional_tracer_keys,
        bool missing_required,
        std::map<std::string, amrex::Real>& budgets) const;
};

/**
 * \class BiologyFactory
 * 
 * \brief Factory for creating biology model instances
 *
 * Static factory that maps BiologyModelType enum to concrete model classes.
 * Responsible for runtime selection and instantiation of biology plugins.
 */
class BiologyFactory {
public:
    /**
     * \brief Create a biology model instance based on type
     *
     * \param[in] model_type    BiologyModelType enum value
     * \param[in] ncons         Number of conservative tracers
     *
     * \return Unique pointer to the created model (or nullptr if none)
     */
    static std::unique_ptr<BiologyModel> CreateModel(
        const std::string& model_name,
        int ncons,
        const TracerIndexMap& tracer_indices
    );
};

} // namespace REMORA_Biology

#endif // REMORA_BIOLOGY_H_
