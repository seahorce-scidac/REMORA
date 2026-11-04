#ifndef REMORA_BIOLOGY_H_
#define REMORA_BIOLOGY_H_

#include <algorithm>
#include <memory>
#include <string>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>

struct SolverChoice;

namespace REMORABiology {

/** \brief Initial-condition source, independent of physical IC selection. */
enum class BiologyICType {
    follow_ic_type,
    analytic,
    netcdf
};

/** \brief Configured inventory and output aliases, frozen before scalar allocation. */
struct BiologyDeclaration {
    amrex::Vector<std::string> tracer_names;
    amrex::Vector<std::string> plot_aliases;

    /** \brief State component of tracer \p name in a block starting at
     *         \p component_start, or -1 if the model does not declare it. */
    int component_index (const std::string& name, int component_start) const {
        auto it = std::find(tracer_names.begin(), tracer_names.end(), name);
        if (it == tracer_names.end()) {
            return -1;
        }
        return component_start +
               static_cast<int>(it - tracer_names.begin());
    }
};

/** \brief What ProblemBase::init_analytic_biology receives.
 *
 * Analytic profiles are problem-specific and live with the problem in
 * Exec/REMORA_Prob.cpp, not in the ecosystem model; they locate their fields
 * with BiologyDeclaration::component_index. Physical fields have already been
 * initialized. Only the declared biology components, starting at
 * component_start, may be changed by the profile.
 */
struct BiologyInitContext {
    int level;
    const amrex::Geometry& geometry;
    const SolverChoice& solver_choice;
    amrex::MultiFab& state;
    int component_start;
};

/** \brief Host views at the existing biology source-update stage.
 *
 * old_state holds concentrations. new_state receives additive Hz-weighted
 * increments; it must not be overwritten with scratch concentrations. Optional
 * forcing views default to null and are validated by the concrete model.
 * Fill it with designated initializers so the pointers cannot be transposed
 * silently. This is a host context: device kernels capture plain
 * views/parameters, not this object.
 */
struct BiologyAdvanceContext {
    const amrex::MultiFab& old_state;
    amrex::MultiFab& new_state;
    const SolverChoice& solver_choice;
    int component_start;
    int top_layer;
    amrex::Real dt;
    double model_time;
    const amrex::MultiFab* hz = nullptr;
    const amrex::MultiFab* z_r = nullptr;
    const amrex::MultiFab* z_w = nullptr;
    const amrex::MultiFab* mask = nullptr;
    const amrex::MultiFab* shortwave = nullptr;
    const amrex::MultiFab* stress_u = nullptr;
    const amrex::MultiFab* stress_v = nullptr;
    const amrex::MultiFab* wind_u = nullptr;
    const amrex::MultiFab* wind_v = nullptr;
};

/** \brief Common host-side ecosystem interface.
 *
 * Configure once before querying the declaration and allocating scalar state.
 * Concrete models own parameters, option validation, and kernels. Analytic
 * initial profiles belong to the problem (ProblemBase::init_analytic_biology).
 * Runtime selection is confined to make_biology; callers do not downcast.
 * No additional clipping or mathematical constraints belong in this interface.
 */
class Biology {
public:
    virtual ~Biology () = default;
    virtual std::string name () const = 0;
    virtual void configure (const std::string& remora_prefix) = 0;
    virtual const BiologyDeclaration& declaration () const noexcept = 0;
    virtual void advance (const BiologyAdvanceContext& context) const = 0;
};

/** \brief Construct an implemented ecosystem; none/off yield a no-op object. */
std::unique_ptr<Biology> make_biology (const std::string& name);
BiologyICType parse_biology_ic_type (const std::string& name);
std::string biology_ic_type_name (BiologyICType type);

}
#endif
