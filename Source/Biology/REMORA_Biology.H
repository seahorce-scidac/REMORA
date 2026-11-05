#ifndef REMORA_BIOLOGY_H_
#define REMORA_BIOLOGY_H_

#include <algorithm>
#include <map>
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
    constant,
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
    const amrex::MultiFab* bottom_stress_u = nullptr;
    const amrex::MultiFab* bottom_stress_v = nullptr;
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
    // Identifies equations/options with persistent state semantics. Rate changes
    // may be allowed at restart; changing inventory meanings requires migration.
    virtual std::string checkpoint_identity () const { return name(); }
    virtual bool accepts_legacy_checkpoint () const { return false; }
    virtual void configure (const std::string& remora_prefix) = 0;
    virtual const BiologyDeclaration& declaration () const noexcept = 0;
    virtual void advance (const BiologyAdvanceContext& context) const = 0;
};

/** \brief Construct an implemented ecosystem; none/off yield a no-op object. */
std::unique_ptr<Biology> make_biology (const std::string& name);
BiologyICType parse_biology_ic_type (const std::string& name);
std::string biology_ic_type_name (BiologyICType type);

/** \brief Constant initial concentration of each biology tracer, by name. */
using BiologyConstantIC = std::map<std::string, amrex::Real>;

/** \brief Read one constant initial concentration per declared tracer, by name.
 *
 * Each value comes from <remora_prefix>.biology_initial_<tracer name>, in the
 * model's units, in the style of the named constant forcing values. There are
 * no defaults: every declared tracer must be given, and the error lists any
 * that are missing.
 */
BiologyConstantIC
read_constant_biology_ic (const BiologyDeclaration& declaration,
                          const std::string& remora_prefix);

/** \brief Set each declared biology component to its constant value, found by
 *         tracer name, ghost cells included. Every declared tracer must have a
 *         value. Physical fields and dyes are left unchanged. */
void initialize_constant_biology (const BiologyDeclaration& declaration,
                                  BiologyInitContext& context,
                                  const BiologyConstantIC& values);

}
#endif
