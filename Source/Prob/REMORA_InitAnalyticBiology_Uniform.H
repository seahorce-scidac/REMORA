#ifndef REMORA_INIT_ANALYTIC_BIOLOGY_UNIFORM_H_
#define REMORA_INIT_ANALYTIC_BIOLOGY_UNIFORM_H_
#include <REMORA_Biology.H>

namespace REMORABiology {

// Explicit application profile, valid for any declared ecosystem inventory.
// There are no scientific defaults; values are in declaration order/units.
inline void
initialize_uniform_biology(
  const Biology& biology,
  BiologyInitContext& context,
  const amrex::Vector<amrex::Real>& values)
{
    const int count =
      static_cast<int>(biology.declaration().tracer_names.size());
    if (static_cast<int>(values.size()) != count) {
        amrex::Abort("remora.biology_initial must contain one concentration "
                     "per declared biology tracer");
    }
    for (int n = 0; n < count; ++n) {
        context.state.setVal(values[n], context.component_start + n, 1, 0);
    }
}

} // namespace REMORABiology
#endif
