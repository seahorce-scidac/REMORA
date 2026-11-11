/** Closed-form flux identities and neutral-log-layer GLS scaling. */
#include <REMORA_BottomExchange.H>
#include <AMReX.H>
#include <AMReX_BLassert.H>
#include <AMReX_Print.H>

#include <algorithm>
#include <cmath>
#include <limits>

using namespace amrex;

/** Compare dimensionless arithmetic or fluxes in a chosen consistent unit system.
 * @param[in   ] actual helper result
 * @param[in   ] expected independent value
 */
void
require_close (Real actual, Real expected)
{
    const Real tolerance = Real(32) * std::numeric_limits<Real>::epsilon();
    AMREX_ALWAYS_ASSERT(std::abs(actual-expected) <=
        tolerance * std::max(std::abs(expected), Real(1.e-12)));
}

/** Check donor fluxes and GLS ratios without creating a model or a mesh.
 * @param[in   ] argc argument count
 * @param[in   ] argv input file and parameter overrides
 */
int
main (int argc, char** argv)
{
    amrex::Initialize(argc, argv);
    require_close(BottomExchange::tracer_flux(Real(2), Real(14), Real(3)), Real(14));
    require_close(BottomExchange::tracer_flux(Real(-2), Real(0), Real(3)), Real(-6));
    require_close(BottomExchange::tracer_flux(Real(0), Real(0), Real(3)), Real(0));
    const Real friction = Real(0.01);
    const Real speed = Real(1.e-5);
    const Real length = Real(1);
    const Real k = friction * friction / Real(0.3);
    const Real mixing = Real(0.4) * friction * length;
    const Real dissipation = friction * friction * friction / (Real(0.4) * length);
    const auto result = BottomExchange::scales(speed, length, mixing, mixing,
        mixing / Real(1.3), k, dissipation, friction * friction, Real(0.1));
    require_close(result.transport, Real(0.00325));
    require_close(result.adjustment, Real(1) / Real(750));
    require_close(result.wall, Real(0.01));
    Print() << "Bottom exchange arithmetic passed\n";
    amrex::Finalize();
}
