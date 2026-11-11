/** Exercise bottom-exchange coupling APIs through a complete REMORA run.
 * The inputs supply native AMReX geometry, as an embedding application does.
 */
#include <REMORA.H>
#include <AMReX.H>
#include <AMReX_ParmParse.H>

#include <cmath>
#include <limits>

using namespace amrex;
std::string inputs_name;

/** Initialize the model, submit owned snapshots and check exported bottom state.
 * @param[in   ] argc argument count
 * @param[in   ] argv input file and parameter overrides
 */
int
main (int argc, char** argv)
{
    inputs_name = argc > 1 ? argv[1] : "";
    amrex::Initialize(argc, argv);
    {
        ParmParse pp("remora"), coupling("coupler"), check("check");
        REMORA ocean;
        ocean.InitData();
        Real q = Real(0);
        Real valid_end = Real(100);
        Real reverse_time = Real(0);
        pp.get("bottom_q", q);
        coupling.query("valid_end", valid_end);
        coupling.query("reverse_time", reverse_time);
        std::string type;
        pp.get("bottom_exchange_type", type);
        BoxArray layout;
        DistributionMapping distribution;
        ocean.GetBottomExchangeLayout(layout, distribution);
        MultiFab water(layout, distribution, 1, 0);
        const auto& names = ocean.GetTracerNames();
        MultiFab incoming(layout, distribution, names.size(), 0);
        if (type == "coupler") {
            water.setVal(q);
            if (reverse_time > Real(0)) {
                ocean.PackBottomState(0, incoming);
            } else {
                for (int n = 0; n < names.size(); ++n) {
                    Real value = std::numeric_limits<Real>::quiet_NaN();
                    pp.query(("bottom_" + names[n]).c_str(), value);
                    incoming.setVal(value, n, 1, incoming.nGrowVect());
                }
            }
            ocean.SetBottomExchangeData(water, incoming, ocean.get_t_new(0), valid_end);
            // The producer can immediately reuse its storage.
            water.setVal(Real(123));
            incoming.setVal(Real(123));
        }
        Real ring_scale = Real(1);
        check.query("rho_ring_scale", ring_scale);
        if (type == "netcdf" && ring_scale != Real(1)) {
            Vector<Real> lo, hi;
            Vector<int> cells;
            pp.getarr("prob_lo", lo);
            pp.getarr("prob_hi", hi);
            pp.getarr("n_cell", cells);
            const Real area = (hi[0]-lo[0]) * (hi[1]-lo[1]) / Real(cells[0]*cells[1]);
            const auto* transport = ocean.GetBottomExchangeTransport(0);
            AMREX_ALWAYS_ASSERT(transport);
            for (MFIter mfi(*transport); mfi.isValid(); ++mfi) {
                const Box edge = makeSlab(mfi.fabbox(), 1, -1) & mfi.fabbox();
                if (!edge.ok()) { continue; }
                const Real expected = ring_scale * q * area;
                const Real tolerance = Real(32) * std::numeric_limits<Real>::epsilon() *
                                       std::abs(expected);
                const Real minimum = (*transport)[mfi].min<RunOn::Device>(edge, 0);
                const Real maximum = (*transport)[mfi].max<RunOn::Device>(edge, 0);
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(std::abs(minimum-expected) <= tolerance &&
                    std::abs(maximum-expected) <= tolerance, "Supplied NetCDF rho ring changed");
            }
        }
        int steps = 0;
        pp.get("max_step", steps);
        Real time = ocean.get_t_new(0);
        bool reversed = false;
        for (int step = 0; step < steps; ++step) {
            if (reverse_time > Real(0) && time >= reverse_time && !reversed) {
                water.setVal(-q);
                ocean.PackBottomState(0, incoming);
                ocean.SetBottomExchangeData(water, incoming, time, valid_end);
                reversed = true;
            }
            const Real advanced = ocean.EvolveOneStep(time, Real(0));
            if (advanced <= Real(0)) { break; }
            time += advanced;
        }
        ocean.WriteAtFinalTime();
        const int temp = ocean.component_index("temp");
        Real expected_temp = Real(0);
        pp.get("T0", expected_temp);
        for (int lev = 0; lev <= ocean.finestLevel(); ++lev) {
            MultiFab bottom(ocean.vec_mskr[lev]->boxArray(), ocean.cons_new[lev]->DistributionMap(),
                            names.size(), 0);
            ocean.PackBottomState(lev, bottom);
            const Real tolerance = Real(32) * std::numeric_limits<Real>::epsilon() *
                                   std::abs(expected_temp);
            AMREX_ALWAYS_ASSERT(std::abs(bottom.min(temp)-expected_temp) <= tolerance);
            AMREX_ALWAYS_ASSERT(std::abs(bottom.max(temp)-expected_temp) <= tolerance);
            const auto* transport = ocean.GetBottomExchangeTransport(lev);
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(transport && !transport->contains_nan(),
                                             "Uninitialized bottom transport halo");
        }
        Print() << "Bottom exchange APIs passed\n";
    }
    amrex::Finalize();
}
