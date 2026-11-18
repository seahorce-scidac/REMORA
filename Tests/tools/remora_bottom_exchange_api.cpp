/** Exercise bottom-exchange coupling APIs through a complete REMORA run.
 * The inputs supply native AMReX geometry, as an embedding application does.
 */
#include <REMORA.H>
#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_Reduce.H>

#include <cmath>
#include <limits>

using namespace amrex;
std::string inputs_name;

/** Check water and tracer transport sums over each refined parent, and dry children.
 * @param[in   ] ocean initialized model with prepared bottom exchange
 */
void
check_coast (REMORA& ocean)
{
    int partial = 0;
    for (int lev = 1; lev <= ocean.finestLevel(); ++lev) {
        const auto& fine = *ocean.GetBottomExchangeTransport(lev);
        const auto& parent = *ocean.GetBottomExchangeTransport(lev-1);
        const IntVect ratio = ocean.refRatio(lev-1);
        const Real children = Real(ratio[0]*ratio[1]);
        MultiFab averaged(amrex::coarsen(fine.boxArray(), ratio), fine.DistributionMap(),
                          fine.nComp(), 0);
        MultiFab expected(averaged.boxArray(), averaged.DistributionMap(), fine.nComp(), 0);
        average_down(fine, averaged, 0, fine.nComp(), ratio);
        averaged.mult(children);
        expected.ParallelCopy(parent);
        MultiFab::Subtract(averaged, expected, 0, 0, fine.nComp(), 0);
        for (int n = 0; n < fine.nComp(); ++n) {
            const Real tolerance = Real(128)*std::numeric_limits<Real>::epsilon()*
                                   amrex::max(Real(1), expected.norm0(n));
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(averaged.norm0(n) <= tolerance,
                                             "Refined bottom transport sum changed");
        }
        MultiFab wet(averaged.boxArray(), averaged.DistributionMap(), 1, 0);
        average_down(*ocean.vec_mskr[lev], wet, 0, 1, ratio);
        ReduceOps<ReduceOpSum> op;
        ReduceData<int> data(op);
        for (MFIter mfi(wet, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const auto mask = wet.const_array(mfi);
            op.eval(mfi.tilebox(), data, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                -> GpuTuple<int> {
                return {mask(i,j,k) > Real(0) && mask(i,j,k) < Real(1) ? 1 : 0};
            });
        }
        partial += get<0>(data.value(op));
        MultiFab dry(fine.boxArray(), fine.DistributionMap(), fine.nComp(), 0);
        for (MFIter mfi(dry, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const auto out = dry.array(mfi);
            const auto in = fine.const_array(mfi);
            const auto mask = ocean.vec_mskr[lev]->const_array(mfi);
            ParallelFor(mfi.tilebox(), fine.nComp(),
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                out(i,j,k,n) = mask(i,j,k) == Real(0) ? in(i,j,k,n) : Real(0);
            });
        }
        for (int n = 0; n < fine.nComp(); ++n) {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(dry.norm0(n) == Real(0),
                                             "Bottom transport reached a dry child");
        }
    }
    ParallelDescriptor::ReduceIntSum(partial);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(partial > 0, "Coast test has no partly wet parents");
    Print() << "Bottom transport sums and dry children passed (" << partial
            << " partly wet parents)\n";
}

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
            for (int lev = 0; lev <= ocean.finestLevel(); ++lev) {
                const auto dx = ocean.Geom(lev).CellSizeArray();
                const Real area = dx[0]*dx[1];
                const auto* transport = ocean.GetBottomExchangeTransport(lev);
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
        }
        bool coast = false;
        check.query("coast", coast);
        if (coast) { check_coast(ocean); }
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
            if (coast) { check_coast(ocean); }
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
            for (MFIter mfi(bottom, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                const auto out = bottom.array(mfi);
                const auto mask = ocean.vec_mskr[lev]->const_array(mfi);
                ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    out(i,j,k,temp) -= expected_temp * mask(i,j,k);
                });
            }
            if (!coast) { AMREX_ALWAYS_ASSERT(bottom.norm0(temp) <= tolerance); }
            AMREX_ALWAYS_ASSERT(!bottom.contains_nan());
            const auto* transport = ocean.GetBottomExchangeTransport(lev);
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(transport && !transport->contains_nan(),
                                             "Uninitialized bottom transport halo");
        }
        Print() << "Bottom exchange APIs passed\n";
    }
    amrex::Finalize();
}
