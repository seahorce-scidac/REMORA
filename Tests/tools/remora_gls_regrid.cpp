#include <REMORA.H>
#include <AMReX_MultiFabUtil.H>

using namespace amrex;

std::string inputs_name;

namespace {
void test_defaults ()
{
    ParmParse amr("amr"), remora("remora"), geometry("geometry");
    Vector<int> cells, periodic;
    Vector<Real> lo, hi;
    remora.getarr("n_cell", cells);
    remora.getarr("is_periodic", periodic);
    remora.getarr("prob_lo", lo);
    remora.getarr("prob_hi", hi);
    amr.addarr("n_cell", cells);
    // Keep the initial patch in a single box and allow the small hand-built layouts below;
    // main.cpp's production defaults would round or split these test boxes.
    amr.add("max_grid_size", 2048);
    amr.add("max_grid_size_z", cells[2]);
    amr.add("blocking_factor", 1);
    amr.add("n_error_buf", 0);
    amr.add("n_proper", 2);
    amr.add("refine_grid_layout_x", 1);
    amr.add("refine_grid_layout_y", 1);
    amr.add("refine_grid_layout_z", 0);
    geometry.addarr("prob_lo", lo);
    geometry.addarr("prob_hi", hi);
    geometry.addarr("is_periodic", periodic);
#ifdef REMORA_USE_GPU
    Vector<int> tile_size = {1024000, 1024000, 1024000};
#else
    Vector<int> tile_size = {16, 16, 1024000};
#endif
    ParmParse("fabarray").addarr("mfiter_tile_size", tile_size);
}

Vector<MultiFab*> fields (REMORA& model, int lev)
{
    Vector<MultiFab*> result;
    for (const auto& field : model.gls_fields()) {
        result.push_back((*field.levels)[lev].get());
    }
    return result;
}

void check_periodic_ghosts (const MultiFab& mf, const Geometry& geom)
{
    AMREX_ALWAYS_ASSERT(!mf.contains_nan(0, mf.nComp(), mf.nGrowVect()));
    AMREX_ALWAYS_ASSERT(!mf.contains_inf(0, mf.nComp(), mf.nGrowVect()));
    MultiFab exchanged(mf.boxArray(), mf.DistributionMap(), mf.nComp(), mf.nGrowVect());
    MultiFab::Copy(exchanged, mf, 0, 0, mf.nComp(), mf.nGrowVect());
    exchanged.FillBoundary(geom.periodicity());
    MultiFab::Subtract(exchanged, mf, 0, 0, mf.nComp(), mf.nGrowVect());
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        exchanged.norm0(0, mf.nComp(), mf.nGrowVect()) == Real(0),
        "New GLS level has ghost cells that differ from their periodic images");
}

// Predictor endpoints enter vertical advection without being rebuilt by gls_prestep.
void set_predictor_endpoints (MultiFab& mf, const Geometry& geom,
                             Real bottom_value, Real surface_value)
{
    const int bottom = geom.Domain().smallEnd(2);
    const int surface = geom.Domain().bigEnd(2) + 1;
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
        auto values = mf.array(mfi);
        ParallelFor(mfi.growntilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            if (k == bottom) { values(i,j,k,2) = bottom_value; }
            if (k == surface) { values(i,j,k,2) = surface_value; }
        });
    }
}

void check_equal (const MultiFab& actual, const MultiFab& expected, int ncomp)
{
    AMREX_ALWAYS_ASSERT(!actual.contains_nan(0, actual.nComp(), 0));
    AMREX_ALWAYS_ASSERT(!actual.contains_inf(0, actual.nComp(), 0));
    MultiFab difference(actual.boxArray(), actual.DistributionMap(), ncomp, 0);
    MultiFab::Copy(difference, actual, 0, 0, ncomp, 0);
    MultiFab::Subtract(difference, expected, 0, 0, ncomp, 0);
    for (int comp = 0; comp < ncomp; ++comp) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(difference.norm0(comp) == Real(0),
                                        "GLS regrid changed retained state or coarse slot mapping");
    }
}
}

int main (int argc, char* argv[])
{
#ifdef AMREX_USE_MPI
    MPI_Init(&argc, &argv);
#endif
    amrex::Initialize(argc, argv, true, MPI_COMM_WORLD, test_defaults);
    inputs_name = argv[1];
    {
        REMORA model;
        model.InitData();
        const auto descriptors = model.gls_fields();
        {
            Real time = model.EvolveOneStep(Real(0), Real(1));
            auto coarse = fields(model, 0);
            auto fine = fields(model, 1);
            Vector<std::unique_ptr<MultiFab>> saved(descriptors.size());
            for (int f = 0; f < descriptors.size(); ++f) {
                // After one coarse step, coarse slot 1 is latest. Subcycling has fine parity 0.
                for (int n = 0; n < coarse[f]->nComp(); ++n) {
                    coarse[f]->setVal(Real(0.02)*(n+1), n, 1, coarse[f]->nGrowVect());
                    fine[f]->setVal(Real(0.003)*(n+1), n, 1, fine[f]->nGrowVect());
                }
                if (descriptors[f].has_history) {
                    set_predictor_endpoints(*coarse[f], model.Geom(0), Real(0.07), Real(0.08));
                    set_predictor_endpoints(*fine[f], model.Geom(1), Real(0.011), Real(0.013));
                }
                saved[f] = std::make_unique<MultiFab>(fine[f]->boxArray(), fine[f]->DistributionMap(),
                                                     fine[f]->nComp(), 0);
                MultiFab::Copy(*saved[f], *fine[f], 0, 0, fine[f]->nComp(), 0);
            }
            // A remake of the same layout must retain all three GLS components.
            model.RemakeLevel(1, time, model.boxArray(1), model.DistributionMap(1));
            fine = fields(model, 1);
            for (int f = 0; f < descriptors.size(); ++f) {
                check_equal(*fine[f], *saved[f], fine[f]->nComp());
            }
            // Enlarge and move the patch: retain some cells, remove others, and add new cells.
            Box box = model.Geom(1).Domain();
            box.setSmall(0, 4);
            box.setBig(0, 27);
            box.setSmall(1, 4);
            box.setBig(1, 35);
            BoxArray changed(box);
            changed.maxSize(IntVect(12, 16, box.length(2)));
            DistributionMapping mapping(changed);
            AMREX_ALWAYS_ASSERT(changed != model.boxArray(1));
            model.RemakeLevel(1, time, changed, mapping);
            model.SetBoxArray(1, changed);
            model.SetDistributionMap(1, mapping);
            fine = fields(model, 1);
            for (int f = 0; f < descriptors.size(); ++f) {
                int count = fine[f]->nComp();
                MultiFab expected(fine[f]->boxArray(), fine[f]->DistributionMap(), count, 0);
                for (int n = 0; n < count; ++n) {
                    Real value = descriptors[f].has_history && n < 2 ? Real(0.04) : Real(0.02)*(n+1);
                    expected.setVal(value, n, 1, 0);
                }
                if (descriptors[f].has_history) {
                    set_predictor_endpoints(expected, model.Geom(1), Real(0.07), Real(0.08));
                }
                expected.ParallelCopy(*saved[f], 0, 0, count, 0, 0, model.Geom(1).periodicity());
                check_equal(*fine[f], expected, count);
            }
            // Level creation must use the coarse latest slot even when fine parity differs.
            model.ClearLevel(1);
            model.MakeNewLevelFromCoarse(1, time, changed, mapping);
            fine = fields(model, 1);
            for (int f = 0; f < descriptors.size(); ++f) {
                int count = fine[f]->nComp();
                MultiFab expected(fine[f]->boxArray(), fine[f]->DistributionMap(), count, 0);
                for (int n = 0; n < count; ++n) {
                    Real value = descriptors[f].has_history && n < 2 ? Real(0.04) : Real(0.02)*(n+1);
                    expected.setVal(value, n, 1, 0);
                }
                if (descriptors[f].has_history) {
                    set_predictor_endpoints(expected, model.Geom(1), Real(0.07), Real(0.08));
                }
                check_equal(*fine[f], expected, count);
            }
            // Coarse interpolation is bounded, while old fine values below the bounds survive.
            for (int f = 0; f < descriptors.size(); ++f) {
                coarse[f]->setVal(Real(-1));
                fine[f]->setVal(Real(-2));
            }
            model.RemakeLevel(1, time, changed, mapping);
            fine = fields(model, 1);
            for (int f = 0; f < descriptors.size(); ++f) {
                for (int n = 0; n < fine[f]->nComp(); ++n) {
                    AMREX_ALWAYS_ASSERT(fine[f]->min(n) == Real(-2));
                    AMREX_ALWAYS_ASSERT(fine[f]->max(n) == Real(-2));
                }
            }
            model.ClearLevel(1);
            model.MakeNewLevelFromCoarse(1, time, changed, mapping);
            fine = fields(model, 1);
            for (int f = 0; f < descriptors.size(); ++f) {
                AMREX_ALWAYS_ASSERT(!fine[f]->contains_nan(0, fine[f]->nComp(), 0));
                for (int n = 0; n < fine[f]->nComp(); ++n) {
                    const auto& bounds = descriptors[f].lower;
                    Real bound = bounds.size() == 1 ? bounds[0] : bounds[n];
                    AMREX_ALWAYS_ASSERT(fine[f]->min(n) >= bound);
                }
            }
            // Smooth data exposes tangential nodal interpolation errors at a periodic seam.
            const int nx = model.Geom(0).Domain().length(0);
            for (auto* mf : coarse) {
                for (MFIter mfi(*mf); mfi.isValid(); ++mfi) {
                    auto values = mf->array(mfi);
                    ParallelFor(mfi.growntilebox(), mf->nComp(),
                        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
                    {
                        values(i,j,k,n) = Real(2) + std::sin(Real(2)*Real(3.141592653589793)*
                            (Real(i)+Real(0.5))/nx) + Real(0.01)*(j+k+n);
                    });
                }
            }
            Box seam_box = model.Geom(1).Domain();
            seam_box.setSmall(1, 4);
            seam_box.setBig(1, 35);
            BoxArray seam(seam_box);
            DistributionMapping seam_mapping(seam);
            model.ClearLevel(1);
            model.MakeNewLevelFromCoarse(1, time, seam, seam_mapping);
            model.SetBoxArray(1, seam);
            model.SetDistributionMap(1, seam_mapping);
            fine = fields(model, 1);
            for (auto* mf : fine) {
                check_periodic_ghosts(*mf, model.Geom(1));
            }
            amrex::Print() << "GLS changing-layout transfer, history slots and bounds passed\n";
        }
        model.Evolve();
    }
    amrex::Finalize();
#ifdef AMREX_USE_MPI
    MPI_Finalize();
#endif
    return 0;
}
