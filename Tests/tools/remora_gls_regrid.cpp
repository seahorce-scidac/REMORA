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

std::array<MultiFab*, 7> fields (REMORA& model, int lev)
{
    return {model.vec_tke[lev].get(), model.vec_gls[lev].get(),
            model.vec_Lscale[lev].get(), model.vec_Akv[lev].get(),
            model.vec_Akt[lev].get(), model.vec_Akk[lev].get(), model.vec_Akp[lev].get()};
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
        bool skip_transfer = false;
        ParmParse test("test");
        test.query("skip_transfer", skip_transfer);
        if (!skip_transfer) {
            Real time = model.EvolveOneStep(Real(0), Real(1));
            auto coarse = fields(model, 0);
            auto fine = fields(model, 1);
            std::array<std::unique_ptr<MultiFab>, 7> saved;
            for (int f = 0; f < 7; ++f) {
                // After one coarse step, coarse slot 1 is latest. Subcycling has fine parity 0.
                for (int n = 0; n < coarse[f]->nComp(); ++n) {
                    coarse[f]->setVal(Real(0.02)*(n+1), n, 1, coarse[f]->nGrowVect());
                    fine[f]->setVal(Real(0.003)*(n+1), n, 1, fine[f]->nGrowVect());
                }
                saved[f] = std::make_unique<MultiFab>(fine[f]->boxArray(), fine[f]->DistributionMap(),
                                                     fine[f]->nComp(), 0);
                MultiFab::Copy(*saved[f], *fine[f], 0, 0, fine[f]->nComp(), 0);
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
            for (int f = 0; f < 7; ++f) {
                int count = f < 2 ? 2 : fine[f]->nComp();
                MultiFab expected(fine[f]->boxArray(), fine[f]->DistributionMap(), count, 0);
                for (int n = 0; n < count; ++n) {
                    expected.setVal(f < 2 ? Real(0.04) : Real(0.02)*(n+1), n, 1, 0);
                }
                expected.ParallelCopy(*saved[f], 0, 0, count, 0, 0, model.Geom(1).periodicity());
                check_equal(*fine[f], expected, count);
                if (f < 2) {
                    int subcycle = 0;
                    ParmParse("remora").query("do_substep", subcycle);
                    MultiFab latest(*fine[f], make_alias, subcycle ? 0 : 1, 1);
                    MultiFab scratch(*fine[f], make_alias, 2, 1);
                    check_equal(scratch, latest, 1);
                }
            }
            // Level creation must use the coarse latest slot even when fine parity differs.
            model.ClearLevel(1);
            model.MakeNewLevelFromCoarse(1, time, changed, mapping);
            fine = fields(model, 1);
            for (int f = 0; f < 7; ++f) {
                int count = f < 2 ? 2 : fine[f]->nComp();
                MultiFab expected(fine[f]->boxArray(), fine[f]->DistributionMap(), count, 0);
                for (int n = 0; n < count; ++n) {
                    expected.setVal(f < 2 ? Real(0.04) : Real(0.02)*(n+1), n, 1, 0);
                }
                check_equal(*fine[f], expected, count);
                if (f < 2) {
                    int subcycle = 0;
                    ParmParse("remora").query("do_substep", subcycle);
                    MultiFab latest(*fine[f], make_alias, subcycle ? 0 : 1, 1);
                    MultiFab scratch(*fine[f], make_alias, 2, 1);
                    check_equal(scratch, latest, 1);
                }
            }
            // Coarse interpolation is bounded, while old fine values below the bounds survive.
            const std::array<const char*, 7> keys = {"gls_Kmin", "gls_Pmin", "", "Akv_bak",
                                                   "Akt_bak", "Akk_bak", "Akp_bak"};
            std::array<Real, 7> lower = {Real(7.6e-6), Real(1.e-12), Real(0), Real(1.e-5),
                                         Real(1.e-6), Real(5.e-6), Real(5.e-6)};
            ParmParse remora("remora");
            for (int f = 0; f < 7; ++f) {
                if (f != 2) { remora.query(keys[f], lower[f]); }
                coarse[f]->setVal(Real(-1));
                fine[f]->setVal(Real(-2));
            }
            model.RemakeLevel(1, time, changed, mapping);
            fine = fields(model, 1);
            for (int f = 0; f < 7; ++f) {
                for (int n = 0; n < fine[f]->nComp(); ++n) {
                    AMREX_ALWAYS_ASSERT(fine[f]->min(n) == Real(-2));
                    AMREX_ALWAYS_ASSERT(fine[f]->max(n) == Real(-2));
                }
            }
            model.ClearLevel(1);
            model.MakeNewLevelFromCoarse(1, time, changed, mapping);
            fine = fields(model, 1);
            for (int f = 0; f < 7; ++f) {
                AMREX_ALWAYS_ASSERT(!fine[f]->contains_nan(0, fine[f]->nComp(), 0));
                for (int n = 0; n < fine[f]->nComp(); ++n) {
                    Real bound = lower[f];
                    if (f == 4) {
                        remora.query(n == 0 ? "Akt_bak_temp" : "Akt_bak_salt", bound);
                    }
                    AMREX_ALWAYS_ASSERT(fine[f]->min(n) >= bound);
                }
            }
            for (int f = 0; f < 7; ++f) {
                coarse[f]->setVal(lower[f]);
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
