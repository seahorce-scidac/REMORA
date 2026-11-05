#include <REMORA_Biology.H>
#include <REMORA_ECBBiology.H>
#include <REMORA_DataStruct.H>

#include <AMReX_Gpu.H>
#include <AMReX_MFIter.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <cmath>
#include <fstream>
#include <limits>

using namespace amrex;
using namespace REMORABiology;

namespace {

// Fixtures are produced by the pinned ROMS kernel, not the native port.
// generate_reference.py runs every case with this density and ROMS Cp.
constexpr Real reference_rho0 = Real(1025.0);
static_assert(Cp == Real(3985.0), "generate_reference.py uses Cp = 3985");

// With no remora.ecb.* parameters set, ECB configures from the ROMS
// bio_ECB.in defaults and passes its own parameter checks.
void check_ecb_defaults ()
{
    auto model = make_biology("ecb");
    model->configure("ecb_defaults");
    AMREX_ALWAYS_ASSERT(model->name() == "ecb");
    AMREX_ALWAYS_ASSERT(model->declaration().tracer_names.size() == 18);

    ECBParameters parameters;
    parameters.init_params("ecb_defaults");
    AMREX_ALWAYS_ASSERT(parameters.BioIter == 1);
    AMREX_ALWAYS_ASSERT(parameters.Vp0 == Real(0.55));
    AMREX_ALWAYS_ASSERT(parameters.wLDet == Real(5.0));
    AMREX_ALWAYS_ASSERT(parameters.pCO2air == Real(370.0));
}

}

void check_ecb_reference ()
{
    check_ecb_defaults();

    std::ifstream input(REMORA_ECB_REFERENCE_FILE);
    AMREX_ALWAYS_ASSERT(input.good());
    int count;
    input >> count;
    Vector<std::string> names(count);
    Vector<Real> values(count);
    for (auto& name : names) {
        input >> name;
    }
    for (auto& value : values) {
        input >> value;
    }
    int cases;
    input >> cases;
    Real worst = 0;
    for (int test = 0; test < cases; ++test) {
        int bulk, nz, iterations;
        Real seconds, light, stress, bottom, wind, mask_value;
        input >> bulk >> nz >> iterations >> seconds >> light >> stress >> bottom >> wind >> mask_value;
        Vector<Real> thickness(nz), interfaces(nz+1), old_values(nz*20), incoming(nz*20), expected(nz*20);
        for (auto* array : {&thickness, &interfaces, &old_values, &incoming, &expected}) {
            for (auto& value : *array) {
                input >> value;
            }
        }
        AMREX_ALWAYS_ASSERT(input.good());

        const std::string prefix = "ecb_reference_" + std::to_string(test);
        ParmParse pp(prefix + ".ecb");
        pp.add("BioIter", iterations);
        for (int n = 0; n < count; ++n) {
            pp.add(names[n].c_str(), values[n]);
        }
        Vector<Real> initial(18);
        for (int n = 0; n < 18; ++n) {
            initial[n] = old_values[(n+2)*nz];
        }
        auto model = make_biology("ECB");
        model->configure(prefix);
        AMREX_ALWAYS_ASSERT(model->name() == "ecb");
        AMREX_ALWAYS_ASSERT(model->declaration().tracer_names.size() == 18);
        AMREX_ALWAYS_ASSERT(model->declaration().tracer_names[7] == "DON");
        AMREX_ALWAYS_ASSERT(model->declaration().tracer_names[15] == "Rdon");

        const Box domain(IntVect(0), IntVect(0,0,nz-1));
        const RealBox physical({0.,0.,-50.}, {1.,1.,0.});
        const int periodic[3] = {0,0,0};
        Geometry geometry(domain, &physical, 0, periodic);
        BoxArray boxes(domain);
        DistributionMapping mapping(boxes);
        const int dyes = test % 4;
        const int base = 2 + dyes;

        // Managed buffers also make this test executable in GPU builds.
        const MFInfo managed = MFInfo().SetArena(The_Managed_Arena());
        MultiFab old(boxes, mapping, base+18, 0, managed);
        MultiFab next(boxes, mapping, base+18, 0, managed);
        MultiFab hz(boxes, mapping, 1, 0, managed);
        BoxArray wboxes = amrex::convert(boxes, IntVect(0,0,1));
        MultiFab zw(wboxes, mapping, 1, 0, managed);
        BoxArray horizontal(Box(IntVect(0), IntVect(0)));
        MultiFab shortwave(horizontal, mapping, 1, 1, managed);
        MultiFab mask(horizontal, mapping, 1, 1, managed);
        MultiFab su(horizontal, mapping, 1, 1, managed);
        MultiFab sv(horizontal, mapping, 1, 1, managed);
        MultiFab bstr_u(horizontal, mapping, 1, 1, managed);
        MultiFab bstr_v(horizontal, mapping, 1, 1, managed);
        MultiFab wu(horizontal, mapping, 1, 1, managed);
        MultiFab wv(horizontal, mapping, 1, 1, managed);
        old.setVal(-17.);
        next.setVal(-29.);
        // columns.txt stores the kinematic srflx (degC m/s) that ROMS consumed;
        // REMORA's shortwave is in W/m2, so convert with the same rho0 and Cp.
        shortwave.setVal(light*reference_rho0*Cp);
        mask.setVal(mask_value);
        su.setVal(stress);
        sv.setVal(stress*Real(0.5));
        bstr_u.setVal(bottom);
        bstr_v.setVal(bottom*Real(0.5));
        wu.setVal(wind);
        wv.setVal(0.);
        Gpu::streamSynchronize();
        for (MFIter mfi(old); mfi.isValid(); ++mfi) {
            auto a = old.array(mfi), b = next.array(mfi), h = hz.array(mfi), z = zw.array(mfi);
            for (int k = 0; k < nz; ++k) {
                h(0,0,k) = thickness[k];
                for (int c = 0; c < 20; ++c) {
                    const int component = c < 2 ? c : base + c - 2;
                    a(0,0,k,component) = old_values[c*nz+k];
                    b(0,0,k,component) = incoming[c*nz+k];
                }
            }
            for (int k = 0; k <= nz; ++k) {
                z(0,0,k) = interfaces[k];
            }
        }

        SolverChoice solver;
        solver.bulk_fluxes = bulk != 0;
        solver.rho0 = reference_rho0;
        BiologyInitContext init{.level = 0,
                                .geometry = geometry,
                                .solver_choice = solver,
                                .state = next,
                                .component_start = base};

        // Set the constant initial values by tracer name, as an inputs file
        // would with remora.biology_ic_type = constant.
        ParmParse ic_pp(prefix);
        for (int n = 0; n < 18; ++n) {
            ic_pp.add(("biology_initial_" + model->declaration().tracer_names[n]).c_str(), initial[n]);
        }
        const BiologyConstantIC constant_ic = read_constant_biology_ic(model->declaration(), prefix);
        AMREX_ALWAYS_ASSERT(constant_ic.size() == 18);
        for (int n = 0; n < 18; ++n) {
            AMREX_ALWAYS_ASSERT(constant_ic.at(model->declaration().tracer_names[n]) == initial[n]);
        }
        initialize_constant_biology(model->declaration(), init, constant_ic);
        if (test == 0) {
            // Ghost cells are set too, and physical fields and dyes are not.
            MultiFab grown(boxes, mapping, base+18, 2, managed);
            grown.setVal(Real(-5.0));
            BiologyInitContext grown_init{.level = 0,
                                          .geometry = geometry,
                                          .solver_choice = solver,
                                          .state = grown,
                                          .component_start = base};
            initialize_constant_biology(model->declaration(), grown_init, constant_ic);
            for (int c = 0; c < 18; ++c) {
                AMREX_ALWAYS_ASSERT(grown.min(base+c, 2) == initial[c]);
                AMREX_ALWAYS_ASSERT(grown.max(base+c, 2) == initial[c]);
            }
            for (int c = 0; c < base; ++c) {
                AMREX_ALWAYS_ASSERT(grown.min(c, 2) == Real(-5.0));
                AMREX_ALWAYS_ASSERT(grown.max(c, 2) == Real(-5.0));
            }
        }
        Gpu::streamSynchronize();
        for (MFIter mfi(next); mfi.isValid(); ++mfi) {
            auto b = next.array(mfi);
            for (int k = 0; k < nz; ++k) {
                for (int c = 0; c < 18; ++c) {
                    AMREX_ALWAYS_ASSERT(b(0,0,k,base+c) == initial[c]);
                    b(0,0,k,base+c) = incoming[(c+2)*nz+k];
                }
            }
        }

        const BiologyAdvanceContext context{.old_state = old,
                                            .new_state = next,
                                            .solver_choice = solver,
                                            .component_start = base,
                                            .top_layer = nz - 1,
                                            .dt = seconds,
                                            .model_time = 0.,
                                            .hz = &hz,
                                            .z_w = &zw,
                                            .mask = &mask,
                                            .shortwave = &shortwave,
                                            .stress_u = &su,
                                            .stress_v = &sv,
                                            .wind_u = &wu,
                                            .wind_v = &wv,
                                            .bottom_stress_u = &bstr_u,
                                            .bottom_stress_v = &bstr_v};
        model->advance(context);
        Gpu::streamSynchronize();

        for (MFIter mfi(next); mfi.isValid(); ++mfi) {
            auto a = old.const_array(mfi), b = next.const_array(mfi);
            for (int k = 0; k < nz; ++k) {
                for (int c = 0; c < 20; ++c) {
                    const int component = c < 2 ? c : base + c - 2;
                    const Real actual = b(0,0,k,component);
                    const Real reference = expected[c*nz+k];
                    const Real error = std::abs(actual - reference) / amrex::max(Real(1.), std::abs(reference));
                    worst = amrex::max(worst, error);
                    if (!std::isfinite(actual) || error > Real(128.)*std::numeric_limits<Real>::epsilon()) {
                        Print() << "ECB case " << test << " component " << c << " layer " << k
                                << " actual " << actual << " reference " << reference
                                << " scaled error " << error << "\n";
                        amrex::Abort("ECB differs from clean ROMS reference");
                    }
                    AMREX_ALWAYS_ASSERT(a(0,0,k,component) == old_values[c*nz+k]);
                    if (c < 2 || c >= 17 || iterations == 0 || mask_value == 0) {
                        AMREX_ALWAYS_ASSERT(actual == incoming[c*nz+k]);
                    }
                }
                for (int dye = 2; dye < base; ++dye) {
                    AMREX_ALWAYS_ASSERT(b(0,0,k,dye) == Real(-29.));
                }
            }
        }
    }
    Print() << "ECB: " << cases << " clean-ROMS column fixtures; maximum scaled error " << worst << "\n";
}
