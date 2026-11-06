#include <REMORA_Biology.H>
#include <REMORA_FennelBiology.H>
#include <REMORA_ECBBiology.H>
#include <REMORA_DataStruct.H>

#include <AMReX.H>
#include <AMReX_MFIter.H>
#include <AMReX_ParmParse.H>

#include <cmath>

using namespace amrex;
using namespace REMORABiology;

namespace {

// Apply the BioToy/Upwelling profile exactly as Problem::init_analytic_biology
// does: the fragment sees only the abstract model and the init context.
void apply_biotoy_profile (Biology const& biology, BiologyInitContext& context)
{
#include "Prob/REMORA_InitAnalyticBiology_BioToy.H"
}

// Every tracer has a state component base + its position in the declaration,
// and component_index finds it by name; an undeclared name gives -1.
void check_component_index (BiologyDeclaration const& declaration)
{
    for (int dyes : {0, 1, 3}) {
        const int base = Tracer_comp + dyes;
        for (int n = 0; n < int(declaration.tracer_names.size()); ++n) {
            AMREX_ALWAYS_ASSERT(declaration.component_index(declaration.tracer_names[n], base) ==
                                base + n);
        }
        AMREX_ALWAYS_ASSERT(declaration.component_index("not_a_tracer", base) == -1);
    }
}

// Fennel declares the seven nitrogen tracers, then each optional pool its
// options enable, in the ROMS order, for every combination of options.
void check_fennel_layouts ()
{
    for (int options = 0; options < 32; ++options) {
        const bool river_don = (options & 1) != 0;
        const bool po4 = (options & 2) != 0;
        const bool carbon = (options & 4) != 0;
        const bool oxygen = (options & 8) != 0;
        const bool odu = (options & 16) != 0;
        const std::string prefix = "layout_" + std::to_string(options);
        ParmParse pp(prefix + ".fennel");
        pp.add("river_don", int(river_don));
        pp.add("po4", int(po4));
        pp.add("carbon", int(carbon));
        pp.add("oxygen", int(oxygen));
        pp.add("odu", int(odu));
        auto model = make_biology("FENNEL");
        model->configure(prefix);
        auto const& declaration = model->declaration();
        AMREX_ALWAYS_ASSERT(model->name() == "fennel");
        AMREX_ALWAYS_ASSERT(declaration.plot_aliases == Vector<std::string>{"fennel"});

        Vector<std::string> expected{"NO3", "NH4", "chlorophyll", "phytoplankton",
                                     "zooplankton", "LdetritusN", "SdetritusN"};
        if (river_don) { expected.push_back("RdetritusN"); }
        if (po4) { expected.push_back("PO4"); }
        if (carbon) {
            for (const auto* name : {"LdetritusC", "SdetritusC", "TIC", "alkalinity"}) {
                expected.push_back(name);
            }
        }
        if (carbon && river_don) { expected.push_back("RdetritusC"); }
        if (oxygen) { expected.push_back("oxygen"); }
        if (odu) { expected.push_back("ODU"); }
        AMREX_ALWAYS_ASSERT(declaration.tracer_names == expected);
        check_component_index(declaration);
    }
}

// ECB always declares the same 18 tracers.
void check_ecb_layout ()
{
    auto model = make_biology("ECB");
    model->configure("ecb_layout");
    auto const& declaration = model->declaration();
    AMREX_ALWAYS_ASSERT(model->name() == "ecb");
    AMREX_ALWAYS_ASSERT(declaration.plot_aliases == Vector<std::string>{"ecb"});
    const Vector<std::string> expected{
        "NO3", "NH4", "chlorophyll", "phytoplankton", "zooplankton", "LdetritusN",
        "SdetritusN", "DON", "PO4", "LdetritusC", "SdetritusC", "TIC", "alkalinity",
        "DOC", "oxygen", "Rdon", "Rdoc", "calcium"};
    AMREX_ALWAYS_ASSERT(declaration.tracer_names == expected);
    check_component_index(declaration);
}

// With no parameters set, each model configures from its ROMS defaults and
// passes its own parameter checks.
void check_defaults ()
{
    FennelParameters fennel;
    fennel.init_params("fennel_defaults");
    AMREX_ALWAYS_ASSERT(fennel.BioIter == 1);
    AMREX_ALWAYS_ASSERT(fennel.Vp0 == Real(1.0));
    AMREX_ALWAYS_ASSERT(fennel.wLDet == Real(1.0));
    AMREX_ALWAYS_ASSERT(fennel.pCO2air == Real(370.0));

    ECBParameters ecb;
    ecb.init_params("ecb_defaults");
    AMREX_ALWAYS_ASSERT(ecb.BioIter == 1);
    AMREX_ALWAYS_ASSERT(ecb.Vp0 == Real(0.55));
    AMREX_ALWAYS_ASSERT(ecb.wLDet == Real(5.0));
    AMREX_ALWAYS_ASSERT(ecb.pCO2air == Real(370.0));
}

// Constant initial values are read by tracer name and written to each
// tracer's component, ghost cells included, leaving physical fields and dyes.
void check_constant_ic ()
{
    const Box domain(IntVect(0), IntVect(3));
    const RealBox physical({0., 0., -10.}, {10., 10., 0.});
    const int periodic[3] = {0, 0, 0};
    Geometry geometry(domain, &physical, 0, periodic);
    BoxArray boxes(domain);
    DistributionMapping mapping(boxes);
    SolverChoice solver;
    const int base = Tracer_comp + 2;

    for (const auto* model_name : {"fennel", "ecb"}) {
        const std::string prefix = std::string("constant_") + model_name;
        auto model = make_biology(model_name);
        model->configure(prefix);
        auto const& declaration = model->declaration();
        const int nbio = int(declaration.tracer_names.size());

        ParmParse pp(prefix);
        for (int n = 0; n < nbio; ++n) {
            pp.add(("biology_initial_" + declaration.tracer_names[n]).c_str(), Real(n + 1));
        }
        const BiologyConstantIC values = read_constant_biology_ic(declaration, prefix);
        AMREX_ALWAYS_ASSERT(int(values.size()) == nbio);

        MultiFab state(boxes, mapping, base + nbio, 2);
        state.setVal(Real(-5.0));
        BiologyInitContext context{.level = 0,
                                   .geometry = geometry,
                                   .solver_choice = solver,
                                   .state = state,
                                   .component_start = base};
        initialize_constant_biology(declaration, context, values);
        for (int n = 0; n < nbio; ++n) {
            const int comp = declaration.component_index(declaration.tracer_names[n], base);
            AMREX_ALWAYS_ASSERT(state.min(comp, 2) == Real(n + 1));
            AMREX_ALWAYS_ASSERT(state.max(comp, 2) == Real(n + 1));
        }
        for (int comp = 0; comp < base; ++comp) {
            AMREX_ALWAYS_ASSERT(state.min(comp, 2) == Real(-5.0));
            AMREX_ALWAYS_ASSERT(state.max(comp, 2) == Real(-5.0));
        }
    }
}

// The Fennel analytic profile writes the declared tracers, by name, and
// nothing else.
void check_fennel_profile ()
{
    const Box domain(IntVect(0), IntVect(3));
    const RealBox physical({0., 0., -10.}, {10., 10., 0.});
    const int periodic[3] = {0, 0, 0};
    Geometry geometry(domain, &physical, 0, periodic);
    BoxArray boxes(domain);
    DistributionMapping mapping(boxes);
    SolverChoice solver;
    const int base = Tracer_comp + 3;
    MultiFab state(boxes, mapping, base + 16, 0);
    state.setVal(Real(-7.0));
    state.setVal(Real(12.0), Temp_comp, 1, 0);
    state.setVal(Real(30.0), Salt_comp, 1, 0);
    BiologyInitContext context{.level = 0,
                               .geometry = geometry,
                               .solver_choice = solver,
                               .state = state,
                               .component_start = base};

    // Nitrogen-only layout: the profile writes the seven base tracers and
    // leaves every optional slot alone, since their names are not declared.
    // The options are set explicitly so a change to a FennelParameters default
    // cannot move the optional slots under this check.
    ParmParse nitrogen_pp("nitrogen.fennel");
    for (const auto* option : {"carbon", "oxygen", "po4", "odu", "river_don"}) {
        nitrogen_pp.add(option, 0);
    }
    auto nitrogen = make_biology("fennel");
    nitrogen->configure("nitrogen");
    apply_biotoy_profile(*nitrogen, context);
    // Temperature is 12, so the silicate proxy is 10 - 4*(12 - 11) = 6.
    const Real sio4 = Real(6.0);
    const Real no3 = Real(1.67) + Real(0.5873)*sio4 + Real(0.0144)*sio4*sio4 +
                     Real(0.0003099)*sio4*sio4*sio4;
    // Each base tracer in Fennel order, with the value ana_biology.h gives it.
    const Vector<Real> nitrogen_values{no3, Real(0.1), Real(0.02), Real(0.08),
                                       Real(0.06), Real(0.02), Real(0.04)};
    for (int n = 0; n < Fennel::NumBaseTracers; ++n) {
        const Real tol = Real(1.0e-12) * std::abs(nitrogen_values[n]);
        const Real expected = nitrogen_values[n];
        AMREX_ALWAYS_ASSERT(std::abs(state.min(base + n) - expected) <= tol);
        AMREX_ALWAYS_ASSERT(std::abs(state.max(base + n) - expected) <= tol);
    }
    for (int icomp = base + Fennel::NumBaseTracers; icomp < state.nComp(); ++icomp) {
        AMREX_ALWAYS_ASSERT(state.min(icomp) == Real(-7.0));
        AMREX_ALWAYS_ASSERT(state.max(icomp) == Real(-7.0));
    }

    ParmParse parameters("profile.fennel");
    for (const auto* option : {"carbon", "oxygen", "po4", "odu", "river_don"}) {
        parameters.add(option, 1);
    }
    auto model = make_biology("fennel");
    model->configure("profile");
    apply_biotoy_profile(*model, context);
    AMREX_ALWAYS_ASSERT(state.min(Temp_comp) == Real(12.0));
    AMREX_ALWAYS_ASSERT(state.min(Salt_comp) == Real(30.0));
    for (int dye = Tracer_comp; dye < base; ++dye) {
        AMREX_ALWAYS_ASSERT(state.min(dye) == Real(-7.0));
        AMREX_ALWAYS_ASSERT(state.max(dye) == Real(-7.0));
    }
    // The final full-layout fields verify that optional profiles do not stop
    // at the old nitrogen-only inventory or overwrite the dye components.
    AMREX_ALWAYS_ASSERT(state.min(base + 8) == Real(0.5));
    AMREX_ALWAYS_ASSERT(state.min(base + 11) == Real(2100.0));
    AMREX_ALWAYS_ASSERT(state.min(base + 15) == Real(1.0));
}

// With no model, or with BioIter = 0, advance returns before requiring any
// forcing and leaves the state untouched.
void check_noops ()
{
    const Box domain(IntVect(0), IntVect(3));
    BoxArray boxes(domain);
    DistributionMapping mapping(boxes);
    SolverChoice solver;
    const int base = Tracer_comp + 1;
    MultiFab state(boxes, mapping, base + 18, 0);
    state.setVal(Real(-7.0));
    MultiFab old_state(boxes, mapping, state.nComp(), 0);
    old_state.setVal(Real(5.0));
    // Forcing pointers are left null by omission.
    const BiologyAdvanceContext advance{.old_state = old_state,
                                        .new_state = state,
                                        .solver_choice = solver,
                                        .component_start = base,
                                        .top_layer = 3,
                                        .dt = Real(10.0),
                                        .model_time = 0.0};
    for (const auto* name : {"none", "OFF"}) {
        auto none = make_biology(name);
        none->configure("unused");
        AMREX_ALWAYS_ASSERT(none->name() == "none");
        AMREX_ALWAYS_ASSERT(none->declaration().tracer_names.empty());
        AMREX_ALWAYS_ASSERT(none->declaration().plot_aliases.empty());
        none->advance(advance);
    }
    for (const auto* name : {"fennel", "ecb"}) {
        const std::string prefix = std::string("disabled_") + name;
        ParmParse pp(prefix + "." + name);
        pp.add("BioIter", 0);
        auto model = make_biology(name);
        model->configure(prefix);
        model->advance(advance);
    }
    for (int comp = 0; comp < state.nComp(); ++comp) {
        AMREX_ALWAYS_ASSERT(state.min(comp) == Real(-7.0));
        AMREX_ALWAYS_ASSERT(state.max(comp) == Real(-7.0));
    }
}

}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        check_fennel_layouts();
        check_ecb_layout();
        check_defaults();
        check_constant_ic();
        check_fennel_profile();
        check_noops();
    }
    amrex::Finalize();
}
