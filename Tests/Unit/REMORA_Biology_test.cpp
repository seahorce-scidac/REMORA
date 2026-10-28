#include <REMORA_Biology.H>
#include <REMORA_FennelBiology.H>
#include <REMORA_prob_common.H>

#include <AMReX.H>
#include <AMReX_ParmParse.H>

using namespace amrex;
using namespace REMORABiology;

namespace {

// An application profile must be able to use the abstract interface without
// knowing Fennel's parameter type or testing an ecosystem name.
class ProfileProbe final : public Biology {
public:
    std::string name () const override { return "profile_probe"; }
    void configure (const std::string&) override {}
    const BiologyDeclaration& declaration () const noexcept override { return m_declaration; }
    void initialize_analytic (BiologyInitContext& context) const override {
        context.state.setVal(Real(123.0), context.component_start, 1, 0);
    }
    void advance (const BiologyAdvanceContext&) const override {}
private:
    BiologyDeclaration m_declaration{{"probe"}, {}};
};

void check_layouts ()
{
    const Vector<std::string> all_names{
        "NO3", "NH4", "chlorophyll", "phytoplankton", "zooplankton",
        "LdetritusN", "SdetritusN", "RdetritusN", "PO4", "LdetritusC",
        "SdetritusC", "TIC", "alkalinity", "RdetritusC", "oxygen", "ODU"};
    for (int options = 0; options < 32; ++options) {
        FennelParameters parameters;
        parameters.river_don = (options & 1) != 0;
        parameters.po4 = (options & 2) != 0;
        parameters.carbon = (options & 4) != 0;
        parameters.oxygen = (options & 8) != 0;
        parameters.odu = (options & 16) != 0;
        const std::string prefix = "layout_" + std::to_string(options);
        ParmParse pp(prefix + ".fennel");
        pp.add("river_don", int(parameters.river_don));
        pp.add("po4", int(parameters.po4));
        pp.add("carbon", int(parameters.carbon));
        pp.add("oxygen", int(parameters.oxygen));
        pp.add("odu", int(parameters.odu));
        auto model = make_biology("FENNEL");
        model->configure(prefix);
        auto const& declaration = model->declaration();
        AMREX_ALWAYS_ASSERT(model->name() == "fennel");
        AMREX_ALWAYS_ASSERT(declaration.plot_aliases == Vector<std::string>{"fennel"});
        AMREX_ALWAYS_ASSERT(int(declaration.tracer_names.size()) == Fennel::num_tracers(parameters));

        for (int dyes : {0, 1, 3}) {
            const int base = Tracer_comp + dyes;
            const auto c = Fennel::components(parameters, base);
            const Vector<int> indices{
                c.no3, c.nh4, c.chlo, c.phyt, c.zoop, c.lden, c.sden,
                c.rden, c.po4, c.ldec, c.sdec, c.tic, c.talk, c.rdec, c.oxyg, c.odu};
            int next = base;
            for (int field = 0; field < int(indices.size()); ++field) {
                if (indices[field] >= 0) {
                    AMREX_ALWAYS_ASSERT(indices[field] == next);
                    AMREX_ALWAYS_ASSERT(declaration.tracer_names[next - base] == all_names[field]);
                    ++next;
                }
            }
            AMREX_ALWAYS_ASSERT(next - base == int(declaration.tracer_names.size()));
        }
    }
}

void check_profiles_and_noops ()
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
    BiologyInitContext context{0, geometry, solver, state, base};

    ProblemBase problem;
    ProfileProbe probe;
    problem.init_analytic_biology(probe, context);
    AMREX_ALWAYS_ASSERT(state.min(base) == Real(123.0));
    AMREX_ALWAYS_ASSERT(state.max(Tracer_comp) == Real(-7.0));

    ParmParse remora_pp("remora");
    remora_pp.add("prob_name", std::string("Upwelling"));
    ParmParse parameters("profile.fennel");
    for (const auto* option : {"carbon", "oxygen", "po4", "odu", "river_don"}) {
        parameters.add(option, 1);
    }
    auto model = make_biology("fennel");
    model->configure("profile");
    problem.init_analytic_biology(*model, context);
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

    MultiFab old_state(boxes, mapping, state.nComp(), 0);
    old_state.setVal(Real(5.0));
    const BiologyAdvanceContext advance{
        old_state, state, solver, base, 3, Real(10.0), 0.0,
        nullptr, nullptr, nullptr, nullptr, nullptr,
        nullptr, nullptr, nullptr, nullptr};
    for (const auto* name : {"none", "OFF"}) {
        auto none = make_biology(name);
        none->configure("unused");
        AMREX_ALWAYS_ASSERT(none->name() == "none");
        AMREX_ALWAYS_ASSERT(none->declaration().tracer_names.empty());
        AMREX_ALWAYS_ASSERT(none->declaration().plot_aliases.empty());
        none->initialize_analytic(context);
        none->advance(advance);
    }
    ParmParse disabled("disabled.fennel");
    disabled.add("BioIter", 0);
    model = make_biology("fennel");
    model->configure("disabled");
    // No forcing pointers are supplied: a disabled update must return before
    // requiring forcing, and retain both physical fields and optional pools.
    model->advance(advance);
    AMREX_ALWAYS_ASSERT(state.min(Temp_comp) == Real(12.0));
    AMREX_ALWAYS_ASSERT(state.min(base + 11) == Real(2100.0));
    AMREX_ALWAYS_ASSERT(state.min(base + 15) == Real(1.0));
}

}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        check_layouts();
        check_profiles_and_noops();
    }
    amrex::Finalize();
}
