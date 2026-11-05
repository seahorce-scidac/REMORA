#include <REMORA_Biology.H>
#include <REMORA_FennelBiology.H>
#include <REMORA_ECBBiology.H>
#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Utility.H>

namespace REMORABiology {
namespace {

class NoBiology final : public Biology {
public:
    std::string name () const override { return "none"; }
    bool accepts_legacy_checkpoint () const override { return true; }
    void configure (const std::string&) override {}
    const BiologyDeclaration& declaration () const noexcept override { return m_declaration; }
    void advance (const BiologyAdvanceContext&) const override {}
private:
    BiologyDeclaration m_declaration;
};

}

std::unique_ptr<Biology>
make_biology (const std::string& name)
{
    const auto lower = amrex::toLower(name);
    if (lower == "none" || lower == "off") {
        return std::make_unique<NoBiology>();
    }
    if (lower == "fennel") {
        return std::make_unique<FennelBiology>();
    }
    if (lower == "ecb") {
        return std::make_unique<ECBBiology>();
    }
    amrex::Abort("Unknown remora.biology_model: " + name);
    return {};
}

BiologyICType
parse_biology_ic_type (const std::string& name)
{
    std::string lower = amrex::toLower(name);
    if (lower == "follow" || lower == "follow_ic_type" || lower == "default") {
        return BiologyICType::follow_ic_type;
    } else if (lower == "constant" || lower == "const") {
        return BiologyICType::constant;
    } else if (lower == "analytic") {
        return BiologyICType::analytic;
    } else if (lower == "netcdf") {
        return BiologyICType::netcdf;
    }
    amrex::Abort("remora.biology_ic_type must be one of: follow, constant, analytic, netcdf");
    return BiologyICType::follow_ic_type;
}

std::string
biology_ic_type_name (BiologyICType type)
{
    switch (type) {
    case BiologyICType::follow_ic_type:
        return "follow";
    case BiologyICType::constant:
        return "constant";
    case BiologyICType::analytic:
        return "analytic";
    case BiologyICType::netcdf:
        return "netcdf";
    }
    amrex::Abort("Invalid biology IC type");
    return "follow";
}

using amrex::Abort;
using amrex::ParmParse;
using amrex::Real;
using amrex::Vector;

Vector<Real>
read_constant_biology_ic (const BiologyDeclaration& declaration,
                          const std::string& remora_prefix)
{
    ParmParse pp(remora_prefix);
    Vector<Real> values(declaration.tracer_names.size());
    std::string missing;
    for (int n = 0; n < static_cast<int>(declaration.tracer_names.size()); ++n) {
        const std::string& name = declaration.tracer_names[n];
        if (!pp.query(("biology_initial_" + name).c_str(), values[n])) {
            missing += " " + remora_prefix + ".biology_initial_" + name;
        }
    }
    if (!missing.empty()) {
        Abort(remora_prefix + ".biology_ic_type = constant needs a value for every"
              " biology tracer; missing:" + missing);
    }
    return values;
}

void
initialize_constant_biology (const BiologyDeclaration& declaration,
                             BiologyInitContext& context,
                             const Vector<Real>& values)
{
    const int count = static_cast<int>(declaration.tracer_names.size());
    AMREX_ALWAYS_ASSERT(static_cast<int>(values.size()) == count);
    for (int n = 0; n < count; ++n) {
        context.state.setVal(values[n], context.component_start + n, 1,
                             context.state.nGrowVect());
    }
}

}
