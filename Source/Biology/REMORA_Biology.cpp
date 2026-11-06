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

void
require_parameter (bool ok, const std::string& prefix, const char* name,
                   const char* condition, Real value)
{
    if (!ok) {
        Abort(prefix + "." + name + " must be " + condition + "; got " + std::to_string(value));
    }
}

BiologyConstantIC
read_constant_biology_ic (const BiologyDeclaration& declaration,
                          const std::string& remora_prefix)
{
    ParmParse pp(remora_prefix);
    BiologyConstantIC values;
    std::string missing;
    for (const auto& name : declaration.tracer_names) {
        Real value;
        if (pp.query(("biology_initial_" + name).c_str(), value)) {
            values[name] = value;
        } else {
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
                             const BiologyConstantIC& values)
{
    for (const auto& name : declaration.tracer_names) {
        const auto value = values.find(name);
        if (value == values.end()) {
            Abort("No constant biology initial value for tracer " + name);
        }
        const int comp = declaration.component_index(name, context.component_start);
        context.state.setVal(value->second, comp, 1, context.state.nGrowVect());
    }
}

}
