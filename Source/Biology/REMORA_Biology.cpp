#include <REMORA_Biology.H>
#include <REMORA_FennelBiology.H>
#include <AMReX.H>
#include <AMReX_Utility.H>

namespace REMORABiology {
namespace {

class NoBiology final : public Biology {
public:
    std::string name () const override { return "none"; }
    void configure (const std::string&) override {}
    const BiologyDeclaration& declaration () const noexcept override { return m_declaration; }
    void initialize_analytic (BiologyInitContext&) const override {}
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
    amrex::Abort("Unknown remora.biology_model: " + name);
    return {};
}

BiologyICType
parse_biology_ic_type (const std::string& name)
{
    std::string lower = amrex::toLower(name);
    if (lower == "follow" || lower == "follow_ic_type" || lower == "default") {
        return BiologyICType::follow_ic_type;
    } else if (lower == "analytic") {
        return BiologyICType::analytic;
    } else if (lower == "netcdf") {
        return BiologyICType::netcdf;
    }
    amrex::Abort("remora.biology_ic_type must be one of: follow, analytic, netcdf");
    return BiologyICType::follow_ic_type;
}

std::string
biology_ic_type_name (BiologyICType type)
{
    switch (type) {
    case BiologyICType::follow_ic_type:
        return "follow";
    case BiologyICType::analytic:
        return "analytic";
    case BiologyICType::netcdf:
        return "netcdf";
    }
    amrex::Abort("Invalid biology IC type");
    return "follow";
}

}
