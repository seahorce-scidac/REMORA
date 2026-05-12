#include "REMORA_biology.H"
#include "Biology/REMORA_bio_fennel.H"
#include "Biology/REMORA_bio_placeholder.H"
#include "REMORA_IndexDefines.H"
#include <AMReX_Print.H>
#include <algorithm>
#include <cctype>

namespace REMORA_Biology {

namespace {

enum class PluginKind {
    fennel,
    placeholder
};

struct ModelDispatchEntry {
    const char* canonical_name;
    const char* display_label;
    PluginKind kind;
    std::vector<std::string> aliases;
};

const std::vector<ModelDispatchEntry>& GetModelDispatchRegistry()
{
    static const std::vector<ModelDispatchEntry> registry = {
        {"fennel", "Fennel", PluginKind::fennel, {}},
        {"ecb", "ECB", PluginKind::placeholder, {}},
        {"ecosim", "ECOSIM", PluginKind::placeholder, {}},
        {"hypoxia_srm", "Hypoxia SRM", PluginKind::placeholder, {"hypoxia-srm"}},
        {"nemuro", "NEMURO", PluginKind::placeholder, {}},
        {"npzd_franks", "NPZD-Franks", PluginKind::placeholder, {"npzd-franks"}},
        {"npzd_powell", "NPZD-Powell", PluginKind::placeholder, {"npzd-powell"}},
        {"npzd_iron", "NPZD-Iron", PluginKind::placeholder, {"npzd-iron"}},
        {"oyster_floats", "Oyster Floats", PluginKind::placeholder, {"oyster-floats"}},
        {"red_tide", "Red Tide", PluginKind::placeholder, {"red-tide"}},
        {"zaecb", "ZAECB", PluginKind::placeholder, {}}
    };
    return registry;
}

const ModelDispatchEntry* FindModelDispatchEntry(const std::string& normalized)
{
    for (const auto& entry : GetModelDispatchRegistry()) {
        if (normalized == entry.canonical_name) {
            return &entry;
        }
        for (const auto& alias : entry.aliases) {
            if (normalized == alias) {
                return &entry;
            }
        }
    }
    return nullptr;
}

int CountAvailableTracerKeys(const TracerIndexMap& tracer_indices,
                            const std::vector<std::string>& keys,
                            int ncons)
{
    int count = 0;
    for (const auto& key : keys) {
        auto it = tracer_indices.find(key);
        if (it != tracer_indices.end() && it->second >= 0 && it->second < ncons) {
            ++count;
        }
    }
    return count;
}

int PackedDonOffsetForKey(const std::string& key)
{
    if (key == "NO3") return BiologyTracers::NO3;
    if (key == "NH4") return BiologyTracers::NH4;
    if (key == "Phyt") return BiologyTracers::Phyt;
    if (key == "Chlo") return BiologyTracers::Chlo;
    if (key == "ZoopS") return BiologyTracers::ZoopS;
    if (key == "ZoopL") return BiologyTracers::ZoopL;
    if (key == "SDet_N") return BiologyTracers::SDet_N;
    if (key == "LDet_N") return BiologyTracers::LDet_N;
    if (key == "DON") return BiologyTracers::DON;
    if (key == "Oxygen") return BiologyTracers::Oxygen;
    if (key == "TIC") return BiologyTracers::TIC;
    if (key == "TAlk") return BiologyTracers::TAlk;
    return -1;
}

} // namespace

std::string NormalizeBiologyModelName(const std::string& model_name)
{
    std::string normalized = model_name;
    std::transform(normalized.begin(), normalized.end(), normalized.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return normalized;
}

ModelTracerRequirements GetModelTracerRequirements(const std::string& model_name)
{
    ModelTracerRequirements req;
    req.normalized_model_name = NormalizeBiologyModelName(model_name);

    const std::string& normalized = req.normalized_model_name;

    if (normalized == "fennel" || normalized == "ecb" || normalized == "ecosim" ||
        normalized == "hypoxia_srm" || normalized == "hypoxia-srm" || normalized == "zaecb") {
        req.required_tracer_keys = {"NO3", "NH4", "Phyt", "Chlo", "ZoopS", "ZoopL", "SDet_N", "LDet_N", "DON"};
        req.optional_tracer_keys = {"Oxygen", "TIC", "TAlk"};
        req.uses_packed_don_layout = true;
        req.packed_min_nscalar = BiologyTracers::DON + 1;
        return req;
    }

    if (normalized == "nemuro" || normalized.rfind("nemuro", 0) == 0) {
        req.required_tracer_keys = {"NO3", "NH4", "Phyt", "ZoopS", "SDet_N"};
        return req;
    }

    if (normalized == "npzd_franks" || normalized == "npzd-franks" ||
        normalized == "npzd_powell" || normalized == "npzd-powell" ||
        normalized == "npzd_iron" || normalized == "npzd-iron" ||
        normalized.rfind("npzd_", 0) == 0 || normalized.rfind("npzd-", 0) == 0) {
        req.required_tracer_keys = {"NO3", "Phyt", "ZoopS", "SDet_N"};
        return req;
    }

    if (normalized == "oyster_floats" || normalized == "oyster-floats") {
        req.required_tracer_keys = {"NO3", "Phyt"};
        return req;
    }

    if (normalized == "red_tide" || normalized == "red-tide") {
        req.required_tracer_keys = {"NO3", "NH4", "Phyt"};
        return req;
    }

    return req;
}

TracerIndexMap BuildDefaultTracerMapForModel(const std::string& model_name, int base_comp)
{
    TracerIndexMap tracer_map;
    const auto req = GetModelTracerRequirements(model_name);
    if (!req.uses_packed_don_layout) {
        return tracer_map;
    }

    auto add_key = [&](const std::string& key) {
        const int offset = PackedDonOffsetForKey(key);
        if (offset >= 0) {
            tracer_map[key] = base_comp + offset;
        }
    };

    for (const auto& key : req.required_tracer_keys) {
        add_key(key);
    }
    for (const auto& key : req.optional_tracer_keys) {
        add_key(key);
    }
    return tracer_map;
}

bool MissingAnyRequiredTracers(const TracerIndexMap& tracer_indices,
                               int ncons,
                               const std::vector<std::string>& required_keys)
{
    for (const auto& key : required_keys) {
        auto it = tracer_indices.find(key);
        if (it == tracer_indices.end() || it->second < 0 || it->second >= ncons) {
            return true;
        }
    }
    return false;
}

void BiologyModel::AddScaffoldBudgetFields(
    const std::string& model_name,
    int ncons,
    const TracerIndexMap& tracer_indices,
    long long calls,
    const std::vector<std::string>& required_tracer_keys,
    const std::vector<std::string>& optional_tracer_keys,
    bool missing_required,
    std::map<std::string, amrex::Real>& budgets) const
{
    const int required_total = static_cast<int>(required_tracer_keys.size());
    const int required_present = CountAvailableTracerKeys(tracer_indices, required_tracer_keys, ncons);
    const int required_missing = required_total - required_present;
    const int optional_total = static_cast<int>(optional_tracer_keys.size());
    const int optional_present = CountAvailableTracerKeys(tracer_indices, optional_tracer_keys, ncons);

    budgets["biology_scaffold_calls"] = static_cast<amrex::Real>(calls);
    budgets["biology_scaffold_ncons"] = static_cast<amrex::Real>(ncons);
    budgets["biology_scaffold_tracer_map_size"] = static_cast<amrex::Real>(tracer_indices.size());
    budgets["biology_scaffold_required_total"] = static_cast<amrex::Real>(required_total);
    budgets["biology_scaffold_required_present"] = static_cast<amrex::Real>(required_present);
    budgets["biology_scaffold_required_missing"] = static_cast<amrex::Real>(required_missing);
    budgets["biology_scaffold_optional_total"] = static_cast<amrex::Real>(optional_total);
    budgets["biology_scaffold_optional_present"] = static_cast<amrex::Real>(optional_present);
    budgets["biology_scaffold_missing_required_flag"] =
        missing_required ? amrex::Real(1.0) : amrex::Real(0.0);

    amrex::ignore_unused(model_name);
}

/**
 * \brief Factory implementation: create biology model by name
 *
 * Maps model names (from BIOLOGY_MODEL parameter) to concrete plugin classes.
 * Currently supports: none, fennel
 * Future: ecb, ecosim, hypoxia_srm, nemuro, npzd variants, etc.
 */
std::unique_ptr<BiologyModel> BiologyFactory::CreateModel(
    const std::string& model_name,
    int ncons,
    const TracerIndexMap& tracer_indices
)
{
    amrex::Print() << "[REMORA] Biology tracer registry entries: " << tracer_indices.size() << "\n";

    const auto normalized = NormalizeBiologyModelName(model_name);

    if (normalized.empty() || normalized == "none" || normalized == "off" ||
        normalized == "false" || normalized == "disabled") {
        amrex::Print() << "[REMORA] Biology disabled (model_name = \"" << model_name << "\")\n";
        return nullptr;
    }

    if (const auto* entry = FindModelDispatchEntry(normalized)) {
        amrex::Print() << "[REMORA] Creating " << entry->display_label << " biology model plugin\n";
        if (entry->kind == PluginKind::fennel) {
            return std::make_unique<FennelBiology>(ncons, tracer_indices);
        }
        return std::make_unique<PlaceholderBiology>(entry->canonical_name, ncons, tracer_indices);
    }

    // Future-facing family hooks so new names get a targeted message.
    if (normalized.rfind("npzd_", 0) == 0 || normalized.rfind("npzd-", 0) == 0) {
        amrex::Print() << "[REMORA] Creating NPZD-family placeholder for: " << model_name << "\n";
        return std::make_unique<PlaceholderBiology>(model_name, ncons, tracer_indices);
    }

    if (normalized.rfind("nemuro", 0) == 0) {
        amrex::Print() << "[REMORA] Creating NEMURO-family placeholder for: " << model_name << "\n";
        return std::make_unique<PlaceholderBiology>(model_name, ncons, tracer_indices);
    }

    if (normalized.rfind("ecosim", 0) == 0 || normalized.rfind("ecb", 0) == 0) {
        amrex::Print() << "[REMORA] Creating ecosystem-family placeholder for: " << model_name << "\n";
        return std::make_unique<PlaceholderBiology>(model_name, ncons, tracer_indices);
    }

    amrex::Print() << "[REMORA] Creating generic placeholder biology model for unknown type: "
                   << model_name << "\n";
    return std::make_unique<PlaceholderBiology>(model_name, ncons, tracer_indices);
}

} // namespace REMORA_Biology
