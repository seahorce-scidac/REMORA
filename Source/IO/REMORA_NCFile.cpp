#include <string>

#include "REMORA_NCFile.H"
#include "REMORA_NCInterface.H"
#include "REMORA_TimeUnits.H"


/**
 * @param fname Name of NetCDF file
 * @param var_name Name of variable
 * @param attr_name Name of attribute to read
 * @returns attribute value
 */
std::string ReadNetCDFVarAttrStr (const std::string& fname,
                                  const std::string& var_name,
                                  const std::string& attr_name)
{
    std::string attr_val;
    auto ncf = ncutils::NCFile::open(fname, NC_NOCLOBBER);
    ncmpi_begin_indep_data(ncf.ncid);
    if (amrex::ParallelDescriptor::IOProcessor())
    {
        if (!ncf.has_var(var_name)) {
            amrex::Print() << "Trying to read attribute " << attr_name << " from variable " << var_name << " that does not exist!" << std::endl;
        }
        if (!ncf.var(var_name).has_attr(attr_name)) {
            amrex::Print() << "Trying to read attribute " << attr_name << " that does not exist from variable " << var_name << "!" << std::endl;
        }
        attr_val = ncf.var(var_name).get_attr(attr_name);
    }
    ncf.close();
    return attr_val;
}

/**
 * @param fname Name of NetCDF file
 * @param var_name Name of variable
 * @param attr_name Name of attribute to read
 * @returns whether attribute was found
 */
bool QueryNetCDFVarAttrStr (const std::string& fname,
                            const std::string& var_name,
                            const std::string& attr_name)
{
    int has_var = 0;
    auto ncf = ncutils::NCFile::open(fname, NC_NOCLOBBER);
    ncmpi_begin_indep_data(ncf.ncid);
    if (amrex::ParallelDescriptor::IOProcessor())
    {
        if (!ncf.has_var(var_name)) {
            amrex::Print() << "Trying to check for attribute " << attr_name << " from variable " << var_name << " that does not exist!" << std::endl;
        }
        has_var = ncf.var(var_name).has_attr(attr_name) ? 1 : 0;
    }
    ncf.close();
    int ioproc = amrex::ParallelDescriptor::IOProcessorNumber();
    amrex::ParallelDescriptor::Bcast(&has_var, 1, ioproc);
    return (has_var != 0);
}

double ReadNetCDFTimeOffsetDays (const std::string& fname,
                                 const std::string& time_name,
                                 double time_ref,
                                 bool units_required,
                                 bool announce)
{
    double offset_days = 0.0;
    auto ncf = ncutils::NCFile::open(fname, NC_NOCLOBBER);
    ncmpi_begin_indep_data(ncf.ncid);
    if (amrex::ParallelDescriptor::IOProcessor())
    {
        if (!ncf.has_var(time_name)) {
            amrex::Abort("Time variable '" + time_name + "' not found in " + fname);
        }
        auto time_var = ncf.var(time_name);
        const std::string where = "time variable '" + time_name + "' in " + fname;
        const std::string ref_date = remora_ref_date_string(time_ref);

        if (!time_var.has_attr("units")) {
            if (units_required) {
                amrex::Abort("No units attribute on " + where +
                             ". Expected 'days since yyyy-mm-dd hh:mm:ss'.");
            }
            amrex::Warning("Units attribute not found on time variable " + time_name +
                           ". Assuming days on the model clock");
        } else {
            const std::string units = time_var.get_attr("units");
            REMORATimeUnits parsed;
            std::string error;
            if (!remora_parse_time_units(units, parsed, error)) {
                amrex::Abort("Could not read units '" + units + "' of " + where + ": " + error);
            }
            if (parsed.unit != "days" && parsed.unit != "day") {
                amrex::Print() << "Units of " << time_name << " given as: " << units << std::endl;
                amrex::Abort("Units must be in days.");
            }

            // A file on another calendar counts days differently, and no shift
            // would put it on this clock.
            if (time_var.has_attr("calendar")) {
                const std::string calendar = time_var.get_attr("calendar");
                if (!remora_calendar_matches_time_ref(calendar, time_ref)) {
                    amrex::Abort("The calendar of " + where + " is '" + calendar +
                                 "', but remora.time_ref selects the " +
                                 remora_ref_calendar(time_ref) + " calendar. Convert the " +
                                 "file, or choose a remora.time_ref on its calendar.");
                }
            }

            if (parsed.has_reference) {
                offset_days = remora_time_units_offset_days(parsed, time_ref);
                if (announce) {
                    amrex::Print() << "Time variable '" << time_name << "' in " << fname
                                   << " counts days since "
                                   << remora_time_units_date_string(parsed);
                    if (offset_days == 0.0) {
                        amrex::Print() << ", the remora.time_ref reference date" << std::endl;
                    } else {
                        amrex::Print() << "; shifting by " << offset_days
                                       << " days to count from the remora.time_ref reference date "
                                       << ref_date << std::endl;
                    }
                }
            } else if (announce) {
                amrex::Print() << "Time variable '" << time_name << "' in " << fname
                               << " has units '" << units << "' with no reference date;"
                               << " taking its times as days since the remora.time_ref"
                               << " reference date " << ref_date << std::endl;
            }
        }
    }
    ncf.close();

    // Every rank places the times, so every rank needs the same offset.
    const int ioproc = amrex::ParallelDescriptor::IOProcessorNumber();
    amrex::ParallelDescriptor::Bcast(&offset_days, 1, ioproc);
    return offset_days;
}

/**
 * @param fname Name of NetCDF file
 * @param var_names Names of the variables to test for
 * @returns whether every named variable is present in the file
 *
 * Use this before calling BuildFABsFromNetCDFFile() on optional variables:
 * ReadNetCDFFile() calls ncf.var(name) unconditionally, so a missing variable
 * is a hard failure rather than a recoverable one. The answer is broadcast so
 * every rank agrees before any collective read is attempted.
 */
bool QueryNetCDFHasVars (const std::string& fname,
                         const amrex::Vector<std::string>& var_names)
{
    int has_all = 1;
    auto ncf = ncutils::NCFile::open(fname, NC_NOCLOBBER);
    ncmpi_begin_indep_data(ncf.ncid);
    if (amrex::ParallelDescriptor::IOProcessor())
    {
        for (const auto& var_name : var_names) {
            if (!ncf.has_var(var_name)) { has_all = 0; break; }
        }
    }
    ncf.close();
    amrex::ParallelDescriptor::Bcast(&has_all, 1, amrex::ParallelDescriptor::IOProcessorNumber());
    return (has_all != 0);
}
