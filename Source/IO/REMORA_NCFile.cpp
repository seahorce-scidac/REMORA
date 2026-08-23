#include <string>

#include "REMORA_NCFile.H"
#include "REMORA_NCInterface.H"


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
