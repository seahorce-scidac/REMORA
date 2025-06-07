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
        attr_val = ncf.var(var_name).get_attr(attr_name);
    }
    ncf.close();
    return attr_val;
}
