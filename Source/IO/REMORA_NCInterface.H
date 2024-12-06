/** \file nc_interface.H
 *
 *  Interface to NetCDF library
 *
 *  Defines convenience wrappers to interact with a NetCDF file in a more
 *  OOP-like manner.
 */

#ifndef REMORA_NC_INTERFACE_H
#define REMORA_NC_INTERFACE_H

#ifdef REMORA_USE_NETCDF
#include <string>
#include <unordered_map>
#include <vector>

#include <pnetcdf.h>

namespace ncutils {

//! Wrapper around NetCDF data types
struct NCDType
{
    static constexpr nc_type Int = NC_INT;
#ifdef AMREX_USE_FLOAT
    static constexpr nc_type Real = NC_FLOAT;
    using RType = float;
#else
    static constexpr nc_type Real = NC_DOUBLE;
    using RType = double;
#endif
};

//! Representation of NetCDF dimension
struct NCDim
{
    //! File/Group Identifier
    const int ncid;

    //! Dimension ID used with NetCDF API
    const int dimid;

    //! Name of this dimension
    std::string name() const;

    //! Length of this dimension
    MPI_Offset len() const;
};

//! Representation of a NetCDF variable
struct NCVar
{
    //! File/Group identifier
    const int ncid;

    //! Variable ID used with NetCDF API
    const int varid;

    //! Name of this variable
    std::string name() const;

    //! Number of array dimensions for this variable
    int ndim() const;

    //! Shape of the array (size in each array dimension)
    std::vector<MPI_Offset> shape() const;

    //! Write out the entire variable
    void put(const double*) const;
    void put(const float*) const;
    void put(const int*) const;

    //! Write out a slice of data
    void
    put(const double*,
        const std::vector<MPI_Offset>&,
        const std::vector<MPI_Offset>&) const;

    //! Write out a slice of data, collective
    void
    put_all(const double*,
        const std::vector<MPI_Offset>&,
        const std::vector<MPI_Offset>&) const;
    //! Write out a slice of data with with strides (see hyperslab definition in
    //! NetCDF)
    void iput(
            const double* dptr,
            const std::vector<MPI_Offset>& start,
            const std::vector<MPI_Offset>& count,
            int * request) const ;

    void
    put(const double*,
        const std::vector<MPI_Offset>&,
        const std::vector<MPI_Offset>&,
        const std::vector<MPI_Offset>&) const;
    //! Write out a slice of data with with strides , collective
     void
     put_all(const double*,
         const std::vector<MPI_Offset>&,
         const std::vector<MPI_Offset>&,
         const std::vector<MPI_Offset>&) const;
    //! Write out a slice of data
    void
    put(const float*,
        const std::vector<MPI_Offset>&,
        const std::vector<MPI_Offset>&) const;
    //! Write out a slice of data, collective
    void
    put_all(const float*,
        const std::vector<MPI_Offset>&,
        const std::vector<MPI_Offset>&) const;
    //! Write out a slice of data with with strides (see hyperslab definition in
    //! NetCDF)
    void
    put(const float*,
        const std::vector<MPI_Offset>&,
        const std::vector<MPI_Offset>&,
        const std::vector<MPI_Offset>&) const;
    //! Write out a slice of data with with strides, collective
    void
    put_all(const float*,
        const std::vector<MPI_Offset>&,
        const std::vector<MPI_Offset>&,
        const std::vector<MPI_Offset>&) const;


    void put(const int*, const std::vector<MPI_Offset>&, const std::vector<MPI_Offset>&) const;
    void put_all(const int*, const std::vector<MPI_Offset>&, const std::vector<MPI_Offset>&) const;

    void
    put(const int*,
        const std::vector<MPI_Offset>&,
        const std::vector<MPI_Offset>&,
        const std::vector<MPI_Offset>&) const;
    void
    put_all(const int*,
        const std::vector<MPI_Offset>&,
        const std::vector<MPI_Offset>&,
        const std::vector<MPI_Offset>&) const;


    void put(const char**, const std::vector<MPI_Offset>&, const std::vector<MPI_Offset>&) const;
    void
    put(
    const char** dptr,
    const std::vector<MPI_Offset>& start,
    const std::vector<MPI_Offset>& count,
    const std::vector<MPI_Offset>& stride) const;

    //! Read the entire variable from file
    void get(double*) const;
    void get(float*) const;
    void get(int*) const;

    //! Read a chunk of data from the file
    void
    get(double*, const std::vector<MPI_Offset>&, const std::vector<MPI_Offset>&) const;

    //! Read a chunk of data from the file, collective
    void
    get_all(double*, const std::vector<MPI_Offset>&, const std::vector<MPI_Offset>&) const;

    //! Read a chunk of data with strides
    void
    get(double*,
        const std::vector<MPI_Offset>&,
        const std::vector<MPI_Offset>&,
        const std::vector<MPI_Offset>&) const;

    //! Read a chunk of data with strides, collective
    void
    get_all(double*,
        const std::vector<MPI_Offset>&,
        const std::vector<MPI_Offset>&,
        const std::vector<MPI_Offset>&) const;

    //! Read a chunk of data from the file
    void
    get(float*, const std::vector<MPI_Offset>&, const std::vector<MPI_Offset>&) const;
    // collective
    void get_all(float*, const std::vector<MPI_Offset>&, const std::vector<MPI_Offset>&) const;

    //! Read a chunk of data with strides
    void
    get(float*,
        const std::vector<MPI_Offset>&,
        const std::vector<MPI_Offset>&,
        const std::vector<MPI_Offset>&) const;

    //! Read a chunk of data with strides, collective
    void
    get_all(float*,
        const std::vector<MPI_Offset>&,
        const std::vector<MPI_Offset>&,
        const std::vector<MPI_Offset>&) const;

    void
    get(int*, const std::vector<MPI_Offset>&, const std::vector<MPI_Offset>&) const;
    void
    get_all(int*, const std::vector<MPI_Offset>&, const std::vector<MPI_Offset>&) const;

    void
    get(int*,
        const std::vector<MPI_Offset>&,
        const std::vector<MPI_Offset>&,
        const std::vector<MPI_Offset>&) const;
   // collective
    void
    get_all(int*,
        const std::vector<MPI_Offset>&,
        const std::vector<MPI_Offset>&,
        const std::vector<MPI_Offset>&) const;


    //! Read a chunk of data from the file
    void
    get(char*, const std::vector<MPI_Offset>&, const std::vector<MPI_Offset>&) const;

    //! Read a chunk of data with strides
    void
    get(char*,
        const std::vector<MPI_Offset>&,
        const std::vector<MPI_Offset>&,
        const std::vector<MPI_Offset>&) const;

    bool has_attr(const std::string& name) const;
    void put_attr(const std::string& name, const std::string& value) const;
    void
    put_attr(const std::string& name, const std::vector<double>& value) const;
    void
    put_attr(const std::string& name, const std::vector<float>& value) const;
    void put_attr(const std::string& name, const std::vector<int>& value) const;

    std::string get_attr(const std::string& name) const;
    void get_attr(const std::string& name, std::vector<double>& value) const;
    void get_attr(const std::string& name, std::vector<float>& value) const;
    void get_attr(const std::string& name, std::vector<int>& value) const;
};


/** Representation of a NetCDF file
 *
 *  Provide wrappes to create and open file
 */
class NCFile
{
public:

    static NCFile create(
        const std::string& name,
        const int cmode = NC_CLOBBER | NC_MPIIO,
        MPI_Comm comm = MPI_COMM_WORLD,
        MPI_Info info = MPI_INFO_NULL);

    static NCFile open(
        const std::string& name,
        const int cmode = NC_NOWRITE,
        MPI_Comm comm = MPI_COMM_WORLD,
        MPI_Info info = MPI_INFO_NULL);

    ~NCFile();

    void close();

    //! Number of dimensions
    int num_dimensions() const;

    //! Number of variables
    int num_variables() const;

    //! Number of attributes
    int num_attributes() const;

    //! Check if a dimension exists by name
    bool has_dim(const std::string&) const;

    //! Check if a variable exists by name
    bool has_var(const std::string&) const;

    //! Check if an attribute exists
    bool has_attr(const std::string&) const;

    //! Get the dimension instance by name
    NCDim dim(const std::string&) const;

    //! Get the variable instance by name
    NCVar var(const std::string&) const;

    //! Define new dimension
    NCDim def_dim(const std::string&, const size_t len) const;

    //! Define a scalar variable, i.e., 0-dimensional array
    NCVar def_scalar(const std::string& name, const nc_type dtype) const;

    //! Define an array
    NCVar def_array(
        const std::string& name,
        const nc_type dtype,
        const std::vector<std::string>&) const;

    //! Define an array with a fill value
    NCVar def_array_fill(
        const std::string& name,
        const nc_type dtype,
        const std::vector<std::string>&,
        const void* fill_val) const;

    //! Define a variable (wrapper for def_array)
    NCVar def_var(
        const std::string& name,
        const nc_type dtype,
        const std::vector<std::string>& dnames) const
    {
        return def_array(name, dtype, dnames);
    }

    //! Define a variable (wrapper for def_array)
    NCVar def_var_fill(
        const std::string& name,
        const nc_type dtype,
        const std::vector<std::string>& dnames,
        const void* fill_val) const
    {
        return def_array_fill(name, dtype, dnames, fill_val);
    }

    void put_attr(const std::string& name, const std::string& value) const;
    void put_attr(const std::string& name, const std::vector<double>& value) const;
    void put_attr(const std::string& name, const std::vector<float>& value) const;
    void put_attr(const std::string& name, const std::vector<int>& value) const;

    std::string get_attr(const std::string& name) const;
    void get_attr(const std::string& name, std::vector<double>& value) const;
    void get_attr(const std::string& name, std::vector<float>& value) const;
    void get_attr(const std::string& name, std::vector<int>& value) const;

    // for non-blocking calls
    void wait_all( int num_requests, int * requests);

    //! Return a list of all dimensions defined in this group
    std::vector<NCDim> all_dims() const;

    //! Return a list of all variables defined in this group
    std::vector<NCVar> all_vars() const;

    //! Enter definition mode (not needed for NetCDF4 format)
    void enter_def_mode() const;

    //! Exit definition mode
    void exit_def_mode() const;

    NCFile(const int id) : ncid(id), is_open{true} {}

    const int ncid{-1};
    bool is_open{false};
};

} // namespace ncutils

#else
// why do we need this if netcdf is not defined ?

namespace ncutils {

struct NCDim {
    const int ncid { -1 };
    const int dimid { -1 };
};

struct NCVar {
    const int ncid { -1 };
    const int varid { -1 };
};

class NCFile {
public:
    const int ncid { -1 };
};

} // namespace ncutils

#endif

#endif /* NC_INTERFACE_H */
