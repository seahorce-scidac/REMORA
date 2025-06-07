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

    //! Write out the entire double variable
    void put(const double* ptr) const;
    //! Write out the entire float variable
    void put(const float* ptr) const;
    //! Write out the entire int variable
    void put(const int* ptr) const;

    //! Write out a slice of data
    void
    put(const double* dptr,
        const std::vector<MPI_Offset>& start,
        const std::vector<MPI_Offset>& count) const;

    //! Write out a slice of data, collective
    void
    put_all(const double* dptr,
        const std::vector<MPI_Offset>& start,
        const std::vector<MPI_Offset>& count) const;
    //! Write out a slice of data with with strides (see hyperslab definition in NetCDF)
    void iput(
            const double* dptr,
            const std::vector<MPI_Offset>& start,
            const std::vector<MPI_Offset>& count,
            int * request) const ;

    //! Write out a slice of data with with strides
    void
    put(const double* dptr,
        const std::vector<MPI_Offset>& start,
        const std::vector<MPI_Offset>& count,
        const std::vector<MPI_Offset>& stride) const;
    //! Write out a slice of data with with strides , collective
     void
     put_all(const double* dptr,
         const std::vector<MPI_Offset>& start,
         const std::vector<MPI_Offset>& count,
         const std::vector<MPI_Offset>& stride) const;
    //! Write out a slice of data
    void
    put(const float* dptr,
        const std::vector<MPI_Offset>& start,
        const std::vector<MPI_Offset>& count) const;
    //! Write out a slice of data, collective
    void
    put_all(const float* dptr,
        const std::vector<MPI_Offset>& start,
        const std::vector<MPI_Offset>& count) const;
    //! Write out a slice of data with with strides (see hyperslab definition in NetCDF)
    void
    put(const float* dptr,
        const std::vector<MPI_Offset>& start,
        const std::vector<MPI_Offset>& count,
        const std::vector<MPI_Offset>& stride) const;
    //! Write out a slice of data with with strides, collective
    void
    put_all(const float* dptr,
        const std::vector<MPI_Offset>& start,
        const std::vector<MPI_Offset>& count,
        const std::vector<MPI_Offset>& stride) const;


    //! Write out a slice of data
    void put(const int* dptr, const std::vector<MPI_Offset>& start, const std::vector<MPI_Offset>& count) const;
    //! Write out a slice of data, collective
    void put_all(const int* dptr, const std::vector<MPI_Offset>& start, const std::vector<MPI_Offset>& count) const;

    //! Write out a slice of data with with strides (see hyperslab definition in NetCDF)
    void
    put(const int* dptr,
        const std::vector<MPI_Offset>& start,
        const std::vector<MPI_Offset>& count,
        const std::vector<MPI_Offset>& stride) const;
    //! Write out a slice of data with with strides, collective
    void
    put_all(const int* dptr,
        const std::vector<MPI_Offset>& start,
        const std::vector<MPI_Offset>& count,
        const std::vector<MPI_Offset>& stride) const;


    //! Write out a slice of string data
    void put(const char** dptr, const std::vector<MPI_Offset>& start, const std::vector<MPI_Offset>& count) const;
    //! Write out a slice of string data with stride
    void
    put(
    const char** dptr,
    const std::vector<MPI_Offset>& start,
    const std::vector<MPI_Offset>& count,
    const std::vector<MPI_Offset>& stride) const;

    //! Read the entire variable from file
    void get(double* ptr) const;
    //! Read the entire variable from file
    void get(float* ptr) const;
    //! Read the entire variable from file
    void get(int* ptr) const;

    //! Read a chunk of data from the file
    void
    get(double* dptr, const std::vector<MPI_Offset>& start, const std::vector<MPI_Offset>& count) const;

    //! Read a chunk of data from the file, collective
    void
    get_all(double* dptr, const std::vector<MPI_Offset>& start, const std::vector<MPI_Offset>& count) const;

    //! Read a chunk of data with strides
    void
    get(double* dptr,
        const std::vector<MPI_Offset>& start,
        const std::vector<MPI_Offset>& count,
        const std::vector<MPI_Offset>& stride) const;

    //! Read a chunk of data with strides, collective
    void
    get_all(double* dptr,
        const std::vector<MPI_Offset>& start,
        const std::vector<MPI_Offset>& count,
        const std::vector<MPI_Offset>& stride) const;

    //! Read a chunk of data from the file
    void
    get(float* dptr, const std::vector<MPI_Offset>& start, const std::vector<MPI_Offset>& count) const;
    //! Read a chunk of data from file, collective
    void get_all(float* dptr, const std::vector<MPI_Offset>& start, const std::vector<MPI_Offset>& count) const;

    //! Read a chunk of data with strides
    void
    get(float* dptr,
        const std::vector<MPI_Offset>& start,
        const std::vector<MPI_Offset>& count,
        const std::vector<MPI_Offset>& stride) const;

    //! Read a chunk of data with strides, collective
    void
    get_all(float* dptr,
        const std::vector<MPI_Offset>& start,
        const std::vector<MPI_Offset>& count,
        const std::vector<MPI_Offset>& stride) const;

    //! Read a chunk of data from the file
    void
    get(int* dptr, const std::vector<MPI_Offset>& start, const std::vector<MPI_Offset>& count) const;
    //! Read a chunk of data from the file, collective
    void
    get_all(int* dptr, const std::vector<MPI_Offset>& start, const std::vector<MPI_Offset>& count) const;

    //! Read a chunk of data with strides
    void
    get(int* dptr,
        const std::vector<MPI_Offset>& start,
        const std::vector<MPI_Offset>& count,
        const std::vector<MPI_Offset>& stride) const;
    //! Read a chunk of data with strides, collective
    void
    get_all(int* dptr,
        const std::vector<MPI_Offset>& start,
        const std::vector<MPI_Offset>& count,
        const std::vector<MPI_Offset>& stride) const;


    //! Read a chunk of data from the file
    void
    get(char* dptr, const std::vector<MPI_Offset>& start, const std::vector<MPI_Offset>& count) const;

    //! Read a chunk of data with strides
    void
    get(char* dptr,
        const std::vector<MPI_Offset>& start,
        const std::vector<MPI_Offset>& count,
        const std::vector<MPI_Offset>& stride) const;

    //! Whether a variable has an attribute with name
    bool has_attr(const std::string& name) const;
    //! Set attribute "name" to "value"
    void put_attr(const std::string& name, const std::string& value) const;
    //! Set attribute "name" to "value"
    void
    put_attr(const std::string& name, const std::vector<double>& value) const;
    //! Set attribute "name" to "value"
    void
    put_attr(const std::string& name, const std::vector<float>& value) const;
    //! Set attribute "name" to "value"
    void put_attr(const std::string& name, const std::vector<int>& value) const;

    //! Read attribute from file
    std::string get_attr(const std::string& name) const;
    //! Read attribute from file
    void get_attr(const std::string& name, std::vector<double>& value) const;
    //! Read attribute from file
    void get_attr(const std::string& name, std::vector<float>& value) const;
    //! Read attribute from file
    void get_attr(const std::string& name, std::vector<int>& value) const;
};


/** Representation of a NetCDF file
 *
 *  Provide wrappers to create and open file
 */
class NCFile
{
public:

    //! Create a file
    static NCFile create(
        const std::string& name,
        const int cmode = NC_CLOBBER | NC_MPIIO,
        MPI_Comm comm = MPI_COMM_WORLD,
        MPI_Info info = MPI_INFO_NULL);

    //! Open an existing file
    static NCFile open(
        const std::string& name,
        const int cmode = NC_NOWRITE,
        MPI_Comm comm = MPI_COMM_WORLD,
        MPI_Info info = MPI_INFO_NULL);

    ~NCFile();

    //! Close file object
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
        const std::vector<std::string>& dnames,
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

    //! Set file attribute to value
    void put_attr(const std::string& name, const std::string& value) const;
    //! Set file attribute to value
    void put_attr(const std::string& name, const std::vector<double>& value) const;
    //! Set file attribute to value
    void put_attr(const std::string& name, const std::vector<float>& value) const;
    //! Set file attribute to value
    void put_attr(const std::string& name, const std::vector<int>& value) const;

    //! Read file attribute from file
    std::string get_attr(const std::string& name) const;
    //! Read file attribute from file
    void get_attr(const std::string& name, std::vector<double>& value) const;
    //! Read file attribute from file
    void get_attr(const std::string& name, std::vector<float>& value) const;
    //! Read file attribute from file
    void get_attr(const std::string& name, std::vector<int>& value) const;

    //! wait for non-blocking calls to finish
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
