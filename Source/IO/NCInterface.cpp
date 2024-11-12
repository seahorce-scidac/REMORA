#include <cstdio>

#include "NCInterface.H"
#include <AMReX.H>
#include <AMReX_Print.H>

#define abort_func amrex::Abort

namespace ncutils {

namespace {

char recname[NC_MAX_NAME + 1];

void check_ncmpi_error(int ierr) {
    if (ierr != NC_NOERR) {
        printf("\n%s\n\n", ncmpi_strerror(ierr));
        abort_func("Encountered NetCDF error; aborting");
    }
}
} // namespace

std::string NCDim::name() const {
    check_ncmpi_error(ncmpi_inq_dimname(ncid, dimid, recname));
    return std::string(recname);
}

MPI_Offset NCDim::len() const {
    MPI_Offset dlen;
    check_ncmpi_error(ncmpi_inq_dimlen(ncid, dimid, &dlen));
    return dlen;
}

std::string NCVar::name() const {
    check_ncmpi_error(ncmpi_inq_varname(ncid, varid, recname));
    return std::string(recname);
}

int NCVar::ndim() const {
    int ndims;
    check_ncmpi_error(ncmpi_inq_varndims(ncid, varid, &ndims));
    return ndims;
}

std::vector<MPI_Offset> NCVar::shape() const {
    int ndims = ndim();
    std::vector<int> dimids(ndims);
    std::vector<MPI_Offset> vshape(ndims);

    check_ncmpi_error(ncmpi_inq_vardimid(ncid, varid, dimids.data()));

    for (int i = 0; i < ndims; ++i) {
        check_ncmpi_error(ncmpi_inq_dimlen(ncid, dimids[i], &vshape[i]));
    }

    return vshape;
}

void NCVar::put(const double *ptr) const {
    check_ncmpi_error(ncmpi_put_var_double(ncid, varid, ptr));
}

void NCVar::put(const float *ptr) const {
    check_ncmpi_error(ncmpi_put_var_float(ncid, varid, ptr));
}

void NCVar::put(const int *ptr) const {
    check_ncmpi_error(ncmpi_put_var_int(ncid, varid, ptr));
}

void NCVar::put(const double *dptr, const std::vector<MPI_Offset> &start, const std::vector<MPI_Offset> &count) const {
    check_ncmpi_error(ncmpi_put_vara_double(ncid, varid, start.data(), count.data(), dptr));
}

//! Write out a slice of data, collective
void NCVar::put_all(const double *dptr, const std::vector<MPI_Offset> &start, const std::vector<MPI_Offset> &count) const {
    check_ncmpi_error(ncmpi_put_vara_double_all(ncid, varid, start.data(), count.data(), dptr));
}

//! Write out a slice of data, non-blocking
void NCVar::iput(const double *dptr, const std::vector<MPI_Offset> &start, const std::vector<MPI_Offset> &count,
        int *request) const {
    check_ncmpi_error(ncmpi_iput_vara_double(ncid, varid, start.data(), count.data(), dptr, request));
}

void NCVar::put(const double *dptr, const std::vector<MPI_Offset> &start, const std::vector<MPI_Offset> &count,
        const std::vector<MPI_Offset> &stride) const {
    check_ncmpi_error(ncmpi_put_vars_double(ncid, varid, start.data(), count.data(), stride.data(), dptr));
}

void NCVar::put_all(const double *dptr, const std::vector<MPI_Offset> &start, const std::vector<MPI_Offset> &count,
        const std::vector<MPI_Offset> &stride) const {
    check_ncmpi_error(ncmpi_put_vars_double_all(ncid, varid, start.data(), count.data(), stride.data(), dptr));
}

void NCVar::put(const float *dptr, const std::vector<MPI_Offset> &start, const std::vector<MPI_Offset> &count) const {
    check_ncmpi_error(ncmpi_put_vara_float(ncid, varid, start.data(), count.data(), dptr));
}

void NCVar::put_all(const float *dptr, const std::vector<MPI_Offset> &start, const std::vector<MPI_Offset> &count) const {
    check_ncmpi_error(ncmpi_put_vara_float_all(ncid, varid, start.data(), count.data(), dptr));
}

void NCVar::put(const float *dptr, const std::vector<MPI_Offset> &start, const std::vector<MPI_Offset> &count,
        const std::vector<MPI_Offset> &stride) const {
    check_ncmpi_error(ncmpi_put_vars_float(ncid, varid, start.data(), count.data(), stride.data(), dptr));
}

void NCVar::put_all(const float *dptr, const std::vector<MPI_Offset> &start, const std::vector<MPI_Offset> &count,
        const std::vector<MPI_Offset> &stride) const {
    check_ncmpi_error(ncmpi_put_vars_float_all(ncid, varid, start.data(), count.data(), stride.data(), dptr));
}

void NCVar::put(const int *dptr, const std::vector<MPI_Offset> &start, const std::vector<MPI_Offset> &count) const {
    check_ncmpi_error(ncmpi_put_vara_int(ncid, varid, start.data(), count.data(), dptr));
}

void NCVar::put_all(const int *dptr, const std::vector<MPI_Offset> &start, const std::vector<MPI_Offset> &count) const {
    check_ncmpi_error(ncmpi_put_vara_int_all(ncid, varid, start.data(), count.data(), dptr));
}

void NCVar::put(const int *dptr, const std::vector<MPI_Offset> &start, const std::vector<MPI_Offset> &count,
        const std::vector<MPI_Offset> &stride) const {
    check_ncmpi_error(ncmpi_put_vars_int(ncid, varid, start.data(), count.data(), stride.data(), dptr));
}

void NCVar::put_all(const int *dptr, const std::vector<MPI_Offset> &start, const std::vector<MPI_Offset> &count,
        const std::vector<MPI_Offset> &stride) const {
    check_ncmpi_error(ncmpi_put_vars_int_all(ncid, varid, start.data(), count.data(), stride.data(), dptr));
}

void NCVar::put(const char **dptr, const std::vector<MPI_Offset> &start, const std::vector<MPI_Offset> &count) const {
    check_ncmpi_error(ncmpi_put_vara_text(ncid, varid, start.data(), count.data(), *dptr));
}

void NCVar::put(const char **dptr, const std::vector<MPI_Offset> &start, const std::vector<MPI_Offset> &count,
        const std::vector<MPI_Offset> &stride) const {
    check_ncmpi_error(ncmpi_put_vars_text(ncid, varid, start.data(), count.data(), stride.data(), *dptr));
}

void NCVar::get(double *ptr) const {
    check_ncmpi_error(ncmpi_get_var_double(ncid, varid, ptr));
}

void NCVar::get(float *ptr) const {
    check_ncmpi_error(ncmpi_get_var_float(ncid, varid, ptr));
}

void NCVar::get(int *ptr) const {
    check_ncmpi_error(ncmpi_get_var_int(ncid, varid, ptr));
}

void NCVar::get(double *dptr, const std::vector<MPI_Offset> &start, const std::vector<MPI_Offset> &count) const {
    check_ncmpi_error(ncmpi_get_vara_double(ncid, varid, start.data(), count.data(), dptr));
}

void NCVar::get_all(double *dptr, const std::vector<MPI_Offset> &start, const std::vector<MPI_Offset> &count) const {
    check_ncmpi_error(ncmpi_get_vara_double_all(ncid, varid, start.data(), count.data(), dptr));
}

void NCVar::get(double *dptr, const std::vector<MPI_Offset> &start, const std::vector<MPI_Offset> &count,
        const std::vector<MPI_Offset> &stride) const {
    check_ncmpi_error(ncmpi_get_vars_double(ncid, varid, start.data(), count.data(), stride.data(), dptr));
}

void NCVar::get_all(double *dptr, const std::vector<MPI_Offset> &start, const std::vector<MPI_Offset> &count,
        const std::vector<MPI_Offset> &stride) const {
    check_ncmpi_error(ncmpi_get_vars_double_all(ncid, varid, start.data(), count.data(), stride.data(), dptr));
}

void NCVar::get(float *dptr, const std::vector<MPI_Offset> &start, const std::vector<MPI_Offset> &count) const {
    check_ncmpi_error(ncmpi_get_vara_float(ncid, varid, start.data(), count.data(), dptr));
}

void NCVar::get_all(float *dptr, const std::vector<MPI_Offset> &start, const std::vector<MPI_Offset> &count) const {
    check_ncmpi_error(ncmpi_get_vara_float_all(ncid, varid, start.data(), count.data(), dptr));
}

void NCVar::get(float *dptr, const std::vector<MPI_Offset> &start, const std::vector<MPI_Offset> &count,
        const std::vector<MPI_Offset> &stride) const {
    check_ncmpi_error(ncmpi_get_vars_float(ncid, varid, start.data(), count.data(), stride.data(), dptr));
}

void NCVar::get_all(float *dptr, const std::vector<MPI_Offset> &start, const std::vector<MPI_Offset> &count,
        const std::vector<MPI_Offset> &stride) const {
    check_ncmpi_error(ncmpi_get_vars_float_all(ncid, varid, start.data(), count.data(), stride.data(), dptr));
}

void NCVar::get(int *dptr, const std::vector<MPI_Offset> &start, const std::vector<MPI_Offset> &count) const {
    check_ncmpi_error(ncmpi_get_vara_int(ncid, varid, start.data(), count.data(), dptr));
}

void NCVar::get_all(int *dptr, const std::vector<MPI_Offset> &start, const std::vector<MPI_Offset> &count) const {
    check_ncmpi_error(ncmpi_get_vara_int(ncid, varid, start.data(), count.data(), dptr));
}

void NCVar::get(int *dptr, const std::vector<MPI_Offset> &start, const std::vector<MPI_Offset> &count,
        const std::vector<MPI_Offset> &stride) const {
    check_ncmpi_error(ncmpi_get_vars_int(ncid, varid, start.data(), count.data(), stride.data(), dptr));
}

void NCVar::get_all(int *dptr, const std::vector<MPI_Offset> &start, const std::vector<MPI_Offset> &count,
        const std::vector<MPI_Offset> &stride) const {
    check_ncmpi_error(ncmpi_get_vars_int_all(ncid, varid, start.data(), count.data(), stride.data(), dptr));
}

void NCVar::get(char *dptr, const std::vector<MPI_Offset> &start, const std::vector<MPI_Offset> &count) const {
    check_ncmpi_error(ncmpi_get_vara_text(ncid, varid, start.data(), count.data(), dptr));
}

void NCVar::get(char *dptr, const std::vector<MPI_Offset> &start, const std::vector<MPI_Offset> &count,
        const std::vector<MPI_Offset> &stride) const {
    check_ncmpi_error(ncmpi_get_vars_text(ncid, varid, start.data(), count.data(), stride.data(), dptr));
}

bool NCVar::has_attr(const std::string &name) const {
    int ierr;
    MPI_Offset lenp;
    ierr = ncmpi_inq_att(ncid, varid, name.data(), NULL, &lenp);
    return (ierr == NC_NOERR);
}

void NCVar::put_attr(const std::string &name, const std::string &value) const {
    check_ncmpi_error(ncmpi_put_att_text(ncid, varid, name.data(), value.size(), value.data()));
}

void NCVar::put_attr(const std::string &name, const std::vector<double> &value) const {
    check_ncmpi_error(ncmpi_put_att_double(ncid, varid, name.data(), NC_DOUBLE, value.size(), value.data()));
}

void NCVar::put_attr(const std::string &name, const std::vector<float> &value) const {
    check_ncmpi_error(ncmpi_put_att_float(ncid, varid, name.data(), NC_FLOAT, value.size(), value.data()));
}

void NCVar::put_attr(const std::string &name, const std::vector<int> &value) const {
    check_ncmpi_error(ncmpi_put_att_int(ncid, varid, name.data(), NC_INT, value.size(), value.data()));
}

std::string NCVar::get_attr(const std::string &name) const {
    MPI_Offset lenp;
    std::vector<char> aval;
    check_ncmpi_error(ncmpi_inq_attlen(ncid, varid, name.data(), &lenp));
    aval.resize(lenp);
    check_ncmpi_error(ncmpi_get_att_text(ncid, varid, name.data(), aval.data()));
    return std::string { aval.begin(), aval.end() };
}

void NCVar::get_attr(const std::string &name, std::vector<double> &values) const {
    MPI_Offset lenp;
    check_ncmpi_error(ncmpi_inq_attlen(ncid, varid, name.data(), &lenp));
    values.resize(lenp);
    check_ncmpi_error(ncmpi_get_att_double(ncid, varid, name.data(), values.data()));
}

void NCVar::get_attr(const std::string &name, std::vector<float> &values) const {
    MPI_Offset lenp;
    check_ncmpi_error(ncmpi_inq_attlen(ncid, varid, name.data(), &lenp));
    values.resize(lenp);
    check_ncmpi_error(ncmpi_get_att_float(ncid, varid, name.data(), values.data()));
}

void NCVar::get_attr(const std::string &name, std::vector<int> &values) const {
    MPI_Offset lenp;
    check_ncmpi_error(ncmpi_inq_attlen(ncid, varid, name.data(), &lenp));
    values.resize(lenp);
    check_ncmpi_error(ncmpi_get_att_int(ncid, varid, name.data(), values.data()));
}

NCDim NCFile::dim(const std::string &name) const {
    int newid;
    check_ncmpi_error(ncmpi_inq_dimid(ncid, name.data(), &newid));
    return NCDim { ncid, newid };
}

NCDim NCFile::def_dim(const std::string &name, const size_t len) const {
    int newid;
    check_ncmpi_error(ncmpi_def_dim(ncid, name.data(), (MPI_Offset) len, &newid));
    return NCDim { ncid, newid };
}

NCVar NCFile::def_scalar(const std::string &name, const nc_type dtype) const {
    int newid;
    check_ncmpi_error(ncmpi_def_var(ncid, name.data(), dtype, 0, NULL, &newid));
    return NCVar { ncid, newid };
}

NCVar NCFile::def_array(const std::string &name, const nc_type dtype, const std::vector<std::string> &dnames) const {
    int newid;
    int ndims = dnames.size();
    std::vector<int> dimids(ndims);
    for (int i = 0; i < ndims; ++i) {
        dimids[i] = dim(dnames[i]).dimid;
    }

    check_ncmpi_error(ncmpi_def_var(ncid, name.data(), dtype, ndims, dimids.data(), &newid));
    return NCVar { ncid, newid };
}

NCVar NCFile::def_array_fill(const std::string &name, const nc_type dtype, const std::vector<std::string> &dnames,
        const void *fill_val) const {
    int newid;
    int ndims = dnames.size();
    std::vector<int> dimids(ndims);
    for (int i = 0; i < ndims; ++i) {
        dimids[i] = dim(dnames[i]).dimid;
    }

    check_ncmpi_error(ncmpi_def_var(ncid, name.data(), dtype, ndims, dimids.data(), &newid));
    check_ncmpi_error(ncmpi_def_var_fill(ncid, newid, NC_FILL, fill_val));
    return NCVar { ncid, newid };
}

NCVar NCFile::var(const std::string &name) const {
    int varid;
    check_ncmpi_error(ncmpi_inq_varid(ncid, name.data(), &varid));
    return NCVar { ncid, varid };
}

int NCFile::num_dimensions() const {
    int ndims;
    check_ncmpi_error(ncmpi_inq(ncid, &ndims, NULL, NULL, NULL));
    return ndims;
}

int NCFile::num_attributes() const {
    int nattrs;
    check_ncmpi_error(ncmpi_inq(ncid, NULL, NULL, &nattrs, NULL));
    return nattrs;
}

int NCFile::num_variables() const {
    int nvars;
    check_ncmpi_error(ncmpi_inq(ncid, NULL, &nvars, NULL, NULL));
    return nvars;
}

bool NCFile::has_dim(const std::string &name) const {
    int ierr = ncmpi_inq_dimid(ncid, name.data(), NULL);
    return (ierr == NC_NOERR);
}

bool NCFile::has_var(const std::string &name) const {
    int ierr = ncmpi_inq_varid(ncid, name.data(), NULL);
    return (ierr == NC_NOERR);
}

bool NCFile::has_attr(const std::string &name) const {
    int ierr;
    MPI_Offset lenp;
    ierr = ncmpi_inq_att(ncid, NC_GLOBAL, name.data(), NULL, &lenp);
    return (ierr == NC_NOERR);
}

void NCFile::put_attr(const std::string &name, const std::string &value) const {
    check_ncmpi_error(ncmpi_put_att_text(ncid, NC_GLOBAL, name.data(), value.size(), value.data()));
}

void NCFile::put_attr(const std::string &name, const std::vector<double> &value) const {
    check_ncmpi_error(ncmpi_put_att_double(ncid, NC_GLOBAL, name.data(), NC_DOUBLE, value.size(), value.data()));
}

void NCFile::put_attr(const std::string &name, const std::vector<float> &value) const {
    check_ncmpi_error(ncmpi_put_att_float(ncid, NC_GLOBAL, name.data(), NC_FLOAT, value.size(), value.data()));
}

void NCFile::put_attr(const std::string &name, const std::vector<int> &value) const {
    check_ncmpi_error(ncmpi_put_att_int(ncid, NC_GLOBAL, name.data(), NC_INT, value.size(), value.data()));
}

std::string NCFile::get_attr(const std::string &name) const {
    MPI_Offset lenp;
    std::vector<char> aval;
    check_ncmpi_error(ncmpi_inq_attlen(ncid, NC_GLOBAL, name.data(), &lenp));
    aval.resize(lenp);
    check_ncmpi_error(ncmpi_get_att_text(ncid, NC_GLOBAL, name.data(), aval.data()));
    return std::string { aval.begin(), aval.end() };
}

void NCFile::get_attr(const std::string &name, std::vector<double> &values) const {
    MPI_Offset lenp;
    check_ncmpi_error(ncmpi_inq_attlen(ncid, NC_GLOBAL, name.data(), &lenp));
    values.resize(lenp);
    check_ncmpi_error(ncmpi_get_att_double(ncid, NC_GLOBAL, name.data(), values.data()));
}

void NCFile::get_attr(const std::string &name, std::vector<float> &values) const {
    MPI_Offset lenp;
    check_ncmpi_error(ncmpi_inq_attlen(ncid, NC_GLOBAL, name.data(), &lenp));
    values.resize(lenp);
    check_ncmpi_error(ncmpi_get_att_float(ncid, NC_GLOBAL, name.data(), values.data()));
}

void NCFile::get_attr(const std::string &name, std::vector<int> &values) const {
    MPI_Offset lenp;
    check_ncmpi_error(ncmpi_inq_attlen(ncid, NC_GLOBAL, name.data(), &lenp));
    values.resize(lenp);
    check_ncmpi_error(ncmpi_get_att_int(ncid, NC_GLOBAL, name.data(), values.data()));
}

std::vector<NCDim> NCFile::all_dims() const {
    std::vector<NCDim> adims;
    int ndims = num_dimensions();
    adims.reserve(ndims);
    for (int i = 0; i < ndims; ++i) {
        adims.emplace_back(NCDim { ncid, i });
    }
    return adims;
}

std::vector<NCVar> NCFile::all_vars() const {
    std::vector<NCVar> avars;
    int nvars = num_variables();
    avars.reserve(nvars);
    for (int i = 0; i < nvars; ++i) {
        avars.emplace_back(NCVar { ncid, i });
    }
    return avars;
}

void NCFile::enter_def_mode() const {
    int ierr;
    ierr = ncmpi_redef(ncid);

    // Ignore already in define mode error
    if (ierr == NC_EINDEFINE)
        return;
    // Handle all other errors
    check_ncmpi_error(ierr);
}

void NCFile::exit_def_mode() const {
    check_ncmpi_error(ncmpi_enddef(ncid));
}

//Uncomment for parallel NetCDF
NCFile NCFile::create(const std::string &name, const int cmode, MPI_Comm comm, MPI_Info info) {
    int ncid;
    check_ncmpi_error(ncmpi_create(comm, name.data(), cmode, info, &ncid));
    return NCFile(ncid);
}

//Uncomment for parallel NetCDF
NCFile NCFile::open(const std::string &name, const int cmode, MPI_Comm comm, MPI_Info info) {
    int ncid;
    check_ncmpi_error(ncmpi_open(comm, name.data(), cmode, info, &ncid));
    return NCFile(ncid);
}

void NCFile::wait_all(int num_requests, int *requests) {
    std::vector<int> statuses(num_requests);
    ncmpi_wait_all(ncid, num_requests, requests, &statuses[0]);
}

NCFile::~NCFile() {
    if (is_open)
        check_ncmpi_error(ncmpi_close(ncid));
}

void NCFile::close() {
    is_open = false;
    check_ncmpi_error(ncmpi_close(ncid));
}
} // namespace ncutils
