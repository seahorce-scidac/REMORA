/**
 * \file REMORA_NCFormats_test.cpp
 *
 * Reads one small file, written in every format the linked PnetCDF can create, through the
 * readers REMORA's NetCDF inputs go through: ReadNetCDFFile, the attribute and variable
 * queries, and the time-units offset.
 *
 * CDF-5 is always tested. NetCDF-4 and NetCDF-4 classic model are tested when PnetCDF has
 * its NetCDF-4 driver (PNETCDF_DRIVER_NETCDF4), and so is a compressed NetCDF-4 fixture,
 * Tests/test_files/NCFormats/nc_format_test_deflate.nc, passed as the first argument. They
 * reach netCDF-C and parallel HDF5 underneath, so this is also the check that the driver
 * was built and linked correctly.
 *
 * Run it on more than one rank. The readers open a file on every rank and then read on the
 * I/O rank alone while the others wait at the close. netCDF-C turns on HDF5's collective
 * metadata reads for a parallel open, so a read that needed metadata on one rank only would
 * hang here rather than pass.
 *
 * The expected values are the ones written, so every comparison is exact: the readers
 * convert to double, which holds every short, float and double written here.
 */

#include "REMORA_NCFile.H"

#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace {

int failures = 0;
int checks = 0;

void check (const std::string& what, bool ok)
{
    ++checks;
    if (!ok) { ++failures; std::printf("FAIL %s\n", what.c_str()); }
}

void check_near (const std::string& what, double got, double expect)
{
    ++checks;
    if (got != expect) {
        ++failures;
        std::printf("FAIL %s: got %.17g, expected %.17g\n", what.c_str(), got, expect);
    }
}

constexpr int nt = 2;
constexpr int ny = 3;
constexpr int nx = 4;

double zeta_value (int t, int j, int i) { return 0.5 * t + 0.25 * j - 0.125 * i; }
float  h_value    (int j, int i)        { return 10.0f + 2.0f * j + 0.5f * i; }
short  mask_value (int j, int i)        { return static_cast<short>((i + j) % 2); }

/**
 * Write the test file collectively. The I/O rank writes every value and the other ranks
 * take part with empty requests, which is how a collective write on a record variable has
 * to go in NetCDF-4: extending the unlimited dimension is collective in HDF5.
 */
void write_file (const std::string& fname, int cmode)
{
    auto ncf = ncutils::NCFile::create(fname, cmode, amrex::ParallelDescriptor::Communicator());
    ncf.def_dim("ocean_time", NC_UNLIMITED);
    ncf.def_dim("eta_rho", ny);
    ncf.def_dim("xi_rho", nx);

    auto time = ncf.def_var("ocean_time", NC_DOUBLE, {"ocean_time"});
    time.put_attr("units", "days since 2001-01-01 00:00:00");
    time.put_attr("calendar", "proleptic_gregorian");
    auto zeta = ncf.def_var("zeta", NC_DOUBLE, {"ocean_time", "eta_rho", "xi_rho"});
    zeta.put_attr("units", "meter");
    auto h = ncf.def_var("h", NC_FLOAT, {"eta_rho", "xi_rho"});
    h.put_attr("long_name", "bathymetry at RHO-points");
    ncf.def_var("mask_rho", NC_SHORT, {"eta_rho", "xi_rho"});
    ncf.put_attr("title", "REMORA NetCDF format test");
    ncf.exit_def_mode();

    const bool io = amrex::ParallelDescriptor::IOProcessor();
    const MPI_Offset have = io ? 1 : 0;

    std::vector<double> times = {3.0, 4.5};
    std::vector<double> zetas;
    for (int t = 0; t < nt; ++t) {
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) { zetas.push_back(zeta_value(t, j, i)); }
        }
    }
    std::vector<float> hs;
    std::vector<short> masks;
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            hs.push_back(h_value(j, i));
            masks.push_back(mask_value(j, i));
        }
    }

    time.put_all(times.data(), {0}, {have * nt});
    zeta.put_all(zetas.data(), {0, 0, 0}, {have * nt, have * ny, have * nx});
    h.put_all(hs.data(), {0, 0}, {have * ny, have * nx});
    // NCVar has no short overload; this is the only short variable REMORA writes.
    const MPI_Offset start[2] = {0, 0};
    const MPI_Offset count[2] = {have * ny, have * nx};
    const int ierr = ncmpi_put_vara_short_all(ncf.ncid, ncf.var("mask_rho").varid,
                                              start, count, masks.data());
    check(fname + ": write mask_rho", ierr == NC_NOERR);
    ncf.close();
}

void read_and_check (const std::string& fname, bool expect_netcdf4)
{
    const bool io = amrex::ParallelDescriptor::IOProcessor();

    {
        auto ncf = ncutils::NCFile::open(fname, NC_NOWRITE, amrex::ParallelDescriptor::Communicator());
        check(fname + ": is_netcdf4", ncf.is_netcdf4() == expect_netcdf4);
        ncf.close();
    }

    // Every variable, whole.
    amrex::Vector<std::string> names = {"zeta", "h", "mask_rho", "ocean_time"};
    amrex::Vector<NDArray<double>> arrays(names.size());
    ReadNetCDFFile(fname, names, arrays);
    if (io) {
        check(fname + ": zeta shape", arrays[0].get_vshape() == std::vector<MPI_Offset>{nt, ny, nx});
        check(fname + ": h shape", arrays[1].get_vshape() == std::vector<MPI_Offset>{ny, nx});
        const double* zeta = arrays[0].get_data();
        const double* h    = arrays[1].get_data();
        const double* mask = arrays[2].get_data();
        const double* time = arrays[3].get_data();
        for (int t = 0; t < nt; ++t) {
            for (int j = 0; j < ny; ++j) {
                for (int i = 0; i < nx; ++i) {
                    check_near(fname + ": zeta", zeta[(t * ny + j) * nx + i], zeta_value(t, j, i));
                }
            }
        }
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                check_near(fname + ": h (float)", h[j * nx + i], h_value(j, i));
                check_near(fname + ": mask_rho (short)", mask[j * nx + i], mask_value(j, i));
            }
        }
        check_near(fname + ": ocean_time[1]", time[1], 4.5);
    }

    // One record, as the time-series readers take it.
    amrex::Vector<std::string> one_name = {"zeta"};
    amrex::Vector<NDArray<double>> one(1);
    ReadNetCDFFile(fname, one_name, one, true, 1);
    if (io) {
        check(fname + ": one-time shape", one[0].get_vshape() == std::vector<MPI_Offset>{1, ny, nx});
        check_near(fname + ": one-time zeta", one[0].get_data()[2 * nx + 3], zeta_value(1, 2, 3));
    }

    // The queries broadcast their answers, so every rank checks them.
    check(fname + ": QueryNetCDFHasVars",
          QueryNetCDFHasVars(fname, {"zeta", "h"}) && !QueryNetCDFHasVars(fname, {"zeta", "salt"}));
    check(fname + ": QueryNetCDFVarAttrStr",
          QueryNetCDFVarAttrStr(fname, "h", "long_name") && !QueryNetCDFVarAttrStr(fname, "h", "units"));
    if (io) {
        check(fname + ": ReadNetCDFVarAttrStr",
              ReadNetCDFVarAttrStr(fname, "zeta", "units") == "meter");
    } else {
        ReadNetCDFVarAttrStr(fname, "zeta", "units");
    }

    // 2001-01-01 is 366 days after a 2000-01-01 time_ref.
    check_near(fname + ": time offset",
               ReadNetCDFTimeOffsetDays(fname, "ocean_time", 20000101.0, true, false), 366.0);
}

} // namespace

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv, false);
    {
        struct Format { const char* fname; int cmode; bool netcdf4; };
        std::vector<Format> formats = {
            {"nc_format_test_cdf5.nc", NC_CLOBBER | NC_64BIT_DATA, false},
#if PNETCDF_DRIVER_NETCDF4
            {"nc_format_test_netcdf4.nc", NC_CLOBBER | NC_NETCDF4, true},
            {"nc_format_test_netcdf4_classic.nc", NC_CLOBBER | NC_NETCDF4 | NC_CLASSIC_MODEL, true},
#endif
        };
#if !PNETCDF_DRIVER_NETCDF4
        amrex::Print() << "PnetCDF has no NetCDF-4 driver; testing CDF-5 only" << std::endl;
#endif

        for (const auto& f : formats) {
            write_file(f.fname, f.cmode);
            read_and_check(f.fname, f.netcdf4);
        }

        // PnetCDF cannot create a compressed variable, so the compressed case is a
        // fixture: write_file's CDF-5 output put through 'nccopy -k nc4 -d 4 -s'.
        // netCDF-C refuses independent access to a compressed variable, so this is
        // the file that fails if a reader reads in independent mode.
#if PNETCDF_DRIVER_NETCDF4
        if (argc > 1) {
            read_and_check(argv[1], true);
        } else {
            check("compressed NetCDF-4 fixture passed as the first argument", false);
        }
#endif

        int local_failures = failures;
        amrex::ParallelDescriptor::ReduceIntSum(local_failures);
        amrex::Print() << checks << " checks on rank 0, " << local_failures
                       << " failures over all ranks, " << formats.size() << " formats" << std::endl;
        failures = local_failures;
    }
    amrex::Finalize();
    return (failures == 0) ? 0 : 1;
}
