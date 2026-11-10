/** Exercise the real NetCDF reader against independently evaluated synthetic harmonics. */
#include "REMORA_NCTimeSeriesBoundary.H"

#include <AMReX_BoxIterator.H>
#include <AMReX_ParmParse.H>

#include <cmath>
#include <limits>

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    int failures = 0;
    {
        using namespace amrex;
        std::string directory;
        ParmParse pp;
        pp.get("data_dir", directory);
        const RealBox physical({0.0, 0.0, -50.0}, {8000.0, 10000.0, 0.0});
        const int periodic[3] = {0, 0, 0};
        Vector<Geometry> geom(2);
        geom[0].define(Box(IntVect(0,0,0), IntVect(7,9,3)), &physical, 0, periodic);
        geom[1].define(Box(IntVect(0,0,0), IntVect(15,29,3)), &physical, 0, periodic);
        const GpuArray<bool, 6> all_sides = {true, true, false, true, true, false};
        const Array<int, 6> boundary_number = {0, 1, -1, 2, 3, -1};
        const Array<std::string, 3> fields = {"zeta", "ubar", "vbar"};
        const Array<IntVect, 3> index_types = {IntVect(0,0,0), IntVect(1,0,0), IntVect(0,1,0)};
        const double pi = std::acos(-1.0);
        const double tolerance = 100 * std::numeric_limits<Real>::epsilon();
        for (int lev = 0; lev < 2; ++lev) {
            const int rx = lev == 0 ? 1 : 2;
            const int ry = lev == 0 ? 1 : 3;
            for (int q = 0; q < 3; ++q) {
                for (bool enabled : {false, true}) {
                    NCTimeSeriesBoundary series(lev, geom, {directory + "/subtidal.nc"},
                        fields[q], "ocean_time", index_types[q], all_sides, true, rx, ry, 20000101.0);
                    series.Initialize();
                    GpuArray<bool, 6> tidal_sides = all_sides;
                    if (q == 1) { tidal_sides[1] = tidal_sides[4] = false; }
                    if (q == 2) { tidal_sides[0] = tidal_sides[3] = false; }
                    if (enabled) { series.InitializeTides(directory + "/tides.nc", tidal_sides); }
                    Array<FArrayBox*, 6> result = {&series.xlo_dat_interp, &series.ylo_dat_interp,
                        nullptr, &series.xhi_dat_interp, &series.yhi_dat_interp, nullptr};
                    // Repeated and decreasing requests must rebuild from the subtidal records.
                    for (double elapsed : {0.0, 60.0, 60.0, 30.0, 300.0, 0.0}) {
                        series.update_interpolated_to_time(86400.0 + elapsed);
                        Gpu::streamSynchronize();
                        for (int ori : {0, 1, 3, 4}) {
                            const int along = (ori == 0 || ori == 3) ? 1 : 0;
                            const int ratio = along == 0 ? rx : ry;
                            const auto actual = result[ori]->const_array();
                            for (BoxIterator it(result[ori]->box()); it.ok(); ++it) {
                                const IntVect point = it();
                                const int coordinate = point[along];
                                const double position = coordinate < 0 ? -1.0
                                    : static_cast<double>(coordinate) / ratio;
                                const double p = position + 1 - index_types[q][along];
                                double expected = 0.02 + 0.002*boundary_number[ori]
                                    + 0.0001*p + 0.003*elapsed/86400.0;
                                if (enabled && tidal_sides[ori]) {
                                    const int lower = static_cast<int>(std::floor(p));
                                    const double fraction = p - lower;
                                    auto harmonic = [&] (int index) {
                                        const double scale = q == 0 ? 0.1 : (q == 1 ? 0.01 : 0.02);
                                        const double phase0 = index % 2 == 0 ? 359.0 : 1.0;
                                        const double phase1 = 120.0 + 10.0*boundary_number[ori];
                                        return scale*(1+0.01*index) * (
                                            std::cos(2*pi*elapsed/240.0 - phase0*pi/180.0)
                                            + 0.3*std::cos(2*pi*elapsed/600.0 - phase1*pi/180.0));
                                    };
                                    expected += (1-fraction)*harmonic(lower) + fraction*harmonic(lower+1);
                                }
                                const double error = std::abs(actual(point) - expected);
                                if (!(error < tolerance)) {
                                    if (failures < 5) {
                                        Print() << "Tidal boundary mismatch: level " << lev << " "
                                                << fields[q] << " side " << ori << " point " << point
                                                << " elapsed " << elapsed << " error " << error << '\n';
                                    }
                                    ++failures;
                                }
                            }
                        }
                    }
                }
            }
        }
        ParallelDescriptor::ReduceIntSum(failures);
        Print() << "Tidal boundary reader: " << failures << " failures\n";
    }
    amrex::Finalize();
    return failures != 0;
}
