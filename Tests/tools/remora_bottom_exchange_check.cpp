/** Check flat-column exchange budgets from REMORA checkpoint output.
 * The test decks use a uniform horizontal grid and linear vertical stretching.
 * No REMORA object or solver initialization is used by this reader.
 */
#include <AMReX.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Reduce.H>
#include <AMReX_VisMF.H>
#ifdef REMORA_USE_NETCDF
#include <pnetcdf.h>
#endif

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>

using namespace amrex;

namespace {

/** Check an independent expected value, including rounding in Real state.
 * @param[in   ] name quantity and units being checked
 * @param[in   ] actual measured value
 * @param[in   ] expected closed-form value
 * @param[in   ] scale absolute scale for a value near zero
 */
void
require_close (const std::string& name, double actual, double expected, double scale = 1)
{
    const double tolerance = std::max(2.e-10, 16.0 * std::numeric_limits<Real>::epsilon());
    Print() << name << ": " << actual << " expected " << expected << '\n';
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(std::isfinite(actual) &&
        std::abs(actual - expected) <= tolerance * std::max(scale, std::abs(expected)), name);
}

/** Checkpoint fields needed for the flat-column budget assertion. */
struct State {
    double time = 0;
    double dt = 0;
    Vector<std::string> names;
    Vector<std::unique_ptr<MultiFab>> cell, depth, surface, u, w;

    /** Resolve a state component using the checkpoint registry.
     * @param[in   ] name registered component name
     */
    int component_index (const std::string& name) const
    {
        const auto entry = std::find(names.begin(), names.end(), name);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(entry != names.end(), "Missing checkpoint field: " + name);
        return static_cast<int>(entry - names.begin());
    }
};

/** Read a checkpoint field onto the cell field's distribution.
 * @param[in   ] path serialized MultiFab path
 * @param[in   ] cells cell-centered layout and distribution
 * @param[in   ] slab whether the field has one cell in z
 * @param[in   ] index staggered index type
 */
std::unique_ptr<MultiFab>
read_field (const std::string& path, const MultiFab& cells, bool slab, const IntVect& index)
{
    MultiFab input;
    VisMF::Read(input, path);
    BoxList boxes = cells.boxArray().boxList();
    if (slab) {
        for (auto& box : boxes) { box.setRange(2, 0); }
    }
    BoxArray layout(std::move(boxes));
    layout.convert(index);
    auto result = std::make_unique<MultiFab>(layout, cells.DistributionMap(), input.nComp(), 0);
    result->ParallelCopy(input);
    return result;
}

/** Read time, component names and physical fields from a checkpoint.
 * @param[in   ] path checkpoint directory
 */
State
read_state (const std::string& path)
{
    State state;
    std::ifstream header(path + "/Header");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(header.good(), "Missing checkpoint: " + path);
    std::string line;
    std::getline(header, line);
    int finest = 0;
    int count = 0;
    header >> finest;
    for (int n = 0; n < 8; ++n) { header >> count; }
    std::getline(header, line);
    std::getline(header, line); // Step counters.
    std::getline(header, line);
    std::istringstream(line) >> state.dt;
    std::getline(header, line); // Absolute times; elapsed time follows the layouts.
    while (std::getline(header, line)) {
        std::istringstream fields(line);
        std::string key;
        fields >> key;
        if (key == "elapsed_time") { fields >> state.time; }
        if (key == "cons_names") {
            fields >> count;
            state.names.resize(count);
            for (auto& name : state.names) { fields >> name; }
        }
    }
    for (int lev = 0; lev <= finest; ++lev) {
        const std::string prefix = path + "/Level_" + std::to_string(lev) + "/";
        auto cells = std::make_unique<MultiFab>();
        VisMF::Read(*cells, prefix + "Cell");
        state.depth.push_back(read_field(prefix + "h", *cells, true, IntVect(0)));
        state.surface.push_back(read_field(prefix + "Zt_avg1", *cells, true, IntVect(0)));
        state.u.push_back(read_field(prefix + "XFace", *cells, false, IntVect(1,0,0)));
        state.w.push_back(read_field(prefix + "ZFace", *cells, false, IntVect(0,0,1)));
        state.cell.push_back(std::move(cells));
    }
    return state;
}

/** Integrate volume, temperature, salt, dye and u momentum without covered coarse cells.
 * @param[in   ] state checkpoint data
 * @param[in   ] lo model lower corner (m)
 * @param[in   ] hi model upper corner (m)
 * @param[in   ] cells level-0 cell counts
 * @param[in   ] ratio horizontal refinement ratio
 */
std::array<double,5>
budget (const State& state, const Vector<Real>& lo, const Vector<Real>& hi,
        const Vector<int>& cells, const IntVect& ratio)
{
    std::array<double,5> total{};
    const int temp = state.component_index("temp");
    const int salt = state.component_index("salt");
    const int dye = state.component_index("tracer");
    const int nz = cells[2];
    double area = double(hi[0] - lo[0]) * double(hi[1] - lo[1]) / (cells[0] * cells[1]);
    for (int lev = 0; lev < state.cell.size(); ++lev) {
        const auto& field = *state.cell[lev];
        iMultiFab mask;
        if (lev + 1 < state.cell.size()) {
            mask = makeFineMask(field.boxArray(), field.DistributionMap(),
                                state.cell[lev+1]->boxArray(), ratio, 1, 0);
        }
        ReduceOps<ReduceOpSum,ReduceOpSum,ReduceOpSum,ReduceOpSum,ReduceOpSum> op;
        ReduceData<double,double,double,double,double> data(op);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(field, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const auto c = field.const_array(mfi);
            const auto h = state.depth[lev]->const_array(mfi);
            const auto z = state.surface[lev]->const_array(mfi);
            const auto u = state.u[lev]->const_array(mfi);
            const auto fine = mask.isDefined() ? mask.const_array(mfi) : Array4<int const>();
            op.eval(mfi.tilebox(), data, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                -> GpuTuple<double,double,double,double,double>
            {
                const double volume = fine && !fine(i,j,k) ? 0
                    : area * (double(h(i,j,0)) + double(z(i,j,0))) / nz;
                return {volume, volume * double(c(i,j,k,temp)), volume * double(c(i,j,k,salt)),
                        volume * double(c(i,j,k,dye)),
                        volume * 0.5 * (double(u(i,j,k)) + double(u(i+1,j,k)))};
            });
        }
        const auto sum = data.value(op);
        total[0] += get<0>(sum);
        total[1] += get<1>(sum);
        total[2] += get<2>(sum);
        total[3] += get<3>(sum);
        total[4] += get<4>(sum);
        area /= ratio[0] * ratio[1];
    }
    ParallelDescriptor::ReduceRealSum(total.data(), total.size());
    return total;
}

#ifdef REMORA_USE_NETCDF
/** Write analytic timed forcing as test input for the ordinary model executable.
 * @param[in   ] pp model inputs
 */
void
write_forcing (ParmParse& pp)
{
    Vector<std::string> files;
    Vector<int> cells;
    pp.getarr("bottom_exchange_file", files);
    pp.getarr("n_cell", cells);
    const char* names[] = {"bottom_q", "bottom_temp", "bottom_salt", "bottom_tracer"};
    auto check = [] (int error) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(error == NC_NOERR, ncmpi_strerror(error));
    };
    Real ramp = 0;
    ParmParse("remora.prob").query("bottom_ramp_time", ramp);
    for (int file = 0; file < files.size(); ++file) {
        const int records = files.size() == 1 ? 2 : 1;
        int nc, eta, xi, time;
        std::array<int,4> vars{};
        check(ncmpi_create(ParallelDescriptor::Communicator(), files[file].c_str(),
                           NC_CLOBBER, MPI_INFO_NULL, &nc));
        check(ncmpi_def_dim(nc, "eta_rho", cells[1]+2, &eta));
        check(ncmpi_def_dim(nc, "xi_rho", cells[0]+2, &xi));
        check(ncmpi_def_dim(nc, "bottom_time", records, &time));
        int clock;
        check(ncmpi_def_var(nc, "bottom_time", NC_DOUBLE, 1, &time, &clock));
        const std::string units = "days since 2000-01-01 00:00:00";
        check(ncmpi_put_att_text(nc, clock, "units", units.size(), units.c_str()));
        const int dimensions[] = {time, eta, xi};
        for (int n = 0; n < 4; ++n) {
            check(ncmpi_def_var(nc, names[n], NC_DOUBLE, 3, dimensions, &vars[n]));
        }
        check(ncmpi_enddef(nc));
        check(ncmpi_begin_indep_data(nc));
        if (ParallelDescriptor::IOProcessor()) {
            const double seconds[] = {double(file) * 16, 16};
            double times[] = {seconds[0]/86400, seconds[1]/86400};
            check(ncmpi_put_var_double(nc, clock, times));
            const int size = (cells[0]+2) * (cells[1]+2);
            for (int n = 0; n < 4; ++n) {
                Real value = 0;
                Real ring_scale = Real(1);
                ParmParse("check").query("rho_ring_scale", ring_scale);
                pp.get(names[n], value);
                Vector<double> values(size * records, double(value));
                for (int record = 0; record < records; ++record) {
                    const double factor = n == 0 && ramp > 0 ? 1 + seconds[record]/ramp : 1;
                    for (int i = 0; i < size; ++i) {
                        const int row = i / (cells[0]+2);
                        const double ring = n == 0 && (row == 0 || row == cells[1]+1)
                                          ? double(ring_scale) : 1;
                        values[record*size+i] *= factor * ring;
                    }
                }
                check(ncmpi_put_var_double(nc, vars[n], values.data()));
            }
        }
        check(ncmpi_end_indep_data(nc));
        check(ncmpi_close(nc));
    }
}
#endif

} // namespace

/** Check a normal model run, or prepare its analytic NetCDF input.
 * @param[in   ] argc argument count
 * @param[in   ] argv input file and parameter overrides
 */
int
main (int argc, char** argv)
{
    amrex::Initialize(argc, argv);
    {
        ParmParse pp("remora"), check("check"), prob("remora.prob");
        bool make_forcing = false;
        check.query("write_forcing", make_forcing);
        if (make_forcing) {
#ifdef REMORA_USE_NETCDF
            write_forcing(pp);
#else
            amrex::Abort("Forcing input creation requires PNetCDF");
#endif
        } else {
            const auto initial = read_state("chk00000");
            std::string last;
            for (const auto& item : std::filesystem::directory_iterator(".")) {
                const auto name = item.path().filename().string();
                if (item.is_directory() && name.rfind("chk", 0) == 0 &&
                    name.size() == 8 && name > last) { last = name; }
            }
            const auto final = read_state(last);
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(final.time > initial.time, "No evolved checkpoint");
            Vector<Real> lo, hi;
            Vector<int> cells, refinements{2,2,1};
            pp.getarr("prob_lo", lo);
            pp.getarr("prob_hi", hi);
            pp.getarr("n_cell", cells);
            ParmParse("amr").queryarr("ref_ratio_vect", refinements);
            const IntVect ratio(refinements[0], refinements[1], refinements[2]);
            const auto before = budget(initial, lo, hi, cells, ratio);
            const auto after = budget(final, lo, hi, cells, ratio);
            std::string type;
            pp.get("bottom_exchange_type", type);
            Real q = 0, salt = 0, temp = 0, dye = 0;
            pp.get("bottom_q", q);
            pp.get("bottom_temp", temp);
            pp.get("bottom_salt", salt);
            pp.get("bottom_tracer", dye);
            Real ramp = 0, q_amplitude = 0, c_amplitude = 0, reverse_time = 0;
            prob.query("bottom_ramp_time", ramp);
            prob.query("bottom_q_amplitude", q_amplitude);
            prob.query("bottom_tracer_amplitude", c_amplitude);
            ParmParse("coupler").query("reverse_time", reverse_time);
            const int dye_comp = initial.component_index("tracer");
            const Real initial_max = initial.cell[0]->max(dye_comp);
            Print() << "Initial passive tracer maximum: " << initial_max << '\n';
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(initial_max > Real(0.1), "Degenerate initial dye");
            const double duration = final.time - initial.time;
            const double integrated = duration + (ramp > 0
                ? (final.time*final.time - initial.time*initial.time)/(2*double(ramp)) : 0);
            const double area = double(hi[0]-lo[0]) * double(hi[1]-lo[1]);
            const double added = type == "none" ? 0 : reverse_time > 0
                ? double(q)*area*(2*double(reverse_time)-duration) : double(q)*area*integrated;
            require_close("volume (m^3)", after[0], before[0] + added);
            const double values[] = {temp, salt, dye};
            for (int n = 1; n <= 3; ++n) {
                const double covariance = n == 3 ? 1 + double(q_amplitude)*c_amplitude/4 : 1;
                const double expected = reverse_time > 0 ? before[n]
                    : q >= 0 ? before[n] + added*values[n-1]*covariance
                    : before[n] * after[0]/before[0];
                require_close(initial.names[initial.component_index(n == 1 ? "temp" :
                    n == 2 ? "salt" : "tracer")] + " inventory", after[n], expected);
            }
            require_close("u momentum (m^4/s)", after[4],
                          q >= 0 ? before[4] : before[4]*after[0]/before[0],
                          before[0]*std::abs(double(q)));
            if (q_amplitude == 0) {
                const double expected = type == "none" ? 0 : reverse_time > 0 ? -q
                    : double(q)*(ramp > 0 ? 1+(final.time-final.dt/2)/ramp : 1);
                for (const auto& w : final.w) {
                    require_close("vertical velocity min (m/s)", w->min(0), expected, 0.001);
                    require_close("vertical velocity max (m/s)", w->max(0), expected, 0.001);
                }
            }
            if (q_amplitude == 0 && type != "coupler") {
                const double expected = type == "none" ? 0 : double(q);
                for (const auto& w : initial.w) {
                    require_close("initial vertical velocity min (m/s)", w->min(0), expected, 0.001);
                    require_close("initial vertical velocity max (m/s)", w->max(0), expected, 0.001);
                }
            }
            std::string mixing;
            pp.query("vertical_mixing_type", mixing);
            if (mixing == "GLS") {
                MultiFab tke;
                VisMF::Read(tke, last + "/Level_0/tke");
                AMREX_ALWAYS_ASSERT(tke.min(0) > Real(0) && !tke.contains_nan());
                const auto& c = *final.cell[0];
                const int salt_comp = final.component_index("salt");
                const Real bottom = c.min(salt_comp);
                const Real top = c.max(salt_comp);
                Print() << "Salinity range (psu): " << bottom << " " << top << '\n';
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(bottom < top-Real(1.e-3),
                    "Freshwater did not destabilize the bottom layer");
            }
            Print() << "Bottom exchange budgets passed\n";
        }
    }
    amrex::Finalize();
}
