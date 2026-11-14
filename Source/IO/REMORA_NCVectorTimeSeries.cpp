#include "REMORA_NCVectorTimeSeries.H"

#ifdef REMORA_USE_NETCDF

#include <AMReX_Print.H>

using namespace amrex;

NCVectorTimeSeries::NCVectorTimeSeries (const Vector<std::string>& file_names,
                                        const std::string& x_name, const std::string& y_name,
                                        const std::string& time_name, const Box& domain,
                                        MultiFab* x_model, MultiFab* y_model,
                                        MultiFab* x_rho, MultiFab* y_rho,
                                        double time_ref, const LatLonRegridOptions& options,
                                        const MultiFab* lon, const MultiFab* lat)
    : m_file_names(file_names), m_name{x_name, y_name}, m_time_name(time_name),
      m_domain(domain), m_model{x_model, y_model}, m_rho{x_rho, y_rho},
      m_time_ref(time_ref), m_options(&options), m_lon(lon), m_lat(lat)
{}

void
NCVectorTimeSeries::Initialize ()
{
    using Kind = NCTimeSeries::SourceGridKind;
    Kind kind[2];
    for (int dir = 0; dir < 2; ++dir) {
        kind[dir] = NCTimeSeries::source_grid(m_file_names, m_name[dir], m_domain,
                                              m_model[dir]->ixType(), true);
    }
    const bool any_latlon = m_options->enabled && (kind[0] == Kind::latlon || kind[1] == Kind::latlon);
    const bool any_model = kind[0] == Kind::model || kind[1] == Kind::model;
    // Model-grid points that are rho points already, as for the winds
    const bool model_on_rho = m_model[0]->ixType().cellCentered() &&
                              m_model[1]->ixType().cellCentered();
    const std::string pair = m_name[0] + " and " + m_name[1];

    if (any_latlon && any_model) {
        if (!model_on_rho) {
            amrex::Abort("One of " + pair + " in the forcing files is on the model grid's u or v "
                         "points and the other on a lat-lon grid; give both on one kind of grid.");
        }
        // Both read on rho points: they pair up unless one has to be rotated with the other
        if (m_options->rotate) {
            amrex::Abort("One of " + pair + " comes from a lat-lon forcing grid and the other from "
                         "the model grid, so they cannot be rotated onto the model grid together. "
                         "Give both the same source, or set remora.frc_rotate = off if the file's "
                         "components are already along the grid axes.");
        }
    }
    m_on_rho = any_latlon;
    m_earth_frame = any_latlon && m_options->rotate;

    for (int dir = 0; dir < 2; ++dir) {
        MultiFab* target = m_on_rho ? m_rho[dir] : m_model[dir];
        m_series[dir] = std::make_unique<NCTimeSeries>(m_file_names, m_name[dir], m_time_name,
                                                       m_domain, target, true, false, m_time_ref,
                                                       m_options, m_lon, m_lat);
        m_series[dir]->Initialize();
    }

    // A component read on rho points in place of u or v points must have come from a lat-lon
    // grid (or be uniform). One that matched the rho shape instead is a model-grid file with
    // the component on the wrong points.
    if (m_on_rho && !model_on_rho) {
        for (int dir = 0; dir < 2; ++dir) {
            if (!m_series[dir]->regridded() && !m_series[dir]->spatially_uniform()) {
                const Box& d = m_domain;
                amrex::Abort(pair + " in the forcing files have the shape of the model grid's "
                             "rho points. On the model grid, " + m_name[0] + " is on u points, "
                             "(eta_rho, xi_u) = (" + std::to_string(d.length(1)+2) + ", " +
                             std::to_string(d.length(0)+1) + "), and " + m_name[1] + " on v points, "
                             "(eta_v, xi_rho) = (" + std::to_string(d.length(1)+1) + ", " +
                             std::to_string(d.length(0)+2) + "); a lat-lon file needs longitude "
                             "and latitude axes.");
            }
        }
    }
}

void
NCVectorTimeSeries::update_interpolated_to_time (double time, int lev, MultiFab* x, MultiFab* y,
                                                 const Vector<Geometry>& geom,
                                                 const Vector<IntVect>& ref_ratio)
{
    m_series[0]->update_interpolated_to_time(time, lev, x, geom, ref_ratio);
    m_series[1]->update_interpolated_to_time(time, lev, y, geom, ref_ratio);
}

NCTimeSeries&
NCVectorTimeSeries::component (int dir)
{
    AMREX_ALWAYS_ASSERT(dir == 0 || dir == 1);
    return *m_series[dir];
}

#endif // REMORA_USE_NETCDF
