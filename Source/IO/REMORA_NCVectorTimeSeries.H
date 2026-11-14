#ifndef REMORA_NCVECTORTIMESERIES_H_
#define REMORA_NCVECTORTIMESERIES_H_

#ifdef REMORA_USE_NETCDF

#include <memory>
#include <string>

#include <AMReX_MultiFab.H>

#include "REMORA_LatLonRegrid.H"
#include "REMORA_NCTimeSeries.H"

/** \brief The two components of a vector read from NetCDF files, such as Uwind and Vwind or
 *         sustr and svstr
 *
 * On the model grid a file gives the components along the grid's axes, at the points of the
 * model targets: rho points for the winds, u and v points for the stress. On a lat-lon grid it
 * gives them east and north. Both are then interpolated to rho points, into the rho targets,
 * and the caller turns them onto the grid's axes there and, for the stress, averages them to
 * the faces (REMORA::rho_vector_to_grid_axes). A component given as one value per time pairs
 * with either kind of grid.
 *
 * Collective: construct, initialize and update on every rank.
 */
class NCVectorTimeSeries
{
public:
    /**
     * @param[in] file_names   forcing files
     * @param[in] x_name       name of the xi (or east) component in the files
     * @param[in] y_name       name of the eta (or north) component
     * @param[in] time_name    name of the time variable, or empty for the default
     * @param[in] domain       level-0 domain
     * @param[in] x_model      level-0 target of the xi component on the model grid
     * @param[in] y_model      level-0 target of the eta component on the model grid
     * @param[in] x_rho        level-0 rho-point target of the east component from a lat-lon grid
     * @param[in] y_rho        level-0 rho-point target of the north component
     * @param[in] time_ref     remora.time_ref
     * @param[in] options      how a lat-lon field is interpolated and rotated (remora.frc_*)
     * @param[in] lon          level-0 rho-point longitude
     * @param[in] lat          level-0 rho-point latitude
     */
    NCVectorTimeSeries (const amrex::Vector<std::string>& file_names,
                        const std::string& x_name, const std::string& y_name,
                        const std::string& time_name, const amrex::Box& domain,
                        amrex::MultiFab* x_model, amrex::MultiFab* y_model,
                        amrex::MultiFab* x_rho, amrex::MultiFab* y_rho,
                        double time_ref, const LatLonRegridOptions& options,
                        const amrex::MultiFab* lon, const amrex::MultiFab* lat);

    /** \brief Decide where each component comes from, check that the two pair up, and read
     *         the times. Aborts on a pair that cannot be put on the grid's axes. */
    void Initialize ();

    /** \brief Whether the components are read on rho points, into the rho targets: true when
     *         one comes from a lat-lon grid */
    bool on_rho () const { return m_on_rho; }

    /** \brief Whether the components are east and north, to be turned onto the grid's axes */
    bool earth_frame () const { return m_earth_frame; }

    /** \brief Interpolate both components to time for level lev, into x and y: the rho
     *         targets' level-lev arrays when on_rho, else the model targets' */
    void update_interpolated_to_time (double time, int lev,
                                      amrex::MultiFab* x, amrex::MultiFab* y,
                                      const amrex::Vector<amrex::Geometry>& geom,
                                      const amrex::Vector<amrex::IntVect>& ref_ratio);

    /** \brief One component's series, for reading it alone when the coupled driver supplies
     *         the other, which is only right when it needs no rotation */
    NCTimeSeries& component (int dir);

private:
    amrex::Vector<std::string> m_file_names;
    std::string m_name[2];
    std::string m_time_name;
    amrex::Box m_domain;
    amrex::MultiFab* m_model[2];
    amrex::MultiFab* m_rho[2];
    double m_time_ref;
    const LatLonRegridOptions* m_options;
    const amrex::MultiFab* m_lon;
    const amrex::MultiFab* m_lat;

    std::unique_ptr<NCTimeSeries> m_series[2];
    bool m_on_rho = false;
    bool m_earth_frame = false;
};

#endif // REMORA_USE_NETCDF
#endif // REMORA_NCVECTORTIMESERIES_H_
