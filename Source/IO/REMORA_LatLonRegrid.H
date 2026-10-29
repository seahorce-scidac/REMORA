#ifndef REMORA_LATLONREGRID_H_
#define REMORA_LATLONREGRID_H_

#ifdef REMORA_USE_NETCDF

#include <string>
#include <vector>

#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>

/** \brief How a NetCDF forcing field on a lat-lon grid is interpolated onto the model grid
 *
 * Set from remora.frc_regrid, remora.frc_lon_varname, remora.frc_lat_varname and
 * remora.frc_fill_max_rings.
 */
struct LatLonRegridOptions
{
    /** Interpolate a field whose shape is not the model grid's (remora.frc_regrid = auto);
     *  when false (off) such a field is an error */
    bool enabled = true;
    /** Name of the longitude variable; empty to take it from the field's coordinates attribute */
    std::string lon_name;
    /** Name of the latitude variable; empty to take it from the field's coordinates attribute */
    std::string lat_name;
    /** How many source cells out to look for a valid value when every corner of a target
     *  point's cell is missing data */
    int fill_max_rings = 3;
};

/**
 * \brief Bilinear interpolation of a 2D field from a rectilinear lat-lon grid onto
 *        REMORA's rho points, the way ROMS's regrid does for forcing.
 *
 * The source grid is given by 1D longitude and latitude axes, each strictly monotone in
 * either direction; the longitude may be in [-180, 180) or [0, 360), may cross the branch
 * cut, and may cover the globe, in which case it is treated as periodic. The target points
 * are the model's lon_rho and lat_rho on level 0 over the domain and its one-cell ring, the
 * cells a file on the model grid fills.
 *
 * The interpolation weights are computed once, in define(), and stored per target cell in
 * the layout the coupling remap uses: up to four weights, and the (i, j, k) source index of
 * each. apply() reads one record of the field, sets aside source values equal to the
 * variable's _FillValue or missing_value (or not finite), renormalizes the remaining weights,
 * and gives a target point whose four corners are all missing the nearest valid source value
 * within LatLonRegridOptions::fill_max_rings cells.
 */
class LatLonRegrid
{
public:
    /**
     * \brief Find the source grid of a variable and build the weights. Collective.
     *
     * @param[in] file_name   file holding the variable
     * @param[in] field_name  the variable, dimensioned (time, lat, lon)
     * @param[in] field_shape its shape, from the IO rank
     * @param[in] target_lon  rho-point longitude on level 0, degrees east
     * @param[in] target_lat  rho-point latitude on level 0, degrees north
     * @param[in] like        the forcing MultiFab the field fills, for the layout
     * @param[in] domain      level 0 domain
     * @param[in] options     remora.frc_* settings
     */
    void define (const std::string& file_name, const std::string& field_name,
                 const std::vector<long long>& field_shape,
                 const amrex::MultiFab& target_lon, const amrex::MultiFab& target_lat,
                 const amrex::MultiFab& like, const amrex::Box& domain,
                 const LatLonRegridOptions& options);

    /** \brief Whether another file's variable lies on the same source grid. Collective. */
    bool same_source_grid (const std::string& file_name, const std::string& field_name,
                           const std::vector<long long>& field_shape) const;

    /** \brief Read record itime of the variable and interpolate it into dst, over the domain
     *         and its one-cell ring. Collective. */
    void apply (const std::string& file_name, const std::string& field_name, int itime,
                amrex::MultiFab& dst) const;

    /** \brief One line describing the source grid, for the log */
    std::string summary () const;

    /** \brief The cells filled: the valid box of mfi grown by the ring, inside the ring */
    amrex::Box target_box (const amrex::MFIter& mfi) const;

    // Kernels. Public because CUDA does not allow a device lambda in a private member.
    void build_weights (const amrex::MultiFab& lon, const amrex::MultiFab& lat);
    void interpolate (const amrex::FArrayBox& src, amrex::MultiFab& dst) const;

    /** \brief Number of stencil points per target cell */
    static constexpr int stencil_size = 4;

private:
    /** Axes as read and the source shape, for comparing files */
    std::vector<double> lon_file, lat_file;
    int nlon = 0, nlat = 0;

    /** Axes put in increasing order and, for longitude, made continuous */
    std::vector<double> lon_inc, lat_inc;
    bool flip_lon = false, flip_lat = false;
    bool periodic_lon = false;

    /** Missing-data markers */
    bool has_fill = false, has_missing = false;
    double fill_value = 0.0, missing_value = 0.0;

    std::string lon_name, lat_name;
    int fill_max_rings = 3;
    amrex::Box domain;

    /** Weights and the (i, j, k) source index of each, per target cell */
    amrex::MultiFab weights;
    amrex::iMultiFab indices;
    /** Fractional source index (i, j) of each target point, for the nearest-valid search */
    amrex::MultiFab position;
};

#endif // REMORA_USE_NETCDF
#endif // REMORA_LATLONREGRID_H_
