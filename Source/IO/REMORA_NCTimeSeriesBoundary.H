#ifndef _REMORA_NCTIMESERIESBOUNDARY_H_
#define _REMORA_NCTIMESERIESBOUNDARY_H_

#ifdef REMORA_USE_NETCDF

#include <string>

#include <AMReX_AmrCore.H>

#include <REMORA.H>

/** \brief A class to hold and interpolate time series data read from a NetCDF file
 *
 * The class only ever holds the two time steps necessary to interpoplate to the current
 * simulation time.
 */
class NCTimeSeriesBoundary
{
    public:
        /** \brief Constructor */
        NCTimeSeriesBoundary (int a_lev,
                              const amrex::Vector<amrex::Geometry> a_geom,
                              const amrex::Vector<std::string>& a_file_name,
                              const std::string a_field_name,
                              const std::string a_time_name,
                              const amrex::IntVect a_index_type,
                              const amrex::GpuArray<bool, AMREX_SPACEDIM*2>* a_var_need_data,
                              bool a_is2d, int rx, int ry);

        /** \brief Read in time array from file and allocate data arrays */
        void Initialize ();

        /** \brief Calculate interpolated values at time, reading in data as necessary */
        void update_interpolated_to_time (amrex::Real time);

        /// FArrayBox to store data at inteprolated time at xlo boundary */
        amrex::FArrayBox xlo_dat_interp;

        /// FArrayBox to store data at inteprolated time at ylo boundary */
        amrex::FArrayBox ylo_dat_interp;

        /// FArrayBox to store data at inteprolated time at xhi boundary */
        amrex::FArrayBox xhi_dat_interp;

        /// FArrayBox to store data at inteprolated time at yhi boundary */
        amrex::FArrayBox yhi_dat_interp;

    private:
        /** \brief Read in data from file at time index itime and fill into mf */
        void read_in_at_time (amrex::FArrayBox& fab_xlo,
                              amrex::FArrayBox& fab_xhi,
                              amrex::FArrayBox& fab_ylo,
                              amrex::FArrayBox& fab_yhi,
                              int itime);

        /** \brief spatially interpolate boundary fabs from level 0 data (read in) to current level */
        void interp_fab (amrex::FArrayBox& dat_crse,
                         amrex::FArrayBox& dat_fine);

        /// File name to read from
        amrex::Vector<std::string> file_names;
        /// Field name in netcdf file
        std::string field_name;
        /// Variable names that will be read from file
        amrex::Vector<std::string> nc_var_names;
        /// Field name for time series in netcdf file
        std::string time_name;
        /// Domain
        amrex::Box domain;
        /// Index type for field to fill
        amrex::IntVect index_type;

        /// Level at which we are holding the boundary data
        int m_lev;

        /// Geometry at all levels
        amrex::Vector<amrex::Geometry> m_geom;

        /// Refinement ratios relative to level 0
        int m_rx; int m_ry;

        /// Whether the field we're reading in is 2d
        bool is2d;

        /// Time points in netcdf file
        amrex::Vector<amrex::Real> bry_times;
        /// File index to access a particular time
        amrex::Vector<int> file_for_time;
        /// Offset to access a particular time within its file
        amrex::Vector<int> file_itime_offset;
        /// Time index immediately before the last time interpolated to
        int i_time_before;
        /// Time in ocean_times immediately before the last time interpolated to
        amrex::Real time_before;
        /// Time in ocean_times immediately after the last time interpolated to
        amrex::Real time_after;

        /// FArrayBox to store data at time_before at xlo boundary
        amrex::FArrayBox xlo_dat_before;
        /// FArrayBox to store data at time_after at xlo boundary
        amrex::FArrayBox xlo_dat_after;

        /// FArrayBox to store data at time_before at ylo boundary
        amrex::FArrayBox ylo_dat_before;
        /// FArrayBox to store data at time_after at ylo boundary
        amrex::FArrayBox ylo_dat_after;

        /// FArrayBox to store data at time_before at xhi boundary
        amrex::FArrayBox xhi_dat_before;
        /// FArrayBox to store data at time_after at xhi boundary
        amrex::FArrayBox xhi_dat_after;

        /// FArrayBox to store data at time_before at yhi boundary
        amrex::FArrayBox yhi_dat_before;
        /// FArrayBox to store data at time_after at yhi boundary
        amrex::FArrayBox yhi_dat_after;

        /// Array over boundaries indicating whether we need physical data for this variable
        amrex::GpuArray<bool, AMREX_SPACEDIM*2> var_need_data;
};

#endif // REMORA_USE_NETCDF
#endif //_REMORA_NCTIMESERIESBOUNDARY_H_
