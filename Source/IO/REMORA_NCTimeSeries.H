#ifndef _REMORA_NCTIMESEIRES_H_
#define _REMORA_NCTIMESEIRES_H_

#ifdef REMORA_USE_NETCDF

#include <string>

#include <AMReX_AmrCore.H>

#include <REMORA.H>

/** \brief A class to hold and interpolate time series data read from a NetCDF file
 *
 * The class only ever holds the two time steps necessary to interpoplate to the current
 * simulation time.
 */
class NCTimeSeries
{
    public:
        /** \brief Constructor */
        NCTimeSeries (const std::string a_file_name,
                      const std::string a_field_name,
                      const std::string a_time_name,
                      const amrex::Box& a_domain,
                      amrex::MultiFab* a_mf_var,
                      bool a_is2d, bool a_save_interpolated);

        /** \brief Read in time array from file and allocate data arrays */
        void Initialize ();

        /** \brief Calculate interpolated values at time, reading in data as necessary */
        void update_interpolated_to_time (amrex::Real time);

        /** Container for interpolated data; Only used if save_interpolated == true */
        amrex::MultiFab* mf_interpolated;
    private:
        /** \brief Read in data from file at time index itime and fill into mf */
        void read_in_at_time (amrex::MultiFab* mf, int itime);

        /// File name to read from
        std::string file_name;
        /// Field name in netcdf file
        std::string field_name;
        /// Field name for time series in netcdf file
        std::string time_name;
        /// Domain
        amrex::Box domain;

        /// Whether the field we're reading in is 2d
        bool is2d;
        /** Whether to save interpolated results in mf_interpolated. If false,
         * data will be calculated to mf_var */
        bool save_interpolated;

        /// Time points in netcdf file
        amrex::Vector<amrex::Real> ocean_times;
        /// Time index immediately before the last time interpolated to
        int i_time_before;
        /// Time in ocean_times immediately before the last time interpolated to
        amrex::Real time_before;
        /// Time in ocean_times immediately after the last time interpolated to
        amrex::Real time_after;

        /// Multifab to store data at time_before
        amrex::MultiFab* mf_before;
        /// Multifab to store data at time_after
        amrex::MultiFab* mf_after;
        /** Pointer to REMORA data that corresponds to the variable being interpolated.
         * Filled by update_interpolated_to_time if save_interpolated==false. Otherwise,
         * just used for getting box array, nodality, and distribution mapping */
        amrex::MultiFab* mf_var;
};

#endif // REMORA_USE_NETCDF
#endif //_REMORA_NCTIMESEIRES_H_
