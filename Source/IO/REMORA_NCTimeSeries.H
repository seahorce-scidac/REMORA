#ifndef _REMORA_NCTIMESEIRES_H_
#define _REMORA_NCTIMESEIRES_H_

#ifdef REMORA_USE_NETCDF

#include <string>
#include <memory>

#include <AMReX_AmrCore.H>

/** \brief A class to hold and interpolate time series data read from a NetCDF file
 *
 * The class only ever holds the two time steps necessary to interpoplate to the current
 * simulation time.
 */
class NCTimeSeries
{
    public:
        /** \brief Constructor */
        NCTimeSeries (const amrex::Vector<std::string>& a_file_names,
                      const std::string a_field_name,
                      const std::string a_time_name,
                      const amrex::Box& a_domain,
                      amrex::MultiFab* a_mf_var,
                      bool a_is2d, bool a_save_interpolated);

        /** \brief Read in time array from file and allocate data arrays */
        void Initialize ();

        /** \brief Calculate interpolated values at time and fill data for level lev */
        void update_interpolated_to_time (amrex::Real time, int lev,
                                          amrex::MultiFab* mf_lev,
                                          const amrex::Vector<amrex::Geometry>& geom,
                                          const amrex::Vector<amrex::IntVect>& ref_ratio);

        /** \brief Access interpolated data saved for a specific level */
        const amrex::MultiFab* get_interpolated_mf (int lev) const;
    private:
        /** \brief Read in data from file at time index itime and fill into mf */
        void read_in_at_time (amrex::MultiFab* mf, int itime);

        /** \brief Build cumulative refinement ratio from level 0 to lev */
        amrex::IntVect cumulative_ref_ratio (int lev,
                                             const amrex::Vector<amrex::IntVect>& ref_ratio) const;

        /// File names to read from
        amrex::Vector<std::string> file_names;
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

        /// Multifab to store data at time_before
        amrex::MultiFab* mf_before;
        /// Multifab to store data at time_after
        amrex::MultiFab* mf_after;
        /// Multifab storing temporally interpolated data on level 0
        amrex::MultiFab* mf_interp_lev0;
        /// Interpolated data on each requested AMR level if save_interpolated=true
        amrex::Vector<std::unique_ptr<amrex::MultiFab>> mf_interpolated_lev;
        /** Pointer to REMORA data that corresponds to the variable being interpolated.
         * Filled by update_interpolated_to_time if save_interpolated==false. Otherwise,
         * just used for getting box array, nodality, and distribution mapping */
        amrex::MultiFab* mf_var;
};

#endif // REMORA_USE_NETCDF
#endif //_REMORA_NCTIMESEIRES_H_
