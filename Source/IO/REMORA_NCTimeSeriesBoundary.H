#ifndef _REMORA_NCTIMESERIESBOUNDARY_H_
#define _REMORA_NCTIMESERIESBOUNDARY_H_

#ifdef REMORA_USE_NETCDF

#include <array>
#include <string>

#include <AMReX_AmrCore.H>
#include <AMReX_Orientation.H>

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
                              const amrex::GpuArray<bool, AMREX_SPACEDIM*2>& a_var_need_data,
                              bool a_is2d, int rx, int ry, double a_time_ref, int a_verbose);

        /** \brief Read in time array from file and allocate data arrays */
        void Initialize ();

        /** \brief Calculate interpolated values at time, in seconds on the model clock,
         *         reading in data as necessary. A time equal to the last one leaves the
         *         interpolated data as they are. */
        void update_interpolated_to_time (double time);

        /** \brief Data interpolated to the last time given to update_interpolated_to_time, on
         *         the strip one point thick outside horizontal side ori. A side that takes no
         *         boundary data holds 0. */
        const amrex::FArrayBox& dat_interp (amrex::Orientation ori) const
        {
            return m_dat_interp[ori];
        }

        /// Refinement ratios relative to level 0
        int m_rx; int m_ry;

        /** \brief spatially interpolate the level-0 strip of side ori (read in) to this level's */
        void interp_fab (const amrex::FArrayBox& dat_crse, amrex::FArrayBox& dat_fine,
                         amrex::Orientation ori);

    private:
        /** \brief Read in data from file at time index itime and fill into mf */
        void read_in_at_time (std::array<amrex::FArrayBox, AMREX_SPACEDIM*2>& fabs, int itime);

        /** \brief The strip one point thick just outside horizontal side ori of dom, at this
         *         field's index type */
        amrex::Box strip_box (const amrex::Box& dom, amrex::Orientation ori) const;

        /** \brief The file variables read, one per side in m_sides and in that order */
        amrex::Vector<std::string> var_names () const;

        /// File name to read from
        amrex::Vector<std::string> file_names;
        /// Field name in netcdf file
        std::string field_name;
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

        /// Whether the field we're reading in is 2d
        bool is2d;

        /** remora.time_ref: the reference date the model clock counts from. The
         * file's times are placed on that clock from the reference date in their
         * units attribute. Double, like every reference date. */
        double time_ref;

        /// Time points in netcdf file, in seconds on the model clock
        amrex::Vector<double> bry_times;
        /// File index to access a particular time
        amrex::Vector<int> file_for_time;
        /// Offset to access a particular time within its file
        amrex::Vector<int> file_itime_offset;
        /// Time index immediately before the last time interpolated to
        int i_time_before;
        /// Time in ocean_times immediately before the last time interpolated to
        double time_before;
        /// Time in ocean_times immediately after the last time interpolated to
        double time_after;
        /// Whether m_dat_interp holds data interpolated to interp_time
        bool have_interp = false;
        /// Time m_dat_interp was last interpolated to
        double interp_time = 0.0;

        /// Data at time_before, at time_after, and interpolated to the current time, on each
        /// horizontal side's strip, indexed by Orientation; the z entries stay empty
        std::array<amrex::FArrayBox, AMREX_SPACEDIM*2> m_dat_before;
        std::array<amrex::FArrayBox, AMREX_SPACEDIM*2> m_dat_after;
        std::array<amrex::FArrayBox, AMREX_SPACEDIM*2> m_dat_interp;

        /// The sides that take boundary data. Every read, broadcast and interpolation loops over
        /// these and only these, and the file variables read are derived from them.
        amrex::Vector<amrex::Orientation> m_sides;

        /// REMORA's verbosity: at 1 or more, print a line for each record read
        int m_verbose = 0;

        /// Array over boundaries indicating whether we need physical data for this variable
        amrex::GpuArray<bool, AMREX_SPACEDIM*2> var_need_data;
};

#endif // REMORA_USE_NETCDF
#endif //_REMORA_NCTIMESERIESBOUNDARY_H_
