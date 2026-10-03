#ifndef _REMORA_NCTIMESEIRESRIVER_H_
#define _REMORA_NCTIMESEIRESRIVER_H_

#ifdef REMORA_USE_NETCDF

#include <string>

#include <AMReX_AmrCore.H>

#include <REMORA.H>

class NCTimeSeriesRiver
{
    public:
        NCTimeSeriesRiver (const amrex::Vector<std::string>& a_file_names,
                           const std::string a_field_name,
                           const std::string a_time_name,
                           const int a_nz, const int a_use_vert_integ=0,
                           const int a_is_transport=0);

        void Initialize ();

        /// Calculate interpolated values at time, in seconds on the model clock
        void update_interpolated_to_time (double time);

        /// Container for interpolated data; Only used if save_interpolated == true
        amrex::FArrayBox* fab_interp;

        /// Whether the field is specified in the z-dimension
        int has_z;
        /// Whether to use the depth-integrated value read from file; Used
        /// when !has_z
        int use_vert_integ;
        /// Whether this field is the river transport, as opposed to a tracer
        /// concentration. Only the transport is distributed in the vertical by
        /// river_Vshape (ROMS Qsrc = Qbar*Qshape)
        int is_transport;

        /// FABs to pointers of river data
        amrex::FArrayBox* fab_before;
        amrex::FArrayBox* fab_after;
        /// Vshape data if needed
        amrex::FArrayBox* fab_vshape = nullptr;

        /**
         * Read in data from file at time index itime and fill into FAB
         */
        void read_in_at_time (amrex::FArrayBox* vec, int itime);
    private:

        /// File names to read from
        amrex::Vector<std::string> file_names;
        /// Field name in netcdf file
        std::string field_name;
        /// Field name for time series in netcdf file
        std::string time_name;

        /// Number of rivers
        int nriv;
        /// Number of vertical points
        int nz;
        /// How many vertical cells there are in the data fabs
        int nzbox;

        /// Time points in netcdf file, in seconds on the model clock
        amrex::Vector<double> river_times;
        /// File index to access a particular time
        amrex::Vector<int> file_for_time;
        /// Offset to access a particular time within its file
        amrex::Vector<int> file_itime_offset;
        /// Time index immediately before the last time interpolated to
        int i_time_before;
        /// Time in ocean_time immediately before the last time interpolated to
        double time_before;
        /// Time in ocean_time immediately after the last time interpolated to
        double time_after;

};

#endif // REMORA_USE_NETCDF
#endif //_REMORA_NCTIMESEIRESRIVER_H_
