#ifndef _REMORA_NCTIMESEIRESRIVER_H_
#define _REMORA_NCTIMESEIRESRIVER_H_

#ifdef REMORA_USE_NETCDF

#include <string>

#include <AMReX_AmrCore.H>

#include <REMORA.H>

class NCTimeSeriesRiver
{
    public:
        NCTimeSeriesRiver (const std::string a_file_name,
                           const std::string a_field_name,
                           const std::string a_time_name,
                           const int a_nz, const int a_use_vert_integ=0);

        void Initialize ();

        /// Calculate
        void update_interpolated_to_time (amrex::Real time);

        /// Container for interpolated data; Only used if save_interpolated == true
        amrex::FArrayBox* fab_interp;

        /// Whether the field is specified in the z-dimension
        int has_z;
        /// Whether to use the depth-integrated value read from file; Used
        /// when !has_z
        int use_vert_integ;

        /// FABs to pointers of river data
        amrex::FArrayBox* fab_before;
        amrex::FArrayBox* fab_after;
        /// Vshape data if needed
        amrex::FArrayBox* fab_vshape;

        /**
         * Read in data from file at time index itime and fill into FAB
         */
        void read_in_at_time (amrex::FArrayBox* vec, int itime);
    private:

        /// File name to read from
        std::string file_name;
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

        /// Time points in netcdf file
        amrex::Vector<amrex::Real> river_times;
        /// Time index immediately before the last time interpolated to
        int i_time_before;
        /// Time in ocean_time immediately before the last time interpolated to
        amrex::Real time_before;
        /// Time in ocean_time immediately after the last time interpolated to
        amrex::Real time_after;

};

#endif // REMORA_USE_NETCDF
#endif //_REMORA_NCTIMESEIRESRIVER_H_
