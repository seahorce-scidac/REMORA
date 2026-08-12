#include <REMORA.H>

using namespace amrex;

/**
 * @param[in] time   Time at which to get updated river variables
 */
void
REMORA::update_rivers (
#ifdef REMORA_USE_NETCDF
        Real time)
#else
        Real /*time*/)
#endif
{
#ifdef REMORA_USE_NETCDF
    river_source_transport->update_interpolated_to_time(time);
    river_source_transportbar->update_interpolated_to_time(time);
    for (int i_comp=0; i_comp < ncons; i_comp++) {
        if (solverChoice.do_rivers_cons[i_comp]) {
            river_source_cons[i_comp]->update_interpolated_to_time(time);
        }
    }
#endif
}
