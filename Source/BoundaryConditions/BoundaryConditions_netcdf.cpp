#include "ROMSX.H"

using namespace amrex;

#ifdef ROMSX_USE_NETCDF
/*
 * Impose boundary conditions using data read in from netcdf boundary files
 *
 * @param[out] mfs  Vector of MultiFabs to be filled
 * @param[in] time  time at which the data should be filled
 */

void
ROMSX::fill_from_bdyfiles (MultiFab& mf_to_fill, const Real time, const int bdy_var_type)
{
    int lev = 0;

    // amrex::Print() << "TIME  " << time << std::endl;

    // Time interpolation
    Real dT = bdy_time_interval;
    // amrex::Print() << "DT    " << dT << std::endl;
    // amrex::Print() << "START " << start_bdy_time << std::endl;

    Real time_since_start = time - start_bdy_time;
    int n_time = static_cast<int>( time_since_start /  dT);

    amrex::Real alpha = (time_since_start - n_time * dT) / dT;
    AMREX_ALWAYS_ASSERT( alpha >= 0. && alpha <= 1.0);
    amrex::Real oma   = 1.0 - alpha;

    // Which variable are we filling
    int ivar = bdy_var_type;

    //
    // Note that "domain" is mapped onto the type of box the data is in
    //
    Box domain = geom[lev].Domain();
    domain.convert(mf_to_fill.boxArray().ixType());

    int ncomp;

    // If we are doing the scalars then do salt as well as temp
    if (ivar == BdyVars::t) {
        ncomp = 2;
    } else {
        ncomp = 1;
    }

    // This must be true for the logic below to work
    AMREX_ALWAYS_ASSERT(Temp_comp == 0);
    AMREX_ALWAYS_ASSERT(Salt_comp == 1);

    for (int icomp = 0; icomp < ncomp; icomp++) // This is to do both temp and salt if doing scalars
    {
        // We have data at fixed time intervals we will call dT
        // Then to interpolate, given time, we can define n = (time/dT)
        // and alpha = (time - n*dT) / dT, then we define the data at time
        // as  alpha * (data at time n+1) + (1 - alpha) * (data at time n)
        const auto& bdatxlo_n   = bdy_data_xlo[n_time  ][ivar+icomp].const_array();
        const auto& bdatxlo_np1 = bdy_data_xlo[n_time+1][ivar+icomp].const_array();
        const auto& bdatxhi_n   = bdy_data_xhi[n_time  ][ivar+icomp].const_array();
        const auto& bdatxhi_np1 = bdy_data_xhi[n_time+1][ivar+icomp].const_array();
        const auto& bdatylo_n   = bdy_data_ylo[n_time  ][ivar+icomp].const_array();
        const auto& bdatylo_np1 = bdy_data_ylo[n_time+1][ivar+icomp].const_array();
        const auto& bdatyhi_n   = bdy_data_yhi[n_time  ][ivar+icomp].const_array();
        const auto& bdatyhi_np1 = bdy_data_yhi[n_time+1][ivar+icomp].const_array();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        // Currently no tiling in order to get the logic right
        for (MFIter mfi(mf_to_fill,false); mfi.isValid(); ++mfi)
        {
            Box mf_box(mf_to_fill[mfi.index()].box());

            // We only proceed if the ghost cells of this FAB are not contained by domain
            if (!domain.contains(mf_box)) {

                // Compute intersections of the FAB to be filled and the bdry data boxes
                Box xlo = bdy_data_xlo[n_time][ivar].box() & mf_box;
                Box xhi = bdy_data_xhi[n_time][ivar].box() & mf_box;
                Box ylo = bdy_data_ylo[n_time][ivar].box() & mf_box;
                Box yhi = bdy_data_yhi[n_time][ivar].box() & mf_box;

                const Array4<Real>& dest_arr = mf_to_fill.array(mfi,icomp);

                if (!xlo.isEmpty()) {
                    ParallelFor(xlo, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        dest_arr(i,j,k,icomp) = oma   * bdatxlo_n  (i,j,k,0)
                                              + alpha * bdatxlo_np1(i,j,k,0);
                    });
                }

                if (!xhi.isEmpty()) {
                    ParallelFor(xhi, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        dest_arr(i,j,k,icomp) = oma   * bdatxhi_n  (i,j,k,0)
                                              + alpha * bdatxhi_np1(i,j,k,0);
                    });
                }

                if (!ylo.isEmpty()) {
                    ParallelFor(ylo, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        dest_arr(i,j,k,icomp) = oma   * bdatylo_n  (i,j,k,0)
                                              + alpha * bdatylo_np1(i,j,k,0);
                    });
                }

                if (!yhi.isEmpty()) {
                    ParallelFor(yhi, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        dest_arr(i,j,k,icomp) = oma    * bdatyhi_n  (i,j,k,0)
                                               + alpha * bdatyhi_np1(i,j,k,0);
                    });
                }

            } // FAB not contained in domain
        } // mfi
    } // icomp

}
#endif
