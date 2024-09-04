#include "REMORA.H"

using namespace amrex;

#ifdef REMORA_USE_NETCDF
/*
 * Impose boundary conditions using data read in from netcdf boundary files
 *
 * @param[out] mfs  Vector of MultiFabs to be filled
 * @param[in] time  time at which the data should be filled
 */

void
REMORA::fill_from_bdyfiles (MultiFab& mf_to_fill, const MultiFab& mf_mask, const Real time, const int bccomp, const int bdy_var_type, const int icomp_to_fill)
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
    AMREX_ALWAYS_ASSERT( alpha >= 0. && alpha <= 1.0_rt);
    amrex::Real oma   = 1.0_rt - alpha;

    // Which variable are we filling
    int ivar = bdy_var_type;

    //
    // Note that "domain" is mapped onto the type of box the data is in
    //
    Box domain = geom[lev].Domain();

    const auto& mf_index_type = mf_to_fill.boxArray().ixType();
    domain.convert(mf_index_type);

    const auto& dom_lo = amrex::lbound(domain);
    const auto& dom_hi = amrex::ubound(domain);

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

    // Make sure we can interpolate in time
    AMREX_ALWAYS_ASSERT(n_time + 1 < bdy_data_xlo.size());

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

        const bool apply_west = domain_bcs_type[bccomp+icomp].lo(0) == REMORABCType::clamped;
        const bool apply_east = domain_bcs_type[bccomp+icomp].hi(0) == REMORABCType::clamped;
        const bool apply_south = domain_bcs_type[bccomp+icomp].lo(1) == REMORABCType::clamped;
        const bool apply_north = domain_bcs_type[bccomp+icomp].hi(1) == REMORABCType::clamped;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        // Currently no tiling in order to get the logic right
        for (MFIter mfi(mf_to_fill,false); mfi.isValid(); ++mfi)
        {
            Box mf_box(mf_to_fill[mfi.index()].box());

            // Compute intersections of the FAB to be filled and the bdry data boxes
            Box xlo = bdy_data_xlo[n_time][ivar].box() & mf_box;
            Box xhi = bdy_data_xhi[n_time][ivar].box() & mf_box;
            Box ylo = bdy_data_ylo[n_time][ivar].box() & mf_box;
            Box yhi = bdy_data_yhi[n_time][ivar].box() & mf_box;

            xlo.setSmall(0,lbound(mf_box).x);
            xhi.setBig  (0,ubound(mf_box).x);
            ylo.setSmall(1,lbound(mf_box).y);
            yhi.setBig  (1,ubound(mf_box).y);

            Box xlo_ylo = xlo & ylo;
            Box xlo_yhi = xlo & yhi;
            Box xhi_ylo = xhi & ylo;
            Box xhi_yhi = xhi & yhi;

            const Array4<Real>& dest_arr = mf_to_fill.array(mfi);
            const Array4<const Real>& mask_arr = mf_mask.array(mfi);

            if (!xlo.isEmpty() && apply_west) {
                ParallelFor(xlo, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    dest_arr(i,j,k,icomp+icomp_to_fill) = (oma   * bdatxlo_n  (ubound(xlo).x,j,k,0)
                                          + alpha * bdatxlo_np1(ubound(xlo).x,j,k,0)) * mask_arr(i,j,0);
                });
            }

            if (!xhi.isEmpty() && apply_east) {
                ParallelFor(xhi, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    dest_arr(i,j,k,icomp+icomp_to_fill) = (oma   * bdatxhi_n  (lbound(xhi).x,j,k,0)
                                          + alpha * bdatxhi_np1(lbound(xhi).x,j,k,0)) * mask_arr(i,j,0);
                });
            }

            if (!ylo.isEmpty() && apply_south) {
                ParallelFor(ylo, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    dest_arr(i,j,k,icomp+icomp_to_fill) = (oma   * bdatylo_n  (i,ubound(ylo).y,k,0)
                                          + alpha * bdatylo_np1(i,ubound(ylo).y,k,0)) * mask_arr(i,j,0);
                });
            }

            if (!yhi.isEmpty() && apply_north) {
                ParallelFor(yhi, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    dest_arr(i,j,k,icomp+icomp_to_fill) = (oma    * bdatyhi_n  (i,lbound(yhi).y,k,0)
                                           + alpha * bdatyhi_np1(i,lbound(yhi).y,k,0)) * mask_arr(i,j,0);
                });
            }

            // If we've applied boundary conditions to either side, update the corner
            if (!xlo_ylo.isEmpty() && (apply_west || apply_south)) {
                ParallelFor(xlo_ylo, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    dest_arr(i,j,k,icomp+icomp_to_fill) = 0.5 * (dest_arr(i,dom_lo.y+mf_index_type[1],k,icomp+icomp_to_fill) + dest_arr(dom_lo.x+mf_index_type[0],j,k,icomp+icomp_to_fill));
                });
            }
            if (!xlo_yhi.isEmpty() && (apply_west || apply_north)) {
                ParallelFor(xlo_yhi, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    dest_arr(i,j,k,icomp+icomp_to_fill) = 0.5 * (dest_arr(i,dom_hi.y-mf_index_type[1],k,icomp+icomp_to_fill) + dest_arr(dom_lo.x+mf_index_type[0],j,k,icomp+icomp_to_fill));
                });
            }
            if (!xhi_ylo.isEmpty() && (apply_east || apply_south)) {
                ParallelFor(xhi_ylo, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    dest_arr(i,j,k,icomp+icomp_to_fill) = 0.5 * (dest_arr(i,dom_lo.y+mf_index_type[1],k,icomp+icomp_to_fill) + dest_arr(dom_hi.x-mf_index_type[0],j,k,icomp+icomp_to_fill));
                });
            }
            if (!xhi_yhi.isEmpty() && (apply_east || apply_north)) {
                ParallelFor(xhi_yhi, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    dest_arr(i,j,k,icomp+icomp_to_fill) = 0.5 * (dest_arr(i,dom_hi.y-mf_index_type[1],k,icomp+icomp_to_fill) + dest_arr(dom_hi.x-mf_index_type[0],j,k,icomp+icomp_to_fill));
                });
            }
        } // mfi
    } // icomp
}
#endif
