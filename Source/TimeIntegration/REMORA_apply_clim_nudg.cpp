#include <REMORA.H>

using namespace amrex;

/**
 * @param[in   ] bx             box to apply climatology on
 * @param[in   ] ioff           offset in x-direction
 * @param[in   ] joff           offset in y-direction
 * @param[inout] var            variable to update
 * @param[in   ] var_old        variable to compare against for nudging
 * @param[in   ] var_clim       climatology value to nudge towards
 * @param[in   ] clim_coeff     nudging time scale (1/s)
 * @param[in   ] Hz             vertical cell height
 * @param[in   ] pm             1/dx
 * @param[in   ] pn             1/dy
 * @param[in   ] dt_lev         time step
 */
void
REMORA::apply_clim_nudg (const Box& bx,
                         int ioff, int joff,
                         const Array4<Real      >& var,
                         const Array4<Real const>& var_old,
                         const Array4<Real const>& var_clim,
                         const Array4<Real const>& clim_coeff,
                         const Array4<Real const>& Hz,
                         const Array4<Real const>& pm,
                         const Array4<Real const>& pn,
                         const Real dt_lev)
{
    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        Real cff = 0.5_rt * (clim_coeff(i-ioff,j-joff,k) + clim_coeff(i,j,k));
        if (ioff==1 || joff==1) {
            Real om = 2.0_rt / (pm(i-ioff,j-joff,0)+pm(i,j,0));
            Real on = 2.0_rt / (pn(i-ioff,j-joff,0)+pn(i,j,0));
            cff *= 0.5_rt * (Hz(i-ioff,j-joff,k) + Hz(i,j,k)) * om * on;
        } else {
            cff *= dt_lev;
        }
        var(i,j,k) += cff * (var_clim(i,j,k) - var_old(i,j,k));
    });
}
