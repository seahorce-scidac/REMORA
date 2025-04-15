#include <REMORA.H>

using namespace amrex;

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
                         const Real dt)
{
    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        Real cff = 0.5_rt * (clim_coeff(i-ioff,j-joff,k) + clim_coeff(i,j,k));
        if (ioff==1 || joff==1) {
            Real om = 2.0_rt / (pm(i-ioff,j-joff,0)+pm(i,j,0));
            Real on = 2.0_rt / (pn(i-ioff,j-joff,0)+pn(i,j,0));
            cff *= 0.5_rt * (Hz(i-ioff,j-joff,k) + Hz(i,j,k)) * om * on;
        } else {
            cff *= dt;
        }
        var(i,j,k) += cff * (var_clim(i,j,k) - var_old(i,j,k));
    });
}
