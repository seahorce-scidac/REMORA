#include <ROMSX.H>

using namespace amrex;

/**
 * set_massflux_3d
 *
 * @param[in ] u    velocity in x-direction
 * @param[out] Hu   H-weighted velocity in x-direction
 * @param[in ] on_u weighting for u
 * @param[in ] v    velocity in y-direction
 * @param[out] Hv   H-weighted velocity in y-direction
 * @param[in ] om_v weighting for v
 */

void
ROMSX::set_massflux_3d (const Array4<Real const>& u,
                        const Array4<Real      >& Hu,
                        const Array4<Real const>& on_u,
                        const Array4<Real const>& v,
                        const Array4<Real      >& Hv,
                        const Array4<Real const>& om_v,
                        const Array4<Real const>& Hz)
{
    //
    //-----------------------------------------------------------------------
    //  Compute horizontal mass fluxes, Hz*u/n and Hz*v/m.
    //-----------------------------------------------------------------------
    //
    ParallelFor(Box(Hu), [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        Hu(i,j,k)=0.5*(Hz(i,j,k)+Hz(i-1,j,k))*u(i,j,k)* on_u(i,j,0);
    });

    ParallelFor(Box(Hv), [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        Hv(i,j,k)=0.5*(Hz(i,j,k)+Hz(i,j-1,k))*v(i,j,k)* om_v(i,j,0);
    });
}
