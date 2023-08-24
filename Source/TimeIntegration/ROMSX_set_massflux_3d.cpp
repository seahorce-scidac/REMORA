#include <ROMSX.H>
#include <Utils.H>

using namespace amrex;

void
ROMSX::set_massflux_3d (const Box& phi_bx, const int ioff, const int joff,
                        Array4<Real> phi,
                        Array4<Real> Hphi,
                        Array4<Real> Hz, Array4<Real> om_v_or_on_u,
                        const int nrhs)
{
    //
    //-----------------------------------------------------------------------
    //  Compute horizontal mass fluxes, Hz*u/n and Hz*v/m.
    //-----------------------------------------------------------------------
    //
    amrex::ParallelFor(phi_bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        if (Box(Hz).contains(i-ioff,j-joff,k))
            {
                Hphi(i,j,k)=0.5*(Hz(i,j,k)+Hz(i-ioff,j-joff,k))*phi(i,j,k,nrhs)* om_v_or_on_u(i,j,0);
            } else {
                Hphi(i,j,k)=(Hz(i,j,k))*phi(i,j,k,nrhs)* om_v_or_on_u(i,j,0);
            }
    });
}
