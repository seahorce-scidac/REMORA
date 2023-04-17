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
    const int Mn = Geom(0).Domain().size()[0];
    const int Mm = Geom(0).Domain().size()[1];
    //
    //-----------------------------------------------------------------------
    //  Compute horizontal mass fluxes, Hz*u/n and Hz*v/m.
    //-----------------------------------------------------------------------
    //
    amrex::ParallelFor(phi_bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        if ((ioff!=0&&i-1>=-2)||(joff!=0&&j-1>=-2))
            {
                Hphi(i,j,k)=0.5*(Hz(i,j,k)+Hz(i-ioff,j-joff,k))*phi(i,j,k,nrhs)* om_v_or_on_u(i,j,0);
            } else {
                Hphi(i,j,k)=(Hz(i,j,k))*phi(i,j,k,nrhs)* om_v_or_on_u(i,j,0);
            }
    });
    //Print()<<FArrayBox(Hz)<<std::endl;
    //Print()<<FArrayBox(om_v_or_on_u)<<std::endl;
    //Print()<<FArrayBox(phi)<<std::endl;
    //Print()<<FArrayBox(Hphi)<<std::endl;
}
