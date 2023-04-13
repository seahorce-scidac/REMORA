#include <ROMSX.H>
#include <Utils.H>

using namespace amrex;

void
ROMSX::update_massflux_3d (const Box& phi_bx, const int ioff, const int joff,
                     Array4<Real> phi_arr,
                     Array4<Real> Hphi,
                     Array4<Real> DC_arr, const int nnew)
{

    amrex::ParallelFor(phi_bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        Hphi(i,j,k) = 0.5 * (Hphi(i,j,k)+phi_arr(i,j,k,nnew)*DC_arr(i,j,k));
    });

}
