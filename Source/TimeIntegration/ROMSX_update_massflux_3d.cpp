#include <ROMSX.H>
#include <Utils.H>

using namespace amrex;

void
ROMSX::update_massflux_3d (const Box& phi_bx, const int ioff, const int joff,
                     Array4<Real> phi,
                     Array4<Real> Hphi,
                     Array4<Real> Hz, Array4<Real> om_v_or_on_u,
                     Array4<Real> Dphi_avg1,
                     Array4<Real> Dphi_avg2, Array4<Real> DC,
                     Array4<Real> FC, Array4<Real> CF, const int nnew)
{
    const int Mn = Geom(0).Domain().size()[0];
    const int Mm = Geom(0).Domain().size()[1];
    auto phi_bxD=phi_bx;
    phi_bxD.makeSlab(2,0);
    //Copied depth of water column calculation from DepthStretchTransform
    amrex::ParallelFor(Box(DC),
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                DC(i,j,k)=0.0;
                CF(i,j,k)=0.0;
            });
    amrex::ParallelFor(Box(FC),
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                FC(i,j,k)=0.0;
            });
    amrex::LoopOnCpu(phi_bx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                DC(i,j,k)=0.5*om_v_or_on_u(i,j,0)*(Hz(i,j,k)+Hz(i-ioff,j-joff,k));
                DC(i,j,-1)=DC(i,j,-1)+DC(i,j,k);
                CF(i,j,-1)=CF(i,j,-1)+DC(i,j,k)*phi(i,j,k,nnew);
            });

    amrex::ParallelFor(phi_bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        Real cff1=DC(i,j,-1);
        if(k==0) {
        DC(i,j,-1)=1.0/DC(i,j,-1);
        CF(i,j,-1)=DC(i,j,-1)*(CF(i,j,-1)-Dphi_avg1(i,j,0));
        }
        //Vertical mean correction on boundary points
#if 0
        if((i<0)||(j<0)||(i>=Mn+1)||(j>=Mm+1))
            phi(i,j,k) -= CF(i,j,-1);
#endif
        Hphi(i,j,k) = 0.5 * (Hphi(i,j,k)+phi(i,j,k,nnew)*DC(i,j,k));
        FC(i,j,0) = FC(i,j,0)+Hphi(i,j,k); //recursive
    });

    amrex::ParallelFor(phi_bxD,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        FC(i,j,0) = DC(i,j,-1)*(FC(i,j,0)-Dphi_avg2(i,j,0)); //recursive
    });
    amrex::ParallelFor(phi_bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        Hphi(i,j,k) = Hphi(i,j,k)-DC(i,j,k)*FC(i,j,0);
    });

}
