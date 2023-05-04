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
    Box validbx(IntVect(AMREX_D_DECL(phi_bx.smallEnd(0)+2,phi_bx.smallEnd(1)+2,phi_bx.smallEnd(2))),
                       IntVect(AMREX_D_DECL(phi_bx.bigEnd(0)-2,phi_bx.bigEnd(1)-2,phi_bx.bigEnd(2))));
    //Copied depth of water column calculation from DepthStretchTransform
    //Compute thicknesses of U-boxes DC(i,j,0:N-1), total depth of the water column DC(i,j,-1), and
    // incorrect vertical mean CF(i,j,-1)
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
            });
    amrex::LoopOnCpu(validbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                DC(i,j,-1)=DC(i,j,-1)+DC(i,j,k);
                CF(i,j,-1)=CF(i,j,-1)+DC(i,j,k)*phi(i,j,k,nnew);
            });

    auto geomdata = Geom(0).data();
    bool NSPeriodic = geomdata.isPeriodic(1);
    bool EWPeriodic = geomdata.isPeriodic(0);
    // Note this loop is in the opposite direction in k in ROMS but does not
    // appear to affect results
    amrex::ParallelFor(validbx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        Real cff1=DC(i,j,-1);
        if(k==0) {
            DC(i,j,-1)=1.0/DC(i,j,-1);
            CF(i,j,-1)=DC(i,j,-1)*(CF(i,j,-1)-Dphi_avg1(i,j,0));
        }
        //Vertical mean correction on boundary points
        // Maybe wrong? Will it just get obliterated on the FillBoundary
        if(!(NSPeriodic&&EWPeriodic)) {
            if((((i<0)||(i>=Mn+1))&&!EWPeriodic)||(((j<0)||(j>=Mm+1))&&!NSPeriodic)) {
                phi(i,j,k) -= CF(i,j,-1);
                Abort("Untested vertical mean");
            }
        }

        //Compute correct mass flux, Hz*v/m
        Hphi(i,j,k) = 0.5 * (Hphi(i,j,k)+phi(i,j,k,nnew)*DC(i,j,k));
        FC(i,j,0) = FC(i,j,0)+Hphi(i,j,k); //recursive
    });

    amrex::LoopOnCpu(phi_bxD,
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
