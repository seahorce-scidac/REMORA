#include <ROMSX.H>

using namespace amrex;

void
ROMSX::update_massflux_3d (const Box& phi_bx, const int ioff, const int joff,
                           Array4<Real> phi, Array4<Real> Hphi,
                           Array4<Real const> Hz, Array4<Real> om_v_or_on_u,
                           Array4<Real const> Dphi_avg1,
                           Array4<Real const> Dphi_avg2, Array4<Real> DC,
                           Array4<Real> FC, Array4<Real> CF, const int nnew)
{
    const int Mn = Geom(0).Domain().size()[0];
    const int Mm = Geom(0).Domain().size()[1];
    auto N = Geom(0).Domain().size()[2]-1; // Number of vertical "levs" aka, NZ
    auto phi_bxD=phi_bx;
    auto phi_bx_g1z=phi_bx;
    phi_bxD.makeSlab(2,0);
    phi_bx_g1z.grow(IntVect(0,0,1));

    auto geomdata = Geom(0).data();
    bool NSPeriodic = geomdata.isPeriodic(1);
    bool EWPeriodic = geomdata.isPeriodic(0);

    //Copied depth of water column calculation from DepthStretchTransform
    //Compute thicknesses of U-boxes DC(i,j,0:N-1), total depth of the water column DC(i,j,-1), and
    // incorrect vertical mean CF(i,j,-1)

    ParallelFor(phi_bx_g1z, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            DC(i,j,k)=0.0;
            CF(i,j,k)=0.0;
        });

    ParallelFor(phi_bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            FC(i,j,k)=0.0;
        });

    Gpu::streamSynchronize();

    //This takes advantage of Hz being an extra grow cell size
    ParallelFor(phi_bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            DC(i,j,k)=0.5*om_v_or_on_u(i,j,0)*(Hz(i,j,k)+Hz(i-ioff,j-joff,k));
        });

    ParallelFor(phi_bxD, [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
            for(int k=0; k<=N; k++) {
            DC(i,j,-1)=DC(i,j,-1)+DC(i,j,k);
            CF(i,j,-1)=CF(i,j,-1)+DC(i,j,k)*phi(i,j,k,nnew);
            }
        });

    // Note this loop is in the opposite direction in k in ROMS but does not
    // appear to affect results

    ParallelFor(phi_bxD, [=] AMREX_GPU_DEVICE (int i, int j, int )
    {
        for (int k=0; k<=N; k++) {
            if(k==0) {
                DC(i,j,-1)=1.0/DC(i,j,-1);
                CF(i,j,-1)=DC(i,j,-1)*(CF(i,j,-1)-Dphi_avg1(i,j,0));
            }

            if (!(NSPeriodic&&EWPeriodic))
            {
                if((((i<0)||(i>=Mn+1))&&!EWPeriodic)||(((j<0)||(j>=Mm+1))&&!NSPeriodic)) {
                    phi(i,j,k) -= CF(i,j,-1);
                    //                Abort("Untested vertical mean");
                }
            }

            Hphi(i,j,k) = 0.5 * (Hphi(i,j,k)+phi(i,j,k,nnew)*DC(i,j,k));
            FC(i,j,0) = FC(i,j,0)+Hphi(i,j,k); //recursive
        } // k
    });

    ParallelFor(phi_bxD, [=] AMREX_GPU_DEVICE (int i, int j, int )
    {
        FC(i,j,0) = DC(i,j,-1)*(FC(i,j,0)-Dphi_avg2(i,j,0)); //recursive
    });

    ParallelFor(phi_bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        Hphi(i,j,k) = Hphi(i,j,k)-DC(i,j,k)*FC(i,j,0);
    });

}
