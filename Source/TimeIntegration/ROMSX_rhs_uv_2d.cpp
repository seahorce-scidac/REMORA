#include <ROMSX.H>
#include <Utils.H>

using namespace amrex;

//
// rhs_2d
//

void
ROMSX::rhs_2d (const Box& bx,
               Array4<Real> uold  , Array4<Real> vold,
               Array4<Real> ru, Array4<Real> rv,
               /*Array4<Real> rufrc, Array4<Real> rvfrc,
                 Array4<Real> sustr, Array4<Real> svstr,*/
               Array4<Real> Huon, Array4<Real> Hvom,
               /*
               Array4<Real> on_u, Array4<Real> om_v,
               Array4<Real> om_u, Array4<Real> on_v,
               Array4<Real> W   , Array4<Real> FC,
               */
               int nrhs, int N)
{
    //copy the tilebox
    Box gbx1 = bx;
    Box gbx2 = bx;

    Box ubx = surroundingNodes(bx,0);
    Box vbx = surroundingNodes(bx,1);

    //make only gbx be grown to match multifabs
    gbx2.grow(IntVect(NGROW,NGROW,0));
    gbx1.grow(IntVect(NGROW-1,NGROW-1,0));

    //
    // Scratch space
    //
    FArrayBox fab_Huee(gbx2,1,amrex::The_Async_Arena()); //fab_Huee.setVal(0.0);
    FArrayBox fab_uee(gbx2,1,amrex::The_Async_Arena()); //fab_uee.setVal(0.0);

    FArrayBox fab_Hvee(gbx2,1,amrex::The_Async_Arena()); //fab_Hvee.setVal(0.0);
    FArrayBox fab_vee(gbx2,1,amrex::The_Async_Arena()); //fab_vee.setVal(0.0);

    FArrayBox fab_Hvxx(gbx2,1,amrex::The_Async_Arena()); //fab_Hvxx.setVal(0.0);
    FArrayBox fab_uxx(gbx2,1,amrex::The_Async_Arena()); //fab_uxx.setVal(0.0);

    FArrayBox fab_Huxx(gbx2,1,amrex::The_Async_Arena()); //fab_Huxx.setVal(0.0);
    FArrayBox fab_vxx(gbx2,1,amrex::The_Async_Arena()); //fab_vxx.setVal(0.0);

    FArrayBox fab_UFx(gbx2,1,amrex::The_Async_Arena()); //fab_UFx.setVal(0.0);
    FArrayBox fab_UFe(gbx2,1,amrex::The_Async_Arena()); //fab_UFe.setVal(0.0);
    FArrayBox fab_VFx(gbx2,1,amrex::The_Async_Arena()); //fab_VFx.setVal(0.0);
    FArrayBox fab_VFe(gbx2,1,amrex::The_Async_Arena()); //fab_VFe.setVal(0.0);

    auto Huxx=fab_Huxx.array();
    auto Hvxx=fab_Hvxx.array();
    auto Huee=fab_Huee.array();
    auto Hvee=fab_Hvee.array();
    auto uxx=fab_uxx.array();
    auto uee=fab_uee.array();
    auto vxx=fab_vxx.array();
    auto vee=fab_vee.array();
    auto UFx=fab_UFx.array();
    auto UFe=fab_UFe.array();
    auto VFx=fab_VFx.array();
    auto VFe=fab_VFe.array();

    //check this////////////
    const Real Gadv = 1.0;
    //uxx, uold AKA grad, ubar
    //Huxx, Huon AKA Dgrad, Duon
    //nrhs AKA krhs
    amrex::ParallelFor(gbx2,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        Huee(i,j,k)=0.0;
        uee(i,j,k)=0.0;
        Hvee(i,j,k)=0.0;
        vee(i,j,k)=0.0;

        Huxx(i,j,k)=0.0;
        uxx(i,j,k)=0.0;
        Hvxx(i,j,k)=0.0;
        vxx(i,j,k)=0.0;

        UFx(i,j,k)=0.0;
        UFe(i,j,k)=0.0;
        VFx(i,j,k)=0.0;
        VFe(i,j,k)=0.0;
    });

    amrex::ParallelFor(gbx1,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        //should not include grow cells
        uxx(i,j,k)=uold(i-1,j,k,nrhs)-2.0*uold(i,j,k,nrhs)+uold(i+1,j,k,nrhs);

        //neglecting terms about periodicity since testing only periodic for now
        Huxx(i,j,k)=Huon(i-1,j,k)-2.0*Huon(i,j,k)+Huon(i+1,j,k);
    });
    amrex::ParallelFor(gbx1,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        Real cff=1.0/6.0;
        Real cff1=uold(i  ,j,k,nrhs)+uold(i+1,j,k,nrhs);

        Real cff3=uxx(i,j,k)+uxx(i+1,j,k);

        UFx(i,j,k)=0.25*(cff1-cff*cff3) * (Huon(i,j,k)+ Huon(i+1,j,k)-cff*(Huxx(i,j,k)+ Huxx(i+1,j,k)));

        //should not include grow cells
        uee(i,j,k)=uold(i,j-1,k,nrhs)-2.0*uold(i,j,k,nrhs)+uold(i,j+1,k,nrhs);
    });

    amrex::ParallelFor(gbx1,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        /////////////MIGHT NEED NEW LOOP HERE
        //neglecting terms about periodicity since testing only periodic for now
        Hvxx(i,j,k)=Hvom(i-1,j,k)-2.0*Hvom(i,j,k)+Hvom(i+1,j,k);
    });

    amrex::ParallelFor(gbx1,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
            Real cff=1.0/6.0;
            Real cff1=uold(i,j  ,k,nrhs)+uold(i,j-1,k,nrhs);
            Real cff2=Hvom(i,j,k)+Hvom(i-1,j,k);
            Real cff3=uee(i,j-1,k)+uee(i,j,k);

            UFe(i,j,k)=0.25*(cff1-cff3*cff)*
              (cff2-cff*(Hvxx(i  ,j,k)+Hvxx(i-1,j,k)));

            vxx(i,j,k)=vold(i-1,j,k,nrhs)-2.0*vold(i,j,k,nrhs)+
              vold(i+1,j,k,nrhs);
            //neglecting terms about periodicity since testing only periodic for now
            Huee(i,j,k)=Huon(i,j-1,k)-2.0*Huon(i,j,k)+Huon(i,j+1,k);
            cff1=vold(i  ,j,k,nrhs)+vold(i-1,j,k,nrhs);
            cff2=Huon(i,j,k)+Huon(i,j-1,k);
            auto vxx_im1 = (i == gbx1.smallEnd(0)) ? vxx(i-1,j,k) :
                (vold(i-2,j,k,nrhs)-2.0*vold(i-1,j,k,nrhs)+vold(i,j,k,nrhs));
            cff3=vxx_im1+vxx(i,j,k);

            auto Huee_jm1 = (j == gbx1.smallEnd(1)) ? Huee(i,j-1,k) :
                (Huon(i,j-2,k)-2.0*Huon(i,j-1,k)+Huon(i,j,k));
            VFx(i,j,k)=0.25*(cff1-cff3*cff)* (cff2-cff*(Huee(i,j,k)+ Huee_jm1));
            vee(i,j,k)=vold(i,j-1,k,nrhs)-2.0*vold(i,j,k,nrhs)+
              vold(i,j+1,k,nrhs);
            Hvee(i,j,k)=Hvom(i,j-1,k)-2.0*Hvom(i,j,k)+Hvom(i,j+1,k);
        });

        amrex::ParallelFor(gbx1,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            //neglecting terms about periodicity since testing only periodic for now
            Real cff=1.0/6.0;
            Real cff1=vold(i,j  ,k,nrhs)+vold(i,j+1,k,nrhs);
            Real cff3=vee(i,j,k)+vee(i,j+1,k);

            VFe(i,j,k) = 0.25 * (cff1-cff3*cff) * (Hvom(i,j  ,k)+ Hvom(i,j+1,k) -
                                           cff  * (Hvee(i,j  ,k)+ Hvee(i,j+1,k)));
        });

        amrex::ParallelFor(gbx1,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
              //
              //  Add in horizontal advection.
              //
              Real cff1=UFx(i,j  ,k)-UFx(i-1,j,k);
              Real cff2=UFe(i,j+1,k)-UFe(i  ,j,k);
              Real cff=cff1+cff2;
              //passing in rhs_ubar so need 0 index here
              ru(i,j,k) -= cff;

              cff1=VFx(i+1,j,k)-VFx(i  ,j,k);
              cff2=VFe(i  ,j,k)-VFe(i,j-1,k);
              cff=cff1+cff2;
              //passing in rhs_vbar so need 0 index here
              rv(i,j,k) -= cff;
        });

}
