#include <ROMSX.H>

using namespace amrex;

//
// rhs_2d
//
/**
 * rhs_2d
 *
 * @param[in] bx
 * @param[in] uold
 * @param[in] vold
 * @param[in] ru
 * @param[in] rv
 * @param[in] mf_ru
 * @param[in] mf_rv
 * @param[inout] mf_rufrc
 * @param[inout] mf_rvfrc
 * @param[inout] mf_Zt_avg1
 * @param[inout] mf_DU_avg1
 * @param[inout] mf_DU_avg2
 * @param[inout] mf_DV_avg1
 * @param[inout] mf_DV_avg2
 * @param[inout] mf_rubar
 * @param[inout] mf_rvbar
 * @param[inout] mf_rbar
 * @param[inout] mf_rbar
 * @param[inout] mf_zeta
 * @param[inout] mf_h
 * @param[inout] mf_visc2_p
 * @param[inout] mf_visc2_f
 * @param[in   ] ncomp
 * @param[in   ] dtfast_lev
 * @param[in   ] predictor_2d_step
 * @param[in   ] first_2d_step
 * @param[in   ] my_iif
 * @param[in   ] next_indx1
 */

void
ROMSX::rhs_2d (const Box& bx,
               Array4<Real> uold  , Array4<Real> vold,
               Array4<Real> ru, Array4<Real> rv,
               Array4<Real> Huon, Array4<Real> Hvom,
               int nrhs)
{
    //copy the tilebox
    Box gbx1 = bx;
    Box gbx2 = bx;

    //make only gbx be grown to match multifabs
    gbx2.grow(IntVect(NGROW,NGROW,0));
    gbx1.grow(IntVect(NGROW-1,NGROW-1,0));

    //
    // Scratch space
    //
    FArrayBox fab_Huee(gbx2,1,amrex::The_Async_Arena());
    FArrayBox fab_uee(gbx2,1,amrex::The_Async_Arena());

    FArrayBox fab_Hvee(gbx2,1,amrex::The_Async_Arena());
    FArrayBox fab_vee(gbx2,1,amrex::The_Async_Arena());

    FArrayBox fab_Hvxx(gbx2,1,amrex::The_Async_Arena());
    FArrayBox fab_uxx(gbx2,1,amrex::The_Async_Arena());

    FArrayBox fab_Huxx(gbx2,1,amrex::The_Async_Arena());
    FArrayBox fab_vxx(gbx2,1,amrex::The_Async_Arena());

    FArrayBox fab_UFx(gbx2,1,amrex::The_Async_Arena());
    FArrayBox fab_UFe(gbx2,1,amrex::The_Async_Arena());
    FArrayBox fab_VFx(gbx2,1,amrex::The_Async_Arena());
    FArrayBox fab_VFe(gbx2,1,amrex::The_Async_Arena());

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

    fab_Huee.template setVal<RunOn::Device>(0.);
    fab_Huxx.template setVal<RunOn::Device>(0.);
    fab_uee.template setVal<RunOn::Device>(0.);
    fab_uxx.template setVal<RunOn::Device>(0.);
    fab_UFx.template setVal<RunOn::Device>(0.);
    fab_UFe.template setVal<RunOn::Device>(0.);

    fab_Hvee.template setVal<RunOn::Device>(0.);
    fab_Hvxx.template setVal<RunOn::Device>(0.);
    fab_vee.template setVal<RunOn::Device>(0.);
    fab_vxx.template setVal<RunOn::Device>(0.);
    fab_VFx.template setVal<RunOn::Device>(0.);
    fab_VFe.template setVal<RunOn::Device>(0.);

    ParallelFor(gbx1,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        //should not include grow cells
        uxx(i,j,k)=uold(i-1,j,k,nrhs)-2.0*uold(i,j,k,nrhs)+uold(i+1,j,k,nrhs);

        //neglecting terms about periodicity since testing only periodic for now
        Huxx(i,j,k)=Huon(i-1,j,k)-2.0*Huon(i,j,k)+Huon(i+1,j,k);
    });
    ParallelFor(gbx1,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        Real cff=1.0/6.0;
        Real cff1=uold(i  ,j,k,nrhs)+uold(i+1,j,k,nrhs);

        Real cff3=uxx(i,j,k)+uxx(i+1,j,k);

        UFx(i,j,k)=0.25*(cff1-cff*cff3) * (Huon(i,j,k)+ Huon(i+1,j,k)-cff*(Huxx(i,j,k)+ Huxx(i+1,j,k)));

        //should not include grow cells
        uee(i,j,k)=uold(i,j-1,k,nrhs)-2.0*uold(i,j,k,nrhs)+uold(i,j+1,k,nrhs);
    });

    ParallelFor(gbx1,
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

              ru(i,j,k) -= cff;

              cff1=VFx(i+1,j,k)-VFx(i  ,j,k);
              cff2=VFe(i  ,j,k)-VFe(i,j-1,k);
              cff=cff1+cff2;

              rv(i,j,k) -= cff;
        });

}
