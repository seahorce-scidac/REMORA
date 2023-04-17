#include <ROMSX.H>
#include <Utils.H>

using namespace amrex;

//
// prestep_uv_3d
//

void
ROMSX::prestep_uv_3d (const Box& bx,
                      Array4<Real> uold  , Array4<Real> vold,
                      Array4<Real> u_arr , Array4<Real> v_arr,
                      Array4<Real> ru_arr, Array4<Real> rv_arr,
                      Array4<Real> Hz_arr, Array4<Real> Akv_arr,
                      Array4<Real> on_u, Array4<Real> om_v,
                      Array4<Real> Huon, Array4<Real> Hvom,
                      Array4<Real> pm_arr, Array4<Real> pn_arr,
                      Array4<Real> W   , Array4<Real> DC_arr,
                      Array4<Real> FC_arr  , Array4<Real> z_r_arr,
                      Array4<Real> sustr_arr, Array4<Real> svstr_arr,
                      int iic, int ntfirst, int nnew, int nstp, int nrhs, int N,
                      Real lambda, Real dt_lev)
{
    //copy the tilebox
    Box gbx1 = bx;
    Box gbx11 = bx;
    Box gbx2 = bx;

    Box ubx = surroundingNodes(bx,0);
    Box vbx = surroundingNodes(bx,1);

    Box gbx3uneven(IntVect(AMREX_D_DECL(bx.smallEnd(0)-3,bx.smallEnd(1)-3,bx.smallEnd(2))),
                   IntVect(AMREX_D_DECL(bx.bigEnd(0)+2,bx.bigEnd(1)+2,bx.bigEnd(2))));
    Box gbx2uneven(IntVect(AMREX_D_DECL(bx.smallEnd(0)-2,bx.smallEnd(1)-2,bx.smallEnd(2))),
                   IntVect(AMREX_D_DECL(bx.bigEnd(0)+1,bx.bigEnd(1)+1,bx.bigEnd(2))));
    //make only gbx be grown to match multifabs
    gbx2.grow(IntVect(2,2,0));
    gbx1.grow(IntVect(1,1,0));
    gbx11.grow(IntVect(1,1,1));

    amrex::AllPrint() << "Box(Huon) " << Box(Huon) << std::endl;
    amrex::AllPrint() << "Box(Hvom) " << Box(Hvom) << std::endl;

    //
    //------------------------------------------------------------------------
    //  Vertically integrate horizontal mass flux divergence.
    //------------------------------------------------------------------------
    //
    //Should really use gbx3uneven
    amrex::ParallelFor(gbx2uneven,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        //  Starting with zero vertical velocity at the bottom, integrate
        //  from the bottom (k=0) to the free-surface (k=N).  The w(:,:,N(ng))
        //  contains the vertical velocity at the free-surface, d(zeta)/d(t).
        //  Notice that barotropic mass flux divergence is not used directly.
        //
        if(k==0) {
            W(i,j,k)=0.0;
        } else {
            W(i,j,k) = W(i,j,k-1)- (Huon(i+1,j,k)-Huon(i,j,k)+ Hvom(i,j+1,k)-Hvom(i,j,k));
        }
    });

    //Need to include pre_step3d.F terms

    update_vel_3d(ubx, 1, 0, u_arr, uold, ru_arr, Hz_arr, Akv_arr, DC_arr, FC_arr,
                  sustr_arr, z_r_arr, pm_arr, pn_arr, iic, ntfirst, nnew, nstp, nrhs, N, lambda, dt_lev);

    update_vel_3d(vbx, 0, 1, v_arr, vold, rv_arr, Hz_arr, Akv_arr, DC_arr, FC_arr,
                  svstr_arr, z_r_arr, pm_arr, pn_arr, iic, ntfirst, nnew, nstp, nrhs, N, lambda, dt_lev);
    //    Print()<<FArrayBox(uold)<<std::endl;
    //    Print()<<FArrayBox(u_arr)<<std::endl;
}
