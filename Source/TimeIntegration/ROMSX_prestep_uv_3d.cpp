#include <ROMSX.H>
#include <Utils.H>

using namespace amrex;

//
// prestep_uv_3d
//

void
ROMSX::prestep_uv_3d (const Box& bx,
                      Array4<Real> uold  , Array4<Real> vold,
                      Array4<Real> u , Array4<Real> v,
                      Array4<Real> ru, Array4<Real> rv,
                      Array4<Real> Hz, Array4<Real> Akv,
                      Array4<Real> on_u, Array4<Real> om_v,
                      Array4<Real> Huon, Array4<Real> Hvom,
                      Array4<Real> pm, Array4<Real> pn,
                      Array4<Real> W   , Array4<Real> DC,
                      Array4<Real> FC  , Array4<Real> z_r,
                      Array4<Real> sustr, Array4<Real> svstr,
                      Array4<Real> bustr, Array4<Real> bvstr,
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
    gbx2.grow(IntVect(NGROW,NGROW,0));
    gbx1.grow(IntVect(NGROW-1,NGROW-1,0));
    gbx11.grow(IntVect(NGROW-1,NGROW-1,NGROW-1));

    amrex::AllPrint() << "Box(Huon) " << Box(Huon) << std::endl;
    amrex::AllPrint() << "Box(Hvom) " << Box(Hvom) << std::endl;

    //Need to include pre_step3d.F terms

    update_vel_3d(ubx, 1, 0, u, uold, ru, Hz, Akv, DC, FC,
                  sustr, bustr, z_r, pm, pn, iic, ntfirst, nnew, nstp, nrhs, N, lambda, dt_lev);

    update_vel_3d(vbx, 0, 1, v, vold, rv, Hz, Akv, DC, FC,
                  svstr, bvstr, z_r, pm, pn, iic, ntfirst, nnew, nstp, nrhs, N, lambda, dt_lev);
    //    //Print()<<FArrayBox(uold)<<std::endl;
    //    //Print()<<FArrayBox(u)<<std::endl;
}
