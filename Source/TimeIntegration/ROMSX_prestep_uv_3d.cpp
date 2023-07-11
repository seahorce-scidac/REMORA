#include <ROMSX.H>
#include <Utils.H>

using namespace amrex;

//
// prestep_uv_3d
//

void
ROMSX::prestep_uv_3d (const Box& tbx, const Box& gbx,
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
    Box gbx1 = tbx;
    Box gbx2 = tbx;
    Box tbxp1 = tbx;
    Box tbxp2 = tbx;

    tbxp1.grow(IntVect(NGROW-1,NGROW-1,0));
    tbxp2.grow(IntVect(NGROW,NGROW,0));

    //make only gbx be grown to match multifabs
    gbx2.grow(IntVect(NGROW,NGROW,0));
    gbx1.grow(IntVect(NGROW-1,NGROW-1,0));

    if (verbose > 0) {
        amrex::AllPrint() << "Box(Huon) " << Box(Huon) << std::endl;
        amrex::AllPrint() << "Box(Hvom) " << Box(Hvom) << std::endl;
    }

    //Need to include pre_step3d.F terms

    update_vel_3d(tbxp1, gbx, 1, 0, u, uold, ru, Hz, Akv, DC, FC,
                  sustr, bustr, z_r, pm, pn, iic, ntfirst, nnew, nstp, nrhs, N, lambda, dt_lev);

    update_vel_3d(tbxp1, gbx, 0, 1, v, vold, rv, Hz, Akv, DC, FC,
                  svstr, bvstr, z_r, pm, pn, iic, ntfirst, nnew, nstp, nrhs, N, lambda, dt_lev);
    //    //Print()<<FArrayBox(uold)<<std::endl;
    //    //Print()<<FArrayBox(u)<<std::endl;
}
