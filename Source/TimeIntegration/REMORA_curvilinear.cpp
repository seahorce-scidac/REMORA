#include <REMORA.H>

using namespace amrex;

/**
 * @param[in   ] bx       box to update on
 * @param[in   ] xbx      nodal-x box to update on
 * @param[in   ] ybx      nodal-y box to update on
 * @param[in   ] uold     u-direction velocity
 * @param[in   ] vold     v-direction velocity
 * @param[inout] ru       u-direction velocity RHS
 * @param[inout] rv       v-direction velocity RHS
 * @param[in   ] Hz       vertical cell height
 * @param[in   ] dndx     d(1/n)/d(xi)
 * @param[in   ] dmde     d(1/m)/d(eta)
 * @param[in   ] nrhs     which velocity component to use
 * @param[in   ] nr       which RHS component to update
 */

void
REMORA::curvilinear (const Box& bx,
                 const Box& xbx,
                 const Box& ybx,
                 const Array4<Real const>& uold,
                 const Array4<Real const>& vold,
                 const Array4<Real      >& ru,
                 const Array4<Real      >& rv,
                 const Array4<Real const>& Hz,
                 const Array4<Real const>& dndx,
                 const Array4<Real const>& dmde,
                 int nrhs, int nr)
{
    BL_PROFILE("REMORA::curvilinear()");
    //
    //-----------------------------------------------------------------------
    //  Add in curvilinear transformation terms.
    //-----------------------------------------------------------------------
    //

    int ncomp = 0;
    int UFx_comp = ncomp++;
    int VFe_comp = ncomp++;

    FArrayBox fab(grow(bx,IntVect(1,1,0)),ncomp,The_Async_Arena());

    auto UFx = fab.array(UFx_comp);
    auto VFe = fab.array(VFe_comp);

    ParallelFor(grow(bx,IntVect(1,1,0)),
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        Real cff1 = 0.5_rt * (vold(i,j,k,nrhs) + vold(i  ,j+1,k,nrhs));
        Real cff2 = 0.5_rt * (uold(i,j,k,nrhs) + uold(i+1,j  ,k,nrhs));
        Real cff3 = cff1 * dndx(i,j,0);
        Real cff4 = cff2 * dmde(i,j,0);
        Real cff = Hz(i,j,k) * (cff3 - cff4);
        UFx(i,j,k) = cff * cff1;
        VFe(i,j,k) = cff * cff2;
    });

    ParallelFor(xbx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        ru(i,j,k,nr) += 0.5_rt * (UFx(i,j,k) + UFx(i-1,j,k));
    });

    ParallelFor(ybx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        rv(i,j,k,nr) -= 0.5_rt * (VFe(i,j,k) + VFe(i,j-1,k));
    });
}
