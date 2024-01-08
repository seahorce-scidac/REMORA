#include "AMReX_PhysBCFunct.H"
#include "IndexDefines.H"
#include <REMORA_PhysBCFunct.H>

using namespace amrex;

//
// mf is the multifab to be filled
// icomp is the index into the MultiFab -- if cell-centered this can be any value
//       from 0 to NCONS-1, if face-centered this must be 0
// ncomp is the number of components -- if cell-centered this can be any value
//       from 1 to NCONS as long as icomp+ncomp <= NCONS-1.  If face-centered this
//       must be 1
// nghost is how many ghost cells to be filled
// time is the time at which the data should be filled
// bccomp is the index into both domain_bcs_type_bcr and bc_extdir_vals for icomp = 0  --
//     so this follows the BCVars enum
//
void REMORAPhysBCFunct::operator() (MultiFab& mf, int icomp, int ncomp, IntVect const& nghost,
                                   Real time, int bccomp)
{
    if (m_geom.isAllPeriodic()) return;

    BL_PROFILE("REMORAPhysBCFunct::()");

    const auto& domain = m_geom.Domain();
    const auto dxInv   = m_geom.InvCellSizeArray();

    // Create a grown domain box containing valid + periodic cells
    Box gdomain = amrex::convert(domain, mf.boxArray().ixType());
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        if (m_geom.isPeriodic(i)) {
            gdomain.grow(i, nghost[i]);
        }
    }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        {
            if (mf.boxArray()[0].ixType() == IndexType(IntVect(0,0,0))) {

                // Cell-centered arrays only
                for (MFIter mfi(mf); mfi.isValid(); ++mfi)
                {
                    const Array4<Real>& dest_arr = mf.array(mfi);
                    Box bx = mfi.validbox(); bx.grow(nghost);

                    if (!gdomain.contains(bx)) {
                        AMREX_ALWAYS_ASSERT(icomp == 0 && icomp+ncomp <= NCONS);
                        impose_cons_bcs(dest_arr,bx,domain,dxInv,icomp,ncomp,time,bccomp);
                    }
                } // mfi

            } else {

                // Face-based arrays only
                for (MFIter mfi(mf); mfi.isValid(); ++mfi)
                {
                    Box bx = mfi.validbox(); bx.grow(nghost);
                    if (!gdomain.contains(bx)) {
                        if(bx.ixType() == IndexType(IntVect(1,0,0))) {
                            const Array4<Real>& dest_arr = mf.array(mfi,0);
                            impose_xvel_bcs(dest_arr,bx,domain,dxInv,time,bccomp);
                        } else if (bx.ixType() == IndexType(IntVect(0,1,0))) {
                            const Array4<Real>& dest_arr = mf.array(mfi,0);
                            impose_yvel_bcs(dest_arr,bx,domain,dxInv,time,bccomp);
                        }
                    }
                } // mfi
            } // box type
        } // OpenMP
} // operator()
