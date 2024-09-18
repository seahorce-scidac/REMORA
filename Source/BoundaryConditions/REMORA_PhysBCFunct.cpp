#include "AMReX_PhysBCFunct.H"
#include "IndexDefines.H"
#include <REMORA_PhysBCFunct.H>

using namespace amrex;

//
// mf is the multifab to be filled
// msk is the land/sea mask for the variable
// icomp is the index into the MultiFab -- if cell-centered this can be any value
//       from 0 to NCONS-1, if face-centered can be any value from 0 to 2 (inclusive)
// ncomp is the number of components -- if cell-centered this can be any value
//       from 1 to NCONS as long as icomp+ncomp <= NCONS-1.  If face-centered this
//       must be 1
// nghost is how many ghost cells to be filled
// time is the time at which the data should be filled
// bccomp is the index into both domain_bcs_type_bcr and bc_extdir_vals for icomp = 0  --
//     so this follows the BCVars enum
// mf_calc is the multifab for the variable used in calculations of boundary, if needed
//
void REMORAPhysBCFunct::operator() (MultiFab& mf, const MultiFab& msk, int icomp, int ncomp, IntVect const& nghost,
                                   Real time, int bccomp,int n_not_fill, const MultiFab& mf_calc,
                                   const MultiFab& mf_msku, const MultiFab& mf_mskv)
{
    if (m_geom.isAllPeriodic()) return;

    BL_PROFILE("REMORAPhysBCFunct::()");

    const auto& domain = m_geom.Domain();
    const auto dxInv   = m_geom.InvCellSizeArray();

    // Create a grown domain box containing valid + periodic cells
    Box gdomain = amrex::convert(domain, mf.boxArray().ixType());
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        if (m_geom.isPeriodic(i) and bccomp != BCVars::foextrap_bc) {
            gdomain.grow(i, nghost[i]);
        }
    }
    const bool null_mf_calc = (mf_calc.ok()) ? false : true;
    const bool null_mf_msku = (mf_msku.ok()) ? false : true;
    const bool null_mf_mskv = (mf_mskv.ok()) ? false : true;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        {
            if (mf.boxArray()[0].ixType() == IndexType(IntVect(0,0,0))) {

                // Cell-centered arrays only
                for (MFIter mfi(mf); mfi.isValid(); ++mfi)
                {
                    const Array4<Real>& dest_arr = mf.array(mfi);
                    const Array4<const Real>& msk_arr = msk.array(mfi);
                    //const Array4<const Real>& calc_arr = mf_calc.array(mfi);
                    const Array4<const Real>& calc_arr = (!null_mf_calc) ? mf_calc.const_array(mfi) : Array4<const Real>();
                    const Array4<const Real>& msku_arr = (!null_mf_msku) ? mf_msku.const_array(mfi) : Array4<const Real>();
                    const Array4<const Real>& mskv_arr = (!null_mf_mskv) ? mf_mskv.const_array(mfi) : Array4<const Real>();
                    Box bx = mfi.validbox(); bx.grow(nghost);

                    if (!gdomain.contains(bx)) {
                        impose_cons_bcs(dest_arr,bx,domain,dxInv,msk_arr,
                                msku_arr,mskv_arr,calc_arr,icomp,ncomp,time,bccomp,n_not_fill);
                    }
                } // mfi

            } else {

                // Face-based arrays only
                for (MFIter mfi(mf); mfi.isValid(); ++mfi)
                {
                    Box bx = mfi.validbox(); bx.grow(nghost);
                    const Array4<const Real>& msk_arr = msk.array(mfi);
                    Array4<const Real> calc_arr = Array4<const Real>();
                    if (!null_mf_calc) {
                         calc_arr = mf_calc.const_array(mfi);
                    }
                    if (!gdomain.contains(bx)) {
                        if(bx.ixType() == IndexType(IntVect(1,0,0))) {
                            const Array4<Real>& dest_arr = mf.array(mfi,icomp);
                            impose_xvel_bcs(dest_arr,bx,domain,dxInv,msk_arr,calc_arr,time,bccomp);
                        } else if (bx.ixType() == IndexType(IntVect(0,1,0))) {
                            const Array4<Real>& dest_arr = mf.array(mfi,icomp);
                            impose_yvel_bcs(dest_arr,bx,domain,dxInv,msk_arr,calc_arr,time,bccomp);
                        } else if (bx.ixType() == IndexType(IntVect(0,0,1))) {
                            const Array4<Real>& dest_arr = mf.array(mfi,icomp);
                            impose_zvel_bcs(dest_arr,bx,domain,dxInv,msk_arr,time,bccomp);
                        }
                    }
                } // mfi
            } // box type
        } // OpenMP
} // operator()
