#include "AMReX_PhysBCFunct.H"
#include "IndexDefines.H"
#include <ROMSX_PhysBCFunct.H>

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
void ROMSXPhysBCFunct::operator() (MultiFab& mf, int icomp, int ncomp, IntVect const& nghost,
                                   Real time, int bccomp)
{
    if (m_geom.isAllPeriodic()) return;

    BL_PROFILE("ROMSXPhysBCFunct::()");

    const auto& domain = m_geom.Domain();
    const auto dx      = m_geom.CellSizeArray();
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
            Vector<BCRec> bcrs(ncomp);

            // Do all BCs except MOST
            for (MFIter mfi(mf); mfi.isValid(); ++mfi)
            {
                const Array4<Real>& dest_arr = mf.array(mfi);
                Box bx = mfi.validbox(); bx.grow(nghost);

                Array4<const Real> velx_arr;
                Array4<const Real> vely_arr;

                //! if there are cells not in the valid + periodic grown box
                //! we need to fill them here
                //!
                if (!gdomain.contains(bx) || (mf[0].box().ixType() == IndexType(IntVect(0,0,1))) )
                {
                    if (mf[0].box().ixType() == IndexType(IntVect(1,0,0)))
                    {
                        AMREX_ALWAYS_ASSERT(ncomp == 1 && icomp == 0);
                        impose_xvel_bcs(dest_arr,bx,domain,
                                        dxInv,time,bccomp);

                    } else if (mf[0].box().ixType() == IndexType(IntVect(0,1,0)))
                    {
                        AMREX_ALWAYS_ASSERT(ncomp == 1 && icomp == 0);
                        impose_yvel_bcs(dest_arr,bx,domain,
                                        dxInv,time,bccomp);

                    } else if (mf[0].box().ixType() == IndexType(IntVect(0,0,1)))
                    {
                        AMREX_ALWAYS_ASSERT(ncomp == 1 && icomp == 0);
                        impose_zvel_bcs(dest_arr,bx,domain,
                                        velx_arr,vely_arr,dx,dxInv,
                                        time,bccomp);

                    } else if (mf[0].box().ixType() == IndexType(IntVect(0,0,0)))
                    {
                        AMREX_ALWAYS_ASSERT(icomp == 0 && icomp+ncomp <= NCONS);
                        impose_cons_bcs(dest_arr,bx,domain,
                                        dxInv,icomp,ncomp,time,bccomp);
                    } else {
                        amrex::Abort("Dont know this box type in ROMSX_PhysBC");
                    }

                    // ****************************************************************************
                    // Based on BCRec for the domain, we need to make BCRec for this Box
                    // bccomp is used as starting index for m_domain_bcs_type
                    //      0 is used as starting index for bcrs
                    // ****************************************************************************
                    amrex::setBC(bx, domain, bccomp, 0, ncomp, m_domain_bcs_type, bcrs);

                    // xlo: ori = 0
                    // ylo: ori = 1
                    // zlo: ori = 2
                    // xhi: ori = 3
                    // yhi: ori = 4
                    // zhi: ori = 5

                    amrex::Gpu::DeviceVector<BCRec> bcrs_d(ncomp);
#ifdef AMREX_USE_GPU
                    Gpu::htod_memcpy_async
                        (bcrs_d.data(), bcrs.data(), sizeof(BCRec)*ncomp);
#else
                    std::memcpy
                        (bcrs_d.data(), bcrs.data(), sizeof(BCRec)*ncomp);
#endif

                        Gpu::streamSynchronize(); // because of bcrs_d
                } // !gdomain.contains(bx)
            } // MFIter
        } // OpenMP
} // operator()
