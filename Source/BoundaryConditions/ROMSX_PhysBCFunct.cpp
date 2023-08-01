#include "AMReX_PhysBCFunct.H"
#include "IndexDefines.H"
#include <ROMSX_PhysBCFunct.H>

using namespace amrex;

//
// mf is the multifab to be filled
// icomp is the index into the MultiFab -- if cell-centered this can be any value
//       from 0 to NVAR-1, if face-centered this must be 0
// ncomp is the number of components -- if cell-centered (var_idx = 0) this can be any value
//       from 1 to NVAR as long as icomp+ncomp <= NVAR-1.  If face-centered this
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

    MultiFab* xvel_ptr   = nullptr;
    MultiFab* yvel_ptr   = nullptr;

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
                if (!gdomain.contains(bx) || (m_var_idx == Vars::zvel))
                {
                    if (m_var_idx == Vars::xvel) {
                        AMREX_ALWAYS_ASSERT(ncomp == 1 && icomp == 0);
                        impose_xvel_bcs(dest_arr,bx,domain,
                                        dxInv,time,bccomp);

                    } else if (m_var_idx == Vars::yvel) {
                        AMREX_ALWAYS_ASSERT(ncomp == 1 && icomp == 0);
                        impose_yvel_bcs(dest_arr,bx,domain,
                                        dxInv,time,bccomp);

                    } else if (m_var_idx == Vars::zvel) {
                        AMREX_ALWAYS_ASSERT(ncomp == 1 && icomp == 0);
                        impose_zvel_bcs(dest_arr,bx,domain,
                                        velx_arr,vely_arr,dx,dxInv,
                                        time,bccomp);

                    } else if (m_var_idx == Vars::cons) {
                        AMREX_ALWAYS_ASSERT(icomp == 0 && icomp+ncomp <= NVAR);
                        impose_cons_bcs(dest_arr,bx,domain,
                                        dxInv,icomp,ncomp,time,bccomp);
                    } else {
                        amrex::Abort("Dont know this var_idx in ROMSX_PhysBC");
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

#ifdef ROMSX_USE_NETCDF
                    const amrex::BCRec* bc_ptr = bcrs_d.data();
                    if (m_init_type == "real") {
                        int icomp_for_wrfbdy, ncomp_for_wrfbdy, bccomp_for_wrfbdy;
                        if (m_var_idx == Vars::cons) {
                            icomp_for_wrfbdy = Temp_comp;
                            bccomp_for_wrfbdy = BCVars::Temp_bc_comp;
                            ncomp_for_wrfbdy = 1; // (Because we are currently only filling U, V, W, T)
                        } else {
                            icomp_for_wrfbdy = icomp;
                            bccomp_for_wrfbdy = bccomp;
                            ncomp_for_wrfbdy = 1; // (Because we are currently only filling U, V, W, T)
                        }
                        fill_from_wrfbdy(m_lev, bx, dest_arr, icomp_for_wrfbdy, bccomp_for_wrfbdy, ncomp_for_wrfbdy,
                                         domain, bc_ptr,
                                         time, m_bdy_time_interval);
                    }
#endif
                        Gpu::streamSynchronize(); // because of bcrs_d
                } // !gdomain.contains(bx)
            } // MFIter
        } // OpenMP
    } // operator()
