#include "AMReX_PhysBCFunct.H"
#include <REMORA_PhysBCFunct.H>

using namespace amrex;

//
// mf is the multifab to be filled
// icomp is the index into the MultiFab -- if cell-centered this can be any value
//       from 0 to NCONS-1, if face-centered this must be 0
// ncomp is the number of components -- if cell-centered (var_idx = 0) this can be any value
//       from 1 to NCONS as long as icomp+ncomp <= NCONS-1.  If face-centered this
//       must be 1
// time is the time at which the data should be filled
// bccomp is the index into both domain_bcs_type_bcr and bc_extdir_vals for icomp = 0  --
//     so this follows the BCVars enum
// n_not_fill perimeter of cells in x and y where BCs are not applied for conditions other than ext_dir.
//     Foextrap is done based on the values at dom_lo-n_not_fill and dom_hi+n_not_fill. Reflecting conditions
//     are untested.
//
void REMORAPhysBCFunct::impose_cons_bcs (const Array4<Real>& dest_arr, const Box& bx, const Box& domain,
                                        const GpuArray<Real,AMREX_SPACEDIM> /*dxInv*/, const Array4<const Real>& mskr,
                                        int icomp, int ncomp, Real /*time*/, int bccomp, int n_not_fill)
{
    BL_PROFILE_VAR("impose_cons_bcs()",impose_cons_bcs);
    const auto& dom_lo = amrex::lbound(domain);
    const auto& dom_hi = amrex::ubound(domain);

    // Based on BCRec for the domain, we need to make BCRec for this Box
    // bccomp is used as starting index for m_domain_bcs_type
    //      0 is used as starting index for bcrs
    Vector<BCRec> bcrs(ncomp);
    amrex::setBC(bx, domain, bccomp, 0, ncomp, m_domain_bcs_type, bcrs);

    // xlo: ori = 0
    // ylo: ori = 1
    // zlo: ori = 2
    // xhi: ori = 3
    // yhi: ori = 4
    // zhi: ori = 5

    amrex::Gpu::DeviceVector<BCRec> bcrs_d(ncomp);
#ifdef AMREX_USE_GPU
    Gpu::htod_memcpy_async(bcrs_d.data(), bcrs.data(), sizeof(BCRec)*ncomp);
#else
    std::memcpy(bcrs_d.data(), bcrs.data(), sizeof(BCRec)*ncomp);
#endif
    const amrex::BCRec* bc_ptr = bcrs_d.data();

    GpuArray<GpuArray<Real, AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NCONS> l_bc_extdir_vals_d;

    for (int i = 0; i < ncomp; i++) {
        for (int ori = 0; ori < 2*AMREX_SPACEDIM; ori++) {
            l_bc_extdir_vals_d[i][ori] = m_bc_extdir_vals[bccomp+i][ori];
        }
    }

    GeometryData const& geomdata = m_geom.data();
    bool is_periodic_in_x = geomdata.isPeriodic(0);
    bool is_periodic_in_y = geomdata.isPeriodic(1);


    // First do all ext_dir bcs
    if (!is_periodic_in_x)
    {
        Box bx_xlo(bx);  bx_xlo.setBig  (0,dom_lo.x-1);
        Box bx_xhi(bx);  bx_xhi.setSmall(0,dom_hi.x+1);
        ParallelFor(
            bx_xlo, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                if (bc_ptr[n].lo(0) == REMORABCType::ext_dir) {
                    dest_arr(i,j,k,icomp+n) = l_bc_extdir_vals_d[n][0] * mskr(i,j,0);
                }
            },
            bx_xhi, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                if (bc_ptr[n].hi(0) == REMORABCType::ext_dir) {
                    dest_arr(i,j,k,icomp+n) = l_bc_extdir_vals_d[n][3] * mskr(i,j,0);
                }
            }
        );
    }

    if (!is_periodic_in_y)
    {
        Box bx_ylo(bx);  bx_ylo.setBig  (1,dom_lo.y-1);
        Box bx_yhi(bx);  bx_yhi.setSmall(1,dom_hi.y+1);
        ParallelFor(
            bx_ylo, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                if (bc_ptr[n].lo(1) == REMORABCType::ext_dir) {
                    dest_arr(i,j,k,icomp+n) = l_bc_extdir_vals_d[n][1] * mskr(i,j,0);
                }
            },
            bx_yhi, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                if (bc_ptr[n].hi(1) == REMORABCType::ext_dir) {
                    dest_arr(i,j,k,icomp+n) = l_bc_extdir_vals_d[n][4] * mskr(i,j,0);
                }
            }
        );
    }

    {
        Box bx_zlo(bx);  bx_zlo.setBig  (2,dom_lo.z-1);
        Box bx_zhi(bx);  bx_zhi.setSmall(2,dom_hi.z+1);
        ParallelFor(
            bx_zlo, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                if (bc_ptr[n].lo(2) == REMORABCType::ext_dir) {
                    dest_arr(i,j,k,icomp+n) = l_bc_extdir_vals_d[n][2] * mskr(i,j,0);
                }
            },
            bx_zhi, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                if (bc_ptr[n].hi(2) == REMORABCType::ext_dir) {
                    dest_arr(i,j,k,icomp+n) = l_bc_extdir_vals_d[n][5] * mskr(i,j,0);
                }
            }
        );
    }

    // Next do ghost cells in x-direction but not reaching out in y
    // The corners we miss here will be covered in the y-loop below or by periodicity
    if (!is_periodic_in_x)
    {
        // Populate ghost cells on lo-x and hi-x domain boundaries
        Box bx_xlo(bx);  bx_xlo.setBig  (0,dom_lo.x-1-n_not_fill);
                         bx_xlo.setSmall(2,std::max(dom_lo.z,bx.smallEnd(2)));
                         bx_xlo.setBig  (2,std::min(dom_hi.z,bx.bigEnd(2)));
        Box bx_xhi(bx);  bx_xhi.setSmall(0,dom_hi.x+1+n_not_fill);
                         bx_xhi.setSmall(2,std::max(dom_lo.z,bx.smallEnd(2)));
                         bx_xhi.setBig  (2,std::min(dom_hi.z,bx.bigEnd(2)));
        ParallelFor(bx_xlo, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                int iflip = dom_lo.x - 1 - i;
                if (bc_ptr[n].lo(0) == REMORABCType::foextrap) {
                    dest_arr(i,j,k,icomp+n) =  dest_arr(dom_lo.x-n_not_fill,j,k,icomp+n);
                } else if (bc_ptr[n].lo(0) == REMORABCType::reflect_even) {
                    dest_arr(i,j,k,icomp+n) =  dest_arr(iflip,j,k,icomp+n);
                } else if (bc_ptr[n].lo(0) == REMORABCType::reflect_odd) {
                    dest_arr(i,j,k,icomp+n) = -dest_arr(iflip,j,k,icomp+n);
                }
            },
            bx_xhi, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                int iflip =  2*dom_hi.x + 1 - i;
                if (bc_ptr[n].hi(0) == REMORABCType::foextrap) {
                    dest_arr(i,j,k,icomp+n) =  dest_arr(dom_hi.x+n_not_fill,j,k,icomp+n);
                } else if (bc_ptr[n].hi(0) == REMORABCType::reflect_even) {
                    dest_arr(i,j,k,icomp+n) =  dest_arr(iflip,j,k,icomp+n);
                } else if (bc_ptr[n].hi(0) == REMORABCType::reflect_odd) {
                    dest_arr(i,j,k,icomp+n) = -dest_arr(iflip,j,k,icomp+n);
                }
            }
        );
    }

    if (!is_periodic_in_y)
    {
        // Populate ghost cells on lo-y and hi-y domain boundaries
        Box bx_ylo(bx);  bx_ylo.setBig  (1,dom_lo.y-1-n_not_fill);
                         bx_ylo.setSmall(2,std::max(dom_lo.z,bx.smallEnd(2)));
                         bx_ylo.setBig  (2,std::min(dom_hi.z,bx.bigEnd(2)));
        Box bx_yhi(bx);  bx_yhi.setSmall(1,dom_hi.y+1+n_not_fill);
                         bx_yhi.setSmall(2,std::max(dom_lo.z,bx.smallEnd(2)));
                         bx_yhi.setBig  (2,std::min(dom_hi.z,bx.bigEnd(2)));
        ParallelFor(
            bx_ylo, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                int jflip = dom_lo.y - 1 - j;
                if (bc_ptr[n].lo(1) == REMORABCType::foextrap) {
                    dest_arr(i,j,k,icomp+n) =  dest_arr(i,dom_lo.y-n_not_fill,k,icomp+n);
                } else if (bc_ptr[n].lo(1) == REMORABCType::reflect_even) {
                    dest_arr(i,j,k,icomp+n) =  dest_arr(i,jflip,k,icomp+n);
                } else if (bc_ptr[n].lo(1) == REMORABCType::reflect_odd) {
                    dest_arr(i,j,k,icomp+n) = -dest_arr(i,jflip,k,icomp+n);
                }
            },
            bx_yhi, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                int jflip =  2*dom_hi.y + 1 - j;
                if (bc_ptr[n].hi(1) == REMORABCType::foextrap) {
                    dest_arr(i,j,k,icomp+n) =  dest_arr(i,dom_hi.y+n_not_fill,k,icomp+n);
                } else if (bc_ptr[n].hi(1) == REMORABCType::reflect_even) {
                    dest_arr(i,j,k,icomp+n) =  dest_arr(i,jflip,k,icomp+n);
                } else if (bc_ptr[n].hi(1) == REMORABCType::reflect_odd) {
                    dest_arr(i,j,k,icomp+n) = -dest_arr(i,jflip,k,icomp+n);
                }
            }
        );
    }

    {
        Box bx_zlo(bx);  bx_zlo.setBig  (2,std::max(dom_lo.z-1,bx.smallEnd(2)));
        Box bx_zhi(bx);  bx_zhi.setSmall(2,std::min(dom_hi.z+1,bx.bigEnd(2)));
        // Populate ghost cells on lo-z and hi-z domain boundaries

        if (bx_zlo.ok()) {
            ParallelFor(bx_zlo, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                int kflip = dom_lo.z - 1 - i;
                if (bc_ptr[n].lo(2) == REMORABCType::foextrap) {
                    dest_arr(i,j,k,icomp+n) =  dest_arr(i,j,dom_lo.z,icomp+n);
                } else if (bc_ptr[n].lo(2) == REMORABCType::reflect_even) {
                    dest_arr(i,j,k,icomp+n) =  dest_arr(i,j,kflip,icomp+n);
                } else if (bc_ptr[n].lo(2) == REMORABCType::reflect_odd) {
                    dest_arr(i,j,k,icomp+n) = -dest_arr(i,j,kflip,icomp+n);
                }
            });
        }

        if (bx_zlo.ok()) {
            ParallelFor(bx_zhi, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                int kflip =  2*dom_hi.z + 1 - i;
                if (bc_ptr[n].hi(2) == REMORABCType::foextrap) {
                    dest_arr(i,j,k,icomp+n) =  dest_arr(i,j,dom_hi.z,icomp+n);
                } else if (bc_ptr[n].hi(2) == REMORABCType::reflect_even) {
                    dest_arr(i,j,k,icomp+n) =  dest_arr(i,j,kflip,icomp+n);
                } else if (bc_ptr[n].hi(2) == REMORABCType::reflect_odd) {
                    dest_arr(i,j,k,icomp+n) = -dest_arr(i,j,kflip,icomp+n);
                }
            });
        }
    }

    Gpu::streamSynchronize();
}
