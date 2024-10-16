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
void REMORAPhysBCFunct::impose_cons_bcs (const Array4<Real>& dest_arr, const Box& bx, const Box& valid_bx, const Box& domain,
                                        const GpuArray<Real,AMREX_SPACEDIM> /*dxInv*/, const Array4<const Real>& mskr,
                                        const Array4<const Real>& msku, const Array4<const Real>& mskv,
                                        const Array4<const Real>& calc_arr,
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

    GpuArray<GpuArray<Real, AMREX_SPACEDIM*2>,AMREX_SPACEDIM+NCONS+8> l_bc_extdir_vals_d;

    for (int i = 0; i < ncomp; i++) {
        for (int ori = 0; ori < 2*AMREX_SPACEDIM; ori++) {
            l_bc_extdir_vals_d[i][ori] = m_bc_extdir_vals[bccomp+i][ori];
        }
    }

    GeometryData const& geomdata = m_geom.data();
    bool is_periodic_in_x = geomdata.isPeriodic(0);
    bool is_periodic_in_y = geomdata.isPeriodic(1);
    const Real eps= 1.0e-20_rt;

    // First do all ext_dir bcs
    if (!is_periodic_in_x)
    {
        Box bx_xlo(bx);  bx_xlo.setBig  (0,dom_lo.x-1);
        bx_xlo.setSmall(1,valid_bx.smallEnd(1));  bx_xlo.setBig(1,valid_bx.bigEnd(1));
        Box bx_xhi(bx);  bx_xhi.setSmall(0,dom_hi.x+1);
        bx_xhi.setSmall(1,valid_bx.smallEnd(1));  bx_xhi.setBig(1,valid_bx.bigEnd(1));
        ParallelFor(
            bx_xlo, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                if (bc_ptr[n].lo(0) == REMORABCType::ext_dir) {
                    dest_arr(i,j,k,icomp+n) = l_bc_extdir_vals_d[n][0] * mskr(i,j,0);
                } else if (bc_ptr[n].lo(0) == REMORABCType::orlanski_rad) {
                    Real grad_lo_im1   = (calc_arr(dom_lo.x-1,j  ,k,icomp+n) - calc_arr(dom_lo.x-1,j-1,k,icomp+n)) * mskv(i,j,0);
                    Real grad_lo       = (calc_arr(dom_lo.x  ,j  ,k,icomp+n) - calc_arr(dom_lo.x  ,j-1,k,icomp+n)) * mskv(i,j,0);
                    Real grad_lo_imjp1 = (calc_arr(dom_lo.x-1,j+1,k,icomp+n) - calc_arr(dom_lo.x-1,j  ,k,icomp+n)) * mskv(i,j,0);
                    Real grad_lo_jp1   = (calc_arr(dom_lo.x  ,j+1,k,icomp+n) - calc_arr(dom_lo.x  ,j  ,k,icomp+n)) * mskv(i,j,0);
                    Real dTdt = calc_arr(dom_lo.x,j,k,icomp+n) - dest_arr(dom_lo.x  ,j,k,icomp+n);
                    Real dTdx = dest_arr(dom_lo.x,j,k,icomp+n) - dest_arr(dom_lo.x+1,j,k,icomp+n);
                    if (dTdt*dTdx < 0.0_rt) dTdt = 0.0_rt;
                    Real dTde = (dTdt * (grad_lo+grad_lo_jp1) > 0.0_rt) ? grad_lo : grad_lo_jp1;
                    Real cff = std::max(dTdx*dTdx+dTde*dTde,eps);
                    Real Cx = dTdt * dTdx;
                    dest_arr(i,j,k,icomp+n) = (cff * calc_arr(dom_lo.x-1,j,k,icomp+n) + Cx * dest_arr(dom_lo.x,j,k,icomp+n)) * mskr(i,j,0) / (cff+Cx);
                }
            },
            bx_xhi, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                if (bc_ptr[n].hi(0) == REMORABCType::ext_dir) {
                    dest_arr(i,j,k,icomp+n) = l_bc_extdir_vals_d[n][3] * mskr(i,j,0);
                } else if (bc_ptr[n].hi(0) == REMORABCType::orlanski_rad) {
                    Real grad_hi      = (calc_arr(dom_hi.x  ,j  ,k,icomp+n) - calc_arr(dom_hi.x  ,j-1,k,icomp+n)) * mskv(i,j,0);
                    Real grad_hi_ip1  = (calc_arr(dom_hi.x+1,j  ,k,icomp+n) - calc_arr(dom_hi.x+1,j-1,k,icomp+n)) * mskv(i,j,0);
                    Real grad_hi_jp1  = (calc_arr(dom_hi.x  ,j+1,k,icomp+n) - calc_arr(dom_hi.x  ,j  ,k,icomp+n)) * mskv(i,j,0);
                    Real grad_hi_ijp1 = (calc_arr(dom_hi.x+1,j+1,k,icomp+n) - calc_arr(dom_hi.x+1,j  ,k,icomp+n)) * mskv(i,j,0);
                    Real dTdt = calc_arr(dom_hi.x,j,k,icomp+n) - dest_arr(dom_hi.x  ,j,k,icomp+n);
                    Real dTdx = dest_arr(dom_hi.x,j,k,icomp+n) - dest_arr(dom_hi.x-1,j,k,icomp+n);
                    if (dTdt * dTdx < 0.0_rt) dTdt = 0.0_rt;
                    Real dTde = (dTdt * (grad_hi + grad_hi_jp1) > 0.0_rt) ? grad_hi : grad_hi_jp1;
                    Real cff = std::max(dTdx*dTdx + dTde*dTde,eps);
                    Real Cx = dTdt * dTdx;
                    dest_arr(i,j,k,icomp+n) = (cff * calc_arr(dom_hi.x+1,j,k,icomp+n) + Cx * dest_arr(dom_hi.x,j,k,icomp+n)) * mskr(i,j,0) / (cff+Cx);
                }
            }
        );
    }

    if (!is_periodic_in_y)
    {
        Box bx_ylo(bx);  bx_ylo.setBig  (1,dom_lo.y-1);
        bx_ylo.setSmall(0,valid_bx.smallEnd(0)); bx_ylo.setBig(0,valid_bx.bigEnd(0));
        Box bx_yhi(bx);  bx_yhi.setSmall(1,dom_hi.y+1);
        bx_yhi.setSmall(0,valid_bx.smallEnd(0)); bx_yhi.setBig(0,valid_bx.bigEnd(0));
        ParallelFor(
            bx_ylo, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                if (bc_ptr[n].lo(1) == REMORABCType::ext_dir) {
                    dest_arr(i,j,k,icomp+n) = l_bc_extdir_vals_d[n][1] * mskr(i,j,0);
                } else if (bc_ptr[n].lo(1) == REMORABCType::orlanski_rad) {
                    Real grad_lo       = (calc_arr(i  ,dom_lo.y,  k,icomp+n) - calc_arr(i-1,dom_lo.y  ,k,icomp+n)) * msku(i,j,0);
                    Real grad_lo_jm1   = (calc_arr(i  ,dom_lo.y-1,k,icomp+n) - calc_arr(i-1,dom_lo.y-1,k,icomp+n)) * msku(i,j,0);
                    Real grad_lo_ip1   = (calc_arr(i+1,dom_lo.y  ,k,icomp+n) - calc_arr(i  ,dom_lo.y  ,k,icomp+n)) * msku(i,j,0);
                    Real grad_lo_ipjm1 = (calc_arr(i+1,dom_lo.y-1,k,icomp+n) - calc_arr(i  ,dom_lo.y-1,k,icomp+n)) * msku(i,j,0);
                    Real dTdt = calc_arr(i,dom_lo.y,k,icomp+n) - dest_arr(i,dom_lo.y  ,k,icomp+n);
                    Real dTde = dest_arr(i,dom_lo.y,k,icomp+n) - dest_arr(i,dom_lo.y+1,k,icomp+n);
                    if (dTdt * dTde < 0.0_rt) dTdt = 0.0_rt;
                    Real dTdx = (dTdt * (grad_lo + grad_lo_ip1) > 0.0_rt) ? grad_lo : grad_lo_ip1;
                    Real cff = std::max(dTdx*dTdx + dTde*dTde, eps);
                    Real Ce = dTdt*dTde;
                    dest_arr(i,j,k,icomp+n) = (cff * calc_arr(i,dom_lo.y-1,k,icomp+n) + Ce * dest_arr(i,dom_lo.y,k,icomp+n)) * mskr(i,j,0) / (cff+Ce);
                }
            },
            bx_yhi, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                if (bc_ptr[n].hi(1) == REMORABCType::ext_dir) {
                    dest_arr(i,j,k,icomp+n) = l_bc_extdir_vals_d[n][4] * mskr(i,j,0);
                } else if (bc_ptr[n].hi(1) == REMORABCType::orlanski_rad) {
                    Real grad_hi      = (calc_arr(i  ,dom_hi.y  ,k,icomp+n) - calc_arr(i-1,dom_hi.y  ,k,icomp+n)) * msku(i,j,0);
                    Real grad_hi_jp1  = (calc_arr(i  ,dom_hi.y+1,k,icomp+n) - calc_arr(i-1,dom_hi.y+1,k,icomp+n)) * msku(i,j,0);
                    Real grad_hi_ip1  = (calc_arr(i+1,dom_hi.y  ,k,icomp+n) - calc_arr(i  ,dom_hi.y  ,k,icomp+n)) * msku(i,j,0);
                    Real grad_hi_ijp1 = (calc_arr(i+1,dom_hi.y+1,k,icomp+n) - calc_arr(i  ,dom_hi.y+1,k,icomp+n)) * msku(i,j,0);
                    Real dTdt = calc_arr(i,dom_hi.y,k,icomp+n) - dest_arr(i,dom_hi.y  ,k,icomp+n);
                    Real dTde = dest_arr(i,dom_hi.y,k,icomp+n) - dest_arr(i,dom_hi.y-1,k,icomp+n);
                    if (dTdt * dTde < 0.0_rt) dTdt = 0.0_rt;
                    Real dTdx = (dTdt * (grad_hi + grad_hi_ip1) > 0.0_rt) ? grad_hi : grad_hi_ip1;
                    Real cff = std::max(dTdx*dTdx + dTde*dTde, eps);
                    Real Ce = dTdt*dTde;
                    dest_arr(i,j,k,icomp+n) = (cff*calc_arr(i,dom_hi.y+1,k,icomp+n) + Ce*dest_arr(i,dom_hi.y,k,icomp+n)) * mskr(i,j,0) / (cff+Ce);
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

    Box bx_xlo(bx);  bx_xlo.setBig  (0,dom_lo.x-1-n_not_fill);
                     bx_xlo.setSmall(2,std::max(dom_lo.z,bx.smallEnd(2)));
                     bx_xlo.setBig  (2,std::min(dom_hi.z,bx.bigEnd(2)));
    Box bx_xhi(bx);  bx_xhi.setSmall(0,dom_hi.x+1+n_not_fill);
                     bx_xhi.setSmall(2,std::max(dom_lo.z,bx.smallEnd(2)));
                     bx_xhi.setBig  (2,std::min(dom_hi.z,bx.bigEnd(2)));
    Box bx_ylo(bx);  bx_ylo.setBig  (1,dom_lo.y-1-n_not_fill);
                     bx_ylo.setSmall(2,std::max(dom_lo.z,bx.smallEnd(2)));
                     bx_ylo.setBig  (2,std::min(dom_hi.z,bx.bigEnd(2)));
    Box bx_yhi(bx);  bx_yhi.setSmall(1,dom_hi.y+1+n_not_fill);
                     bx_yhi.setSmall(2,std::max(dom_lo.z,bx.smallEnd(2)));
                     bx_yhi.setBig  (2,std::min(dom_hi.z,bx.bigEnd(2)));
    // Calculate intersections for corners before adjusting to exclude them
    Box xlo_ylo = bx_xlo & bx_ylo;
    Box xhi_ylo = bx_xhi & bx_ylo;
    Box xlo_yhi = bx_xlo & bx_yhi;
    Box xhi_yhi = bx_xhi & bx_yhi;
//    bx_xlo.setSmall(1,valid_bx.smallEnd(1));  bx_xlo.setBig(1,valid_bx.bigEnd(1));
//    bx_xhi.setSmall(1,valid_bx.smallEnd(1));  bx_xhi.setBig(1,valid_bx.bigEnd(1));
//    bx_ylo.setSmall(0,valid_bx.smallEnd(0)); bx_ylo.setBig(0,valid_bx.bigEnd(0));
//    bx_yhi.setSmall(0,valid_bx.smallEnd(0)); bx_yhi.setBig(0,valid_bx.bigEnd(0));
    // Next do ghost cells in x-direction but not reaching out in y
    // The corners we miss here will be covered in the y-loop below or by periodicity
    if (!is_periodic_in_x or bccomp==BCVars::foextrap_bc)
    {
        // Populate ghost cells on lo-x and hi-x domain boundaries
        ParallelFor(bx_xlo, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                int iflip = dom_lo.x - 1 - i;
                int inner = (bc_ptr[n].lo(0) == REMORABCType::orlanski_rad) ? 1 : 0;
                if (bc_ptr[n].lo(0) == REMORABCType::foextrap || bc_ptr[n].lo(0) == REMORABCType::clamped || bc_ptr[n].lo(0) == REMORABCType::chapman || bc_ptr[n].lo(0) == REMORABCType::orlanski_rad ||
                    bc_ptr[n].lo(0) == REMORABCType::orlanski_rad_nudge) {
                    dest_arr(i,j,k,icomp+n) =  dest_arr(dom_lo.x-n_not_fill-inner,j,k,icomp+n);
                } else if (bc_ptr[n].lo(0) == REMORABCType::reflect_even) {
                    dest_arr(i,j,k,icomp+n) =  dest_arr(iflip,j,k,icomp+n);
                } else if (bc_ptr[n].lo(0) == REMORABCType::reflect_odd) {
                    dest_arr(i,j,k,icomp+n) = -dest_arr(iflip,j,k,icomp+n);
                }
            },
            bx_xhi, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                int iflip =  2*dom_hi.x + 1 - i;
                int inner = (bc_ptr[n].hi(0) == REMORABCType::orlanski_rad) ? 1 : 0;
                if (bc_ptr[n].hi(0) == REMORABCType::foextrap || bc_ptr[n].hi(0) == REMORABCType::clamped || bc_ptr[n].hi(0) == REMORABCType::chapman || bc_ptr[n].hi(0) == REMORABCType::orlanski_rad ||
                    bc_ptr[n].hi(0) == REMORABCType::orlanski_rad_nudge) {
                    dest_arr(i,j,k,icomp+n) =  dest_arr(dom_hi.x+n_not_fill+inner,j,k,icomp+n);
                } else if (bc_ptr[n].hi(0) == REMORABCType::reflect_even) {
                    dest_arr(i,j,k,icomp+n) =  dest_arr(iflip,j,k,icomp+n);
                } else if (bc_ptr[n].hi(0) == REMORABCType::reflect_odd) {
                    dest_arr(i,j,k,icomp+n) = -dest_arr(iflip,j,k,icomp+n);
                }
            }
        );
    }

    if (!is_periodic_in_y or bccomp==BCVars::foextrap_bc)
    {
        // Populate ghost cells on lo-y and hi-y domain boundaries
        ParallelFor(
            bx_ylo, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                int jflip = dom_lo.y - 1 - j;
                int inner = (bc_ptr[n].lo(1) == REMORABCType::orlanski_rad) ? 1 : 0;
                if (bc_ptr[n].lo(1) == REMORABCType::foextrap || bc_ptr[n].lo(1) == REMORABCType::clamped || bc_ptr[n].lo(1) == REMORABCType::chapman || bc_ptr[n].lo(1) == REMORABCType::orlanski_rad ||
                    bc_ptr[n].lo(1) == REMORABCType::orlanski_rad_nudge) {
                    dest_arr(i,j,k,icomp+n) =  dest_arr(i,dom_lo.y-n_not_fill-inner,k,icomp+n);
                } else if (bc_ptr[n].lo(1) == REMORABCType::reflect_even) {
                    dest_arr(i,j,k,icomp+n) =  dest_arr(i,jflip,k,icomp+n);
                } else if (bc_ptr[n].lo(1) == REMORABCType::reflect_odd) {
                    dest_arr(i,j,k,icomp+n) = -dest_arr(i,jflip,k,icomp+n);
                }
            },
            bx_yhi, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                int jflip =  2*dom_hi.y + 1 - j;
                int inner = (bc_ptr[n].hi(1) == REMORABCType::orlanski_rad) ? 1 : 0;
                if (bc_ptr[n].hi(1) == REMORABCType::foextrap || bc_ptr[n].hi(1) == REMORABCType::clamped || bc_ptr[n].hi(1) == REMORABCType::chapman || bc_ptr[n].lo(1) == REMORABCType::orlanski_rad ||
                    bc_ptr[n].hi(1) == REMORABCType::orlanski_rad_nudge) {
                    dest_arr(i,j,k,icomp+n) =  dest_arr(i,dom_hi.y+n_not_fill+inner,k,icomp+n);
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
    if ((!is_periodic_in_x && !is_periodic_in_y) or bccomp==BCVars::foextrap_bc) {
        // If we've applied boundary conditions to either side, update the corner
        if (!xlo_ylo.isEmpty()) {
            ParallelFor(xlo_ylo, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                if (!(bc_ptr[n].lo(0) == REMORABCType::clamped || bc_ptr[n].lo(0) == REMORABCType::flather ||
                      bc_ptr[n].lo(0) == REMORABCType::chapman || bc_ptr[n].lo(0) == REMORABCType::orlanski_rad_nudge)
                 && !(bc_ptr[n].lo(1) == REMORABCType::clamped || bc_ptr[n].lo(1) == REMORABCType::flather ||
                      bc_ptr[n].lo(1) == REMORABCType::chapman || bc_ptr[n].lo(1) == REMORABCType::orlanski_rad_nudge)) {
                    dest_arr(i,j,k,icomp+n) = 0.5 * (dest_arr(i,dom_lo.y,k,icomp+n)
                                                    + dest_arr(dom_lo.x,j,k,icomp+n));
                }
            });
        }
        if (!xlo_yhi.isEmpty()) {
            ParallelFor(xlo_yhi, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                if (!(bc_ptr[n].lo(0) == REMORABCType::clamped || bc_ptr[n].lo(0) == REMORABCType::flather ||
                      bc_ptr[n].lo(0) == REMORABCType::chapman || bc_ptr[n].lo(0) == REMORABCType::orlanski_rad_nudge)
                 && !(bc_ptr[n].hi(1) == REMORABCType::clamped || bc_ptr[n].hi(1) == REMORABCType::flather ||
                      bc_ptr[n].hi(1) == REMORABCType::chapman || bc_ptr[n].hi(1) == REMORABCType::orlanski_rad_nudge)) {
                    dest_arr(i,j,k,icomp+n) = 0.5 * (dest_arr(i,dom_hi.y,k,icomp+n)
                                                    + dest_arr(dom_lo.x,j,k,icomp+n));
                }
            });
        }
        if (!xhi_ylo.isEmpty()) {
            ParallelFor(xhi_ylo, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                if (!(bc_ptr[n].hi(0) == REMORABCType::clamped || bc_ptr[n].hi(0) == REMORABCType::flather ||
                      bc_ptr[n].hi(0) == REMORABCType::chapman || bc_ptr[n].hi(0) == REMORABCType::orlanski_rad_nudge)
                 && !(bc_ptr[n].lo(1) == REMORABCType::clamped || bc_ptr[n].lo(1) == REMORABCType::flather ||
                      bc_ptr[n].lo(1) == REMORABCType::chapman || bc_ptr[n].lo(1) == REMORABCType::orlanski_rad_nudge)) {
                    dest_arr(i,j,k,icomp+n) = 0.5 * (dest_arr(i,dom_lo.y,k,icomp+n)
                                                    + dest_arr(dom_hi.x,j,k,icomp+n));
                }
            });
        }
        if (!xhi_yhi.isEmpty()) {
            ParallelFor(xhi_yhi, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                if (!(bc_ptr[n].hi(0) == REMORABCType::clamped || bc_ptr[n].hi(0) == REMORABCType::flather ||
                      bc_ptr[n].hi(0) == REMORABCType::chapman || bc_ptr[n].hi(0) == REMORABCType::orlanski_rad_nudge)
                 && !(bc_ptr[n].hi(1) == REMORABCType::clamped || bc_ptr[n].hi(1) == REMORABCType::flather ||
                      bc_ptr[n].hi(1) == REMORABCType::chapman || bc_ptr[n].hi(1) == REMORABCType::orlanski_rad_nudge)) {
                    dest_arr(i,j,k,icomp+n) = 0.5 * (dest_arr(i,dom_hi.y,k,icomp+n)
                                                    + dest_arr(dom_hi.x,j,k,icomp+n));
                }
            });
        }
    }

    Gpu::streamSynchronize();
}
