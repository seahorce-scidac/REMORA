#include "AMReX_PhysBCFunct.H"
#include <REMORA_PhysBCFunct.H>

using namespace amrex;

/**
 * @param[inout] dest_arr      data on which to apply BCs
 * @param[in   ] bx            box to update on
 * @param[in   ] valid_bx      valid box
 * @param[in   ] domain        domain box
 * @param[in   ] dxInv         pm or pn
 * @param[in   ] mskr          land-sea mask on rho-points
 * @param[in   ] msku          land-sea mask on u-points
 * @param[in   ] mskv          land-sea mask on v-points
 * @param[in   ] calc_arr      data to use in the RHS of calculations
 * @param[in   ] icomp         component to update
 * @param[in   ] ncomp         number of components to update, starting from icomp
 * @param[in   ] time          current time
 * @param[in   ] bccomp        index into both domain_bcs_type_bcr and bc_extdir_vals for icomp=0
 * @param[in   ] n_not_fill    perimter of cells in x and y where BCs are not applied for non-ext_dir conditions
 */
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
    const auto* bc_extdir_vals_ptr = m_bc_extdir_vals_d.data();

    GeometryData const& geomdata = m_geom.data();
    bool is_periodic_in_x = geomdata.isPeriodic(0);
    bool is_periodic_in_y = geomdata.isPeriodic(1);
    const Real eps= Real(1.0e-20);

    // If we're doing zeta, then calc_arr only has a single component
    // corresponding to the component to be used in calculating the boundary
    // value. If it's another variable, either we aren't using calc_arr
    // or the components correspond to salt, temp, etc and we loop over ncomp
    // so we leave icomp as is
    int icomp_calc = (bccomp == BCVars::zeta_bc(m_ncons)) ? 0 : icomp;

    Box dest_arr_box = Box(dest_arr);
    // First do all ext_dir bcs
    if (!is_periodic_in_x)
    {
        Box bx_xlo(bx);  bx_xlo.setBig  (0,dom_lo.x-1);
        bx_xlo.setSmall(1,std::max(valid_bx.smallEnd(1)-1,dom_lo.y));  bx_xlo.setBig(1,std::min(valid_bx.bigEnd(1)+1,dom_hi.y));
        Box bx_xhi(bx);  bx_xhi.setSmall(0,dom_hi.x+1);
        bx_xhi.setSmall(1,std::max(valid_bx.smallEnd(1)-1,dom_lo.y));  bx_xhi.setBig(1,std::min(valid_bx.bigEnd(1)+1,dom_hi.y));
        ParallelFor(
            bx_xlo & dest_arr_box, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                if (bc_ptr[n].lo(0) == REMORABCType::ext_dir) {
                    dest_arr(i,j,k,icomp+n) = bc_extdir_vals_ptr[bccomp+n][0] * mskr(i,j,0);
                } else if (bc_ptr[n].lo(0) == REMORABCType::orlanski_rad) {
                    Real grad_lo       = (calc_arr(dom_lo.x  ,j  ,k,icomp_calc+n) - calc_arr(dom_lo.x  ,j-1,k,icomp_calc+n)) * mskv(i,j,0);
                    Real grad_lo_jp1   = (calc_arr(dom_lo.x  ,j+1,k,icomp_calc+n) - calc_arr(dom_lo.x  ,j  ,k,icomp_calc+n)) * mskv(i,j,0);
                    Real dTdt = calc_arr(dom_lo.x,j,k,icomp_calc+n) - dest_arr(dom_lo.x  ,j,k,icomp+n);
                    Real dTdx = dest_arr(dom_lo.x,j,k,icomp+n) - dest_arr(dom_lo.x+1,j,k,icomp+n);
                    if (dTdt*dTdx < zero) dTdt = zero;
                    Real dTde = (dTdt * (grad_lo+grad_lo_jp1) > zero) ? grad_lo : grad_lo_jp1;
                    Real cff = std::max(dTdx*dTdx+dTde*dTde,eps);
                    Real Cx = dTdt * dTdx;
                    dest_arr(i,j,k,icomp+n) = (cff * calc_arr(dom_lo.x-1,j,k,icomp_calc+n) + Cx * dest_arr(dom_lo.x,j,k,icomp+n)) * mskr(i,j,0) / (cff+Cx);
                }
            },
            bx_xhi & dest_arr_box, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                if (bc_ptr[n].hi(0) == REMORABCType::ext_dir) {
                    dest_arr(i,j,k,icomp+n) = bc_extdir_vals_ptr[bccomp+n][3] * mskr(i,j,0);
                } else if (bc_ptr[n].hi(0) == REMORABCType::orlanski_rad) {
                    Real grad_hi      = (calc_arr(dom_hi.x  ,j  ,k,icomp_calc+n) - calc_arr(dom_hi.x  ,j-1,k,icomp_calc+n)) * mskv(i,j,0);
                    Real grad_hi_jp1  = (calc_arr(dom_hi.x  ,j+1,k,icomp_calc+n) - calc_arr(dom_hi.x  ,j  ,k,icomp_calc+n)) * mskv(i,j,0);
                    Real dTdt = calc_arr(dom_hi.x,j,k,icomp_calc+n) - dest_arr(dom_hi.x  ,j,k,icomp+n);
                    Real dTdx = dest_arr(dom_hi.x,j,k,icomp+n) - dest_arr(dom_hi.x-1,j,k,icomp+n);
                    if (dTdt * dTdx < zero) dTdt = zero;
                    Real dTde = (dTdt * (grad_hi + grad_hi_jp1) > zero) ? grad_hi : grad_hi_jp1;
                    Real cff = std::max(dTdx*dTdx + dTde*dTde,eps);
                    Real Cx = dTdt * dTdx;
                    dest_arr(i,j,k,icomp+n) = (cff * calc_arr(dom_hi.x+1,j,k,icomp_calc+n) + Cx * dest_arr(dom_hi.x,j,k,icomp+n)) * mskr(i,j,0) / (cff+Cx);
                }
            }
        );
    }

    if (!is_periodic_in_y)
    {
        Box bx_ylo(bx);  bx_ylo.setBig  (1,dom_lo.y-1);
        bx_ylo.setSmall(0,std::max(valid_bx.smallEnd(0)-1,dom_lo.x)); bx_ylo.setBig(0,std::min(valid_bx.bigEnd(0)+1,dom_hi.x));
        Box bx_yhi(bx);  bx_yhi.setSmall(1,dom_hi.y+1);
        bx_yhi.setSmall(0,std::max(valid_bx.smallEnd(0)-1,dom_lo.x)); bx_yhi.setBig(0,std::min(valid_bx.bigEnd(0)+1,dom_hi.x));
        ParallelFor(
            bx_ylo & dest_arr_box, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                if (bc_ptr[n].lo(1) == REMORABCType::ext_dir) {
                    dest_arr(i,j,k,icomp+n) = bc_extdir_vals_ptr[bccomp+n][1] * mskr(i,j,0);
                } else if (bc_ptr[n].lo(1) == REMORABCType::orlanski_rad) {
                    Real grad_lo       = (calc_arr(i  ,dom_lo.y,  k,icomp_calc+n) - calc_arr(i-1,dom_lo.y  ,k,icomp_calc+n)) * msku(i,j,0);
                    Real grad_lo_ip1   = (calc_arr(i+1,dom_lo.y  ,k,icomp_calc+n) - calc_arr(i  ,dom_lo.y  ,k,icomp_calc+n)) * msku(i,j,0);
                    Real dTdt = calc_arr(i,dom_lo.y,k,icomp_calc+n) - dest_arr(i,dom_lo.y  ,k,icomp+n);
                    Real dTde = dest_arr(i,dom_lo.y,k,icomp+n) - dest_arr(i,dom_lo.y+1,k,icomp+n);
                    if (dTdt * dTde < zero) dTdt = zero;
                    Real dTdx = (dTdt * (grad_lo + grad_lo_ip1) > zero) ? grad_lo : grad_lo_ip1;
                    Real cff = std::max(dTdx*dTdx + dTde*dTde, eps);
                    Real Ce = dTdt*dTde;
                    dest_arr(i,j,k,icomp+n) = (cff * calc_arr(i,dom_lo.y-1,k,icomp_calc+n) + Ce * dest_arr(i,dom_lo.y,k,icomp+n)) * mskr(i,j,0) / (cff+Ce);
                }
            },
            bx_yhi & dest_arr_box, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                if (bc_ptr[n].hi(1) == REMORABCType::ext_dir) {
                    dest_arr(i,j,k,icomp+n) = bc_extdir_vals_ptr[bccomp+n][4] * mskr(i,j,0);
                } else if (bc_ptr[n].hi(1) == REMORABCType::orlanski_rad) {
                    Real grad_hi      = (calc_arr(i  ,dom_hi.y  ,k,icomp_calc+n) - calc_arr(i-1,dom_hi.y  ,k,icomp_calc+n)) * msku(i,j,0);
                    Real grad_hi_ip1  = (calc_arr(i+1,dom_hi.y  ,k,icomp_calc+n) - calc_arr(i  ,dom_hi.y  ,k,icomp_calc+n)) * msku(i,j,0);
                    Real dTdt = calc_arr(i,dom_hi.y,k,icomp_calc+n) - dest_arr(i,dom_hi.y  ,k,icomp+n);
                    Real dTde = dest_arr(i,dom_hi.y,k,icomp+n) - dest_arr(i,dom_hi.y-1,k,icomp+n);
                    if (dTdt * dTde < zero) dTdt = zero;
                    Real dTdx = (dTdt * (grad_hi + grad_hi_ip1) > zero) ? grad_hi : grad_hi_ip1;
                    Real cff = std::max(dTdx*dTdx + dTde*dTde, eps);
                    Real Ce = dTdt*dTde;
                    dest_arr(i,j,k,icomp+n) = (cff*calc_arr(i,dom_hi.y+1,k,icomp_calc+n) + Ce*dest_arr(i,dom_hi.y,k,icomp+n)) * mskr(i,j,0) / (cff+Ce);
                }
            }
        );
    }

    {
        Box bx_zlo(bx);  bx_zlo.setBig  (2,dom_lo.z-1);
        Box bx_zhi(bx);  bx_zhi.setSmall(2,dom_hi.z+1);
        ParallelFor(
            bx_zlo & dest_arr_box, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                if (bc_ptr[n].lo(2) == REMORABCType::ext_dir) {
                    dest_arr(i,j,k,icomp+n) = bc_extdir_vals_ptr[bccomp+n][2] * mskr(i,j,0);
                }
            },
            bx_zhi & dest_arr_box, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                if (bc_ptr[n].hi(2) == REMORABCType::ext_dir) {
                    dest_arr(i,j,k,icomp+n) = bc_extdir_vals_ptr[bccomp+n][5] * mskr(i,j,0);
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
    if (!is_periodic_in_x or bccomp == BCVars::foextrap_bc(m_ncons))
    {
        // Populate ghost cells on lo-x and hi-x domain boundaries
        ParallelFor(bx_xlo & dest_arr_box, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
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
            bx_xhi & dest_arr_box, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
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

    if (!is_periodic_in_y or bccomp == BCVars::foextrap_bc(m_ncons))
    {
        // Populate ghost cells on lo-y and hi-y domain boundaries
        ParallelFor(
            bx_ylo & dest_arr_box, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
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
            bx_yhi & dest_arr_box, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                int jflip =  2*dom_hi.y + 1 - j;
                int inner = (bc_ptr[n].hi(1) == REMORABCType::orlanski_rad) ? 1 : 0;
                if (bc_ptr[n].hi(1) == REMORABCType::foextrap || bc_ptr[n].hi(1) == REMORABCType::clamped || bc_ptr[n].hi(1) == REMORABCType::chapman || bc_ptr[n].hi(1) == REMORABCType::orlanski_rad ||
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
            ParallelFor(bx_zlo & dest_arr_box, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                int kflip = dom_lo.z - 1 - k;
                if (bc_ptr[n].lo(2) == REMORABCType::foextrap) {
                    dest_arr(i,j,k,icomp+n) =  dest_arr(i,j,dom_lo.z,icomp+n);
                } else if (bc_ptr[n].lo(2) == REMORABCType::reflect_even) {
                    dest_arr(i,j,k,icomp+n) =  dest_arr(i,j,kflip,icomp+n);
                } else if (bc_ptr[n].lo(2) == REMORABCType::reflect_odd) {
                    dest_arr(i,j,k,icomp+n) = -dest_arr(i,j,kflip,icomp+n);
                }
            });
        }

        if (bx_zhi.ok()) {
            ParallelFor(bx_zhi & dest_arr_box, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                int kflip =  2*dom_hi.z + 1 - k;
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
    if ((!is_periodic_in_x && !is_periodic_in_y) or bccomp == BCVars::foextrap_bc(m_ncons)) {
        // If we've applied boundary conditions to either side, update the corner
        if (!xlo_ylo.isEmpty()) {
            ParallelFor(xlo_ylo & dest_arr_box, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                if (!(bc_ptr[n].lo(0) == REMORABCType::clamped || bc_ptr[n].lo(0) == REMORABCType::flather ||
                      bc_ptr[n].lo(0) == REMORABCType::chapman || bc_ptr[n].lo(0) == REMORABCType::orlanski_rad_nudge)
                 && !(bc_ptr[n].lo(1) == REMORABCType::clamped || bc_ptr[n].lo(1) == REMORABCType::flather ||
                      bc_ptr[n].lo(1) == REMORABCType::chapman || bc_ptr[n].lo(1) == REMORABCType::orlanski_rad_nudge)) {
                    dest_arr(i,j,k,icomp+n) = half * (dest_arr(i,dom_lo.y,k,icomp+n)
                                                    + dest_arr(dom_lo.x,j,k,icomp+n));
                }
            });
        }
        if (!xlo_yhi.isEmpty()) {
            ParallelFor(xlo_yhi & dest_arr_box, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                if (!(bc_ptr[n].lo(0) == REMORABCType::clamped || bc_ptr[n].lo(0) == REMORABCType::flather ||
                      bc_ptr[n].lo(0) == REMORABCType::chapman || bc_ptr[n].lo(0) == REMORABCType::orlanski_rad_nudge)
                 && !(bc_ptr[n].hi(1) == REMORABCType::clamped || bc_ptr[n].hi(1) == REMORABCType::flather ||
                      bc_ptr[n].hi(1) == REMORABCType::chapman || bc_ptr[n].hi(1) == REMORABCType::orlanski_rad_nudge)) {
                    dest_arr(i,j,k,icomp+n) = half * (dest_arr(i,dom_hi.y,k,icomp+n)
                                                    + dest_arr(dom_lo.x,j,k,icomp+n));
                }
            });
        }
        if (!xhi_ylo.isEmpty()) {
            ParallelFor(xhi_ylo & dest_arr_box, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                if (!(bc_ptr[n].hi(0) == REMORABCType::clamped || bc_ptr[n].hi(0) == REMORABCType::flather ||
                      bc_ptr[n].hi(0) == REMORABCType::chapman || bc_ptr[n].hi(0) == REMORABCType::orlanski_rad_nudge)
                 && !(bc_ptr[n].lo(1) == REMORABCType::clamped || bc_ptr[n].lo(1) == REMORABCType::flather ||
                      bc_ptr[n].lo(1) == REMORABCType::chapman || bc_ptr[n].lo(1) == REMORABCType::orlanski_rad_nudge)) {
                    dest_arr(i,j,k,icomp+n) = half * (dest_arr(i,dom_lo.y,k,icomp+n)
                                                     + dest_arr(dom_hi.x,j,k,icomp+n));
                }
            });
        }
        if (!xhi_yhi.isEmpty()) {
            ParallelFor(xhi_yhi & dest_arr_box, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                if (!(bc_ptr[n].hi(0) == REMORABCType::clamped || bc_ptr[n].hi(0) == REMORABCType::flather ||
                      bc_ptr[n].hi(0) == REMORABCType::chapman || bc_ptr[n].hi(0) == REMORABCType::orlanski_rad_nudge)
                 && !(bc_ptr[n].hi(1) == REMORABCType::clamped || bc_ptr[n].hi(1) == REMORABCType::flather ||
                      bc_ptr[n].hi(1) == REMORABCType::chapman || bc_ptr[n].hi(1) == REMORABCType::orlanski_rad_nudge)) {
                    dest_arr(i,j,k,icomp+n) = half * (dest_arr(i,dom_hi.y,k,icomp+n)
                                                    + dest_arr(dom_hi.x,j,k,icomp+n));
                }
            });
        }
    }

    Gpu::streamSynchronize();
}
