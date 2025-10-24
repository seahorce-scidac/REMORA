#include "REMORA.H"

using namespace amrex;

#ifdef REMORA_USE_NETCDF
/*
 * @param[in   ] lev           level to operate on
 * @param[inout] mf_to_fill    data on which to apply BCs
 * @param[in   ] mf_mask       land-sea mask
 * @param[in   ] time          current time
 * @param[in   ] bccomp        index into both domain_bcs_type_bcr and bc_extdir_vals for icomp=0
 * @param[in   ] bdy_var_type  which netcdf boundary data to fill from
 * @param[in   ] icomp_to_fill component to update
 * @param[in   ] icomp_calc    component to reference from on RHS
 * @param[in   ] mf_calc       data for RHS of calculation
 * @param[in   ] dt_calc       time step for the calculation
 */

void
REMORA::fill_from_bdyfiles (int lev, MultiFab& mf_to_fill, const MultiFab& mf_mask, const Real time, const int bccomp,
                            const int bdy_var_type, const int icomp_to_fill, const int icomp_calc, const MultiFab& mf_calc, const Real dt_calc)
{
    // Which variable are we filling
    int ivar = bdy_var_type;

    //
    // Note that "domain" is mapped onto the type of box the data is in
    //
    Box domain = geom[lev].Domain();

    const auto& mf_index_type = mf_to_fill.boxArray().ixType();
    domain.convert(mf_index_type);

    const auto& dom_lo = amrex::lbound(domain);
    const auto& dom_hi = amrex::ubound(domain);

    int ncomp;

    // If we are doing the scalars then do salt as well as temp
    if (ivar == BdyVars::t) {
        ncomp = 2;
    } else {
        ncomp = 1;
    }

    // This must be true for the logic below to work
    AMREX_ALWAYS_ASSERT(Temp_comp == 0);
    AMREX_ALWAYS_ASSERT(Salt_comp == 1);

    const Real eps= 1.0e-20_rt;
    const bool null_mf_calc = (!mf_calc.ok());

    for (int icomp = 0; icomp < ncomp; icomp++) // This is to do both temp and salt if doing scalars
    {
        // If we're doing zeta, ubar, or vbar, then calc_arr only has a single component
        // corresponding to the component to be used in calculating the boundary
        // value. Since we access icomp + icomp_to_fill_calc, we need icomp_to_fill_calc to be zero.
        // If it's another variable, either we aren't using calc_arr
        // or the components correspond to salt, temp, etc so we leave it as is.
        int icomp_to_fill_calc = (bccomp == BCVars::zeta_bc || bccomp == BCVars::ubar_bc ||
                              bccomp == BCVars::vbar_bc) ? 0 : icomp_to_fill;

        boundary_series[lev][ivar+icomp]->update_interpolated_to_time(time);

        const auto& bdatxlo = boundary_series[lev][ivar+icomp]->xlo_dat_interp.const_array();
        const auto& bdatxhi = boundary_series[lev][ivar+icomp]->xhi_dat_interp.const_array();
        const auto& bdatylo = boundary_series[lev][ivar+icomp]->ylo_dat_interp.const_array();
        const auto& bdatyhi = boundary_series[lev][ivar+icomp]->yhi_dat_interp.const_array();

        const auto& bx_bdatxlo = boundary_series[lev][ivar+icomp]->xlo_dat_interp.box();
        const auto& bx_bdatxhi = boundary_series[lev][ivar+icomp]->xhi_dat_interp.box();
        const auto& bx_bdatylo = boundary_series[lev][ivar+icomp]->ylo_dat_interp.box();
        const auto& bx_bdatyhi = boundary_series[lev][ivar+icomp]->yhi_dat_interp.box();

        if (domain_bcs_type[bccomp+icomp].lo(0) == REMORABCType::flather ||
            domain_bcs_type[bccomp+icomp].hi(0) == REMORABCType::flather ||
            domain_bcs_type[bccomp+icomp].lo(1) == REMORABCType::flather ||
            domain_bcs_type[bccomp+icomp].hi(1) == REMORABCType::flather) {
            boundary_series[lev][BdyVars::zeta]->update_interpolated_to_time(time);
        }
        const auto& bdatxlo_zeta = domain_bcs_type[bccomp+icomp].lo(0) == REMORABCType::flather ?
                                   boundary_series[lev][BdyVars::zeta]->xlo_dat_interp.const_array() : Array4<Real>();
        const auto& bdatxhi_zeta = domain_bcs_type[bccomp+icomp].hi(0) == REMORABCType::flather ?
                                   boundary_series[lev][BdyVars::zeta]->xhi_dat_interp.const_array() : Array4<Real>();
        const auto& bdatylo_zeta = domain_bcs_type[bccomp+icomp].lo(1) == REMORABCType::flather ?
                                   boundary_series[lev][BdyVars::zeta]->ylo_dat_interp.const_array() : Array4<Real>();
        const auto& bdatyhi_zeta = domain_bcs_type[bccomp+icomp].hi(1) == REMORABCType::flather ?
                                   boundary_series[lev][BdyVars::zeta]->yhi_dat_interp.const_array() : Array4<Real>();

        const bool apply_west  = (domain_bcs_type[bccomp+icomp].lo(0) == REMORABCType::clamped) ||
                                 (domain_bcs_type[bccomp+icomp].lo(0) == REMORABCType::flather) ||
                                 (domain_bcs_type[bccomp+icomp].lo(0) == REMORABCType::chapman) ||
                                 (domain_bcs_type[bccomp+icomp].lo(0) == REMORABCType::orlanski_rad_nudge);
        const bool apply_east  = (domain_bcs_type[bccomp+icomp].hi(0) == REMORABCType::clamped) ||
                                 (domain_bcs_type[bccomp+icomp].hi(0) == REMORABCType::flather) ||
                                 (domain_bcs_type[bccomp+icomp].hi(0) == REMORABCType::chapman) ||
                                 (domain_bcs_type[bccomp+icomp].hi(0) == REMORABCType::orlanski_rad_nudge);
        const bool apply_north = (domain_bcs_type[bccomp+icomp].lo(1) == REMORABCType::clamped) ||
                                 (domain_bcs_type[bccomp+icomp].lo(1) == REMORABCType::flather) ||
                                 (domain_bcs_type[bccomp+icomp].lo(1) == REMORABCType::chapman) ||
                                 (domain_bcs_type[bccomp+icomp].lo(1) == REMORABCType::orlanski_rad_nudge);
        const bool apply_south = (domain_bcs_type[bccomp+icomp].hi(1) == REMORABCType::clamped) ||
                                 (domain_bcs_type[bccomp+icomp].hi(1) == REMORABCType::flather) ||
                                 (domain_bcs_type[bccomp+icomp].hi(1) == REMORABCType::chapman) ||
                                 (domain_bcs_type[bccomp+icomp].hi(1) == REMORABCType::orlanski_rad_nudge);

        const bool cell_centered = (mf_index_type[0] == 0 and mf_index_type[1] == 0);

        const Real obcfac = solverChoice.obcfac;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        // Currently no tiling in order to get the logic right
        for (MFIter mfi(mf_to_fill,false); mfi.isValid(); ++mfi)
        {
            Box mf_box(mf_to_fill[mfi.index()].box());

            // Compute intersections of the FAB to be filled and the bdry data boxes
            Box xlo = bx_bdatxlo & mf_box;
            Box xhi = bx_bdatxhi & mf_box;
            Box ylo = bx_bdatylo & mf_box;
            Box yhi = bx_bdatyhi & mf_box;

            xlo.setSmall(0,lbound(mf_box).x);
            xhi.setBig  (0,ubound(mf_box).x);
            ylo.setSmall(1,lbound(mf_box).y);
            yhi.setBig  (1,ubound(mf_box).y);

            Box xlo_ylo = xlo & ylo;
            Box xlo_yhi = xlo & yhi;
            Box xhi_ylo = xhi & ylo;
            Box xhi_yhi = xhi & yhi;

            Box xlo_edge = xlo; xlo_edge.setSmall(0,ubound(xlo).x); xlo_edge.setBig(0,ubound(xlo).x);
            Box xhi_edge = xhi; xhi_edge.setSmall(0,lbound(xhi).x); xhi_edge.setBig(0,lbound(xhi).x);
            Box ylo_edge = ylo; ylo_edge.setSmall(1,ubound(ylo).y); ylo_edge.setBig(1,ubound(ylo).y);
            Box yhi_edge = yhi; yhi_edge.setSmall(1,lbound(yhi).y); yhi_edge.setBig(1,lbound(yhi).y);

            Box xlo_ghost = xlo; xlo_ghost.setBig(0,ubound(xlo).x-1);
            Box xhi_ghost = xhi; xhi_ghost.setSmall(0,lbound(xhi).x+1);
            Box ylo_ghost = ylo; ylo_ghost.setBig(1,ubound(ylo).y-1);
            Box yhi_ghost = yhi; yhi_ghost.setSmall(1,lbound(yhi).y+1);

            const Array4<Real>& dest_arr = mf_to_fill.array(mfi);
            const Array4<const Real>& mask_arr = mf_mask.array(mfi);
            const Array4<const Real>& calc_arr = (!null_mf_calc) ? mf_calc.array(mfi) : Array4<amrex::Real>();
            const Array4<const Real>& h_arr = vec_h[lev]->const_array(mfi);
            const Array4<const Real>& zeta_arr = vec_zeta[lev]->const_array(mfi);
            const Array4<const Real>& pm = vec_pm[lev]->const_array(mfi);
            const Array4<const Real>& pn = vec_pn[lev]->const_array(mfi);

            const Array4<const Real>& msku = vec_msku[lev]->const_array(mfi);
            const Array4<const Real>& mskv = vec_mskv[lev]->const_array(mfi);

            const Array4<const Real> nudg_coeff_out = vec_nudg_coeff[bdy_var_type][lev]->const_array(mfi);

            //
            // We are inside a loop over components so we do one at a time here
            //
            Vector<BCRec> bcrs(1);
            amrex::setBC(mf_box, domain, bccomp+icomp, 0, 1, domain_bcs_type, bcrs);

            // xlo: ori = 0
            // ylo: ori = 1
            // zlo: ori = 2
            // xhi: ori = 3
            // yhi: ori = 4
            // zhi: ori = 5

            auto bcr = bcrs[0];

            // Even though we don't loop over xlo itself, this is the right condition to check, since xlo_edge will always be the same for each grid,
            // but if the grid doesn't include the low x-boundary, the xlo box will be invalid and the execution will be skipped.
            if (!xlo.isEmpty() && apply_west) {
                ParallelFor(grow(xlo_edge,IntVect(0,-1,0)), [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    Real bry_val = bdatxlo(ubound(xlo).x,j,k,0);
                    if (bcr.lo(0) == REMORABCType::clamped) {
                        dest_arr(i,j,k,icomp+icomp_to_fill) = bry_val * mask_arr(i,j,0);
                    } else if (bcr.lo(0) == REMORABCType::flather) {
                        Real bry_val_zeta = bdatxlo_zeta(ubound(xlo).x-1,j,k,0);
                        Real cff = 1.0_rt / (0.5_rt * (h_arr(dom_lo.x-1,j,0) + zeta_arr(dom_lo.x-1,j,0,icomp_calc)
                                                     + h_arr(dom_lo.x,j,0) + zeta_arr(dom_lo.x,j,0,icomp_calc)));
                        Real Cx = std::sqrt(g * cff);
                        dest_arr(i,j,k,icomp+icomp_to_fill) = (bry_val
                                - Cx * (0.5_rt * (zeta_arr(dom_lo.x-1,j,0,icomp_calc) + zeta_arr(dom_lo.x,j,0,icomp_calc))
                                    - bry_val_zeta)) * mask_arr(i,j,0);
                    } else if (bcr.lo(0) == REMORABCType::chapman) {
                        Real cff = dt_calc * 0.5_rt * (pm(dom_lo.x,j-mf_index_type[1],0) + pm(dom_lo.x,j,0));
                        Real cff1 = std::sqrt(g * 0.5_rt * (h_arr(dom_lo.x,j-mf_index_type[1],0)
                                    + zeta_arr(dom_lo.x,j-mf_index_type[1],0,icomp_calc) + h_arr(dom_lo.x,j,0)
                                    + zeta_arr(dom_lo.x,j,0,icomp_calc)));
                        Real Cx = cff * cff1;
                        Real cff2 = 1.0_rt / (1.0_rt + Cx);
                        dest_arr(i,j,k,icomp+icomp_to_fill) = cff2 * (dest_arr(dom_lo.x-1,j,k,icomp_calc)
                                + Cx * dest_arr(dom_lo.x,j,k,icomp+icomp_to_fill)) * mask_arr(i,j,0);
                    } else if (bcr.lo(0) == REMORABCType::orlanski_rad_nudge) {
                        Real grad_lo_im1   = (calc_arr(dom_lo.x+mf_index_type[0]-1,j  ,k,icomp+icomp_to_fill_calc) - calc_arr(dom_lo.x-1+mf_index_type[0],j-1,k,icomp+icomp_to_fill_calc));
                        Real grad_lo       = (calc_arr(dom_lo.x+mf_index_type[0]  ,j  ,k,icomp+icomp_to_fill_calc) - calc_arr(dom_lo.x  +mf_index_type[0],j-1,k,icomp+icomp_to_fill_calc));
                        Real grad_lo_imjp1 = (calc_arr(dom_lo.x+mf_index_type[0]-1,j+1,k,icomp+icomp_to_fill_calc) - calc_arr(dom_lo.x-1+mf_index_type[0],j  ,k,icomp+icomp_to_fill_calc));
                        Real grad_lo_jp1   = (calc_arr(dom_lo.x+mf_index_type[0]  ,j+1,k,icomp+icomp_to_fill_calc) - calc_arr(dom_lo.x  +mf_index_type[0],j  ,k,icomp+icomp_to_fill_calc));
                        if (cell_centered) {
                                grad_lo_im1   *= mskv(i,j,0);
                                grad_lo       *= mskv(i,j,0);
                                grad_lo_imjp1 *= mskv(i,j,0);
                                grad_lo_jp1   *= mskv(i,j,0);
                        }
                        Real dTdt = calc_arr(dom_lo.x+mf_index_type[0],j,k,icomp+icomp_to_fill_calc) - dest_arr(dom_lo.x+mf_index_type[0]  ,j,k,icomp+icomp_to_fill);
                        Real dTdx = dest_arr(dom_lo.x+mf_index_type[0],j,k,icomp+icomp_to_fill) - dest_arr(dom_lo.x+mf_index_type[0]+1,j,k,icomp+icomp_to_fill);
                        Real tau;
                        Real nudg_coeff_out_local = (nudg_coeff_out(i-mf_index_type[0],j-mf_index_type[1],k) +
                                                     nudg_coeff_out(i,j,k)) * 0.5_rt;
                        if (dTdt*dTdx < 0.0_rt) {
                            tau = nudg_coeff_out_local * obcfac * dt_calc;
                            dTdt = 0.0_rt;
                        } else {
                            tau = nudg_coeff_out_local * dt_calc;
                        }
                        Real dTde = (dTdt * (grad_lo+grad_lo_jp1) > 0.0_rt) ? grad_lo : grad_lo_jp1;
                        Real cff = std::max(dTdx*dTdx+dTde*dTde,eps);
                        Real Cx = dTdt * dTdx;
                        dest_arr(i,j,k,icomp+icomp_to_fill) = (cff * calc_arr(dom_lo.x-1+mf_index_type[0],j,k,icomp+icomp_to_fill_calc) + Cx * dest_arr(dom_lo.x+mf_index_type[0],j,k,icomp+icomp_to_fill)) / (cff+Cx);
                        dest_arr(i,j,k,icomp+icomp_to_fill) = mask_arr(i,j,0) * (dest_arr(dom_lo.x-1+mf_index_type[0],j,k,icomp+icomp_to_fill) + tau * (bry_val - calc_arr(dom_lo.x-1+mf_index_type[0],j,k,icomp+icomp_to_fill_calc)));
                    }
                });
                ParallelFor(grow(xlo_ghost,IntVect(0,-1,0)), [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    dest_arr(i,j,k,icomp+icomp_to_fill) = dest_arr(ubound(xlo).x,j,k,icomp+icomp_to_fill);
                });
            }

            // See comment on xlo
            if (!xhi.isEmpty() && apply_east) {
                ParallelFor(grow(xhi_edge,IntVect(0,-1,0)), [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    Real bry_val = bdatxhi(lbound(xhi).x,j,k,0);
                    if (bcr.hi(0) == REMORABCType::clamped) {
                        dest_arr(i,j,k,icomp+icomp_to_fill) = bry_val * mask_arr(i,j,0);
                    } else if (bcr.hi(0) == REMORABCType::flather) {
                        Real bry_val_zeta = bdatxhi_zeta(lbound(xhi).x,j,k,0);
                        Real cff = 1.0_rt / (0.5_rt * (h_arr(dom_hi.x-1,j,0) + zeta_arr(dom_hi.x-1,j,0,icomp_calc)
                                                     + h_arr(dom_hi.x,j,0) + zeta_arr(dom_hi.x,j,0,icomp_calc)));
                        Real Cx = std::sqrt(g * cff);
                        dest_arr(i,j,k,icomp+icomp_to_fill) = (bry_val
                                + Cx * (0.5_rt * (zeta_arr(dom_hi.x-1,j,0,icomp_calc) + zeta_arr(dom_hi.x,j,0,icomp_calc))
                                    - bry_val_zeta)) * mask_arr(i,j,0);
                    } else if (bcr.hi(0) == REMORABCType::chapman) {
                        Real cff = dt_calc * 0.5_rt * (pm(dom_hi.x,j-mf_index_type[1],0) + pm(dom_hi.x,j,0));
                        Real cff1 = std::sqrt(g * 0.5_rt * (h_arr(dom_hi.x,j-mf_index_type[1],0)
                                    + zeta_arr(dom_hi.x,j-mf_index_type[1],0,icomp_calc) + h_arr(dom_hi.x,j,0)
                                    + zeta_arr(dom_hi.x,j,0,icomp_calc)));
                        Real Cx = cff * cff1;
                        Real cff2 = 1.0_rt / (1.0_rt + Cx);
                        dest_arr(i,j,k,icomp+icomp_to_fill) = cff2 * (dest_arr(dom_hi.x+1,j,k,icomp_calc)
                                + Cx * dest_arr(dom_hi.x,j,k,icomp+icomp_to_fill)) * mask_arr(i,j,0);
                    } else if (bcr.hi(0) == REMORABCType::orlanski_rad_nudge) {
                        Real grad_hi      = (calc_arr(dom_hi.x-mf_index_type[0]  ,j  ,k,icomp+icomp_to_fill_calc) - calc_arr(dom_hi.x-mf_index_type[0]  ,j-1,k,icomp+icomp_to_fill_calc));
                        Real grad_hi_ip1  = (calc_arr(dom_hi.x-mf_index_type[0]+1,j  ,k,icomp+icomp_to_fill_calc) - calc_arr(dom_hi.x-mf_index_type[0]+1,j-1,k,icomp+icomp_to_fill_calc));
                        Real grad_hi_jp1  = (calc_arr(dom_hi.x-mf_index_type[0]  ,j+1,k,icomp+icomp_to_fill_calc) - calc_arr(dom_hi.x-mf_index_type[0]  ,j  ,k,icomp+icomp_to_fill_calc));
                        Real grad_hi_ijp1 = (calc_arr(dom_hi.x-mf_index_type[0]+1,j+1,k,icomp+icomp_to_fill_calc) - calc_arr(dom_hi.x-mf_index_type[0]+1,j  ,k,icomp+icomp_to_fill_calc));
                        if (cell_centered) {
                            grad_hi      *= mskv(i,j,0);
                            grad_hi_ip1  *= mskv(i,j,0);
                            grad_hi_jp1  *= mskv(i,j,0);
                            grad_hi_ijp1 *= mskv(i,j,0);
                        }
                        Real dTdt = calc_arr(dom_hi.x-mf_index_type[0],j,k,icomp+icomp_to_fill_calc) - dest_arr(dom_hi.x-mf_index_type[0]  ,j,k,icomp+icomp_to_fill);
                        Real dTdx = dest_arr(dom_hi.x-mf_index_type[0],j,k,icomp+icomp_to_fill) - dest_arr(dom_hi.x-mf_index_type[0]-1,j,k,icomp+icomp_to_fill);
                        Real tau;
                        Real nudg_coeff_out_local = (nudg_coeff_out(i-mf_index_type[0],j-mf_index_type[1],k) +
                                                     nudg_coeff_out(i,j,k)) * 0.5_rt;
                        if (dTdt*dTdx < 0.0_rt) {
                            tau = nudg_coeff_out_local * obcfac * dt_calc;
                            dTdt = 0.0_rt;
                        } else {
                            tau = nudg_coeff_out_local * dt_calc;
                        }
                        if (dTdt * dTdx < 0.0_rt) dTdt = 0.0_rt;
                        Real dTde = (dTdt * (grad_hi + grad_hi_jp1) > 0.0_rt) ? grad_hi : grad_hi_jp1;
                        Real cff = std::max(dTdx*dTdx + dTde*dTde,eps);
                        Real Cx = dTdt * dTdx;
                        dest_arr(i,j,k,icomp+icomp_to_fill) = (cff * calc_arr(dom_hi.x+1-mf_index_type[0],j,k,icomp+icomp_to_fill_calc) + Cx * dest_arr(dom_hi.x-mf_index_type[0],j,k,icomp+icomp_to_fill)) * mask_arr(i,j,0) / (cff+Cx);
                        dest_arr(i,j,k,icomp+icomp_to_fill) = mask_arr(i,j,0) * (dest_arr(dom_hi.x+1-mf_index_type[0],j,k,icomp+icomp_to_fill) + tau * (bry_val - calc_arr(dom_hi.x+1-mf_index_type[0],j,k,icomp+icomp_to_fill_calc)));
                    }
                });
                ParallelFor(grow(xhi_ghost,IntVect(0,-1,0)), [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    dest_arr(i,j,k,icomp+icomp_to_fill) = dest_arr(lbound(xhi).x,j,k,icomp+icomp_to_fill);
                });
            }

            // See comment on xlo
            if (!ylo.isEmpty() && apply_south) {
                ParallelFor(grow(ylo_edge,IntVect(-1,0,0)), [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    Real bry_val = bdatylo(i,ubound(ylo).y,k,0);
                    if (bcr.lo(1) == REMORABCType::clamped) {
                        dest_arr(i,j,k,icomp+icomp_to_fill) = bry_val * mask_arr(i,j,0);
                    } else if (bcr.lo(1) == REMORABCType::flather) {
                        Real bry_val_zeta = bdatylo_zeta(i,ubound(ylo).y-1,k,0);
                        Real cff = 1.0_rt / (0.5_rt * (h_arr(i,dom_lo.y-1,0) + zeta_arr(i,dom_lo.y-1,0,icomp_calc)
                                                     + h_arr(i,dom_lo.y,0) + zeta_arr(i,dom_lo.y,0,icomp_calc)));
                        Real Ce = std::sqrt(g * cff);
                        dest_arr(i,j,k,icomp+icomp_to_fill) = (bry_val
                                - Ce * (0.5_rt * (zeta_arr(i,dom_lo.y-1,0,icomp_calc) + zeta_arr(i,dom_lo.y,0,icomp_calc))
                                    - bry_val_zeta)) * mask_arr(i,j,0);
                    } else if (bcr.lo(1) == REMORABCType::chapman) {
                        Real cff = dt_calc * 0.5_rt * (pn(i-mf_index_type[0],dom_lo.y,0) + pn(i,dom_lo.y,0));
                        Real cff1 = std::sqrt(g * 0.5_rt * (h_arr(i-mf_index_type[0],dom_lo.y,0) +
                                    zeta_arr(i-mf_index_type[0],dom_lo.y,0,icomp_calc) + h_arr(i,dom_lo.y,0)
                                    + zeta_arr(i,dom_lo.y,0,icomp_calc)));
                        Real Ce = cff * cff1;
                        Real cff2 = 1.0_rt / (1.0_rt + Ce);
                        dest_arr(i,j,k,icomp+icomp_to_fill) = cff2 * (dest_arr(i,dom_lo.y-1,k,icomp_calc)
                                + Ce * dest_arr(i,dom_lo.y,k,icomp+icomp_to_fill)) * mask_arr(i,j,0);
                    } else if (bcr.lo(1) == REMORABCType::orlanski_rad_nudge) {
                        Real grad_lo       = (calc_arr(i  ,dom_lo.y+mf_index_type[1],  k,icomp+icomp_to_fill_calc) - calc_arr(i-1,dom_lo.y+mf_index_type[1]  ,k,icomp+icomp_to_fill_calc));
                        Real grad_lo_jm1   = (calc_arr(i  ,dom_lo.y+mf_index_type[1]-1,k,icomp+icomp_to_fill_calc) - calc_arr(i-1,dom_lo.y+mf_index_type[1]-1,k,icomp+icomp_to_fill_calc));
                        Real grad_lo_ip1   = (calc_arr(i+1,dom_lo.y+mf_index_type[1]  ,k,icomp+icomp_to_fill_calc) - calc_arr(i  ,dom_lo.y+mf_index_type[1]  ,k,icomp+icomp_to_fill_calc));
                        Real grad_lo_ipjm1 = (calc_arr(i+1,dom_lo.y+mf_index_type[1]-1,k,icomp+icomp_to_fill_calc) - calc_arr(i  ,dom_lo.y+mf_index_type[1]-1,k,icomp+icomp_to_fill_calc));
                        if (cell_centered) {
                            grad_lo       *= msku(i,j,0);
                            grad_lo_jm1   *= msku(i,j,0);
                            grad_lo_ip1   *= msku(i,j,0);
                            grad_lo_ipjm1 *= msku(i,j,0);
                        }
                        Real dTdt = calc_arr(i,dom_lo.y+mf_index_type[1],k,icomp+icomp_to_fill_calc) - dest_arr(i,dom_lo.y  +mf_index_type[1],k,icomp+icomp_to_fill);
                        Real dTde = dest_arr(i,dom_lo.y+mf_index_type[1],k,icomp+icomp_to_fill) - dest_arr(i,dom_lo.y+1+mf_index_type[1],k,icomp+icomp_to_fill);
                        Real tau;
                        Real nudg_coeff_out_local = (nudg_coeff_out(i-mf_index_type[0],j-mf_index_type[1],k) +
                                                     nudg_coeff_out(i,j,k)) * 0.5_rt;
                        if (dTdt*dTde < 0.0_rt) {
                            tau = nudg_coeff_out_local * obcfac * dt_calc;
                            dTdt = 0.0_rt;
                        } else {
                            tau = nudg_coeff_out_local * dt_calc;
                        }
                        if (dTdt * dTde < 0.0_rt) dTdt = 0.0_rt;
                        Real dTdx = (dTdt * (grad_lo + grad_lo_ip1) > 0.0_rt) ? grad_lo : grad_lo_ip1;
                        Real cff = std::max(dTdx*dTdx + dTde*dTde, eps);
                        Real Ce = dTdt*dTde;
                        dest_arr(i,j,k,icomp+icomp_to_fill) = (cff * calc_arr(i,dom_lo.y-1+mf_index_type[1],k,icomp+icomp_to_fill_calc) + Ce * dest_arr(i,dom_lo.y+mf_index_type[1],k,icomp+icomp_to_fill)) / (cff+Ce);
                        dest_arr(i,j,k,icomp+icomp_to_fill) = mask_arr(i,j,0) * (dest_arr(i,dom_lo.y-1+mf_index_type[1],k,icomp+icomp_to_fill) + tau * (bry_val - calc_arr(i,dom_lo.y-1+mf_index_type[1],k,icomp+icomp_to_fill_calc)));
                    }
                });
                ParallelFor(grow(ylo_ghost,IntVect(-1,0,0)), [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    dest_arr(i,j,k,icomp+icomp_to_fill) = dest_arr(i,ubound(ylo).y,k,icomp+icomp_to_fill);
                });
            }

            // See comment on xlo
            if (!yhi.isEmpty() && apply_north) {
                ParallelFor(grow(yhi_edge,IntVect(-1,0,0)), [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    Real bry_val = bdatyhi(i,lbound(yhi).y,k,0);
                    if (bcr.hi(1) == REMORABCType::clamped) {
                        dest_arr(i,j,k,icomp+icomp_to_fill) = bry_val * mask_arr(i,j,0);
                    } else if (bcr.hi(1) == REMORABCType::flather) {
                        Real bry_val_zeta = bdatyhi_zeta(i,lbound(yhi).y,k,0);
                        Real cff = 1.0_rt / (0.5_rt * (h_arr(i,dom_hi.y-1,0) + zeta_arr(i,dom_hi.y-1,0,icomp_calc)
                                                     + h_arr(i,dom_hi.y,0) + zeta_arr(i,dom_hi.y,0,icomp_calc)));
                        Real Ce = std::sqrt(g * cff);
                        dest_arr(i,j,k,icomp+icomp_to_fill) = (bry_val
                                + Ce * (0.5_rt * (zeta_arr(i,dom_hi.y-1,0,icomp_calc) + zeta_arr(i,dom_hi.y,0,icomp_calc))
                                    - bry_val_zeta)) * mask_arr(i,j,0);
                    } else if (bcr.hi(1) == REMORABCType::chapman) {
                        Real cff = dt_calc * 0.5_rt * (pn(i-mf_index_type[0],dom_hi.y,0) + pn(i,dom_hi.y,0));
                        Real cff1 = std::sqrt(g * 0.5_rt * (h_arr(i-mf_index_type[0],dom_hi.y,0)
                                                          + zeta_arr(i-mf_index_type[0],dom_hi.y,0,icomp_calc) +
                                                            h_arr(i,dom_hi.y,0) + zeta_arr(i,dom_hi.y,0,icomp_calc)));
                        Real Ce = cff * cff1;
                        Real cff2 = 1.0_rt / (1.0_rt + Ce);
                        dest_arr(i,j,k,icomp+icomp_to_fill) = cff2 * (dest_arr(i,dom_hi.y+1,k,icomp_calc)
                                + Ce * dest_arr(i,dom_hi.y,k,icomp+icomp_to_fill)) * mask_arr(i,j,0);
                    } else if (bcr.hi(1) == REMORABCType::orlanski_rad_nudge) {
                        Real grad_hi      = calc_arr(i  ,dom_hi.y-mf_index_type[1]  ,k,icomp+icomp_to_fill_calc) - calc_arr(i-1,dom_hi.y-mf_index_type[1]  ,k,icomp+icomp_to_fill_calc);
                        Real grad_hi_jp1  = calc_arr(i  ,dom_hi.y-mf_index_type[1]+1,k,icomp+icomp_to_fill_calc) - calc_arr(i-1,dom_hi.y-mf_index_type[1]+1,k,icomp+icomp_to_fill_calc);
                        Real grad_hi_ip1  = calc_arr(i+1,dom_hi.y-mf_index_type[1]  ,k,icomp+icomp_to_fill_calc) - calc_arr(i  ,dom_hi.y-mf_index_type[1]  ,k,icomp+icomp_to_fill_calc);
                        Real grad_hi_ijp1 = calc_arr(i+1,dom_hi.y-mf_index_type[1]+1,k,icomp+icomp_to_fill_calc) - calc_arr(i  ,dom_hi.y-mf_index_type[1]+1,k,icomp+icomp_to_fill_calc);
                        if (cell_centered) {
                            grad_hi      *= msku(i,j,0);
                            grad_hi_jp1  *= msku(i,j,0);
                            grad_hi_ip1  *= msku(i,j,0);
                            grad_hi_ijp1 *= msku(i,j,0);
                        }
                        Real dTdt = calc_arr(i,dom_hi.y-mf_index_type[1],k,icomp+icomp_to_fill_calc) - dest_arr(i,dom_hi.y  -mf_index_type[1],k,icomp+icomp_to_fill);
                        Real dTde = dest_arr(i,dom_hi.y-mf_index_type[1],k,icomp+icomp_to_fill) - dest_arr(i,dom_hi.y-1-mf_index_type[1],k,icomp+icomp_to_fill);
                        Real tau;
                        Real nudg_coeff_out_local = (nudg_coeff_out(i-mf_index_type[0],j-mf_index_type[1],k) +
                                                     nudg_coeff_out(i,j,k)) * 0.5_rt;
                        if (dTdt*dTde < 0.0_rt) {
                            tau = nudg_coeff_out_local * obcfac * dt_calc;
                            dTdt = 0.0_rt;
                        } else {
                            tau = nudg_coeff_out_local * dt_calc;
                        }
                        if (dTdt * dTde < 0.0_rt) dTdt = 0.0_rt;
                        Real dTdx = (dTdt * (grad_hi + grad_hi_ip1) > 0.0_rt) ? grad_hi : grad_hi_ip1;
                        Real cff = std::max(dTdx*dTdx + dTde*dTde, eps);
                        Real Ce = dTdt*dTde;
                        dest_arr(i,j,k,icomp+icomp_to_fill) = (cff*calc_arr(i,dom_hi.y+1-mf_index_type[1],k,icomp+icomp_to_fill_calc) + Ce*dest_arr(i,dom_hi.y-mf_index_type[1],k,icomp+icomp_to_fill)) * mask_arr(i,j,0) / (cff+Ce);
                        dest_arr(i,j,k,icomp+icomp_to_fill) = mask_arr(i,j,0) * (dest_arr(i,dom_hi.y+1-mf_index_type[1],k,icomp+icomp_to_fill) + tau * (bry_val - calc_arr(i,dom_hi.y+1-mf_index_type[1],k,icomp+icomp_to_fill_calc)));
                    }
                });
                ParallelFor(grow(yhi_ghost,IntVect(-1,0,0)), [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    dest_arr(i,j,k,icomp+icomp_to_fill) = dest_arr(i,lbound(yhi).y,k,icomp+icomp_to_fill);
                });

            }
            // If we've applied boundary conditions to either side, update the corner
            if (!xlo_ylo.isEmpty() && (apply_west || apply_south)) {
                ParallelFor(xlo_ylo, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    dest_arr(i,j,k,icomp+icomp_to_fill) = 0.5 * (dest_arr(i,dom_lo.y+mf_index_type[1],k,icomp+icomp_to_fill)
                                                               + dest_arr(dom_lo.x+mf_index_type[0],j,k,icomp+icomp_to_fill));
                });
            }
            if (!xlo_yhi.isEmpty() && (apply_west || apply_north)) {
                ParallelFor(xlo_yhi, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    dest_arr(i,j,k,icomp+icomp_to_fill) = 0.5 * (dest_arr(i,dom_hi.y-mf_index_type[1],k,icomp+icomp_to_fill)
                                                               + dest_arr(dom_lo.x+mf_index_type[0],j,k,icomp+icomp_to_fill));
                });
            }
            if (!xhi_ylo.isEmpty() && (apply_east || apply_south)) {
                ParallelFor(xhi_ylo, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    dest_arr(i,j,k,icomp+icomp_to_fill) = 0.5 * (dest_arr(i,dom_lo.y+mf_index_type[1],k,icomp+icomp_to_fill)
                                                               + dest_arr(dom_hi.x-mf_index_type[0],j,k,icomp+icomp_to_fill));
                });
            }
            if (!xhi_yhi.isEmpty() && (apply_east || apply_north)) {
                ParallelFor(xhi_yhi, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    dest_arr(i,j,k,icomp+icomp_to_fill) = 0.5 * (dest_arr(i,dom_hi.y-mf_index_type[1],k,icomp+icomp_to_fill)
                                                               + dest_arr(dom_hi.x-mf_index_type[0],j,k,icomp+icomp_to_fill));
                });
            }
        } // mfi
    } // icomp
}
#endif
