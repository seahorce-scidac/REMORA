#include <REMORA.H>
#include "AMReX_Interp_3D_C.H"
#include "AMReX_PlotFileUtil.H"

using namespace amrex;

static PhysBCFunctNoOp null_bc_for_fill;

template<typename V, typename T>
bool containerHasElement(const V& iterable, const T& query) {
    return std::find(iterable.begin(), iterable.end(), query) != iterable.end();
}

// The nodal displacement a viewer adds to the Cartesian node position it builds from the
// Header: nu = (0, 0, z_phys_nd - (prob_lo_z + k*dz)), with dz and prob_lo from the Geometry the
// Header describes this level with. rz > 1 is the expand_plotvars_to_unif_rr case, where the
// Header's dz is the native one over rz and mf_nd carries rz node layers per native layer: z is
// interpolated linearly in k between native node layers, so node rz*k coincides with native
// node k. Fills component 2 only; the caller zeroes the rest.
static void
fill_nodal_z_displacement (MultiFab& mf_nd, const MultiFab& z_phys_nd, const Geometry& g, int rz)
{
    const Real dz  = g.CellSizeArray()[2];
    const Real zlo = g.ProbLoArray()[2];
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(mf_nd, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Array4<Real>       const& nu = mf_nd.array(mfi);
        Array4<Real const> const& zp = z_phys_nd.const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int kf) noexcept
        {
            const int k = kf / rz;
            const int r = kf - k * rz;
            Real z = zp(i,j,k);
            // z_phys_nd has one z ghost but only k = -1 is filled, so k+1 is read only when
            // its weight is nonzero, which keeps it at or below the top node.
            if (r > 0) { z += (Real(r) / Real(rz)) * (zp(i,j,k+1) - zp(i,j,k)); }
            nu(i,j,kf,2) = z - (zlo + Real(kf) * dz);
        });
    }
}

// Refine a face- or cell-centred MultiFab in z by rz onto dst, for the same expand path:
// piecewise constant where the data is cell-centred in z, linear between layers where it is
// nodal in z (the w faces).
static void
refine_in_z (const MultiFab& src, MultiFab& dst, int rz)
{
    const bool nodal_z = src.boxArray().ixType().nodeCentered(2);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(dst, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Array4<Real>       const& d = dst.array(mfi);
        Array4<Real const> const& s = src.const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int kf) noexcept
        {
            const int k = kf / rz;
            if (nodal_z) {
                const int r = kf - k * rz;
                Real v = s(i,j,k);
                if (r > 0) { v += (Real(r) / Real(rz)) * (s(i,j,k+1) - s(i,j,k)); }
                d(i,j,kf) = v;
            } else {
                d(i,j,kf) = s(i,j,k);
            }
        });
    }
}

// Write plotfile to disk
void
REMORA::WritePlotFile (int istep_for_plot)
{
#ifndef REMORA_USE_NETCDF
    amrex::ignore_unused(istep_for_plot);
#endif
    Vector<std::string> varnames_3d;
    varnames_3d.insert(varnames_3d.end(), plot_var_names_3d.begin(), plot_var_names_3d.end());

    Vector<std::string> varnames_2d;
    varnames_2d.insert(varnames_2d.end(), plot_var_names_2d.begin(), plot_var_names_2d.end());

    Vector<std::string> varnames_2d_rho;
    Vector<std::string> varnames_2d_u;
    Vector<std::string> varnames_2d_v;

    const int ncomp_mf_3d = varnames_3d.size();
    const auto ngrow_vars = IntVect(NGROW-1,NGROW-1,0);

    // These are the ncomp for the 2D cell-centered, x-face-based, y-face-based MultiFabs respectively
    int ncomp_mf_2d_rho = 0;
    int ncomp_mf_2d_u   = 0;
    int ncomp_mf_2d_v   = 0;

    // Check to see if we found all the requested variables
    for (auto plot_name : varnames_2d) {
      {
         if (plot_name == "zeta" ) {varnames_2d_rho.push_back(plot_name); ncomp_mf_2d_rho++;}
         if (plot_name == "h"    ) {varnames_2d_rho.push_back(plot_name); ncomp_mf_2d_rho++;}
         if (plot_name == "f"    ) {varnames_2d_rho.push_back(plot_name); ncomp_mf_2d_rho++;}
         if (plot_name == "visc2") {varnames_2d_rho.push_back(plot_name); ncomp_mf_2d_rho++;}
         for (int n = 0; n < ncons; ++n) {
             const std::string diff2_name = std::string("diff2_") + cons_names[n];
             if (plot_name == diff2_name) {
                 varnames_2d_rho.push_back(plot_name);
                 ncomp_mf_2d_rho++;
             }
         }
         for (int n = 0; n < ncons; ++n) {
             const std::string stflux_name = std::string("stflux_") + cons_names[n];
             if (plot_name == stflux_name) {
                 varnames_2d_rho.push_back(plot_name);
                 ncomp_mf_2d_rho++;
             }
         }
         if (plot_name == "lrflux") {varnames_2d_rho.push_back(plot_name); ncomp_mf_2d_rho++;}
         if (plot_name == "lhflux") {varnames_2d_rho.push_back(plot_name); ncomp_mf_2d_rho++;}
         if (plot_name == "srflux") {varnames_2d_rho.push_back(plot_name); ncomp_mf_2d_rho++;}
         if (plot_name == "shflux") {varnames_2d_rho.push_back(plot_name); ncomp_mf_2d_rho++;}
         if (plot_name == "mask_rho") {varnames_2d_rho.push_back(plot_name); ncomp_mf_2d_rho++;}
         if (plot_name == "ubar" ) {varnames_2d_u.push_back(plot_name); ncomp_mf_2d_u++;}
         if (plot_name == "sustr") {varnames_2d_u.push_back(plot_name); ncomp_mf_2d_u++;}
         if (plot_name == "bustr") {varnames_2d_u.push_back(plot_name); ncomp_mf_2d_u++;}
         if (plot_name == "mask_u") {varnames_2d_u.push_back(plot_name); ncomp_mf_2d_u++;}
         if (plot_name == "vbar" ) {varnames_2d_v.push_back(plot_name); ncomp_mf_2d_v++;}
         if (plot_name == "svstr") {varnames_2d_v.push_back(plot_name); ncomp_mf_2d_v++;}
         if (plot_name == "bvstr") {varnames_2d_v.push_back(plot_name); ncomp_mf_2d_v++;}
         if (plot_name == "mask_v") {varnames_2d_v.push_back(plot_name); ncomp_mf_2d_v++;}
      }
    }

    // We fillpatch here because some of the derived quantities require derivatives
    //     which require ghost cells to be filled. Don't fill the boundary, though.
    for (int lev = 0; lev <= finest_level; ++lev) {
        FillPatchNoBC(lev, t_new[lev], *cons_new[lev], cons_new, BdyVars::t,0,true,false);
        FillPatchNoBC(lev, t_new[lev], *xvel_new[lev], xvel_new, BdyVars::u,0,true,false);
        FillPatchNoBC(lev, t_new[lev], *yvel_new[lev], yvel_new, BdyVars::v,0,true,false);
        FillPatchNoBC(lev, t_new[lev], *zvel_new[lev], zvel_new, BdyVars::null,0,true,false);
        FillPatchNoBC(lev, t_new[lev], *vec_visc2_r[lev], GetVecOfPtrs(vec_visc2_r), BdyVars::null,0,true,false);
        FillPatchNoBC(lev, t_new[lev], *vec_diff2[lev],   GetVecOfPtrs(vec_diff2),   BdyVars::null,0,true,false);
    }

    if (plotfile_type == PlotfileType::amrex) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            mask_arrays_for_write(lev, plotfile_fill_value, zero);
        }
    } else if (plotfile_type == PlotfileType::netcdf) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            mask_arrays_for_write(lev, (Real) netcdf_fill_value, zero);
        }
    } else {
        amrex::Abort("Don't know this plotfile type");
    }

    // Array of 3D MultiFabs to hold the plotfile data
    Vector<MultiFab> plotMF(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        plotMF[lev].define(grids[lev], dmap[lev], ncomp_mf_3d, ngrow_vars);
        plotMF[lev].setVal(1.234e20);
    }

    // Array of 2D MultiFabs to hold the plotfile data
    Vector<MultiFab> mf_2d_rho(finest_level+1);
    Vector<MultiFab> mf_2d_u(finest_level+1);
    Vector<MultiFab> mf_2d_v(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        BoxArray ba(grids[lev]);
        BoxList bl2d = ba.boxList();
        for (auto& b : bl2d) {
            b.setRange(2,0);
        }
        BoxArray ba2d(std::move(bl2d));
        mf_2d_rho[lev].define(ba2d, dmap[lev], ncomp_mf_2d_rho, IntVect(0,0,0));
          mf_2d_u[lev].define(ba2d, dmap[lev], ncomp_mf_2d_u  , IntVect(0,0,0));
          mf_2d_v[lev].define(ba2d, dmap[lev], ncomp_mf_2d_v  , IntVect(0,0,0));
    }


    // Array of MultiFabs for nodal data
    Vector<MultiFab> mf_nd(finest_level+1);
    if (plot_nodal_data) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            BoxArray nodal_grids(grids[lev]); nodal_grids.surroundingNodes();
            mf_nd[lev].define(nodal_grids, dmap[lev], AMREX_SPACEDIM, 0);
            mf_nd[lev].setVal(zero);
        }
    }

    // Vector of MultiFabs for face-centered velocity
    Vector<MultiFab> mf_u(finest_level+1);
    Vector<MultiFab> mf_v(finest_level+1);
    Vector<MultiFab> mf_w(finest_level+1);
    if (plot_staggered_vels) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            BoxArray grid_stag_u(grids[lev]); grid_stag_u.surroundingNodes(0);
            BoxArray grid_stag_v(grids[lev]); grid_stag_v.surroundingNodes(1);
            BoxArray grid_stag_w(grids[lev]); grid_stag_w.surroundingNodes(2);
            mf_u[lev].define(grid_stag_u, dmap[lev], 1, 0);
            mf_v[lev].define(grid_stag_v, dmap[lev], 1, 0);
            mf_w[lev].define(grid_stag_w, dmap[lev], 1, 0);
            MultiFab::Copy(mf_u[lev],*xvel_new[lev],0,0,1,0);
            MultiFab::Copy(mf_v[lev],*yvel_new[lev],0,0,1,0);
            MultiFab::Copy(mf_w[lev],*zvel_new[lev],0,0,1,0);
        }
    }

    // Array of MultiFabs for cell-centered velocity
    Vector<MultiFab> mf_cc_vel(finest_level+1);

    if (containerHasElement(plot_var_names_3d, "x_velocity") ||
        containerHasElement(plot_var_names_3d, "y_velocity") ||
        containerHasElement(plot_var_names_3d, "z_velocity") ||
        containerHasElement(plot_var_names_3d, "vorticity") ) {

        for (int lev = 0; lev <= finest_level; ++lev) {
            mf_cc_vel[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, IntVect(1,1,0));
            mf_cc_vel[lev].setVal(zero); // FillBdyCCVels below leaves corners alone
            average_face_to_cellcenter(mf_cc_vel[lev],0,
                                       Array<const MultiFab*,3>{xvel_new[lev],yvel_new[lev],zvel_new[lev]},IntVect(1,1,0));
            mf_cc_vel[lev].FillBoundary(geom[lev].periodicity());
        } // lev

        // Fill level 0 before the interpolation below carries its ghost cells onto
        // the finer levels. Matches the fill the vorticity tagging criterion uses.
        FillBdyCCVels(0, mf_cc_vel[0]);

        // We need ghost cells if computing vorticity
        amrex::Interpolater* mapper = &cell_cons_interp;
        if ( containerHasElement(plot_var_names_3d, "vorticity") ) {
            for (int lev = 1; lev <= finest_level; ++lev) {
                Vector<MultiFab*> fmf = {&(mf_cc_vel[lev]), &(mf_cc_vel[lev])};
                Vector<Real> ftime    = {t_new[lev], t_new[lev]};
                Vector<MultiFab*> cmf = {&mf_cc_vel[lev-1], &mf_cc_vel[lev-1]};
                Vector<Real> ctime    = {t_new[lev], t_new[lev]};

                MultiFab mf_to_fill;
                amrex::FillPatchTwoLevels(mf_cc_vel[lev], mf_cc_vel[lev].nGrowVect(), IntVect(0,0,0),
                                          t_new[lev], cmf, ctime, fmf, ftime,
                                          0, 0, mf_cc_vel[lev].nComp(), geom[lev-1], geom[lev],
                                          refRatio(lev-1), mapper, domain_bcs_type, foextrap_bc());

                // Redo the reflections foextrap just overwrote at the domain boundary
                FillBdyCCVels(lev, mf_cc_vel[lev]);
            } // lev
        } // if
    } // if

    int icomp_rho = 0;
    for (auto plot_name : varnames_2d_rho)
    {
        if (plot_name == "zeta" ) {
            for (int lev = 0; lev <= finest_level; ++lev) { MultiFab::Copy(mf_2d_rho[lev],*vec_Zt_avg1[lev],0,icomp_rho,1,0); }
            icomp_rho++;
        }
        if (plot_name == "h" ) {
            for (int lev = 0; lev <= finest_level; ++lev) { MultiFab::Copy(mf_2d_rho[lev],*vec_h[lev],0,icomp_rho,1,0); }
            icomp_rho++;
        }
        if (plot_name == "f" ) {
            for (int lev = 0; lev <= finest_level; ++lev) { MultiFab::Copy(mf_2d_rho[lev],*vec_fcor[lev],0,icomp_rho,1,0); }
            icomp_rho++;
        }
        if (plot_name == "visc2" ) {
            for (int lev = 0; lev <= finest_level; ++lev) {
                if (vec_visc2_r[lev]->contains_nan(0, 1, 0, true) || vec_visc2_r[lev]->contains_inf(0, 1, 0, true)) {
                    amrex::Abort("Found while writing output: visc2 contains nan or inf");
                }
                for (MFIter mfi(mf_2d_rho[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                    const Box& bx = mfi.validbox();
                    const int K = mfi.index();
                    auto dst = mf_2d_rho[lev].array(mfi, icomp_rho);
                    auto src = vec_visc2_r[lev]->const_array(K);
                    ParallelFor(makeSlab(bx,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept {
                        dst(i,j,0) = src(i,j,0);
                   });
                }
            }
            icomp_rho++;
        }
        for (int n = 0; n < ncons; ++n) {
            const std::string diff2_name = std::string("diff2_") + cons_names[n];
            if (plot_name == diff2_name) {
                for (int lev = 0; lev <= finest_level; ++lev) {
                    for (MFIter mfi(mf_2d_rho[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                        const Box& bx = mfi.validbox();
                        const int K = mfi.index();
                        auto dst = mf_2d_rho[lev].array(mfi, icomp_rho);
                        auto src = vec_diff2[lev]->const_array(K);
                        ParallelFor(makeSlab(bx,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept {
                            dst(i,j,0) = src(i,j,0,n);
                        });
                    }
                }
                icomp_rho++;
            }
        }
        for (int n = 0; n < ncons; ++n) {
            const std::string stflux_name = std::string("stflux_") + cons_names[n];
            if (plot_name == stflux_name) {
                for (int lev = 0; lev <= finest_level; ++lev) {
                    MultiFab::Copy(mf_2d_rho[lev],*vec_stflux[lev],n,icomp_rho,1,0);
                }
                icomp_rho++;
            }
        }
        if (plot_name == "lrflux" ) {
            if (!solverChoice.bulk_fluxes && !solverChoice.atm2ocn_flux_mode) {
                amrex::Abort("Attempting to write longwave radiation flux to plotfile. Variable not allocated when bulk_fluxes turned off");
            }
            for (int lev = 0; lev <= finest_level; ++lev) { MultiFab::Copy(mf_2d_rho[lev],*vec_lrflx[lev],0,icomp_rho,1,0); }
            icomp_rho++;
        }
        if (plot_name == "lhflux" ) {
            if (!solverChoice.bulk_fluxes && !solverChoice.atm2ocn_flux_mode) {
                amrex::Abort("Attempting to write latent heat flux to plotfile. Variable not allocated when bulk_fluxes turned off");
            }
            for (int lev = 0; lev <= finest_level; ++lev) { MultiFab::Copy(mf_2d_rho[lev],*vec_lhflx[lev],0,icomp_rho,1,0); }
            icomp_rho++;
        }
        if (plot_name == "srflux" ) {
            if (!solverChoice.bulk_fluxes && !solverChoice.atm2ocn_flux_mode) {
                amrex::Abort("Attempting to write shortwave radiation flux to plotfile. Variable not allocated when bulk_fluxes turned off");
            }
            for (int lev = 0; lev <= finest_level; ++lev) { MultiFab::Copy(mf_2d_rho[lev],*vec_srflx[lev],0,icomp_rho,1,0); }
            icomp_rho++;
        }
        if (plot_name == "shflux" ) {
            if (!solverChoice.bulk_fluxes && !solverChoice.atm2ocn_flux_mode) {
                amrex::Abort("Attempting to write sensible heat flux to plotfile. Variable not allocated when bulk_fluxes turned off");
            }
            for (int lev = 0; lev <= finest_level; ++lev) { MultiFab::Copy(mf_2d_rho[lev],*vec_shflx[lev],0,icomp_rho,1,0); }
            icomp_rho++;
        }
        if (plot_name == "mask_rho" ) {
            for (int lev = 0; lev <= finest_level; ++lev) { MultiFab::Copy(mf_2d_rho[lev],*vec_mskr[lev],0,icomp_rho,1,0); }
            icomp_rho++;
        }
    }

    int icomp_u   = 0;
    for (auto plot_name : varnames_2d_u)
    {
        if (plot_name == "ubar" ) {
            for (int lev = 0; lev <= finest_level; ++lev) {
                MultiFab::Copy(mf_2d_u[lev],*vec_ubar[lev],0,icomp_u,1,0);
            }
            icomp_u++;
        }
        if (plot_name == "sustr" ) {
            for (int lev = 0; lev <= finest_level; ++lev) { MultiFab::Copy(mf_2d_u[lev],*vec_sustr[lev],0,icomp_u,1,0); }
            icomp_u++;
        }
        if (plot_name == "bustr" ) {
            for (int lev = 0; lev <= finest_level; ++lev) { MultiFab::Copy(mf_2d_u[lev],*vec_bustr[lev],0,icomp_u,1,0); }
            icomp_u++;
        }
        if (plot_name == "mask_u" ) {
            for (int lev = 0; lev <= finest_level; ++lev) { MultiFab::Copy(mf_2d_u[lev],*vec_msku[lev],0,icomp_u,1,0); }
            icomp_u++;
        }
    }

    int icomp_v   = 0;
    for (auto plot_name : varnames_2d_v)
    {
        if (plot_name == "vbar" ) {
            for (int lev = 0; lev <= finest_level; ++lev) { MultiFab::Copy(mf_2d_v[lev],*vec_vbar[lev],0,icomp_v,1,0); }
            icomp_v++;
        }
        if (plot_name == "svstr" ) {
            for (int lev = 0; lev <= finest_level; ++lev) { MultiFab::Copy(mf_2d_v[lev],*vec_svstr[lev],0,icomp_v,1,0); }
            icomp_v++;
        }
        if (plot_name == "bvstr" ) {
            for (int lev = 0; lev <= finest_level; ++lev) { MultiFab::Copy(mf_2d_v[lev],*vec_bvstr[lev],0,icomp_v,1,0); }
            icomp_v++;
        }
        if (plot_name == "mask_v" ) {
            for (int lev = 0; lev <= finest_level; ++lev) { MultiFab::Copy(mf_2d_v[lev],*vec_mskv[lev],0,icomp_v,1,0); }
            icomp_v++;
        }
    }

    for (int lev = 0; lev <= finest_level; ++lev)
    {
        int mf_comp = 0;

        AMREX_ALWAYS_ASSERT(cons_names.size() == ncons);

        // Check each tracer we are about to write, so the abort can name the offender.
        // One pass per component, rather than an all-component scan repeated once for
        // every variable in the plot list. nGrowVect() keeps the ghost-cell coverage the
        // previous no-argument contains_nan() had.
        for (int i = 0; i < ncons; ++i) {
            if (!containerHasElement(plot_var_names_3d, cons_names[i])) { continue; }
            const IntVect ng = cons_new[lev]->nGrowVect();
            if (cons_new[lev]->contains_nan(i,1,ng) || cons_new[lev]->contains_inf(i,1,ng)) {
                amrex::Abort("Found while writing output: " + cons_names[i] +
                             " contains nan or inf");
            }
        }

        // First, copy any of the conserved state variables into the output plotfile
        for (int i = 0; i < ncons; ++i) {
            if (containerHasElement(plot_var_names_3d, cons_names[i])) {
                MultiFab::Copy(plotMF[lev],*cons_new[lev],i,mf_comp,1,ngrow_vars);
                mf_comp++;
            }
        } // ncons

        // Next, check for velocities
        if (containerHasElement(plot_var_names_3d, "x_velocity")) {
            if (mf_cc_vel[lev].contains_nan(0,1) || mf_cc_vel[lev].contains_inf(0,1)) {
                amrex::Abort("Found while writing output: u velocity contains nan or inf");
            }
            MultiFab::Copy(plotMF[lev], mf_cc_vel[lev], 0, mf_comp, 1, 0);
            mf_comp += 1;
        }
        if (containerHasElement(plot_var_names_3d, "y_velocity")) {
            if (mf_cc_vel[lev].contains_nan(1,1) || mf_cc_vel[lev].contains_inf(1,1)) {
                amrex::Abort("Found while writing output: v velocity contains nan or inf");
            }
            MultiFab::Copy(plotMF[lev], mf_cc_vel[lev], 1, mf_comp, 1, 0);
            mf_comp += 1;
        }
        if (containerHasElement(plot_var_names_3d, "z_velocity")) {
            if (mf_cc_vel[lev].contains_nan(2,1) || mf_cc_vel[lev].contains_inf(2,1)) {
                amrex::Abort("Found while writing output: z velocity contains nan or inf");
            }
            MultiFab::Copy(plotMF[lev], mf_cc_vel[lev], 2, mf_comp, 1, 0);
            mf_comp += 1;
        }

        // Fill cell-centered location
        Real dx = Geom()[lev].CellSizeArray()[0];
        Real dy = Geom()[lev].CellSizeArray()[1];

        // Next, check for location names -- if we write one we write all
        // Note: the locations must be filled before the derived variables, to match
        //       the order of the names built in set3DPlotVariables
        if (containerHasElement(plot_var_names_3d, "x_cc") ||
            containerHasElement(plot_var_names_3d, "y_cc") ||
            containerHasElement(plot_var_names_3d, "z_cc"))
        {
            MultiFab dmf(plotMF[lev], make_alias, mf_comp, AMREX_SPACEDIM);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(dmf, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.tilebox();
                const Array4<Real> loc_arr = dmf.array(mfi);
                const Array4<Real const> zp_arr = vec_z_phys_nd[lev]->const_array(mfi);

                const Real xlo = Geom()[lev].ProbLoArray()[0];
                const Real ylo = Geom()[lev].ProbLoArray()[1];
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    loc_arr(i,j,k,0) = xlo + (i+Real(0.5)) * dx;
                    loc_arr(i,j,k,1) = ylo + (j+Real(0.5)) * dy;
                    loc_arr(i,j,k,2) = Real(0.125) * (zp_arr(i,j  ,k  ) + zp_arr(i+1,j  ,k  ) +
                                                   zp_arr(i,j+1,k  ) + zp_arr(i+1,j+1,k  ) +
                                                   zp_arr(i,j  ,k+1) + zp_arr(i+1,j  ,k+1) +
                                                   zp_arr(i,j+1,k+1) + zp_arr(i+1,j+1,k+1) );
                });
            } // mfi
            mf_comp += AMREX_SPACEDIM;
        } // if containerHasElement

        // Define standard process for calling the functions in Derive.cpp
        auto calculate_derived = [&](const std::string& der_name,
                                     decltype(derived::remora_dernull)& der_function)
        {
            if (containerHasElement(plot_var_names_3d, der_name)) {
                MultiFab dmf(plotMF[lev], make_alias, mf_comp, 1);
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                for (MFIter mfi(dmf, TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    const Box& bx = mfi.tilebox();
                    auto& dfab = dmf[mfi];

                    if (der_name == "vorticity") {
                        auto const& sfab = mf_cc_vel[lev][mfi];
                        der_function(bx, dfab, 0, 1, sfab, vec_pm[lev]->const_array(mfi), vec_pn[lev]->const_array(mfi), vec_mskr[lev]->const_array(mfi), Geom(lev), t_new[0], nullptr, lev);
                    } else {
                        auto const& sfab = (*cons_new[lev])[mfi];
                        der_function(bx, dfab, 0, 1, sfab, vec_pm[lev]->const_array(mfi), vec_pn[lev]->const_array(mfi), vec_mskr[lev]->const_array(mfi), Geom(lev), t_new[0], nullptr, lev);
                    }
                }

                mf_comp++;
            }
        };

        // Note: All derived variables must be computed in order of "derived_names" defined in REMORA.H
        calculate_derived("vorticity",  derived::remora_dervort);

#ifdef REMORA_USE_PARTICLES
        const auto& particles_namelist( particleData.getNames() );
        for (ParticlesNamesVector::size_type i = 0; i < particles_namelist.size(); i++) {
            if (containerHasElement(plot_var_names_3d, std::string(particles_namelist[i]+"_count")))
            {
                MultiFab temp_dat(plotMF[lev].boxArray(), plotMF[lev].DistributionMap(), 1, 0);
                temp_dat.setVal(0);
                particleData[particles_namelist[i]]->Increment(temp_dat, lev);
                MultiFab::Copy(plotMF[lev], temp_dat, 0, mf_comp, 1, 0);
                mf_comp += 1;
            }
        }

        Vector<std::string> particle_mesh_plot_names(0);
        particleData.GetMeshPlotVarNames( particle_mesh_plot_names );
        for (int i = 0; i < particle_mesh_plot_names.size(); i++) {
            std::string plot_var_name(particle_mesh_plot_names[i]);
            if (containerHasElement(plot_var_names_3d, plot_var_name) ) {
                MultiFab temp_dat(plotMF[lev].boxArray(), plotMF[lev].DistributionMap(), 1, 1);
                temp_dat.setVal(0);
                particleData.GetMeshPlotVar(plot_var_name, temp_dat, lev);
                MultiFab::Copy(plotMF[lev], temp_dat, 0, mf_comp, 1, 0);
                mf_comp += 1;
            }
        }
#endif

        if (plot_nodal_data) {
            fill_nodal_z_displacement(mf_nd[lev], *vec_z_phys_nd[lev], Geom(lev), 1);
        }
    } // lev

    if (plotfile_type == PlotfileType::amrex)
    {

    std::string plotfilename = Concatenate(plot_file_name, istep[0], file_min_digits);

    if (finest_level == 0)
    {
        if (plotfile_type == PlotfileType::amrex) {
            amrex::Print() << "Writing plotfile " << plotfilename << "\n";
            if (check_plot_z) { check_plot_nodal_z(GetVecOfConstPtrs(mf_nd), Geom()); }
            WriteMultiLevelPlotfileWithBathymetry(plotfilename, finest_level+1,
                                                  GetVecOfConstPtrs(plotMF),
                                                  GetVecOfConstPtrs(mf_nd),
                                                  GetVecOfConstPtrs(mf_u),
                                                  GetVecOfConstPtrs(mf_v),
                                                  GetVecOfConstPtrs(mf_w),
                                                  GetVecOfConstPtrs(mf_2d_rho),
                                                  GetVecOfConstPtrs(mf_2d_u),
                                                  GetVecOfConstPtrs(mf_2d_v),
                                                  varnames_3d, varnames_2d_rho,
                                                  varnames_2d_u, varnames_2d_v,
                                                  Geom(),
                                                  t_new[0], istep, refRatio());
            writeJobInfo(plotfilename);

#ifdef REMORA_USE_PARTICLES
            particleData.Checkpoint(plotfilename);
#endif

        }

    } else { // multilevel
        if (plotfile_type == PlotfileType::amrex) {
            amrex::Print() << "Writing plotfile " << plotfilename << "\n";
            int lev0 = 0;
            [[maybe_unused]] int desired_ratio = std::max(std::max(ref_ratio[lev0][0],ref_ratio[lev0][1]),ref_ratio[lev0][2]);
            bool any_ratio_one = ( ( (ref_ratio[lev0][0] == 1) || (ref_ratio[lev0][1] == 1) ) ||
                                     (ref_ratio[lev0][2] == 1) );
            for (int lev = 1; lev < finest_level; lev++) {
                any_ratio_one = any_ratio_one ||
                                     ( ( (ref_ratio[lev][0] == 1) || (ref_ratio[lev][1] == 1) ) ||
                                         (ref_ratio[lev][2] == 1) );
            }
            if (any_ratio_one && expand_plotvars_to_unif_rr) {
                Vector<IntVect>   r2(finest_level);
                Vector<Geometry>  g2(finest_level+1);
                Vector<MultiFab> mf2(finest_level+1);

                mf2[0].define(grids[0], dmap[0], ncomp_mf_3d, 0);

                // Copy level 0 as is
                MultiFab::Copy(mf2[0],plotMF[0],0,0,plotMF[0].nComp(),0);

                // Define a new multi-level array of Geometry's so that we pass the new "domain" at lev > 0
                Array<int,AMREX_SPACEDIM> periodicity =
                             {Geom()[0].isPeriodic(0),Geom()[0].isPeriodic(1),Geom()[0].isPeriodic(2)};
                g2[0].define(Geom()[0].Domain(),&(Geom()[0].ProbDomain()),0,periodicity.data());

                r2[0] = IntVect(1,1,ref_ratio[0][0]);
                for (int lev = 1; lev <= finest_level; ++lev) {
                    if (lev > 1) {
                        r2[lev-1][0] = 1;
                        r2[lev-1][1] = 1;
                        r2[lev-1][2] = r2[lev-2][2] * ref_ratio[lev-1][0];
                    }

                    mf2[lev].define(refine(grids[lev],r2[lev-1]), dmap[lev], ncomp_mf_3d, 0);

                    // Set the new problem domain
                    Box d2(Geom()[lev].Domain());
                    d2.refine(r2[lev-1]);

                    g2[lev].define(d2,&(Geom()[lev].ProbDomain()),0,periodicity.data());
                }

                // Everything the Header describes with g2 has to live on g2's grid, not only the
                // cell data: a viewer places level lev's nodes with g2[lev]'s dz, so the nodal
                // displacement and the face velocities are refined in z by the same ratio.
                Vector<MultiFab> mf_nd2(finest_level+1);
                Vector<MultiFab> mf_u2(finest_level+1), mf_v2(finest_level+1), mf_w2(finest_level+1);
                Vector<const MultiFab*> nd2_ptrs = GetVecOfConstPtrs(mf_nd);
                Vector<const MultiFab*> u2_ptrs  = GetVecOfConstPtrs(mf_u);
                Vector<const MultiFab*> v2_ptrs  = GetVecOfConstPtrs(mf_v);
                Vector<const MultiFab*> w2_ptrs  = GetVecOfConstPtrs(mf_w);
                for (int lev = 1; lev <= finest_level; ++lev) {
                    const int rz = r2[lev-1][2];
                    if (plot_nodal_data) {
                        BoxArray nodal2(mf2[lev].boxArray()); nodal2.surroundingNodes();
                        mf_nd2[lev].define(nodal2, dmap[lev], AMREX_SPACEDIM, 0);
                        mf_nd2[lev].setVal(zero);
                        fill_nodal_z_displacement(mf_nd2[lev], *vec_z_phys_nd[lev], g2[lev], rz);
                        nd2_ptrs[lev] = &mf_nd2[lev];
                    }
                    if (plot_staggered_vels) {
                        mf_u2[lev].define(refine(mf_u[lev].boxArray(), r2[lev-1]), dmap[lev], 1, 0);
                        mf_v2[lev].define(refine(mf_v[lev].boxArray(), r2[lev-1]), dmap[lev], 1, 0);
                        mf_w2[lev].define(refine(mf_w[lev].boxArray(), r2[lev-1]), dmap[lev], 1, 0);
                        refine_in_z(mf_u[lev], mf_u2[lev], rz);
                        refine_in_z(mf_v[lev], mf_v2[lev], rz);
                        refine_in_z(mf_w[lev], mf_w2[lev], rz);
                        u2_ptrs[lev] = &mf_u2[lev];
                        v2_ptrs[lev] = &mf_v2[lev];
                        w2_ptrs[lev] = &mf_w2[lev];
                    }
                }

                // Make a vector of BCRec with default values so we can use it here -- note the values
                //      aren't actually used because we do PCInterp
                amrex::Vector<amrex::BCRec> null_dom_bcs;
                null_dom_bcs.resize(mf2[0].nComp());
                for (int n = 0; n < mf2[0].nComp(); n++) {
                    for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
                        null_dom_bcs[n].setLo(dir, REMORABCType::int_dir);
                        null_dom_bcs[n].setHi(dir, REMORABCType::int_dir);
                    }
                }

                // Do piecewise interpolation of mf into mf2
                for (int lev = 1; lev <= finest_level; ++lev) {
                    Interpolater* mapper_c = &pc_interp;
                    InterpFromCoarseLevel(mf2[lev], t_new[lev], plotMF[lev],
                                          0, 0, mf2[lev].nComp(),
                                          geom[lev], g2[lev],
                                          null_bc_for_fill, 0, null_bc_for_fill, 0,
                                          r2[lev-1], mapper_c, null_dom_bcs, 0);
                }

                // Define an effective ref_ratio which is isotropic to be passed into WriteMultiLevelPlotfile
                Vector<IntVect> rr(finest_level);
                for (int lev = 0; lev < finest_level; ++lev) {
                    rr[lev] = IntVect(ref_ratio[lev][0],ref_ratio[lev][1],ref_ratio[lev][0]);
                }

                if (check_plot_z) { check_plot_nodal_z(nd2_ptrs, g2); }
                WriteMultiLevelPlotfileWithBathymetry(plotfilename, finest_level+1,
                                                      GetVecOfConstPtrs(mf2),
                                                      nd2_ptrs, u2_ptrs, v2_ptrs, w2_ptrs,
                                                      GetVecOfConstPtrs(mf_2d_rho),
                                                      GetVecOfConstPtrs(mf_2d_u),
                                                      GetVecOfConstPtrs(mf_2d_v),
                                                      varnames_3d, varnames_2d_rho,
                                                      varnames_2d_u, varnames_2d_v,
                                                      g2,
                                                      t_new[0], istep, rr);
                writeJobInfo(plotfilename);

#ifdef REMORA_USE_PARTICLES
                particleData.Checkpoint(plotfilename);
#endif
            } else {
                if (check_plot_z) { check_plot_nodal_z(GetVecOfConstPtrs(mf_nd), Geom()); }
                WriteMultiLevelPlotfileWithBathymetry(plotfilename, finest_level+1,
                                                      GetVecOfConstPtrs(plotMF),
                                                      GetVecOfConstPtrs(mf_nd),
                                                      GetVecOfConstPtrs(mf_u),
                                                      GetVecOfConstPtrs(mf_v),
                                                      GetVecOfConstPtrs(mf_w),
                                                      GetVecOfConstPtrs(mf_2d_rho),
                                                      GetVecOfConstPtrs(mf_2d_u),
                                                      GetVecOfConstPtrs(mf_2d_v),
                                                      varnames_3d, varnames_2d_rho,
                                                      varnames_2d_u, varnames_2d_v,
                                                      Geom(),
                                                      t_new[0], istep, ref_ratio);
                writeJobInfo(plotfilename);
#ifdef REMORA_USE_PARTICLES
                particleData.Checkpoint(plotfilename);
#endif
            }
        }
    } // end multi-level

    }
#ifdef REMORA_USE_NETCDF
    else if (plotfile_type == PlotfileType::netcdf)
    {
        // Each level goes to its own file: level 0 to _d01, a refined level to _d02 and up.
        // A level with more than one box is not handled; only its first subdomain is written.
        for (int lev = 0; lev <= finest_level; ++lev) {
            plotMF[lev].FillBoundary(geom[lev].periodicity());
            WriteNCPlotFile(istep_for_plot,&plotMF[lev],lev);
        }
    } // end if plotfile_type == netcdf
#endif
    if (plotfile_type == PlotfileType::amrex) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            mask_arrays_for_write(lev, zero, plotfile_fill_value);
        }
    } else if (plotfile_type == PlotfileType::netcdf) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            mask_arrays_for_write(lev, zero, netcdf_fill_value);
        }
    } else {
        amrex::Abort("Don't know this plotfile type");
    }
}

/**
 * @param plotfilename    name of plotfile to write to
 * @param nlevels         number of levels to write out
 * @param mf              MultiFab of data to write out
 * @param mf_nd           Multifab of nodal data to write out
 * @param varnames_3d     3D variable names to write out
 * @param varnames_2d_rho 2D cell-centered variable names to write out
 * @param varnames_2d_u   2D x-face-based variable names to write out
 * @param varnames_2d_v   2D y-face-based variable names to write out
 * @param my_geom         geometry to use for writing plotfile
 * @param time            time at which to output
 * @param level_steps     vector over level of iterations
 * @param rr              refinement ratio to use for writing plotfile
 * @param versionName     version string for VisIt
 * @param levelPrefix     string to prepend to level number
 * @param mfPrefix        subdirectory for multifab data
 * @param extra_dirs      additional subdirectories within plotfile
 */
 void
 REMORA::WriteMultiLevelPlotfileWithBathymetry (const std::string& plotfilename, int nlevels,
                                               const Vector<const MultiFab*>& mf,
                                               const Vector<const MultiFab*>& mf_nd,
                                               const Vector<const MultiFab*>& mf_u,
                                               const Vector<const MultiFab*>& mf_v,
                                               const Vector<const MultiFab*>& mf_w,
                                               const Vector<const MultiFab*>& mf_2d_rho,
                                               const Vector<const MultiFab*>& mf_2d_u,
                                               const Vector<const MultiFab*>& mf_2d_v,
                                               const Vector<std::string>& varnames_3d,
                                               const Vector<std::string>& varnames_2d_rho,
                                               const Vector<std::string>& varnames_2d_u,
                                               const Vector<std::string>& varnames_2d_v,
                                               const Vector<Geometry>& my_geom,
                                               Real time,
                                               const Vector<int>& level_steps,
                                               const Vector<IntVect>& rr,
                                               const std::string &versionName,
                                               const std::string &levelPrefix,
                                               const std::string &mfPrefix,
                                               const Vector<std::string>& extra_dirs) const
{
    BL_PROFILE("WriteMultiLevelPlotfileWithBathymetry()");

    AMREX_ASSERT(nlevels <= mf.size());
    AMREX_ASSERT(nlevels <= ref_ratio.size()+1);
    AMREX_ASSERT(nlevels <= level_steps.size());

    AMREX_ASSERT(mf[0]->nComp() == varnames_3d.size());

    // Every extra set must live on the grid the Header describes with my_geom, or a viewer
    // places it wrongly: the nodal set has one more node layer than the level has cells.
    for (int level = 0; level < nlevels; ++level) {
        if (plot_nodal_data) {
            const Box nb = mf_nd[level]->boxArray().minimalBox();
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                mf_nd[level]->boxArray().ixType().nodeCentered() &&
                nb.length(2) == my_geom[level].Domain().length(2) + 1,
                "plotfile nodal set does not match the geometry written to the Header");
        }
        if (plot_staggered_vels) {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                mf_u[level]->boxArray().minimalBox().length(2) == my_geom[level].Domain().length(2) &&
                mf_w[level]->boxArray().minimalBox().length(2) == my_geom[level].Domain().length(2) + 1,
                "plotfile face sets do not match the geometry written to the Header");
        }
    }

    bool callBarrier(false);
    PreBuildDirectorHierarchy(plotfilename, levelPrefix, nlevels, callBarrier);
    if (!extra_dirs.empty()) {
        for (const auto& d : extra_dirs) {
            const std::string ed = plotfilename+"/"+d;
            PreBuildDirectorHierarchy(ed, levelPrefix, nlevels, callBarrier);
        }
    }
    ParallelDescriptor::Barrier();

    if (ParallelDescriptor::MyProc() == ParallelDescriptor::NProcs()-1) {
        Vector<BoxArray> boxArrays(nlevels);
        for(int level(0); level < boxArrays.size(); ++level) {
            boxArrays[level] = mf[level]->boxArray();
        }

        auto f = [this, plotfilename, nlevels, boxArrays, varnames_3d,
                  varnames_2d_rho, varnames_2d_u, varnames_2d_v, my_geom,
                  time, level_steps, rr, versionName, levelPrefix, mfPrefix]() {
            VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
            std::string HeaderFileName(plotfilename + "/Header");
            std::ofstream HeaderFile;
            HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
            HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |
                                                    std::ofstream::trunc |
                                                    std::ofstream::binary);
            if( ! HeaderFile.good()) FileOpenFailed(HeaderFileName);
            WriteGenericPlotfileHeaderWithBathymetry(HeaderFile, nlevels, boxArrays, varnames_3d,
                                                     varnames_2d_rho, varnames_2d_u, varnames_2d_v,
                                                     my_geom, time, level_steps, rr, versionName,
                                                     levelPrefix, mfPrefix);
        };

        if (AsyncOut::UseAsyncOut()) {
            AsyncOut::Submit(std::move(f));
        } else {
            f();
        }
    }

    std::string mf_nodal_prefix = "Nu_nd";
    std::string mf_uface_prefix = "UFace";
    std::string mf_vface_prefix = "VFace";
    std::string mf_wface_prefix = "WFace";
    std::string mf_2d_rho_prefix = "rho2d";
    std::string mf_2d_u_prefix   = "u2d";
    std::string mf_2d_v_prefix   = "v2d";

    for (int level = 0; level <= finest_level; ++level)
    {
        if (AsyncOut::UseAsyncOut()) {
            VisMF::AsyncWrite(*mf[level],
                              MultiFabFileFullPrefix(level, plotfilename, levelPrefix, mfPrefix),
                              true);
            if (plot_nodal_data) {
                VisMF::AsyncWrite(*mf_nd[level],
                                  MultiFabFileFullPrefix(level, plotfilename, levelPrefix, mf_nodal_prefix),
                                  true);
            }
            if (plot_staggered_vels) {
                VisMF::AsyncWrite(*mf_u[level],
                                  MultiFabFileFullPrefix(level, plotfilename, levelPrefix, mf_uface_prefix),
                                  true);
                VisMF::AsyncWrite(*mf_v[level],
                                  MultiFabFileFullPrefix(level, plotfilename, levelPrefix, mf_vface_prefix),
                                  true);
                VisMF::AsyncWrite(*mf_w[level],
                                  MultiFabFileFullPrefix(level, plotfilename, levelPrefix, mf_wface_prefix),
                                  true);
            }
            if (mf_2d_rho[level]->nComp() > 0) {
                VisMF::AsyncWrite(*mf_2d_rho[level],
                                  MultiFabFileFullPrefix(level, plotfilename, levelPrefix, mf_2d_rho_prefix),
                                  true);
            }
            if (mf_2d_u[level]->nComp() > 0) {
                VisMF::AsyncWrite(*mf_2d_u[level],
                                  MultiFabFileFullPrefix(level, plotfilename, levelPrefix, mf_2d_u_prefix),
                                  true);
            }
            if (mf_2d_v[level]->nComp() > 0) {
                VisMF::AsyncWrite(*mf_2d_v[level],
                                  MultiFabFileFullPrefix(level, plotfilename, levelPrefix, mf_2d_v_prefix),
                                  true);
            }
        } else {
            const MultiFab* data;
            std::unique_ptr<MultiFab> mf_tmp;
            if (mf[level]->nGrowVect() != 0) {
                mf_tmp = std::make_unique<MultiFab>(mf[level]->boxArray(),
                                                    mf[level]->DistributionMap(),
                                                    mf[level]->nComp(), 0, MFInfo(),
                                                    mf[level]->Factory());
                MultiFab::Copy(*mf_tmp, *mf[level], 0, 0, mf[level]->nComp(), 0);
                data = mf_tmp.get();
            } else {
                data = mf[level];
            }
            VisMF::Write(*data       , MultiFabFileFullPrefix(level, plotfilename, levelPrefix, mfPrefix));
            if (plot_nodal_data) {
                VisMF::Write(*mf_nd[level], MultiFabFileFullPrefix(level, plotfilename, levelPrefix, mf_nodal_prefix));
            }
            if (plot_staggered_vels) {
                VisMF::Write(*mf_u[level], MultiFabFileFullPrefix(level, plotfilename, levelPrefix, mf_uface_prefix));
                VisMF::Write(*mf_v[level], MultiFabFileFullPrefix(level, plotfilename, levelPrefix, mf_vface_prefix));
                VisMF::Write(*mf_w[level], MultiFabFileFullPrefix(level, plotfilename, levelPrefix, mf_wface_prefix));
            }
            if (mf_2d_rho[level]->nComp() > 0) {
                VisMF::Write(*mf_2d_rho[level], MultiFabFileFullPrefix(level, plotfilename, levelPrefix, mf_2d_rho_prefix));
            }
            if (mf_2d_u[level]->nComp() > 0) {
                VisMF::Write(*mf_2d_u[level], MultiFabFileFullPrefix(level, plotfilename, levelPrefix, mf_2d_u_prefix));
            }
            if (mf_2d_v[level]->nComp() > 0) {
                VisMF::Write(*mf_2d_v[level], MultiFabFileFullPrefix(level, plotfilename, levelPrefix, mf_2d_v_prefix));
            }
        }
    } // level
}

/**
 * @param HeaderFile      output stream for header
 * @param nlevels         number of levels to write out
 * @param bArray          vector over levels of BoxArrays
 * @param varnames_3d     3D variable names to write out
 * @param varnames_2d     2D variable names to write out
 * @param my_geom         geometry to use for writing plotfile
 * @param time            time at which to output
 * @param level_steps     vector over level of iterations
 * @param my_ref_ratio    refinement ratio to use for writing plotfile
 * @param versionName     version string for VisIt
 * @param levelPrefix     string to prepend to level number
 * @param mfPrefix        subdirectory for multifab data
 */
void
REMORA::WriteGenericPlotfileHeaderWithBathymetry (std::ostream &HeaderFile,
                                                 [[maybe_unused]] int nlevels,
                                                 const Vector<BoxArray> &bArray,
                                                 const Vector<std::string> &varnames_3d,
                                                 const Vector<std::string> &varnames_2d_rho,
                                                 const Vector<std::string> &varnames_2d_u,
                                                 const Vector<std::string> &varnames_2d_v,
                                                 const Vector<Geometry>& my_geom,
                                                 Real time,
                                                 const Vector<int> &level_steps,
                                                 const Vector<IntVect>& my_ref_ratio,
                                                 const std::string &versionName,
                                                 const std::string &levelPrefix,
                                                 const std::string &mfPrefix) const
{
    AMREX_ASSERT(nlevels <= bArray.size());
    AMREX_ASSERT(nlevels <= ref_ratio.size()+1);
    AMREX_ASSERT(nlevels <= level_steps.size());

    int num_extra_mfs = plot_nodal_data ? 1 : 0; // for nodal, if it is written
    if (plot_staggered_vels) {
        num_extra_mfs += 3; // for nodal, which is always on
    }

    HeaderFile.precision(17);

    // ---- this is the generic plot file type name
    HeaderFile << versionName << '\n';

    HeaderFile << varnames_3d.size() << '\n';

    for (int ivar = 0; ivar < varnames_3d.size(); ++ivar) {
        HeaderFile << varnames_3d[ivar] << "\n";
    }
    HeaderFile << AMREX_SPACEDIM << '\n';
    HeaderFile << time << '\n';
    HeaderFile << finest_level << '\n';
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        HeaderFile << my_geom[0].ProbLo(i) << ' ';
    }
    HeaderFile << '\n';
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        HeaderFile << my_geom[0].ProbHi(i) << ' ';
    }
    HeaderFile << '\n';
    for (int i = 0; i < finest_level; ++i) {
        HeaderFile << my_ref_ratio[i][0] << ' ';
        }
    HeaderFile << '\n';
    for (int i = 0; i <= finest_level; ++i) {
        HeaderFile << my_geom[i].Domain() << ' ';
    }
    HeaderFile << '\n';
    for (int i = 0; i <= finest_level; ++i) {
            HeaderFile << level_steps[i] << ' ';
    }
    HeaderFile << '\n';
    for (int i = 0; i <= finest_level; ++i) {
        for (int k = 0; k < AMREX_SPACEDIM; ++k) {
            HeaderFile << my_geom[i].CellSize()[k] << ' ';
        }
        HeaderFile << '\n';
    }
    HeaderFile << (int) my_geom[0].Coord() << '\n';
    HeaderFile << "0\n";

    for (int level = 0; level <= finest_level; ++level) {
        HeaderFile << level << ' ' << bArray[level].size() << ' ' << time << '\n';
        HeaderFile << level_steps[level] << '\n';

        const IntVect& domain_lo = my_geom[level].Domain().smallEnd();
        for (int i = 0; i < bArray[level].size(); ++i)
        {
            // Need to shift because the RealBox ctor we call takes the
            // physical location of index (0,0,0).  This does not affect
            // the usual cases where the domain index starts with 0.
            const Box& b = shift(bArray[level][i], -domain_lo);
            RealBox loc = RealBox(b, my_geom[level].CellSize(), my_geom[level].ProbLo());
            for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                HeaderFile << loc.lo(n) << ' ' << loc.hi(n) << '\n';
            }
        }

        HeaderFile << MultiFabHeaderPath(level, levelPrefix, mfPrefix) << '\n';
    }
        HeaderFile << num_extra_mfs << "\n";
        if (plot_nodal_data) {
            HeaderFile << "3" << "\n";
            HeaderFile << "amrexvec_nu_x" << "\n";
            HeaderFile << "amrexvec_nu_y" << "\n";
            HeaderFile << "amrexvec_nu_z" << "\n";
            std::string mf_nodal_prefix = "Nu_nd";
            for (int level = 0; level <= finest_level; ++level) {
                HeaderFile << MultiFabHeaderPath(level, levelPrefix, mf_nodal_prefix) << '\n';
            }
        }
        if (plot_staggered_vels) {
            HeaderFile << "1" << "\n"; // number of components in the multifab
            HeaderFile << "u_vel" << "\n";
            std::string mf_uface_prefix = "UFace";
            for (int level = 0; level <= finest_level; ++level) {
                HeaderFile << MultiFabHeaderPath(level, levelPrefix, mf_uface_prefix) << '\n';
            }
            HeaderFile << "1" << "\n";
            HeaderFile << "v_vel" << "\n";
            std::string mf_vface_prefix = "VFace";
            for (int level = 0; level <= finest_level; ++level) {
                HeaderFile << MultiFabHeaderPath(level, levelPrefix, mf_vface_prefix) << '\n';
            }
            HeaderFile << "1" << "\n";
            HeaderFile << "w_vel" << "\n";
            std::string mf_wface_prefix = "WFace";
            for (int level = 0; level <= finest_level; ++level) {
                HeaderFile << MultiFabHeaderPath(level, levelPrefix, mf_wface_prefix) << '\n';
            }
        }

        if (varnames_2d_rho.size() > 0) {
            HeaderFile << varnames_2d_rho.size() << "\n"; // number of components in the 2D rho multifab
            for (int ivar = 0; ivar < varnames_2d_rho.size(); ++ivar) {
                HeaderFile << varnames_2d_rho[ivar] << "\n";
            }
            std::string mf_2d_rho_prefix = "rho2d";
            for (int level = 0; level <= finest_level; ++level) {
                HeaderFile << MultiFabHeaderPath(level, levelPrefix, mf_2d_rho_prefix) << "\n";
            }
        }

        if (varnames_2d_u.size() > 0) {
            HeaderFile << varnames_2d_u.size() << "\n"; // number of components in the 2D rho multifab
            for (int ivar = 0; ivar < varnames_2d_u.size(); ++ivar) {
                HeaderFile << varnames_2d_u[ivar] << "\n";
            }
            std::string mf_2d_u_prefix = "u2d";
            for (int level = 0; level <= finest_level; ++level) {
                HeaderFile << MultiFabHeaderPath(level, levelPrefix, mf_2d_u_prefix) << "\n";
            }
        }

        if (varnames_2d_v.size() > 0) {
            HeaderFile << varnames_2d_v.size() << "\n"; // number of components in the 2D v multifab
            for (int ivar = 0; ivar < varnames_2d_v.size(); ++ivar) {
                HeaderFile << varnames_2d_v[ivar] << "\n";
            }
            std::string mf_2d_v_prefix = "v2d";
            for (int level = 0; level <= finest_level; ++level) {
                HeaderFile << MultiFabHeaderPath(level, levelPrefix, mf_2d_v_prefix) << "\n";
            }
        }
}

/**
 * Check the nodal z a viewer rebuilds from what WritePlotFile is about to write.
 *
 * nd[lev] is the nodal displacement on the grid g[lev] describes -- the native geometry, or g2
 * on the expand_plotvars_to_unif_rr path -- so a viewer's node is z = ProbLo(2) + k*dz + nu_z.
 * Per level, the lowest and highest node layers must reproduce z_phys_nd's bottom and top to
 * rounding. Per level pair, at every fine node that coincides with a coarse node, |z_f - z_c|
 * is taken over the patch interior and its perimeter separately: a perimeter node averages
 * fine ghost columns the parent interpolated, so it is reported, and only the interior is held
 * to check_plot_z_tol.
 *
 * @param[in] nd  nodal displacement per level, as passed to the writer
 * @param[in] g   geometry per level, as passed to the writer
 */
void
REMORA::check_plot_nodal_z (const Vector<const MultiFab*>& nd, const Vector<Geometry>& g)
{
    if (!plot_nodal_data) { return; }

    Vector<MultiFab> z(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        const MultiFab& nu = *nd[lev];
        z[lev].define(nu.boxArray(), nu.DistributionMap(), 1, 0);
        const Real dz  = g[lev].CellSizeArray()[2];
        const Real zlo = g[lev].ProbLoArray()[2];
        const int  klo = g[lev].Domain().smallEnd(2);
        const int  Nf  = g[lev].Domain().length(2);     // nodes 0..Nf on this grid
        const int  N   = Geom(lev).Domain().length(2);  // nodes 0..N on the native one

        ReduceOps<ReduceOpMax, ReduceOpMax> reduce_op;
        ReduceData<Real, Real> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;
        for (MFIter mfi(z[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            Array4<Real>       const& za = z[lev].array(mfi);
            Array4<Real const> const& na = nu.const_array(mfi);
            Array4<Real const> const& zp = vec_z_phys_nd[lev]->const_array(mfi);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                za(i,j,k) = zlo + Real(k - klo) * dz + na(i,j,k,2);
            });
            reduce_op.eval(bx, reduce_data, [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                const Real eb = (k == klo)      ? std::abs(za(i,j,k) - zp(i,j,0)) : Real(0.0);
                const Real et = (k == klo + Nf) ? std::abs(za(i,j,k) - zp(i,j,N)) : Real(0.0);
                return {eb, et};
            });
        }
        ReduceTuple hv = reduce_data.value(reduce_op);
        Real err_bot = amrex::get<0>(hv);
        Real err_top = amrex::get<1>(hv);
        ParallelDescriptor::ReduceRealMax(err_bot);
        ParallelDescriptor::ReduceRealMax(err_top);
        amrex::Print() << "Plot nodal z, level " << lev << ": |z_bottom - z_phys| "
                       << err_bot << ", |z_top - z_phys| " << err_top << "\n";
        const Real self_tol = Real(1.e-10) * std::max(Real(1.0), std::abs(zlo));
        if (err_bot > self_tol || err_top > self_tol) {
            amrex::Abort("check_plot_z: the plotfile's nodal z does not reproduce z_phys_nd");
        }
    }

    for (int lev = 1; lev <= finest_level; ++lev)
    {
        IntVect r;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            r[d] = g[lev].Domain().length(d) / g[lev-1].Domain().length(d);
        }

        // Sample the fine nodes that coincide with coarse nodes, then bring them onto the
        // coarse level's layout. A coarse node no fine box reaches keeps the sentinel.
        BoxArray cba = z[lev].boxArray(); cba.coarsen(r);
        MultiFab zfc(cba, z[lev].DistributionMap(), 1, 0);
        for (MFIter mfi(zfc, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            Array4<Real>       const& c = zfc.array(mfi);
            Array4<Real const> const& f = z[lev].const_array(mfi);
            const int r0 = r[0], r1 = r[1], r2 = r[2];
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                c(i,j,k) = f(i*r0, j*r1, k*r2);
            });
        }
        const Real sentinel = Real(1.e30);
        MultiFab zcc(z[lev-1].boxArray(), z[lev-1].DistributionMap(), 1, 0);
        zcc.setVal(sentinel);
        zcc.ParallelCopy(zfc);

        // A coarse node is interior to the patch when the fine level covers all four coarse
        // cells around it; refinement spans the whole depth, so one cell row answers for k.
        iMultiFab covered = makeFineMask(grids[lev-1], dmap[lev-1], IntVect(1,1,0), grids[lev],
                                         refRatio(lev-1), geom[lev-1].periodicity(), 0, 1);
        const int kcell = grids[lev-1].minimalBox().smallEnd(2);

        ReduceOps<ReduceOpMax, ReduceOpMax, ReduceOpSum, ReduceOpSum> reduce_op;
        ReduceData<Real, Real, int, int> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;
        for (MFIter mfi(zcc, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            Array4<Real const> const& a = zcc.const_array(mfi);
            Array4<Real const> const& b = z[lev-1].const_array(mfi);
            Array4<int  const> const& m = covered.const_array(mfi);
            reduce_op.eval(bx, reduce_data, [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                if (a(i,j,k) >= sentinel) { return {Real(0.0), Real(0.0), 0, 0}; }
                const bool interior = m(i-1,j-1,kcell) && m(i,j-1,kcell) &&
                                      m(i-1,j  ,kcell) && m(i,j  ,kcell);
                const Real d = std::abs(a(i,j,k) - b(i,j,k));
                return interior ? ReduceTuple{d, Real(0.0), 1, 0} : ReduceTuple{Real(0.0), d, 0, 1};
            });
        }
        ReduceTuple hv = reduce_data.value(reduce_op);
        Real d_int = amrex::get<0>(hv), d_per = amrex::get<1>(hv);
        int  n_int = amrex::get<2>(hv), n_per = amrex::get<3>(hv);
        ParallelDescriptor::ReduceRealMax(d_int);
        ParallelDescriptor::ReduceRealMax(d_per);
        ParallelDescriptor::ReduceIntSum(n_int);
        ParallelDescriptor::ReduceIntSum(n_per);
        amrex::Print() << "Plot nodal z, levels " << lev-1 << "/" << lev
                       << ": coarse-fine max |dz| interior " << d_int << " (" << n_int << " nodes)"
                       << ", perimeter " << d_per << " (" << n_per << " nodes)\n";
        if (check_plot_z_tol >= Real(0.0) && d_int > check_plot_z_tol) {
            amrex::Abort("check_plot_z: coarse and fine nodal z disagree inside the patch");
        }
    }
}

/**
 * @param lev          level to mask
 * @param fill_value   fill value to mask with
 * @param fill_where   value at cells where we will apply the mask. This is necessary because rivers
 */
void
REMORA::mask_arrays_for_write(int lev, Real fill_value, Real fill_where)
{
    for (MFIter mfi(*cons_new[lev],false); mfi.isValid(); ++mfi) {
        Box gbx1 = mfi.growntilebox(IntVect(NGROW+1,NGROW+1,0));
        Box gbx_coeff = mfi.growntilebox(IntVect(NGROW,NGROW,0));
        Box ubx = mfi.grownnodaltilebox(0,IntVect(NGROW,NGROW,0));
        Box vbx = mfi.grownnodaltilebox(1,IntVect(NGROW,NGROW,0));

        Array4<Real> const& Zt_avg1 = vec_Zt_avg1[lev]->array(mfi);
        Array4<Real> const& ubar = vec_ubar[lev]->array(mfi);
        Array4<Real> const& vbar = vec_vbar[lev]->array(mfi);
        Array4<Real> const& xvel = xvel_new[lev]->array(mfi);
        Array4<Real> const& yvel = yvel_new[lev]->array(mfi);
        Array4<Real> const& visc2 = vec_visc2_r[lev]->array(mfi);
        Array4<Real> const& diff2 = vec_diff2[lev]->array(mfi);
        Array4<Real> const& temp = cons_new[lev]->array(mfi,Temp_comp);
        Array4<Real> const& salt = cons_new[lev]->array(mfi,Salt_comp);

        Array4<Real const> const& mskr = vec_mskr[lev]->array(mfi);
        Array4<Real const> const& msku = vec_msku[lev]->array(mfi);
        Array4<Real const> const& mskv = vec_mskv[lev]->array(mfi);
        const int ncons_local = ncons;

        ParallelFor(makeSlab(gbx1,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
            if (mskr(i,j,0) == zero) {  // Explicitly compare to 0.0
                Zt_avg1(i,j,0) = fill_value;
            }
        });
        ParallelFor(gbx1, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            if (mskr(i,j,0) == zero) {  // Explicitly compare to 0.0
                temp(i,j,k) = fill_value;
                salt(i,j,k) = fill_value;
            }
        });
        ParallelFor(makeSlab(gbx_coeff,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
            if (mskr(i,j,0) == zero) {  // Explicitly compare to 0.0
                visc2(i,j,0) = fill_value;
                for (int n = 0; n < ncons_local; ++n) {
                    diff2(i,j,0,n) = fill_value;
                }
            }
        });
        ParallelFor(makeSlab(ubx,2,0), 3, [=] AMREX_GPU_DEVICE (int i, int j, int , int n)
        {
            if (msku(i,j,0) == zero && ubar(i,j,0)==fill_where) {  // Explicitly compare to 0.0
                ubar(i,j,0,n) = fill_value;
            }
        });
        ParallelFor(makeSlab(vbx,2,0), 3, [=] AMREX_GPU_DEVICE (int i, int j, int , int n)
        {
            if (mskv(i,j,0) == zero && vbar(i,j,0)==fill_where) {  // Explicitly compare to 0.0
                vbar(i,j,0,n) = fill_value;
            }
        });
        ParallelFor(ubx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            if (msku(i,j,0) == zero && xvel(i,j,k)==fill_where) {  // Explicitly compare to 0.0
                xvel(i,j,k) = fill_value;
            }
        });
        ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            if (mskv(i,j,0) == zero && yvel(i,j,k)==fill_where) {  // Explicitly compare to 0.0
                yvel(i,j,k) = fill_value;
            }
        });
    } // mfi
    Gpu::streamSynchronize();
}
