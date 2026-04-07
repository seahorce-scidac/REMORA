#include <REMORA.H>
#include "AMReX_Interp_3D_C.H"
#include "AMReX_PlotFileUtil.H"

using namespace amrex;

PhysBCFunctNoOp null_bc_for_fill;

template<typename V, typename T>
bool containerHasElement(const V& iterable, const T& query) {
    return std::find(iterable.begin(), iterable.end(), query) != iterable.end();
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

    // For scaled_to_grid, viscosity/diffusivity coefficients are vertically homogeneous and
    // time-invariant. For AMReX plotfiles, write them to a dedicated coefficient plotfile and
    // omit them from regular time-series plotfiles. For AMR, we write a coefficient plotfile
    // per plot output so it matches the current grid hierarchy.
    Vector<std::string> varnames_2d_coeff;
    const bool write_coeff_plotfile =
        (solverChoice.horiz_mixing_type == HorizMixingType::scaled_to_grid) &&
        (plotfile_type == PlotfileType::amrex);

    if (write_coeff_plotfile) {
        for (auto const& nm : varnames_2d) {
            if (nm == "visc2") {
                varnames_2d_coeff.push_back(nm);
                continue;
            }
            for (int n = 0; n < NCONS; ++n) {
                if (nm == std::string("diff2_") + cons_names[n]) {
                    varnames_2d_coeff.push_back(nm);
                    break;
                }
            }
        }

        varnames_2d.erase(std::remove_if(varnames_2d.begin(), varnames_2d.end(),
                                         [&] (std::string const& nm) {
                                             if (nm == "visc2") { return true; }
                                             for (int n = 0; n < NCONS; ++n) {
                                                 if (nm == std::string("diff2_") + cons_names[n]) { return true; }
                                             }
                                             return false;
                                         }),
                          varnames_2d.end());
    }

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
         if (plot_name == "visc2") {varnames_2d_rho.push_back(plot_name); ncomp_mf_2d_rho++;}
         if (plot_name == "diff2_temp"  ) {varnames_2d_rho.push_back(plot_name); ncomp_mf_2d_rho++;}
         if (plot_name == "diff2_salt"  ) {varnames_2d_rho.push_back(plot_name); ncomp_mf_2d_rho++;}
         if (plot_name == "diff2_tracer") {varnames_2d_rho.push_back(plot_name); ncomp_mf_2d_rho++;}
         if (plot_name == "ubar" ) {varnames_2d_u.push_back(plot_name); ncomp_mf_2d_u++;}
         if (plot_name == "sustr") {varnames_2d_u.push_back(plot_name); ncomp_mf_2d_u++;}
         if (plot_name == "bustr") {varnames_2d_u.push_back(plot_name); ncomp_mf_2d_u++;}
         if (plot_name == "vbar" ) {varnames_2d_v.push_back(plot_name); ncomp_mf_2d_v++;}
         if (plot_name == "svstr") {varnames_2d_v.push_back(plot_name); ncomp_mf_2d_v++;}
         if (plot_name == "bvstr") {varnames_2d_v.push_back(plot_name); ncomp_mf_2d_v++;}
      }
    }

    // We fillpatch here because some of the derived quantities require derivatives
    //     which require ghost cells to be filled. Don't fill the boundary, though.
    for (int lev = 0; lev <= finest_level; ++lev) {
        FillPatchNoBC(lev, t_new[lev], *cons_new[lev], cons_new, BdyVars::t,0,true,false);
        FillPatchNoBC(lev, t_new[lev], *xvel_new[lev], xvel_new, BdyVars::u,0,true,false);
        FillPatchNoBC(lev, t_new[lev], *yvel_new[lev], yvel_new, BdyVars::v,0,true,false);
        FillPatchNoBC(lev, t_new[lev], *zvel_new[lev], zvel_new, BdyVars::null,0,true,false);
        // These are constant-in-time fields for most runs, but we still ensure
        // ghost cells are valid before writing plotfiles.
        FillPatchNoBC(lev, t_new[lev], *vec_visc2_r[lev], GetVecOfPtrs(vec_visc2_r), BdyVars::null,0,true,false);
        FillPatchNoBC(lev, t_new[lev], *vec_diff2[lev],   GetVecOfPtrs(vec_diff2),   BdyVars::null,0,true,false);
    }

    Real fill_value = 0.0_rt;
    for (int lev = 0; lev <= finest_level; ++lev) {
        mask_arrays_for_write(lev, (Real) fill_value, 0.0_rt);
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

    // Coefficient plotfile (2D rho-point fields).
    //
    // NOTE: REMORA/AMReX is built in 3D (AMREX_SPACEDIM=3). For coefficient output we write a
    // k=0 slab (nz=1) so downstream tools see these fields as 2D (or 3D with a singleton z).
    const int ncomp_mf_2d_rho_coeff = static_cast<int>(varnames_2d_coeff.size());
    Vector<MultiFab> plotMF_zero;
    Vector<MultiFab> mf_2d_rho_coeff;
    Vector<MultiFab> mf_2d_u_zero;
    Vector<MultiFab> mf_2d_v_zero;
    Vector<MultiFab> mf_nd_coeff;
    Vector<Geometry> geom_coeff;
    if (write_coeff_plotfile && (ncomp_mf_2d_rho_coeff > 0)) {
        plotMF_zero.resize(finest_level+1);
        mf_2d_rho_coeff.resize(finest_level+1);
        mf_2d_u_zero.resize(finest_level+1);
        mf_2d_v_zero.resize(finest_level+1);
        mf_nd_coeff.resize(finest_level+1);
        geom_coeff.resize(finest_level+1);
        for (int lev = 0; lev <= finest_level; ++lev) {
            // Build a k=0 slab BoxArray for 2D coefficient output.
            BoxArray ba(grids[lev]);
            BoxList bl2d = ba.boxList();
            for (auto& b : bl2d) { b.setRange(2,0); }
            BoxArray ba2d(std::move(bl2d));

            // "3D" plot MF with zero components, defined on the slab BoxArray so the plotfile
            // header describes a 2D (nz=1) hierarchy.
            plotMF_zero[lev].define(ba2d, dmap[lev], 0, IntVect(0,0,0));

            mf_2d_rho_coeff[lev].define(ba2d, dmap[lev], ncomp_mf_2d_rho_coeff, IntVect(0,0,0));
            mf_2d_u_zero[lev].define  (ba2d, dmap[lev], 0, IntVect(0,0,0));
            mf_2d_v_zero[lev].define  (ba2d, dmap[lev], 0, IntVect(0,0,0));

            // Nodal coordinates for the coefficient plotfile, also on the slab hierarchy.
            BoxArray nodal_grids(ba2d);
            nodal_grids.surroundingNodes();
            mf_nd_coeff[lev].define(nodal_grids, dmap[lev], AMREX_SPACEDIM, 0);
            mf_nd_coeff[lev].setVal(0.0_rt);

            // Geometry for the coefficient plotfile: use a z-slab domain (nz=1).
            Box dom2d = Geom()[lev].Domain();
            dom2d.setRange(2,0);
            Array<int,AMREX_SPACEDIM> periodicity =
                {Geom()[lev].isPeriodic(0),Geom()[lev].isPeriodic(1),Geom()[lev].isPeriodic(2)};
            geom_coeff[lev].define(dom2d, &(Geom()[lev].ProbDomain()), Geom()[lev].Coord(), periodicity.data());
        }
    }

    // Array of MultiFabs for nodal data
    Vector<MultiFab> mf_nd(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        BoxArray nodal_grids(grids[lev]); nodal_grids.surroundingNodes();
        mf_nd[lev].define(nodal_grids, dmap[lev], AMREX_SPACEDIM, 0);
        mf_nd[lev].setVal(0.);
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
            mf_cc_vel[lev].setVal(0.0_rt); // zero out velocity in case we have any wall boundaries
            average_face_to_cellcenter(mf_cc_vel[lev],0,
                                       Array<const MultiFab*,3>{xvel_new[lev],yvel_new[lev],zvel_new[lev]},IntVect(1,1,0));
            mf_cc_vel[lev].FillBoundary(geom[lev].periodicity());
        } // lev

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
                                          refRatio(lev-1), mapper, domain_bcs_type, BCVars::foextrap_bc);
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
                     ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                         dst(i,j,k) = src(i,j,0,0);
                     });
                 }
             }
             icomp_rho++;
         }
         if (plot_name == "diff2_temp" ) {
             for (int lev = 0; lev <= finest_level; ++lev) {
                 for (MFIter mfi(mf_2d_rho[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                     const Box& bx = mfi.validbox();
                     const int K = mfi.index();
                     auto dst = mf_2d_rho[lev].array(mfi, icomp_rho);
                     auto src = vec_diff2[lev]->const_array(K);
                     ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                         dst(i,j,k) = src(i,j,0,Temp_comp);
                     });
                 }
             }
             icomp_rho++;
         }
         if (plot_name == "diff2_salt" ) {
             for (int lev = 0; lev <= finest_level; ++lev) {
                 for (MFIter mfi(mf_2d_rho[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                     const Box& bx = mfi.validbox();
                     const int K = mfi.index();
                     auto dst = mf_2d_rho[lev].array(mfi, icomp_rho);
                     auto src = vec_diff2[lev]->const_array(K);
                     ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                         dst(i,j,k) = src(i,j,0,Salt_comp);
                     });
                 }
             }
             icomp_rho++;
         }
         if (plot_name == "diff2_tracer" ) {
             for (int lev = 0; lev <= finest_level; ++lev) {
                 for (MFIter mfi(mf_2d_rho[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                     const Box& bx = mfi.validbox();
                     const int K = mfi.index();
                     auto dst = mf_2d_rho[lev].array(mfi, icomp_rho);
                     auto src = vec_diff2[lev]->const_array(K);
                     ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                         dst(i,j,k) = src(i,j,0,Tracer_comp);
                     });
                 }
             }
             icomp_rho++;
         }
    }

    int icomp_u   = 0;
    for (auto plot_name : varnames_2d_u)
    {
         if (plot_name == "ubar" ) {
             for (int lev = 0; lev <= finest_level; ++lev) {
                 MultiFab::Copy(mf_2d_u[lev],*vec_DU_avg1[lev],0,icomp_u,1,0);
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
    }

    int icomp_v   = 0;
    for (auto plot_name : varnames_2d_v)
    {
         if (plot_name == "vbar" ) {
             for (int lev = 0; lev <= finest_level; ++lev) { MultiFab::Copy(mf_2d_v[lev],*vec_DV_avg1[lev],0,icomp_v,1,0); }
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
    }

    // Fill the time-invariant coefficient plotfile (2D rho points) if enabled
    if (write_coeff_plotfile && (ncomp_mf_2d_rho_coeff > 0)) {
        int icomp_coeff = 0;
        for (auto const& nm : varnames_2d_coeff) {
            if (nm == "visc2") {
                for (int lev = 0; lev <= finest_level; ++lev) {
                    for (MFIter mfi(mf_2d_rho_coeff[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                        const Box& bx = mfi.validbox();
                        const int K = mfi.index();
                        auto dst = mf_2d_rho_coeff[lev].array(mfi, icomp_coeff);
                        auto src = vec_visc2_r[lev]->const_array(K);
                        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                            dst(i,j,k) = src(i,j,0,0);
                        });
                    }
                }
                icomp_coeff++;
                continue;
            }
            if (nm == "diff2_temp") {
                for (int lev = 0; lev <= finest_level; ++lev) {
                    for (MFIter mfi(mf_2d_rho_coeff[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                        const Box& bx = mfi.validbox();
                        const int K = mfi.index();
                        auto dst = mf_2d_rho_coeff[lev].array(mfi, icomp_coeff);
                        auto src = vec_diff2[lev]->const_array(K);
                        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                            dst(i,j,k) = src(i,j,0,Temp_comp);
                        });
                    }
                }
                icomp_coeff++;
                continue;
            }
            if (nm == "diff2_salt") {
                for (int lev = 0; lev <= finest_level; ++lev) {
                    for (MFIter mfi(mf_2d_rho_coeff[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                        const Box& bx = mfi.validbox();
                        const int K = mfi.index();
                        auto dst = mf_2d_rho_coeff[lev].array(mfi, icomp_coeff);
                        auto src = vec_diff2[lev]->const_array(K);
                        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                            dst(i,j,k) = src(i,j,0,Salt_comp);
                        });
                    }
                }
                icomp_coeff++;
                continue;
            }
            if (nm == "diff2_tracer") {
                for (int lev = 0; lev <= finest_level; ++lev) {
                    for (MFIter mfi(mf_2d_rho_coeff[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                        const Box& bx = mfi.validbox();
                        const int K = mfi.index();
                        auto dst = mf_2d_rho_coeff[lev].array(mfi, icomp_coeff);
                        auto src = vec_diff2[lev]->const_array(K);
                        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                            dst(i,j,k) = src(i,j,0,Tracer_comp);
                        });
                    }
                }
                icomp_coeff++;
                continue;
            }
        }
    }

    for (int lev = 0; lev <= finest_level; ++lev)
    {
        int mf_comp = 0;

        // First, copy any of the conserved state variables into the output plotfile
        AMREX_ALWAYS_ASSERT(cons_names.size() == NCONS);
        for (int i = 0; i < NCONS; ++i) {
            if (containerHasElement(plot_var_names_3d, cons_names[i])) {
                if (cons_new[lev]->contains_nan() || cons_new[lev]->contains_inf()) {
                    amrex::Abort("Found while writing output: Cons (salt, temp, or tracer, etc) contains nan or inf");
                }
                MultiFab::Copy(plotMF[lev],*cons_new[lev],i,mf_comp,1,ngrow_vars);
                mf_comp++;
            }
        } // NCONS

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

        // Fill cell-centered location
        Real dx = Geom()[lev].CellSizeArray()[0];
        Real dy = Geom()[lev].CellSizeArray()[1];

        // Next, check for location names -- if we write one we write all
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

                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    loc_arr(i,j,k,0) = (i+0.5_rt) * dx;
                    loc_arr(i,j,k,1) = (j+0.5_rt) * dy;
                    loc_arr(i,j,k,2) = 0.125_rt * (zp_arr(i,j  ,k  ) + zp_arr(i+1,j  ,k  ) +
                                                   zp_arr(i,j+1,k  ) + zp_arr(i+1,j+1,k  ) +
                                                   zp_arr(i,j  ,k+1) + zp_arr(i+1,j  ,k+1) +
                                                   zp_arr(i,j+1,k+1) + zp_arr(i+1,j+1,k+1) );
                });
            } // mfi
            mf_comp += AMREX_SPACEDIM;
        } // if containerHasElement

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

        MultiFab::Copy(mf_nd[lev],*vec_z_phys_nd[lev],0,2,1,0);
        Real dz = Geom()[lev].CellSizeArray()[2];
        int N = Geom()[lev].Domain().size()[2];

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(mf_nd[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            Array4<Real> mf_arr = mf_nd[lev].array(mfi);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                mf_arr(i,j,k,2) = mf_arr(i,j,k,2) + (N-k) * dz;
            });
        } // mfi

        // Fill the coefficient plotfile nodal z coordinate (k=0 slab) if enabled.
        if (write_coeff_plotfile && (ncomp_mf_2d_rho_coeff > 0)) {
            for (MFIter mfi(mf_nd_coeff[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                Array4<Real> mf_arr = mf_nd_coeff[lev].array(mfi);
                Array4<Real const> zp_arr = vec_z_phys_nd[lev]->const_array(mfi);
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    mf_arr(i,j,k,2) = zp_arr(i,j,k) + (N-k) * dz;
                });
            }
        }

    } // lev

    if ( (plotfile_type == PlotfileType::amrex) ||
         (plotfile_type == PlotfileType::hdf5) )
    {

    std::string plotfilename = Concatenate(plot_file_name, istep[0], file_min_digits);

    if (finest_level == 0)
    {
        if (plotfile_type == PlotfileType::amrex) {
            if (write_coeff_plotfile && (ncomp_mf_2d_rho_coeff > 0)) {
                const std::string coeff_plotfilename = plotfilename + "_hmixcoef";
                amrex::Print() << "Writing horizontal mixing coefficient plotfile "
                               << coeff_plotfilename << "\n";
                Vector<std::string> empty_3d;
                Vector<std::string> empty_2d;
                WriteMultiLevelPlotfileWithBathymetry(coeff_plotfilename, finest_level+1,
                                                      GetVecOfConstPtrs(plotMF_zero),
                                                      GetVecOfConstPtrs(mf_nd_coeff),
                                                      GetVecOfConstPtrs(mf_u),
                                                      GetVecOfConstPtrs(mf_v),
                                                      GetVecOfConstPtrs(mf_w),
                                                      GetVecOfConstPtrs(mf_2d_rho_coeff),
                                                      GetVecOfConstPtrs(mf_2d_u_zero),
                                                      GetVecOfConstPtrs(mf_2d_v_zero),
                                                      empty_3d,
                                                      varnames_2d_coeff, empty_2d, empty_2d,
                                                      geom_coeff,
                                                      t_new[0], istep, refRatio());
                writeJobInfo(coeff_plotfilename);
            }
            amrex::Print() << "Writing plotfile " << plotfilename << "\n";
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

#ifdef REMORA_USE_HDF5
        } else if (plotfile_type == PlotfileType::hdf5) {
            amrex::Print() << "Writing plotfile " << plotfilename+"d01.h5" << "\n";
            WriteMultiLevelPlotfileHDF5(plotfilename, finest_level+1,
                                        GetVecOfConstPtrs(plotMF),
                                        varnames_3d,
                                        Geom(), t_new[0], istep, refRatio());
#endif
        }

    } else { // multilevel
        if (plotfile_type == PlotfileType::amrex) {
            if (write_coeff_plotfile && (ncomp_mf_2d_rho_coeff > 0)) {
                const std::string coeff_plotfilename = plotfilename + "_hmixcoef";
                amrex::Print() << "Writing horizontal mixing coefficient plotfile "
                               << coeff_plotfilename << "\n";
                Vector<std::string> empty_3d;
                Vector<std::string> empty_2d;
                WriteMultiLevelPlotfileWithBathymetry(coeff_plotfilename, finest_level+1,
                                                      GetVecOfConstPtrs(plotMF_zero),
                                                      GetVecOfConstPtrs(mf_nd_coeff),
                                                      GetVecOfConstPtrs(mf_u),
                                                      GetVecOfConstPtrs(mf_v),
                                                      GetVecOfConstPtrs(mf_w),
                                                      GetVecOfConstPtrs(mf_2d_rho_coeff),
                                                      GetVecOfConstPtrs(mf_2d_u_zero),
                                                      GetVecOfConstPtrs(mf_2d_v_zero),
                                                      empty_3d,
                                                      varnames_2d_coeff, empty_2d, empty_2d,
                                                      geom_coeff,
                                                      t_new[0], istep, ref_ratio);
                writeJobInfo(coeff_plotfilename);
            }
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

                WriteMultiLevelPlotfileWithBathymetry(plotfilename, finest_level+1,
                                                      GetVecOfConstPtrs(mf2),
                                                      GetVecOfConstPtrs(mf_nd),
                                                      GetVecOfConstPtrs(mf_u),
                                                      GetVecOfConstPtrs(mf_v),
                                                      GetVecOfConstPtrs(mf_w),
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
    for (int lev = 0; lev <= finest_level; ++lev) {
        mask_arrays_for_write(lev, 0.0_rt, (Real) fill_value);
    }

    }
#ifdef REMORA_USE_NETCDF
    else if (plotfile_type == PlotfileType::netcdf)
    {
        // Currently this is hard-coded to plot only level 0
        AMREX_ASSERT(finest_level == 0);
        int lev = 0;
        plotMF[0].FillBoundary(geom[lev].periodicity());
        WriteNCPlotFile(istep_for_plot,&plotMF[lev]);
    } // end if plotfile_type == netcdf
#endif
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

        auto f = [=]() {
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
            VisMF::AsyncWrite(*mf_nd[level],
                              MultiFabFileFullPrefix(level, plotfilename, levelPrefix, mf_nodal_prefix),
                              true);
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
            VisMF::Write(*mf_nd[level], MultiFabFileFullPrefix(level, plotfilename, levelPrefix, mf_nodal_prefix));
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

    int num_extra_mfs = 1; // for nodal, which is always on
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
        HeaderFile << "3" << "\n";
        HeaderFile << "amrexvec_nu_x" << "\n";
        HeaderFile << "amrexvec_nu_y" << "\n";
        HeaderFile << "amrexvec_nu_z" << "\n";
        std::string mf_nodal_prefix = "Nu_nd";
        for (int level = 0; level <= finest_level; ++level) {
            HeaderFile << MultiFabHeaderPath(level, levelPrefix, mf_nodal_prefix) << '\n';
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

        ParallelFor(makeSlab(gbx1,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
            if (mskr(i,j,0) == 0.0) {  // Explicitly compare to 0.0
                Zt_avg1(i,j,0) = fill_value;
            }
        });
        ParallelFor(gbx1, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            if (mskr(i,j,0) == 0.0) {  // Explicitly compare to 0.0
                temp(i,j,k) = fill_value;
                salt(i,j,k) = fill_value;
            }
        });
        ParallelFor(gbx_coeff, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            if (mskr(i,j,0) == 0.0) {  // Explicitly compare to 0.0
                visc2(i,j,k) = fill_value;
                for (int n = 0; n < NCONS; ++n) {
                    diff2(i,j,k,n) = fill_value;
                }
            }
        });
        ParallelFor(makeSlab(ubx,2,0), 3, [=] AMREX_GPU_DEVICE (int i, int j, int , int n)
        {
            if (msku(i,j,0) == 0.0 && ubar(i,j,0)==fill_where) {  // Explicitly compare to 0.0
                ubar(i,j,0,n) = fill_value;
            }
        });
        ParallelFor(makeSlab(vbx,2,0), 3, [=] AMREX_GPU_DEVICE (int i, int j, int , int n)
        {
            if (mskv(i,j,0) == 0.0 && vbar(i,j,0)==fill_where) {  // Explicitly compare to 0.0
                vbar(i,j,0,n) = fill_value;
            }
        });
        ParallelFor(ubx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            if (msku(i,j,0) == 0.0 && xvel(i,j,k)==fill_where) {  // Explicitly compare to 0.0
                xvel(i,j,k) = fill_value;
            }
        });
        ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            if (mskv(i,j,0) == 0.0 && yvel(i,j,k)==fill_where) {  // Explicitly compare to 0.0
                yvel(i,j,k) = fill_value;
            }
        });
    } // mfi
    Gpu::streamSynchronize();
}
