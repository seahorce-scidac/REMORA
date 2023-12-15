#include <EOS.H>
#include <ROMSX.H>
#include "AMReX_Interp_3D_C.H"
#include "AMReX_PlotFileUtil.H"

using namespace amrex;

void
ROMSX::setPlotVariables (const std::string& pp_plot_var_names, Vector<std::string>& plot_var_names)
{
    ParmParse pp(pp_prefix);

    if (pp.contains(pp_plot_var_names.c_str()))
    {
        std::string nm;

        int nPltVars = pp.countval(pp_plot_var_names.c_str());

        for (int i = 0; i < nPltVars; i++)
        {
            pp.get(pp_plot_var_names.c_str(), nm, i);

            // Add the named variable to our list of plot variables
            // if it is not already in the list
            if (!containerHasElement(plot_var_names, nm)) {
                plot_var_names.push_back(nm);
            }
        }
    } else {
        //
        // The default is to add none of the variables to the list
        //
        plot_var_names.clear();
    }

    // Get state variables in the same order as we define them,
    // since they may be in any order in the input list
    Vector<std::string> tmp_plot_names;

    for (int i = 0; i < NCONS; ++i) {
        if ( containerHasElement(plot_var_names, cons_names[i]) ) {
            tmp_plot_names.push_back(cons_names[i]);
        }
    }
    // Check for velocity since it's not in cons_names
    // If we are asked for any velocity component, we will need them all
    if (containerHasElement(plot_var_names, "x_velocity") ||
        containerHasElement(plot_var_names, "y_velocity") ||
        containerHasElement(plot_var_names, "z_velocity")) {
        tmp_plot_names.push_back("x_velocity");
        tmp_plot_names.push_back("y_velocity");
        tmp_plot_names.push_back("z_velocity");
    }

    // If we are asked for any location component, we will provide them all
    if (containerHasElement(plot_var_names, "x_cc") ||
        containerHasElement(plot_var_names, "y_cc") ||
        containerHasElement(plot_var_names, "z_cc")) {
        tmp_plot_names.push_back("x_cc");
        tmp_plot_names.push_back("y_cc");
        tmp_plot_names.push_back("z_cc");
    }

    for (int i = 0; i < derived_names.size(); ++i) {
        if ( containerHasElement(plot_var_names, derived_names[i]) ) {
#ifdef ROMSX_USE_PARTICLES
            if (particleData.use_tracer_particles || (derived_names[i] != "tracer_particle_count")) {
#endif
               tmp_plot_names.push_back(derived_names[i]);
#ifdef ROMSX_USE_PARTICLES
            }
#endif
        } // if
    } // i

    // Check to see if we found all the requested variables
    for (auto plot_name : plot_var_names) {
      if (!containerHasElement(tmp_plot_names, plot_name)) {
           Warning("\nWARNING: Requested to plot variable '" + plot_name + "' but it is not available");
      }
    }
    plot_var_names = tmp_plot_names;
}

// set plotfile variable names
Vector<std::string>
ROMSX::PlotFileVarNames ( Vector<std::string> plot_var_names ) const
{
    Vector<std::string> names;

    names.insert(names.end(), plot_var_names.begin(), plot_var_names.end());

    return names;

}

// Write plotfile to disk
void
ROMSX::WritePlotFile (int which, Vector<std::string> plot_var_names)
{
    const Vector<std::string> varnames = PlotFileVarNames(plot_var_names);
    const int ncomp_mf = varnames.size();
    const auto ngrow_vars = IntVect(NGROW-1,NGROW-1,0);

    // We fillpatch here because some of the derived quantities require derivatives
    //     which require ghost cells to be filled
    for (int lev = 0; lev <= finest_level; ++lev) {
        FillPatch(lev, t_new[lev], *cons_new[lev], cons_new);
        FillPatch(lev, t_new[lev], *xvel_new[lev], xvel_new);
        FillPatch(lev, t_new[lev], *yvel_new[lev], yvel_new);
        FillPatch(lev, t_new[lev], *zvel_new[lev], zvel_new);
    }

    if (ncomp_mf == 0)
        return;

    Vector<MultiFab> mf(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        mf[lev].define(grids[lev], dmap[lev], ncomp_mf, ngrow_vars);
    }

    Vector<MultiFab> mf_nd(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        BoxArray nodal_grids(grids[lev]); nodal_grids.surroundingNodes();
        mf_nd[lev].define(nodal_grids, dmap[lev], AMREX_SPACEDIM, 0);
        mf_nd[lev].setVal(0.);
    }

    for (int lev = 0; lev <= finest_level; ++lev) {
        int mf_comp = 0;

        // First, copy any of the conserved state variables into the output plotfile
        AMREX_ALWAYS_ASSERT(cons_names.size() == NCONS);
        for (int i = 0; i < NCONS; ++i) {
            if (containerHasElement(plot_var_names, cons_names[i])) {
              MultiFab::Copy(mf[lev],*cons_new[lev],i,mf_comp,1,ngrow_vars);
                mf_comp++;
            }
        }

        // Next, check for velocities and if desired, output them
        if (containerHasElement(plot_var_names, "x_velocity") ||
            containerHasElement(plot_var_names, "y_velocity") ||
            containerHasElement(plot_var_names, "z_velocity")) {
            amrex::Print()<<"For now, print faces as if they are at cell centers"<<std::endl;
            //            average_face_to_cellcenter(mf[lev],mf_comp,
            //                Array<const MultiFab*,3>{&xvel_new[lev],&yvel_new[lev],&zvel_new[lev]});
            //
            // Convert the map-factor-scaled-velocities back to velocities
            //
            MultiFab dmf(mf[lev], make_alias, mf_comp, AMREX_SPACEDIM);
            for (MFIter mfi(dmf, TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.growntilebox();
                const Array4<Real> vel_arr = dmf.array(mfi);
                const Array4<const Real> velx_arr = xvel_new[lev]->const_array(mfi);
                const Array4<const Real> vely_arr = yvel_new[lev]->const_array(mfi);
                const Array4<const Real> velz_arr = zvel_new[lev]->const_array(mfi);
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    vel_arr(i,j,k,0) = velx_arr(i,j,k);
                    vel_arr(i,j,k,1) = vely_arr(i,j,k);
                    vel_arr(i,j,k,2) = velz_arr(i,j,k);
                });
            }
            mf_comp += AMREX_SPACEDIM;
        }

        // Fill cell-centered location
        Real dx = Geom()[lev].CellSizeArray()[0];
        Real dy = Geom()[lev].CellSizeArray()[1];

        // Next, check for location names -- if we write one we write all
        if (containerHasElement(plot_var_names, "x_cc") ||
            containerHasElement(plot_var_names, "y_cc") ||
            containerHasElement(plot_var_names, "z_cc"))
        {
            MultiFab dmf(mf[lev], make_alias, mf_comp, AMREX_SPACEDIM);
            for (MFIter mfi(dmf, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.tilebox();
                const Array4<Real> loc_arr = dmf.array(mfi);
                const Array4<Real const> zp_arr = vec_z_phys_nd[lev]->const_array(mfi);

                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    loc_arr(i,j,k,0) = (i+0.5) * dx;
                    loc_arr(i,j,k,1) = (j+0.5) * dy;
                    loc_arr(i,j,k,2) = Real(0.125) * (zp_arr(i,j  ,k  ) + zp_arr(i+1,j  ,k  ) +
                                                      zp_arr(i,j+1,k  ) + zp_arr(i+1,j+1,k  ) +
                                                      zp_arr(i,j  ,k+1) + zp_arr(i+1,j  ,k+1) +
                                                      zp_arr(i,j+1,k+1) + zp_arr(i+1,j+1,k+1) );
                });
            } // mfi
            mf_comp += AMREX_SPACEDIM;
        } // if containerHasElement

#ifdef ROMSX_USE_PARTICLES
        if (containerHasElement(plot_var_names, "tracer_particle_count"))
        {
            MultiFab temp_dat(mf[lev].boxArray(), mf[lev].DistributionMap(), 1, 0);
            temp_dat.setVal(0);
            particleData.tracer_particles->Increment(temp_dat, lev);
            MultiFab::Copy(mf[lev], temp_dat, 0, mf_comp, 1, 0);
            mf_comp += 1;
        }
#endif

        MultiFab::Copy(mf_nd[lev],*vec_z_phys_nd[lev],0,2,1,0);
        Real dz = Geom()[lev].CellSizeArray()[2];
        for (MFIter mfi(mf_nd[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            Array4<Real> mf_arr = mf_nd[lev].array(mfi);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                mf_arr(i,j,k,2) -= k * dz;
            });
        } // mfi

    } // lev

    std::string plotfilename;
    if (which == 1)
       plotfilename = Concatenate(plot_file_1, istep[0], 5);
    else if (which == 2)
       plotfilename = Concatenate(plot_file_2, istep[0], 5);

    if (finest_level == 0)
    {
        if (plotfile_type == "amrex") {
            amrex::Print() << "Writing plotfile " << plotfilename << "\n";
            WriteMultiLevelPlotfileWithBathymetry(plotfilename, finest_level+1,
                                                  GetVecOfConstPtrs(mf),
                                                  GetVecOfConstPtrs(mf_nd),
                                                  varnames,
                                                  t_new[0], istep);
            writeJobInfo(plotfilename);

#ifdef ROMSX_USE_PARTICLES
            particleData.Checkpoint(plotfilename);
#endif

#ifdef ROMSX_USE_HDF5
        } else if (plotfile_type == "hdf5" || plotfile_type == "HDF5") {
            amrex::Print() << "Writing plotfile " << plotfilename+"d01.h5" << "\n";
            WriteMultiLevelPlotfileHDF5(plotfilename, finest_level+1,
                                        GetVecOfConstPtrs(mf),
                                        varnames,
                                        Geom(), t_new[0], istep, refRatio());
#endif
#ifdef ROMSX_USE_NETCDF
        } else if (plotfile_type == "netcdf" || plotfile_type == "NetCDF") {
             int lev   = 0;
             int which = 0;
             writeNCPlotFile(lev, which, plotfilename, GetVecOfConstPtrs(mf), varnames, istep, t_new[0]);
             total_plot_file_step_1 += 1;
#endif
        } else {
            amrex::Print() << "User specified plot_filetype = " << plotfile_type << std::endl;
            amrex::Abort("Dont know this plot_filetype");
        }

    } else { // multilevel

        Vector<IntVect>   r2(finest_level);
        Vector<Geometry>  g2(finest_level+1);
        Vector<MultiFab> mf2(finest_level+1);

        mf2[0].define(grids[0], dmap[0], ncomp_mf, 0);

        // Copy level 0 as is
        MultiFab::Copy(mf2[0],mf[0],0,0,mf[0].nComp(),0);

        // Define a new multi-level array of Geometry's so that we pass the new "domain" at lev > 0
        Array<int,AMREX_SPACEDIM> periodicity =
                     {Geom()[0].isPeriodic(0),Geom()[0].isPeriodic(1),Geom()[0].isPeriodic(2)};
        g2[0].define(Geom()[0].Domain(),&(Geom()[0].ProbDomain()),0,periodicity.data());

        if (plotfile_type == "amrex") {
            r2[0] = IntVect(1,1,ref_ratio[0][0]);
            for (int lev = 1; lev <= finest_level; ++lev) {
                if (lev > 1) {
                    r2[lev-1][0] = 1;
                    r2[lev-1][1] = 1;
                    r2[lev-1][2] = r2[lev-2][2] * ref_ratio[lev-1][0];
                }

                mf2[lev].define(refine(grids[lev],r2[lev-1]), dmap[lev], ncomp_mf, 0);

                // Set the new problem domain
                Box d2(Geom()[lev].Domain());
                d2.refine(r2[lev-1]);

                g2[lev].define(d2,&(Geom()[lev].ProbDomain()),0,periodicity.data());
            }

            // Do piecewise interpolation of mf into mf2
            for (int lev = 1; lev <= finest_level; ++lev) {
                for (MFIter mfi(mf2[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                    const Box& bx = mfi.tilebox();
                    pcinterp_interp(bx,mf2[lev].array(mfi), 0, mf[lev].nComp(), mf[lev].const_array(mfi),0,r2[lev-1]);
                }
            }

            // Define an effective ref_ratio which is isotropic to be passed into WriteMultiLevelPlotfile
            Vector<IntVect> rr(finest_level);
            for (int lev = 0; lev < finest_level; ++lev) {
                rr[lev] = IntVect(ref_ratio[lev][0],ref_ratio[lev][1],ref_ratio[lev][0]);
            }

            WriteMultiLevelPlotfile(plotfilename, finest_level+1, GetVecOfConstPtrs(mf2), varnames,
                                    g2, t_new[0], istep, rr);
            writeJobInfo(plotfilename);

#ifdef ROMSX_USE_PARTICLES
            particleData.Checkpoint(plotfilename);
#endif

#ifdef ROMSX_USE_NETCDF
        } else if (plotfile_type == "netcdf" || plotfile_type == "NetCDF") {
             for (int lev = 0; lev <= finest_level; ++lev) {
                 for (int which = 0; which < num_boxes_at_level[lev]; which++) {
                     writeNCPlotFile(lev, which, plotfilename, GetVecOfConstPtrs(mf), varnames, istep, t_new[0]);
                     total_plot_file_step_1 += 1;
                 }
             }
#endif
        }
    } // end multi-level
}

void
ROMSX::WriteMultiLevelPlotfileWithBathymetry (const std::string& plotfilename, int nlevels,
                                              const Vector<const MultiFab*>& mf,
                                              const Vector<const MultiFab*>& mf_nd,
                                              const Vector<std::string>& varnames,
                                              Real time,
                                              const Vector<int>& level_steps,
                                              const std::string &versionName,
                                              const std::string &levelPrefix,
                                              const std::string &mfPrefix,
                                              const Vector<std::string>& extra_dirs) const
{
    BL_PROFILE("WriteMultiLevelPlotfileWithBathymetry()");

    BL_ASSERT(nlevels <= mf.size());
    BL_ASSERT(nlevels <= ref_ratio.size()+1);
    BL_ASSERT(nlevels <= level_steps.size());
    BL_ASSERT(mf[0]->nComp() == varnames.size());

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
            WriteGenericPlotfileHeaderWithBathymetry(HeaderFile, nlevels, boxArrays, varnames,
                                                     time, level_steps, versionName,
                                                     levelPrefix, mfPrefix);
        };

        if (AsyncOut::UseAsyncOut()) {
            AsyncOut::Submit(std::move(f));
        } else {
            f();
        }
    }

    std::string mf_nodal_prefix = "Nu_nd";
    for (int level = 0; level <= finest_level; ++level)
    {
        if (AsyncOut::UseAsyncOut()) {
            VisMF::AsyncWrite(*mf[level],
                              MultiFabFileFullPrefix(level, plotfilename, levelPrefix, mfPrefix),
                              true);
            VisMF::AsyncWrite(*mf_nd[level],
                              MultiFabFileFullPrefix(level, plotfilename, levelPrefix, mf_nodal_prefix),
                              true);
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
        }
    }
}

void
ROMSX::WriteGenericPlotfileHeaderWithBathymetry (std::ostream &HeaderFile,
                                                 int nlevels,
                                                 const Vector<BoxArray> &bArray,
                                                 const Vector<std::string> &varnames,
                                                 Real time,
                                                 const Vector<int> &level_steps,
                                                 const std::string &versionName,
                                                 const std::string &levelPrefix,
                                                 const std::string &mfPrefix) const
{
        BL_ASSERT(nlevels <= bArray.size());
        BL_ASSERT(nlevels <= ref_ratio.size()+1);
        BL_ASSERT(nlevels <= level_steps.size());

        HeaderFile.precision(17);

        // ---- this is the generic plot file type name
        HeaderFile << versionName << '\n';

        HeaderFile << varnames.size() << '\n';

        for (int ivar = 0; ivar < varnames.size(); ++ivar) {
            HeaderFile << varnames[ivar] << "\n";
        }
        HeaderFile << AMREX_SPACEDIM << '\n';
        HeaderFile << time << '\n';
        HeaderFile << finest_level << '\n';
        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            HeaderFile << geom[0].ProbLo(i) << ' ';
        }
        HeaderFile << '\n';
        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            HeaderFile << geom[0].ProbHi(i) << ' ';
        }
        HeaderFile << '\n';
        for (int i = 0; i < finest_level; ++i) {
            HeaderFile << ref_ratio[i][0] << ' ';
        }
        HeaderFile << '\n';
        for (int i = 0; i <= finest_level; ++i) {
            HeaderFile << geom[i].Domain() << ' ';
        }
        HeaderFile << '\n';
        for (int i = 0; i <= finest_level; ++i) {
            HeaderFile << level_steps[i] << ' ';
        }
        HeaderFile << '\n';
        for (int i = 0; i <= finest_level; ++i) {
            for (int k = 0; k < AMREX_SPACEDIM; ++k) {
                HeaderFile << geom[i].CellSize()[k] << ' ';
            }
            HeaderFile << '\n';
        }
        HeaderFile << (int) geom[0].Coord() << '\n';
        HeaderFile << "0\n";

        for (int level = 0; level <= finest_level; ++level) {
            HeaderFile << level << ' ' << bArray[level].size() << ' ' << time << '\n';
            HeaderFile << level_steps[level] << '\n';

            const IntVect& domain_lo = geom[level].Domain().smallEnd();
            for (int i = 0; i < bArray[level].size(); ++i)
            {
                // Need to shift because the RealBox ctor we call takes the
                // physical location of index (0,0,0).  This does not affect
                // the usual cases where the domain index starts with 0.
                const Box& b = shift(bArray[level][i], -domain_lo);
                RealBox loc = RealBox(b, geom[level].CellSize(), geom[level].ProbLo());
                for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                    HeaderFile << loc.lo(n) << ' ' << loc.hi(n) << '\n';
                }
            }

            HeaderFile << MultiFabHeaderPath(level, levelPrefix, mfPrefix) << '\n';
        }
        HeaderFile << "1" << "\n";
        HeaderFile << "3" << "\n";
        HeaderFile << "amrexvec_nu_x" << "\n";
        HeaderFile << "amrexvec_nu_y" << "\n";
        HeaderFile << "amrexvec_nu_z" << "\n";
        std::string mf_nodal_prefix = "Nu_nd";
        for (int level = 0; level <= finest_level; ++level) {
            HeaderFile << MultiFabHeaderPath(level, levelPrefix, mf_nodal_prefix) << '\n';
        }
}
