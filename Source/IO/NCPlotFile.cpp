#include <iomanip>
#include <iostream>
#include <string>
#include <ctime>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <AMReX_Utility.H>
#include <AMReX_buildInfo.H>
#include <AMReX_ParmParse.H>

#include "REMORA.H"
#include "NCInterface.H"
#include "NCPlotFile.H"
#include "IndexDefines.H"

using namespace amrex;

void
REMORA::WriteNCPlotFile(int which_step) const
{
    // For right now we assume single level -- we will generalize this later to multilevel
    int lev             = 0;
    int which_subdomain = 0;

    // Create filename
    std::string plt_string;
    std::string plotfilename;
    if (REMORA::write_history_file) {
        plotfilename = "plt_his";
    } else {
        plotfilename = Concatenate("plt", which_step, 5);
    }

    // Set the full IO path for NetCDF output
    std::string FullPath = plotfilename;
    if (lev == 0) {
        const std::string& extension = amrex::Concatenate("_d",lev+1,2);
        FullPath += extension + ".nc";
    } else {
        const std::string& extension = amrex::Concatenate("_d",lev+1+which_subdomain,2);
        FullPath += extension + ".nc";
    }

    bool not_empty_file = false;

    if (REMORA::write_history_file)
    {
        not_empty_file = amrex::FileSystem::Exists(FullPath);

        auto ncf = not_empty_file ?
        ncutils::NCFile::open_par(FullPath, NC_WRITE | NC_NETCDF4 | NC_MPIIO,
                                  amrex::ParallelContext::CommunicatorSub(), MPI_INFO_NULL) :
        ncutils::NCFile::create_par(FullPath, NC_NETCDF4 | NC_MPIIO,
                                    amrex::ParallelContext::CommunicatorSub(), MPI_INFO_NULL);

        amrex::Print() << "Writing level " << lev << " NetCDF history file " << FullPath
                       << "\nFor step "<< which_step <<std::endl;

        WriteNCPlotFile_which(which_step, lev, not_empty_file, which_subdomain, ncf);

     } else {

        // Open new netcdf file to write data
        auto ncf = ncutils::NCFile::create_par(FullPath, NC_NETCDF4 | NC_MPIIO,
                                               amrex::ParallelContext::CommunicatorSub(), MPI_INFO_NULL);
        amrex::Print() << "Writing level " << lev << " NetCDF plot file " << FullPath
                       << "\nFor step "<< which_step <<std::endl;

        WriteNCPlotFile_which(which_step, lev, not_empty_file, which_subdomain, ncf);
     }
}

void
REMORA::WriteNCPlotFile_which(int which_step, int lev, int which_subdomain,
                              bool not_empty, ncutils::NCFile& ncf) const
{
    int iproc = amrex::ParallelContext::MyProcAll();
    int nproc = amrex::ParallelDescriptor::NProcs();

    // Number of cells in this "domain" at this level
    std::vector<int> n_cells;

    int nblocks = grids[lev].size();
    // auto dm = mf[lev].DistributionMap();

    // Number of points in each block at this level
    std::vector<int> offset_s;
    std::vector<int> offset_u;
    std::vector<int> offset_v;
    offset_s.reserve(nproc);
    offset_u.reserve(nproc);
    offset_v.reserve(nproc);

    for (auto n = 0; n < nproc; n++) {
        offset_s[n] = 0;
        offset_u[n] = 0;
        offset_v[n] = 0;
    }

    for (auto ib=0; ib<nblocks; ib++) {
       Box tmp_bx_s(grids[lev][ib]); tmp_bx_s.grow(IntVect(1,1,0));
       offset_s[dmap[lev][ib]] += tmp_bx_s.numPts();

       Box tmp_bx_u(grids[lev][ib]); tmp_bx_u.surroundingNodes(0); tmp_bx_u.grow(IntVect(0,1,0));
       offset_u[dmap[lev][ib]] += tmp_bx_u.numPts();

       Box tmp_bx_v(grids[lev][ib]); tmp_bx_v.surroundingNodes(1); tmp_bx_v.grow(IntVect(1,0,0));
       offset_v[dmap[lev][ib]] += tmp_bx_v.numPts();
    }

    // We only do single-level writes when using NetCDF format
    int flev = lev;

     Box subdomain;
     if (lev == 0) {
         subdomain = geom[lev].Domain();
     } else {
         subdomain = boxes_at_level[lev][which_subdomain];
     }

     int nx = subdomain.length(0);
     int ny = subdomain.length(1);
     int nz = subdomain.length(2);

     n_cells.push_back(nx);
     n_cells.push_back(ny);
     n_cells.push_back(nz);

     int num_pts   = nx*ny*nz;
     int num_s_pts = (nx+2)*(ny+2)*nz;
     int num_u_pts = (nx+1)*(ny+2)*nz;
     int num_v_pts = (nx+2)*(ny+1)*nz;

     const std::string nt_name   = "num_time_steps";
     const std::string ndim_name = "num_geo_dimensions";

     const std::string np_name   = "num_points_per_block";
     const std::string np_s_name = "num_points_per_block_s";
     const std::string np_u_name = "num_points_per_block_u";
     const std::string np_v_name = "num_points_per_block_v";

     const std::string nb_name   = "num_blocks";
     const std::string nx_name   = "NX";
     const std::string ny_name   = "NY";
     const std::string nz_name   = "NZ";
     const std::string flev_name = "FINEST_LEVEL";

     amrex::Print() << "NOT EMPTY" << not_empty << std::endl;
     // if (!not_empty)
     {
         ncf.enter_def_mode();
         ncf.put_attr("title", "REMORA data ");
         ncf.def_dim(nt_name,   NC_UNLIMITED);
         ncf.def_dim(ndim_name, AMREX_SPACEDIM);

         ncf.def_dim(np_name  ,   num_pts);
         ncf.def_dim(np_s_name,   num_s_pts);
         ncf.def_dim(np_u_name,   num_u_pts);
         ncf.def_dim(np_v_name,   num_v_pts);

         ncf.def_dim(nb_name,   nblocks);
         ncf.def_dim(flev_name, flev);

         ncf.def_dim(nx_name,   n_cells[0]);
         ncf.def_dim(ny_name,   n_cells[1]);
         ncf.def_dim(nz_name,   n_cells[2]);

         ncf.def_var("probLo"  ,   NC_FLOAT,  {ndim_name});
         ncf.def_var("probHi"  ,   NC_FLOAT,  {ndim_name});

         ncf.def_var("Geom.smallend", NC_INT, {flev_name, ndim_name});
         ncf.def_var("Geom.bigend"  , NC_INT, {flev_name, ndim_name});
         ncf.def_var("CellSize"     , NC_FLOAT, {flev_name, ndim_name});

         ncf.def_var("x_grid", NC_FLOAT, {np_name});
         ncf.def_var("y_grid", NC_FLOAT, {np_name});
         ncf.def_var("z_grid", NC_FLOAT, {np_name});
         ncf.exit_def_mode();
     }

     ncf.enter_def_mode();
     if (REMORA::write_history_file) {
         ncf.def_var("temp", NC_FLOAT, {nt_name,np_s_name});
         ncf.def_var("salt", NC_FLOAT, {nt_name, np_s_name});
         ncf.def_var("u"   , NC_FLOAT, {nt_name, np_u_name});
         ncf.def_var("v"   , NC_FLOAT, {nt_name, np_v_name});
     } else {
         ncf.def_var("temp", NC_FLOAT, {np_s_name});
         ncf.def_var("salt", NC_FLOAT, {np_s_name});
         ncf.def_var("u"   , NC_FLOAT, {np_u_name});
         ncf.def_var("v"   , NC_FLOAT, {np_v_name});
     }
     ncf.exit_def_mode();

     {
     // We are doing single-level writes but it doesn't have to be level 0
     //
     // Write out the header information.
     //

     Real time = 0.;

     // Right now this is hard-wired to {temp, salt, u, v}
     int n_data_items = 4;

     ncf.put_attr("number_variables", std::vector<int>{n_data_items});
     ncf.put_attr("space_dimension", std::vector<int>{AMREX_SPACEDIM});
     ncf.put_attr("current_time", std::vector<double>{time});
     ncf.put_attr("start_time", std::vector<double>{start_bdy_time});
     ncf.put_attr("CurrentLevel", std::vector<int>{flev});

      Real dx[AMREX_SPACEDIM];
      for (int i = 0; i < AMREX_SPACEDIM; i++)
         dx[i] = geom[lev].CellSize()[i];
      const auto *base = geom[lev].ProbLo();
      RealBox rb(subdomain,dx,base);

      amrex::Vector<Real> probLo;
      amrex::Vector<Real> probHi;
      for (int i = 0; i < AMREX_SPACEDIM; i++) {
        probLo.push_back(rb.lo(i));
        probHi.push_back(rb.hi(i));
      }

      auto nc_probLo = ncf.var("probLo");
      nc_probLo.par_access(NC_COLLECTIVE);
      nc_probLo.put(probLo.data(), {0}, {AMREX_SPACEDIM});

      auto nc_probHi = ncf.var("probHi");
      nc_probHi.par_access(NC_COLLECTIVE);
      nc_probHi.put(probHi.data(), {0}, {AMREX_SPACEDIM});

      amrex::Vector<int> smallend;
      amrex::Vector<int> bigend;
      for (int i = lev; i <= flev; i++) {
        smallend.clear(); bigend.clear();
        for (int j = 0; j < AMREX_SPACEDIM; j++) {
           smallend.push_back(subdomain.smallEnd(j));
             bigend.push_back(subdomain.bigEnd(j));
        }
        auto nc_Geom_smallend = ncf.var("Geom.smallend");
        nc_Geom_smallend.par_access(NC_COLLECTIVE);
        nc_Geom_smallend.put(smallend.data(), {static_cast<long unsigned int>(i-lev), 0}, {1, AMREX_SPACEDIM});

        auto nc_Geom_bigend = ncf.var("Geom.bigend");
        nc_Geom_bigend.par_access(NC_COLLECTIVE);
        nc_Geom_bigend.put(bigend.data(), {static_cast<long unsigned int>(i-lev), 0}, {1, AMREX_SPACEDIM});
      }

      amrex::Vector<Real> CellSize;
      for (int i = lev; i <= flev; i++) {
        CellSize.clear();
        for (double & j : dx) {
          CellSize.push_back(j);
        }
        auto nc_CellSize = ncf.var("CellSize");
        nc_CellSize.par_access(NC_COLLECTIVE);
        nc_CellSize.put(CellSize.data(), {static_cast<long unsigned int>(i-lev), 0}, {1, AMREX_SPACEDIM});
      }

      ncf.put_attr("DefaultGeometry", std::vector<int>{amrex::DefaultGeometry().Coord()});
    }

    std::vector<Real> x_grid;
    std::vector<Real> y_grid;
    std::vector<Real> z_grid;
    long unsigned goffset = 0;
    long unsigned glen    = 0;
    for (int i = 0; i < grids[lev].size(); ++i) {
        auto box = grids[lev][i];
        if (subdomain.contains(box)) {
            RealBox gridloc = RealBox(grids[lev][i], geom[lev].CellSize(), geom[lev].ProbLo());

            x_grid.clear(); y_grid.clear(); z_grid.clear();
            for (auto k1 = 0; k1 < grids[lev][i].length(0); ++k1) {
              for (auto k2 = 0; k2 < grids[lev][i].length(1); ++k2) {
                 for (auto k3 = 0; k3 < grids[lev][i].length(2); ++k3) {
                    x_grid.push_back(gridloc.lo(0)+geom[lev].CellSize(0)*static_cast<Real>(k1));
                    y_grid.push_back(gridloc.lo(1)+geom[lev].CellSize(1)*static_cast<Real>(k2));
                    z_grid.push_back(gridloc.lo(2)+geom[lev].CellSize(2)*static_cast<Real>(k3));
                 }
              }
            }

            goffset += glen;
            glen = grids[lev][i].length(0)*grids[lev][i].length(1)*grids[lev][i].length(2);

            auto nc_x_grid = ncf.var("x_grid");
            auto nc_y_grid = ncf.var("y_grid");
            auto nc_z_grid = ncf.var("z_grid");

            nc_x_grid.par_access(NC_INDEPENDENT);
            nc_y_grid.par_access(NC_INDEPENDENT);
            nc_z_grid.par_access(NC_INDEPENDENT);

            nc_x_grid.put(x_grid.data(), {goffset}, {glen});
            nc_y_grid.put(y_grid.data(), {goffset}, {glen});
            nc_z_grid.put(z_grid.data(), {goffset}, {glen});
       }
    }

    size_t nbox_per_proc = 0;
    long unsigned numpts = 0;

    long unsigned diff_s = 0;
    long unsigned npts_s = 0;

    for (MFIter mfi(*cons_new[lev],false); mfi.isValid(); ++mfi) {
        auto box = mfi.validbox();
        if (subdomain.contains(box)) {
            long unsigned diff = nbox_per_proc*numpts;
            for(auto ip = 1; ip <= iproc; ++ip) {diff += offset_s[ip-1];}

            Box tmp_bx(box); tmp_bx.grow(IntVect(1,1,0));
            long unsigned numpts = tmp_bx.numPts();

            {
            amrex::Print() << "WRITING TEMP " << std::endl;
            FArrayBox tmp_temp(tmp_bx,1);
            tmp_temp.template copy<RunOn::Device>((*cons_new[lev])[mfi.index()],Temp_comp,0,1);

            auto nc_plot_var = ncf.var("temp");
            nc_plot_var.par_access(NC_INDEPENDENT);
            nc_plot_var.put(tmp_temp.dataPtr(), {diff}, {numpts});
            }

            {
            amrex::Print() << "WRITING SALT " << std::endl;
            FArrayBox tmp_salt(tmp_bx,1);
            tmp_salt.template copy<RunOn::Device>((*cons_new[lev])[mfi.index()],Salt_comp,0,1);

            auto nc_plot_var = ncf.var("salt");
            nc_plot_var.par_access(NC_INDEPENDENT);
            nc_plot_var.put(tmp_salt.dataPtr(), {diff}, {numpts});
            }

            nbox_per_proc++;
        }
    }

    long unsigned diff_u = 0;
    long unsigned npts_u = 0;

    // Writing u (we loop over cons to get cell-centered box)
    for (MFIter mfi(*cons_new[lev],false); mfi.isValid(); ++mfi)
    {
        Box box = mfi.validbox();
        int idx = mfi.index();

        if (subdomain.contains(box))
        {
            diff_u = npts_u;

            for (auto ip = 0; ip < iproc; ++ip) {
                diff_u += offset_u[ip];
            }

            Box tmp_bx(box); tmp_bx.surroundingNodes(0); tmp_bx.grow(IntVect(0,1,0));

            long unsigned numpts = tmp_bx.numPts();

            FArrayBox tmp(tmp_bx,1);
            tmp.template copy<RunOn::Device>((*xvel_new[lev])[idx],0,0,1);

            auto nc_plot_var = ncf.var("u");
            nc_plot_var.par_access(NC_INDEPENDENT);
            nc_plot_var.put(tmp.dataPtr(), {diff_u}, {numpts});

            npts_u += numpts;

        } // in subdomain
    } // mfi

    long unsigned diff_v = 0;
    long unsigned npts_v = 0;

    // Writing v (we loop over cons to get cell-centered box)
    for (MFIter mfi(*cons_new[lev],false); mfi.isValid(); ++mfi)
    {
       Box box = mfi.validbox();
        int idx = mfi.index();

        if (subdomain.contains(box))
        {
            diff_v = npts_v;

            for (auto ip = 0; ip <= iproc-1; ++ip) {
                diff_v += offset_v[ip];
            }

            Box tmp_bx(box); tmp_bx.surroundingNodes(1); tmp_bx.grow(IntVect(1,0,0));

            long unsigned numpts = tmp_bx.numPts();

            FArrayBox tmp(tmp_bx,1);
            tmp.template copy<RunOn::Device>((*yvel_new[lev])[idx],0,0,1);

            auto nc_plot_var = ncf.var("v");
            nc_plot_var.par_access(NC_INDEPENDENT);
            nc_plot_var.put(tmp.dataPtr(), {diff_v}, {numpts});

            npts_v += numpts;

        } // in subdomain
    } // mfi

   ncf.close();
}

