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

#include "ROMSX.H"
#include "NCInterface.H"
#include "NCPlotFile.H"
#include "IndexDefines.H"

namespace {
 std::string nc_state_filename = "Plot_State_MF.nc";
} // namespace

using namespace amrex;

void
ROMSX::writeNCPlotFile(int lev, int which_subdomain, const std::string& dir,
                     const Vector<const MultiFab*> &plotMF,
                     const Vector<std::string> &plot_var_names,
                     const Vector<int> level_steps, const Real time) const
{
     // get the processor number
     int iproc = amrex::ParallelContext::MyProcAll();

     // number of cells in this "domain" at this level
     std::vector<int> n_cells;

     // set the full IO path for NetCDF output
     std::string FullPath = dir;
     if (lev == 0) {
         const std::string& extension = amrex::Concatenate("_d",lev+1,2);
         FullPath += extension + ".nc";
     } else {
         const std::string& extension = amrex::Concatenate("_d",lev+1+which_subdomain,2);
         FullPath += extension + ".nc";
     }
#ifdef ROMSX_USE_HISTORYFILE
     bool not_empty_file = true;
     if (true) {
         FullPath = plot_file_1 + "_his.nc";
         not_empty_file = amrex::FileSystem::Exists(FullPath);
     }
#endif

     amrex::Print() << "Writing level " << lev << " NetCDF plot file " << FullPath << "\nFor step "<<istep[0]<<std::endl;

#ifdef ROMSX_USE_HISTORYFILE
     // open netcdf file to write data
     auto ncf = not_empty_file ?
       ncutils::NCFile::open_par(FullPath, NC_WRITE | NC_NETCDF4 | NC_MPIIO,
                                 amrex::ParallelContext::CommunicatorSub(), MPI_INFO_NULL) :
       ncutils::NCFile::create_par(FullPath, NC_NETCDF4 | NC_MPIIO,
                                   amrex::ParallelContext::CommunicatorSub(), MPI_INFO_NULL);
#else
     // open netcdf file to write data
     auto ncf = ncutils::NCFile::create_par(FullPath, NC_NETCDF4 | NC_MPIIO,
                                   amrex::ParallelContext::CommunicatorSub(), MPI_INFO_NULL);
#endif

     int nblocks = grids[lev].size();

     // We only do single-level writes when using NetCDF format
     int flev    = lev;

     Box subdomain;
     if (lev == 0) {
         subdomain = geom[lev].Domain();
     } else {
         subdomain = boxes_at_level[lev][which_subdomain];
     }

     int nx = subdomain.length(0);
     int ny = subdomain.length(1);
     int nz = subdomain.length(2);
     n_cells.push_back(nx); n_cells.push_back(ny); n_cells.push_back(nz);
     int num_pts = nx * ny * nz;

     int n_data_items = plotMF[lev]->nComp();

     const std::string nt_name   = "num_time_steps";
     const std::string ndim_name = "num_geo_dimensions";
     const std::string np_name   = "num_points_per_block";
     const std::string nb_name   = "num_blocks";
     const std::string nx_name   = "NX";
     const std::string ny_name   = "NY";
     const std::string nz_name   = "NZ";
     const std::string flev_name = "FINEST_LEVEL";
#ifdef ROMSX_USE_HISTORYFILE
     if(!not_empty_file) {
#endif
     ncf.enter_def_mode();
     ncf.put_attr("title", "ROMSX NetCDF Plot data output");
     ncf.def_dim(nt_name,   NC_UNLIMITED);
     ncf.def_dim(ndim_name, AMREX_SPACEDIM);
     ncf.def_dim(np_name,   num_pts);
     ncf.def_dim(nb_name,   nblocks);
     ncf.def_dim(flev_name, flev);

     ncf.def_dim(nx_name,   n_cells[0]);
     ncf.def_dim(ny_name,   n_cells[1]);
     ncf.def_dim(nz_name,   n_cells[2]);
     
     ncf.def_var("probLo"  ,   NC_FLOAT,  {ndim_name});
     ncf.def_var("probHi"  ,   NC_FLOAT,  {ndim_name});
     ncf.def_var("refRatio",   NC_INT,    {flev_name});
     ncf.def_var("levelSteps", NC_INT,    {flev_name});

     ncf.def_var("Geom.smallend", NC_INT, {flev_name, ndim_name});
     ncf.def_var("Geom.bigend"  , NC_INT, {flev_name, ndim_name});
     ncf.def_var("CellSize"     , NC_FLOAT, {flev_name, ndim_name});

     ncf.def_var("x_grid", NC_FLOAT, {nb_name, nx_name});
     ncf.def_var("y_grid", NC_FLOAT, {nb_name, ny_name});
     ncf.def_var("z_grid", NC_FLOAT, {nb_name, nz_name});

     ncf.def_var("z_r_grid", NC_FLOAT, {nb_name, nz_name});
     ncf.def_var("z_w_grid", NC_FLOAT, {nb_name, nz_name});
#ifdef ROMSX_USE_HISTORYFILE
     for (int i = 0; i < plot_var_names.size(); i++) {
       ncf.def_var(plot_var_names[i], NC_FLOAT, {nt_name, nz_name, ny_name, nx_name});
     }
     /*
     ncf.def_var("z_r", NC_FLOAT, {nt_name, nz_name, ny_name, nx_name});
     ncf.def_var("z_w", NC_FLOAT, {nt_name, nz_name, ny_name, nx_name});
     */
#else
     for (int i = 0; i < plot_var_names.size(); i++) {
       ncf.def_var(plot_var_names[i], NC_FLOAT, {nz_name, ny_name, nx_name});
     }
     /*
     ncf.def_var("z_r", NC_FLOAT, {nz_name, ny_name, nx_name});
     ncf.def_var("z_w", NC_FLOAT, {nz_name, ny_name, nx_name});
     */
#endif

     ncf.exit_def_mode();
#ifdef ROMSX_USE_HISTORYFILE
     }
#endif
     {
      // We are doing single-level writes but it doesn't have to be level 0
      //
      // Write out the netcdf plotfile head information.
      //
      if (n_data_items == 0)
        amrex::Error("Must specify at least one valid data item to plot");

      ncf.put_attr("number_variables", std::vector<int>{n_data_items});
      ncf.put_attr("space_dimension", std::vector<int>{AMREX_SPACEDIM});
      ncf.put_attr("current_time", std::vector<double>{time});
      ncf.put_attr("FinestLevel", std::vector<int>{finest_level});
      ncf.put_attr("CurrentLevel", std::vector<int>{lev});

      Real dx[AMREX_SPACEDIM];
      for (int i = 0; i < AMREX_SPACEDIM; i++)
         dx[i] = geom[lev].CellSize()[i];
      auto base = geom[lev].ProbLo();
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

      amrex::Vector<int> refRatio;
      if (lev < finest_level)
      refRatio.push_back(ref_ratio[lev][0]);
      auto nc_refRatio = ncf.var("refRatio");
      nc_refRatio.par_access(NC_COLLECTIVE);
      nc_refRatio.put(refRatio.data(), {0}, {static_cast<long unsigned int>(flev-lev)});

      amrex::Vector<int> levelSteps;
      for (int i = lev; i <= flev; i++)
        levelSteps.push_back(level_steps[i]);
      auto nc_levelSteps = ncf.var("levelSteps");
      nc_levelSteps.par_access(NC_COLLECTIVE);
      nc_levelSteps.put(levelSteps.data(), {0}, {static_cast<long unsigned int>(flev-lev)});

      amrex::Vector<int> smallend;
      amrex::Vector<int> bigend;
      for (int i = lev; i <= flev; i++)
      {
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
        for (int j = 0; j < AMREX_SPACEDIM; j++) {
          CellSize.push_back(dx[j]);
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
    std::vector<Real> z_r_grid;
    std::vector<Real> z_w_grid;
    for (int i = 0; i < grids[lev].size(); ++i) {
        auto box = grids[lev][i];
        if (subdomain.contains(box)) {
            RealBox gridloc = RealBox(grids[lev][i], geom[lev].CellSize(), geom[lev].ProbLo());

            x_grid.clear(); y_grid.clear(); z_grid.clear();
            for (auto k1 = 0; k1 < grids[lev][i].length(0); ++k1) {
              x_grid.push_back(gridloc.lo(0)+geom[lev].CellSize(0)*static_cast<Real>(k1));
            }
            for (auto k2 = 0; k2 < grids[lev][i].length(1); ++k2) {
              y_grid.push_back(gridloc.lo(1)+geom[lev].CellSize(1)*static_cast<Real>(k2));
            }
            for (auto k3 = 0; k3 < grids[lev][i].length(2); ++k3) {
              z_grid.push_back(gridloc.lo(2)+geom[lev].CellSize(2)*static_cast<Real>(k3));
            }

            auto xlen = static_cast<long unsigned int>(grids[lev][i].length(0));
            auto ylen = static_cast<long unsigned int>(grids[lev][i].length(1));
            auto zlen = static_cast<long unsigned int>(grids[lev][i].length(2));
            auto nc_x_grid = ncf.var("x_grid");
            nc_x_grid.par_access(NC_COLLECTIVE);
            nc_x_grid.put(x_grid.data(), {static_cast<long unsigned int>(i), 0}, {1, xlen});
            auto nc_y_grid = ncf.var("y_grid");
            nc_y_grid.par_access(NC_COLLECTIVE);
            nc_y_grid.put(y_grid.data(), {static_cast<long unsigned int>(i), 0}, {1, ylen});
            auto nc_z_grid = ncf.var("z_grid");
            nc_z_grid.par_access(NC_COLLECTIVE);
            nc_z_grid.put(z_grid.data(), {static_cast<long unsigned int>(i), 0}, {1, zlen});
       }
   }

   MultiFab* mf_z_w = z_w[lev].get();
   MultiFab* mf_z_r = z_r[lev].get();

   for ( MFIter mfi(*mf_z_w, TilingIfNotGPU()); mfi.isValid(); ++mfi )
   {
      auto box = mfi.validbox();
      auto zlen = static_cast<long unsigned int>(box.length(2));
      Array4<Real> const& z_w_arr = (mf_z_w)->array(mfi);
      Array4<Real> const& z_r_arr = (mf_z_r)->array(mfi);
      z_r_grid.clear(); z_w_grid.clear();
      std::vector<size_t> startp(2);
      std::vector<size_t> countp(2);
      startp[0]=0;
      countp[0]=1;
      startp[1]=box.smallEnd(2);
      countp[1]=0;

      for (auto k4 = box.smallEnd(2); k4 < box.bigEnd(2); ++k4) {
        if(box.contains(amrex::IntVect(0,0,k4)))
        {
          z_r_grid.push_back(z_r_arr(0,0,k4));
          z_w_grid.push_back(z_w_arr(0,0,k4));
          countp[1]=box.length(2);
          amrex::Print()<<z_r_arr(0,0,k4)<<"\t"<<startp[0]<<"\t"<<countp[0]<<std::endl;
        }
      }
      if(countp[1]>0)
      {
        auto nc_z_r_grid = ncf.var("z_r_grid");
        nc_z_r_grid.par_access(NC_COLLECTIVE);
        nc_z_r_grid.put(z_r_grid.data(), startp, countp);
        auto nc_z_w_grid = ncf.var("z_w_grid");
        nc_z_w_grid.par_access(NC_COLLECTIVE);
        nc_z_w_grid.put(z_w_grid.data(), startp, countp);
      }
   }
   size_t nfai = 0;
   const int ncomp = plotMF[lev]->nComp();

   for (amrex::MFIter fai(*plotMF[lev]); fai.isValid(); ++fai) {
       auto box             = fai.validbox();
       if (subdomain.contains(box)) {
           long unsigned numpts = box.numPts();
           auto array_version = plotMF[lev]->array(fai);

           #ifdef ROMSX_USE_HISTORYFILE
           int num_var_dims=AMREX_SPACEDIM+1;
           #else
           int num_var_dims=AMREX_SPACEDIM;
           #endif
           std::vector<size_t> startp(num_var_dims);
           std::vector<size_t> countp(num_var_dims);
           std::vector<ptrdiff_t> stride(num_var_dims);

           stride[num_var_dims-1]=(ptrdiff_t) (&(array_version(1,0,0,0))-&(array_version(0,0,0.0)));
           stride[num_var_dims-2]=(ptrdiff_t) (&(array_version(0,1,0,0))-&(array_version(0,0,0.0)))/stride[num_var_dims];
           stride[num_var_dims-3]=(ptrdiff_t) (&(array_version(0,0,1,0))-&(array_version(0,0,0.0)))/stride[num_var_dims-1];

           for (int i=0;i<3;i++) {
               startp[num_var_dims-1-i]=box.smallEnd(i);
               countp[num_var_dims-1-i]=box.length(i);
           }
#ifdef ROMSX_USE_HISTORYFILE
           startp[0]=total_plot_file_step_1+1;
           countp[0]=1;
           stride[0]=1;
#endif
          for (int k(0); k < ncomp; ++k) {
              auto data = plotMF[lev]->get(fai).dataPtr(k);
              auto nc_plot_var = ncf.var(plot_var_names[k]);
              nc_plot_var.par_access(NC_COLLECTIVE);
              nc_plot_var.put(data, startp, countp);
              //              nc_plot_var.put(data, startp, countp, stride);
          }
          /*
          {
          auto z_w_data = z_w[lev]->get(fai).dataPtr();
          auto nc_plot_var = ncf.var("z_w");
          nc_plot_var.par_access(NC_COLLECTIVE);
          nc_plot_var.put(z_w_data, startp, countp);
          }
          {
          auto z_r_data = z_r[lev]->get(fai).dataPtr();
          auto nc_plot_var = ncf.var("z_r");
          nc_plot_var.par_access(NC_COLLECTIVE);
          nc_plot_var.put(z_r_data, startp, countp);
          }*/
          nfai++;
       }
   }

   ncf.close();
}
