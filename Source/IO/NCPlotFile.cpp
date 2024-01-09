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

     // Number of cells in this "domain" at this level
     std::vector<int> n_cells;

#ifdef REMORA_USE_HISTORYFILE
     std::string plotfilename = Concatenate("plt_his", which_step, 5);
#else
     std::string plotfilename = Concatenate("plt", which_step, 5);
#endif

     // Set the full IO path for NetCDF output
     std::string FullPath = plotfilename;
     if (lev == 0) {
         const std::string& extension = amrex::Concatenate("_d",lev+1,2);
         FullPath += extension + ".nc";
     } else {
         const std::string& extension = amrex::Concatenate("_d",lev+1+which_subdomain,2);
         FullPath += extension + ".nc";
     }

#ifdef REMORA_USE_HISTORYFILE
     bool not_empty_file = true;
     if (true) {
         FullPath = plot_file_1 + "_his.nc";
         not_empty_file = amrex::FileSystem::Exists(FullPath);
     }
#endif

     amrex::Print() << "Writing level " << lev << " NetCDF plot file " << FullPath 
                    << "\nFor step "<< which_step <<std::endl;

#ifdef REMORA_USE_HISTORYFILE
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

     Box subdomain;
     if (lev == 0) {
         subdomain = geom[lev].Domain();
     } else {
         subdomain = boxes_at_level[lev][which_subdomain];
     }

     int dom_nx = subdomain.length(0);
     int dom_ny = subdomain.length(1);
     int dom_nz = subdomain.length(2);
     n_cells.push_back(dom_nx); n_cells.push_back(dom_ny); n_cells.push_back(dom_nz);
     int num_pts = dom_nx * dom_ny * dom_nz;

     const std::string nt_name   = "ocean_time";
     const std::string ndim_name = "num_geo_dimensions";
     const std::string np_name   = "num_points_per_block";
     const std::string nb_name   = "num_blocks";
     const std::string x_r_name   = "x_rho";
     const std::string y_r_name   = "y_rho";

     const std::string s_r_name   = "s_rho";
     const std::string s_w_name   = "s_w";

     const std::string xi_rho_name   = "xi_rho";
     const std::string xi_u_name     = "xi_u";
     const std::string xi_v_name     = "xi_v";
     const std::string xi_psi_name   = "xi_psi";
     const std::string xi_vert_name  = "xi_vert";

     const std::string eta_rho_name  = "eta_rho";
     const std::string eta_u_name    = "eta_u";
     const std::string eta_v_name    = "eta_v";
     const std::string eta_psi_name  = "eta_psi";
     const std::string eta_vert_name = "eta_vert";

     const std::string z_r_name   = "z_rho";
     const std::string z_w_name   = "z_w";

#ifdef REMORA_USE_HISTORYFILE
     if (!not_empty_file)
#endif
     {
         ncf.enter_def_mode();
         ncf.put_attr("title", "REMORA NetCDF Plot data output");
         ncf.def_dim(nt_name,   NC_UNLIMITED);
         ncf.def_dim(ndim_name, AMREX_SPACEDIM);
         ncf.def_dim(np_name,   num_pts);
         ncf.def_dim(nb_name,   nblocks);

         ncf.def_dim(xi_rho_name,  dom_nx+2);
         ncf.def_dim(xi_u_name,    dom_nx+1);
         ncf.def_dim(xi_v_name,    dom_nx+2);
         ncf.def_dim(xi_psi_name,  dom_nx+1);
         ncf.def_dim(xi_vert_name, dom_nx+2);

         ncf.def_dim(eta_rho_name,  dom_ny+2);
         ncf.def_dim(eta_u_name,    dom_ny+2);
         ncf.def_dim(eta_v_name,    dom_ny+1);
         ncf.def_dim(eta_psi_name,  dom_ny+1);
         ncf.def_dim(eta_vert_name, dom_ny+2);

         ncf.def_dim(s_r_name,      dom_nz+2);
         ncf.def_dim(s_w_name,      dom_nz+1);

         ncf.def_dim(x_r_name,     (dom_nx+2)*(dom_ny+2));
         ncf.def_dim(y_r_name,     (dom_nx+2)*(dom_ny+2));
         ncf.def_dim(z_r_name,     (dom_nx+2)*(dom_ny+2)*(dom_nz+2));
         ncf.def_dim(z_w_name,     (dom_nx+2)*(dom_ny+2)*(dom_nz+1));

         // ncf.def_var("z_w"     , NC_FLOAT, {s_r_name, eta_vert_name, xi_vert_name});
         // ncf.def_var("s_rho"   , NC_FLOAT, {s_r_name});
         ncf.def_var("x_rho"   , NC_FLOAT, {eta_rho_name, xi_rho_name});
         ncf.def_var("y_rho"   , NC_FLOAT, {eta_rho_name, xi_rho_name});
         ncf.def_var("z_rho"   , NC_FLOAT, {s_r_name, eta_rho_name, xi_rho_name});

         ncf.def_var("x_grid"  , NC_FLOAT, {nb_name, x_r_name});
         ncf.def_var("y_grid"  , NC_FLOAT, {nb_name, y_r_name});
         ncf.def_var("z_grid"  , NC_FLOAT, {nb_name, z_r_name});
         ncf.def_var("z_w_grid", NC_FLOAT, {nb_name, z_w_name});

#ifdef REMORA_USE_HISTORYFILE
         ncf.def_var("temp", NC_FLOAT, {nt_name, s_r_name, eta_rho_name, xi_rho_name});
         ncf.def_var("salt", NC_FLOAT, {nt_name, s_r_name, eta_rho_name, xi_rho_name});
         ncf.def_var("u"   , NC_FLOAT, {nt_name, s_r_name, eta_u_name  , xi_u_name});
         ncf.def_var("v"   , NC_FLOAT, {nt_name, s_r_name, eta_v_name  , xi_v_name});
#else
         ncf.def_var("temp", NC_FLOAT, {s_r_name, eta_rho_name, xi_rho_name});
         ncf.def_var("salt", NC_FLOAT, {s_r_name, eta_rho_name, xi_rho_name});
         ncf.def_var("u"   , NC_FLOAT, {s_r_name, eta_u_name  , xi_u_name});
         ncf.def_var("v"   , NC_FLOAT, {s_r_name, eta_v_name  , xi_v_name});
#endif
         ncf.exit_def_mode();
     }

     // Here we hardwire our output to plot {temp, salt, u, v, ubar, vbar};
     int n_data_items = 6;

     ncf.put_attr("number_variables", std::vector<int>{n_data_items});

     // std::vector<Real> xi_grid;
     // std::vector<Real> eta_grid;
     // std::vector<Real> vert_grid;

     // std::vector<Real> z_r_grid;
     // std::vector<Real> z_w_grid;

     const unsigned long indexOffset = 1;

#ifdef REMORA_USE_HISTORYFILE
     int num_var_dims=AMREX_SPACEDIM+1; // we include time if writing a history file
#else
     int num_var_dims=AMREX_SPACEDIM;
#endif

      std::vector<size_t> startp(num_var_dims);
      std::vector<size_t> countp(num_var_dims);

#ifdef REMORA_USE_HISTORYFILE
      startp[0] = REMORA::total_nc_plot_file_step+1;
      countp[0] = 1;
#endif

     // Number of points in each block at this level
     std::vector<int> offset_s;
     std::vector<int> offset_u;
     std::vector<int> offset_v;

     int iproc = amrex::ParallelContext::MyProcAll();
     int nproc = ParallelDescriptor::NProcs();

     auto dm = cons_new[lev]->DistributionMap();

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
        offset_s[dm[ib]] += tmp_bx_s.numPts();

        Box tmp_bx_u(grids[lev][ib]); tmp_bx_u.surroundingNodes(0); tmp_bx_u.grow(IntVect(0,1,0));
        offset_u[dm[ib]] += tmp_bx_u.numPts();

        Box tmp_bx_v(grids[lev][ib]); tmp_bx_v.surroundingNodes(1); tmp_bx_v.grow(IntVect(1,0,0));
        offset_v[dm[ib]] += tmp_bx_v.numPts();
     }

     long unsigned diff_s = 0;
     long unsigned npts_s = 0;

     // No tiling
     for (MFIter mfi(*cons_new[lev],false); mfi.isValid(); ++mfi)
     {
         Box box = mfi.validbox();
         int idx = mfi.index();

         if (subdomain.contains(box)) 
         {
             diff_s = npts_s;

             for (auto ip = 0; ip < iproc; ++ip) {
                 diff_s += offset_s[ip];
             }

             Box tmp_bx(box); tmp_bx.grow(IntVect(1,1,0));
             amrex::Print() << "TMP_BX FOR S " << tmp_bx << std::endl; 

             long unsigned numpts = tmp_bx.numPts();

             amrex::Print() << " DIFF NMPTS " << diff_s << " " << numpts << std::endl;

             FArrayBox tmp(tmp_bx,1);

             // Writing temperature (1 ghost cell in lateral directions)
             {
             tmp.template copy<RunOn::Device>((*cons_new[lev])[idx],Temp_comp,0,1);
             auto data = tmp.dataPtr();
             auto nc_plot_var = ncf.var("temp");
             nc_plot_var.par_access(NC_INDEPENDENT);
             nc_plot_var.put(data, {diff_s}, {numpts});
             }
  
             // Writing salt (1 ghost cell in lateral directions)
             {
             tmp.template copy<RunOn::Device>((*cons_new[lev])[idx],Salt_comp,0,1);
             auto data = tmp.dataPtr();
             auto nc_plot_var = ncf.var("salt");
             nc_plot_var.par_access(NC_INDEPENDENT);
             nc_plot_var.put(data, {diff_s}, {numpts});
             }

             npts_s += numpts;

         } // in subdomain
     } // mfi

     Print() << "FINISHED WRITING TEMP AND SALT " << std::endl;

     long unsigned diff_u = 0;
     long unsigned npts_u = 0;
   
     // Writing u
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
             amrex::Print() << "TMP_BX FOR U " << tmp_bx << std::endl; 
 
             long unsigned numpts = tmp_bx.numPts();

             FArrayBox tmp(tmp_bx,1);
             tmp.template copy<RunOn::Device>((*xvel_new[lev])[idx],0,0,1);

             auto data = tmp.dataPtr();

             auto nc_plot_var = ncf.var("u");
             nc_plot_var.par_access(NC_INDEPENDENT);
             nc_plot_var.put(data, {diff_u}, {numpts});

             npts_u += numpts;

         } // in subdomain
     } // mfi

     Print() << "FINISHED WRITING U " << std::endl;

     long unsigned diff_v = 0;
     long unsigned npts_v = 0;
   
     // Writing v
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
             amrex::Print() << "TMP_BX FOR V " << tmp_bx << std::endl; 

             long unsigned numpts = tmp_bx.numPts();

             FArrayBox tmp(tmp_bx,1);
             tmp.template copy<RunOn::Device>((*yvel_new[lev])[idx],0,0,1);

             auto data = tmp.dataPtr();

             auto nc_plot_var = ncf.var("v");
             nc_plot_var.par_access(NC_INDEPENDENT);
             nc_plot_var.put(data, {diff_v}, {numpts});

             npts_v += numpts;

         } // in subdomain
     } // mfi

     Print() << "FINISHED WRITING V " << std::endl;

#if 0
     // Writing others
     for (MFIter mfi(*cons_new[lev],false); mfi.isValid(); ++mfi)
     {
         Box box = mfi.validbox();
         int idx = mfi.index();

         if (subdomain.contains(box)) 
         {

             {
             auto data = vec_x_r[lev]->get(mfi).dataPtr();
             auto nc_plot_var = ncf.var(x_r_name);
             nc_plot_var.par_access(NC_COLLECTIVE);
             nc_plot_var.put(data, {(unsigned long) box.smallEnd(1)+indexOffset,
                                    (unsigned long) box.smallEnd(0)+indexOffset},
                                   {(unsigned long) box.length(1),
                                    (unsigned long) box.length(0)});
             }

             {
             auto data = vec_y_r[lev]->get(mfi).dataPtr();
             auto nc_plot_var = ncf.var(y_r_name);
             nc_plot_var.par_access(NC_COLLECTIVE);
             nc_plot_var.put(data, {(unsigned long) box.smallEnd(1)+indexOffset,
                                    (unsigned long) box.smallEnd(0)+indexOffset},
                                   {(unsigned long) box.length(1),
                                    (unsigned long) box.length(0)});
             }

             {
             auto data = vec_z_r[lev]->get(mfi).dataPtr();
             auto nc_plot_var = ncf.var(z_r_name);
             nc_plot_var.par_access(NC_COLLECTIVE);
             nc_plot_var.put(data, {(unsigned long) box.smallEnd(2),
                                    (unsigned long) box.smallEnd(1)+indexOffset,
                                    (unsigned long) box.smallEnd(0)+indexOffset},
                                   {(unsigned long) box.length(2),
                                    (unsigned long) box.length(1),
                                    (unsigned long) box.length(0)});
             }

             {
             auto data = vec_z_w[lev]->get(mfi).dataPtr();
             auto nc_plot_var = ncf.var(z_w_name);
             nc_plot_var.par_access(NC_COLLECTIVE);
             nc_plot_var.put(data, {(unsigned long) box.smallEnd(2),
                                    (unsigned long) box.smallEnd(1)+indexOffset,
                                    (unsigned long) box.smallEnd(0)+indexOffset},
                                   {(unsigned long) box.length(2),
                                    (unsigned long) box.length(1),
                                    (unsigned long) box.length(0)});
             }

             {
             auto data = vec_s_r[lev]->get(mfi).dataPtr();
             auto nc_plot_var = ncf.var(s_r_name);
             nc_plot_var.par_access(NC_COLLECTIVE);
             nc_plot_var.put(data, {(unsigned long) box.smallEnd(2)}, {(unsigned long) box.length(2)});
             }
         } // if in subdomain
     } // mfi
#endif

#ifdef REMORA_USE_HISTORYFILE
   REMORA::total_nc_plot_file_step += 1 ;
#endif

   ncf.close();
}

