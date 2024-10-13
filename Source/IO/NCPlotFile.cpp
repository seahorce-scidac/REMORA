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
REMORA::WriteNCPlotFile(int which_step)
{
    AMREX_ASSERT(max_level==0);
    // For right now we assume single level -- we will generalize this later to multilevel
    int lev             = 0;
    int which_subdomain = 0;

    // Create filename
    std::string plt_string;
    std::string plotfilename;
    if (REMORA::write_history_file) {
        plotfilename = plot_file_name + "_his";
    } else {
        plotfilename = Concatenate(plot_file_name, which_step, 5);
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

    //
    // Check if this file/directory already exists and if so,
    //       have the IOProcessor move the existing
    //       file/directory to filename.old
    //
    if ( (!REMORA::write_history_file) || (which_step == 0) ) {
        if (amrex::ParallelDescriptor::IOProcessor()) {
            if (amrex::FileExists(FullPath)) {
                std::string newoldname(FullPath + ".old." + amrex::UniqueString());
                amrex::Print() << "WriteNCPlotFile:  " << FullPath
                               << " exists.  Renaming to:  " << newoldname << std::endl;
                if (std::rename(FullPath.c_str(), newoldname.c_str())) {
                    amrex::Abort("WriteNCPlotfile:: std::rename failed");
                }
            }
        }
        ParallelDescriptor::Barrier();
    }

    bool is_history;

    if (REMORA::write_history_file)
    {
        is_history = true;
        bool write_header = !(amrex::FileExists(FullPath));

        auto ncf = write_header ?
        ncutils::NCFile::create(FullPath, NC_CLOBBER,
                                    amrex::ParallelContext::CommunicatorSub(), MPI_INFO_NULL) :
        ncutils::NCFile::open(FullPath, NC_WRITE,
                                  amrex::ParallelContext::CommunicatorSub(), MPI_INFO_NULL);

        amrex::Print() << "Writing into level " << lev << " NetCDF history file " << FullPath << std::endl;

        WriteNCPlotFile_which(lev, which_subdomain, write_header, ncf, is_history);

    } else {

        is_history = false;
        bool write_header = true;

        // Open new netcdf file to write data
        auto ncf = ncutils::NCFile::create(FullPath, NC_CLOBBER,
                                               amrex::ParallelContext::CommunicatorSub(), MPI_INFO_NULL);
        amrex::Print() << "Writing level " << lev << " NetCDF plot file " << FullPath << std::endl;

        WriteNCPlotFile_which(lev, which_subdomain, write_header, ncf, is_history);
    }
}

void
REMORA::WriteNCPlotFile_which(int lev, int which_subdomain,
                              bool write_header, ncutils::NCFile& ncf,
                              bool is_history)
{
    // Number of cells in this "domain" at this level
    std::vector<int> n_cells;

    int nblocks = grids[lev].size();

    // We only do single-level writes when using NetCDF format
    int flev = 1; //max_level;

    Box subdomain;
    if (lev == 0) {
        subdomain = geom[lev].Domain();
    } else {
        subdomain = boxes_at_level[lev][which_subdomain];
    }

    int nx = subdomain.length(0);
    int ny = subdomain.length(1);
    int nz = subdomain.length(2);

    // unsigned long int nt= NC_UNLIMITED;
    if (is_history && max_step < 0) amrex::Abort("Need to know max_step if writing history file");
    unsigned long int nt = is_history ? static_cast<unsigned long int>(max_step / std::min(plot_int,max_step))+1 : 1;

    n_cells.push_back(nx);
    n_cells.push_back(ny);
    n_cells.push_back(nz);

    int num_pts   = nx*ny*nz;
    int num_u_pts = (nx+1)*(ny+2)*nz;
    int num_v_pts = (nx+2)*(ny+1)*nz;

    const std::string nt_name   = "ocean_time";
    const std::string ndim_name = "num_geo_dimensions";

    const std::string np_name   = "num_points_per_block";
    const std::string np_u_name = "num_points_per_block_u";
    const std::string np_v_name = "num_points_per_block_v";

    const std::string nb_name   = "num_blocks";
    const std::string flev_name = "FINEST_LEVEL";

    const std::string nx_name   = "NX";
    const std::string ny_name   = "NY";
    const std::string nz_name   = "NZ";

    const std::string nx_s_name   = "xi_rho";
    const std::string ny_s_name   = "eta_rho";
    const std::string nz_s_name   = "s_rho";

    const std::string nx_u_name   = "xi_u";
    const std::string ny_u_name   = "eta_u";

    const std::string nx_v_name   = "xi_v";
    const std::string ny_v_name   = "eta_v";

    const Real fill_value = 1.0e37_rt;

    if (write_header)
    {
        ncf.enter_def_mode();
        ncf.put_attr("title", "REMORA data ");
        ncf.def_dim(nt_name,   nt);
        ncf.def_dim(ndim_name, AMREX_SPACEDIM);

        ncf.def_dim(nx_s_name  , nx+2);
        ncf.def_dim(ny_s_name  , ny+2);
        ncf.def_dim(nz_s_name  , nz);

        ncf.def_dim(nx_u_name  , nx+1);
        ncf.def_dim(ny_u_name  , ny+2);

        ncf.def_dim(nx_v_name  , nx+2);
        ncf.def_dim(ny_v_name  , ny+1);

        ncf.def_dim(np_name  , num_pts);
        ncf.def_dim(np_u_name, num_u_pts);
        ncf.def_dim(np_v_name, num_v_pts);

        ncf.def_dim(nb_name,   nblocks);
        ncf.def_dim(flev_name, flev);

        ncf.def_dim(nx_name,   n_cells[0]);
        ncf.def_dim(ny_name,   n_cells[1]);
        ncf.def_dim(nz_name,   n_cells[2]);

        ncf.def_var("probLo"  ,   ncutils::NCDType::Real,  {ndim_name});
        ncf.def_var("probHi"  ,   ncutils::NCDType::Real,  {ndim_name});

        ncf.def_var("Geom.smallend", NC_INT, {flev_name, ndim_name});
        ncf.def_var("Geom.bigend"  , NC_INT, {flev_name, ndim_name});
        ncf.def_var("CellSize"     , ncutils::NCDType::Real, {flev_name, ndim_name});

        ncf.def_var("x_grid", ncutils::NCDType::Real, {np_name});
        ncf.def_var("y_grid", ncutils::NCDType::Real, {np_name});
        ncf.def_var("z_grid", ncutils::NCDType::Real, {np_name});
        ncf.def_var("ocean_time", ncutils::NCDType::Real, {nt_name});

        ncf.def_var("h"   , ncutils::NCDType::Real, {         ny_s_name, nx_s_name});
        ncf.def_var_fill("zeta", ncutils::NCDType::Real, {nt_name, ny_s_name, nx_s_name}, &fill_value);
        ncf.def_var_fill("temp", ncutils::NCDType::Real, {nt_name, nz_s_name, ny_s_name, nx_s_name}, &fill_value);
        ncf.def_var_fill("salt", ncutils::NCDType::Real, {nt_name, nz_s_name, ny_s_name, nx_s_name}, &fill_value);
        ncf.def_var_fill("u"   , ncutils::NCDType::Real, {nt_name, nz_s_name, ny_u_name, nx_u_name}, &fill_value);
        ncf.def_var_fill("v"   , ncutils::NCDType::Real, {nt_name, nz_s_name, ny_v_name, nx_v_name}, &fill_value);
        ncf.def_var_fill("ubar", ncutils::NCDType::Real, {nt_name, ny_u_name, nx_u_name}, &fill_value);
        ncf.def_var_fill("vbar", ncutils::NCDType::Real, {nt_name, ny_v_name, nx_v_name}, &fill_value);
        ncf.def_var("sustr", ncutils::NCDType::Real, {nt_name, ny_u_name, nx_u_name});
        ncf.def_var("svstr", ncutils::NCDType::Real, {nt_name, ny_v_name, nx_v_name});

        Real time = 0.;

        // Right now this is hard-wired to {temp, salt, u, v}
        int n_data_items = 4;
        ncf.put_attr("number_variables", std::vector<int>{n_data_items});
        ncf.put_attr("space_dimension", std::vector<int>{AMREX_SPACEDIM});
        ncf.put_attr("current_time", std::vector<double>{time});
        ncf.put_attr("start_time", std::vector<double>{start_bdy_time});
        ncf.put_attr("CurrentLevel", std::vector<int>{flev});
        ncf.put_attr("DefaultGeometry", std::vector<int>{amrex::DefaultGeometry().Coord()});

        ncf.exit_def_mode();

        // We are doing single-level writes but it doesn't have to be level 0
        //
        // Write out the header information.
        //

        Real dx[AMREX_SPACEDIM];
        for (int i = 0; i < AMREX_SPACEDIM; i++) {
            dx[i] = geom[lev].CellSize()[i];
        }
        const auto *base = geom[lev].ProbLo();
        RealBox rb(subdomain,dx,base);

        amrex::Vector<Real> probLo;
        amrex::Vector<Real> probHi;
        for (int i = 0; i < AMREX_SPACEDIM; i++) {
            probLo.push_back(rb.lo(i));
            probHi.push_back(rb.hi(i));
        }


        //nc_probLo.par_access(NC_COLLECTIVE);
        // small variable data written by just the master proc
        ncmpi_begin_indep_data(ncf.ncid);
        if (amrex::ParallelDescriptor::IOProcessor()) // only master proc
        {
            auto nc_probLo = ncf.var("probLo");

           nc_probLo.put(probLo.data(), {0}, {AMREX_SPACEDIM});

           auto nc_probHi = ncf.var("probHi");
        //nc_probHi.par_access(NC_COLLECTIVE);
           nc_probHi.put(probHi.data(), {0}, {AMREX_SPACEDIM});

            amrex::Vector<int> smallend;
            amrex::Vector<int> bigend;
            for (int i = lev; i < flev; i++) {
                smallend.clear(); bigend.clear();
                for (int j = 0; j < AMREX_SPACEDIM; j++) {
                    smallend.push_back(subdomain.smallEnd(j));
                    bigend.push_back(subdomain.bigEnd(j));
                }
                auto nc_Geom_smallend = ncf.var("Geom.smallend");
                //nc_Geom_smallend.par_access(NC_COLLECTIVE);
                nc_Geom_smallend.put(smallend.data(), {static_cast<long unsigned int>(i-lev), 0}, {1, AMREX_SPACEDIM});

                auto nc_Geom_bigend = ncf.var("Geom.bigend");
                //nc_Geom_bigend.par_access(NC_COLLECTIVE);
                nc_Geom_bigend.put(bigend.data(), {static_cast<long unsigned int>(i-lev), 0}, {1, AMREX_SPACEDIM});
            }

            amrex::Vector<Real> CellSize;
            for (int i = lev; i < flev; i++) {
                CellSize.clear();
                for (Real & j : dx) {
                    CellSize.push_back(amrex::Real(j));
                }
                auto nc_CellSize = ncf.var("CellSize");
                //nc_CellSize.par_access(NC_COLLECTIVE);
                nc_CellSize.put(CellSize.data(), {static_cast<long unsigned int>(i-lev), 0}, {1, AMREX_SPACEDIM});
            }

        }
        ncmpi_end_indep_data(ncf.ncid);

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

                //nc_x_grid.par_access(NC_INDEPENDENT);
                //nc_y_grid.par_access(NC_INDEPENDENT);
                //nc_z_grid.par_access(NC_INDEPENDENT);

                nc_x_grid.put_all(x_grid.data(), {goffset}, {glen});
                nc_y_grid.put_all(y_grid.data(), {goffset}, {glen});
                nc_z_grid.put_all(z_grid.data(), {goffset}, {glen});
            }
        }

    } // end if write_header

    ncmpi_begin_indep_data(ncf.ncid);
    //
    // We compute the offsets based on location of the box within the domain
    //
    long unsigned local_start_nt = (is_history ? static_cast<long unsigned>(history_count) : static_cast<long unsigned>(0));
    long unsigned local_nt       = 1;  // We write data for only one time

    {
    auto nc_plot_var = ncf.var("ocean_time");
    //nc_plot_var.par_access(NC_COLLECTIVE);
    nc_plot_var.put(&t_new[lev], {local_start_nt}, {local_nt});
    }
    ncmpi_end_indep_data(ncf.ncid);

    // temporary
    //ncf.close();

    cons_new[lev]->FillBoundary(geom[lev].periodicity());

    mask_arrays_for_write(lev, (Real) fill_value);

    for (MFIter mfi(*cons_new[lev],false); mfi.isValid(); ++mfi)
    {
        auto bx = mfi.validbox();
        if (subdomain.contains(bx))
        {
            //
            // We only include one grow cell at subdomain boundaries, not internal grid boundaries
            //
            Box tmp_bx(bx);
            if (tmp_bx.smallEnd()[0] == subdomain.smallEnd()[0]) tmp_bx.growLo(0,1);
            if (tmp_bx.smallEnd()[1] == subdomain.smallEnd()[1]) tmp_bx.growLo(1,1);
            if (tmp_bx.bigEnd()[0]   == subdomain.bigEnd()[0])   tmp_bx.growHi(0,1);
            if (tmp_bx.bigEnd()[1]   == subdomain.bigEnd()[1])   tmp_bx.growHi(1,1);
            // amrex::Print() << "    BX " << bx << std::endl;
            // amrex::Print() << "TMP_BX " << tmp_bx << std::endl;

            Box tmp_bx_2d(tmp_bx);
            tmp_bx_2d.makeSlab(2,0);

            //
            // These are the dimensions of the data we write for only this box
            //
            long unsigned local_nx = tmp_bx.length()[0];
            long unsigned local_ny = tmp_bx.length()[1];
            long unsigned local_nz = tmp_bx.length()[2];

            // We do the "+1" because the offset needs to start at 0
            long unsigned local_start_x  = static_cast<long unsigned>(tmp_bx.smallEnd()[0]+1);
            long unsigned local_start_y  = static_cast<long unsigned>(tmp_bx.smallEnd()[1]+1);
            long unsigned local_start_z  = static_cast<long unsigned>(tmp_bx.smallEnd()[2]);

            if (write_header) {
                FArrayBox tmp_bathy;
                tmp_bathy.resize(tmp_bx_2d,1,amrex::The_Pinned_Arena());

                tmp_bathy.template copy<RunOn::Device>((*vec_hOfTheConfusingName[lev])[mfi.index()],0,0,1);
                Gpu::streamSynchronize();

                auto nc_plot_var = ncf.var("h");
                //nc_plot_var.par_access(NC_INDEPENDENT);
                nc_plot_var.put_all(tmp_bathy.dataPtr(), {local_start_y,local_start_x},
                                                     {local_ny, local_nx});
            }

            {
            FArrayBox tmp_zeta;
            tmp_zeta.resize(tmp_bx_2d,1,amrex::The_Pinned_Arena());
            tmp_zeta.template copy<RunOn::Device>((*vec_Zt_avg1[lev])[mfi.index()],0,0,1);
            Gpu::streamSynchronize();

            auto nc_plot_var = ncf.var("zeta");
            //nc_plot_var.par_access(NC_INDEPENDENT);
            nc_plot_var.put_all(tmp_zeta.dataPtr(), {local_start_nt,local_start_y,local_start_x},
                                                {local_nt, local_ny, local_nx});
            }

            {
            FArrayBox tmp_temp;
            tmp_temp.resize(tmp_bx,1,amrex::The_Pinned_Arena());
            tmp_temp.template copy<RunOn::Device>((*cons_new[lev])[mfi.index()],Temp_comp,0,1);
            Gpu::streamSynchronize();

            auto nc_plot_var = ncf.var("temp");
            //nc_plot_var.par_access(NC_INDEPENDENT);
            nc_plot_var.put_all(tmp_temp.dataPtr(), {local_start_nt,local_start_z,local_start_y,local_start_x},
                                                {local_nt, local_nz, local_ny, local_nx});
            }

            {
            FArrayBox tmp_salt;
            tmp_salt.resize(tmp_bx,1,amrex::The_Pinned_Arena());
            tmp_salt.template copy<RunOn::Device>((*cons_new[lev])[mfi.index()],Salt_comp,0,1);
            Gpu::streamSynchronize();

            auto nc_plot_var = ncf.var("salt");
            //nc_plot_var.par_access(NC_INDEPENDENT);
            nc_plot_var.put_all(tmp_salt.dataPtr(), {local_start_nt,local_start_z,local_start_y,local_start_x},
                                                {local_nt, local_nz, local_ny, local_nx});
            }
        } // subdomain
    } // mfi

    // Writing u (we loop over cons to get cell-centered box)
    for (MFIter mfi(*cons_new[lev],false); mfi.isValid(); ++mfi)
    {
        Box bx = mfi.validbox();

        if (subdomain.contains(bx))
        {
            //
            // We only include one grow cell at subdomain boundaries, not internal grid boundaries
            //
            Box tmp_bx(bx); tmp_bx.surroundingNodes(0);
            if (tmp_bx.smallEnd()[1] == subdomain.smallEnd()[1]) tmp_bx.growLo(1,1);
            if (tmp_bx.bigEnd()[1]   == subdomain.bigEnd()[1])   tmp_bx.growHi(1,1);
            Box tmp_bx_2d(tmp_bx);
            tmp_bx_2d.makeSlab(2,0);

            //
            // These are the dimensions of the data we write for only this box
            //
            long unsigned local_nx = tmp_bx.length()[0];
            long unsigned local_ny = tmp_bx.length()[1];
            long unsigned local_nz = tmp_bx.length()[2];

            // We do the "+1" because the offset needs to start at 0
            long unsigned local_start_x  = static_cast<long unsigned>(tmp_bx.smallEnd()[0]);
            long unsigned local_start_y  = static_cast<long unsigned>(tmp_bx.smallEnd()[1]+1);
            long unsigned local_start_z  = static_cast<long unsigned>(tmp_bx.smallEnd()[2]);

            {
            FArrayBox tmp;
            tmp.resize(tmp_bx,1,amrex::The_Pinned_Arena());
            tmp.template copy<RunOn::Device>((*xvel_new[lev])[mfi.index()],0,0,1);
            Gpu::streamSynchronize();

            auto nc_plot_var = ncf.var("u");
            //nc_plot_var.par_access(NC_INDEPENDENT);
            nc_plot_var.put_all(tmp.dataPtr(), {local_start_nt,local_start_z,local_start_y,local_start_x},
                                           {local_nt, local_nz, local_ny, local_nx});
            }

            {
            FArrayBox tmp;
            tmp.resize(tmp_bx_2d,1,amrex::The_Pinned_Arena());
            tmp.template copy<RunOn::Device>((*vec_ubar[lev])[mfi.index()],0,0,1);
            Gpu::streamSynchronize();

            auto nc_plot_var = ncf.var("ubar");
            //nc_plot_var.par_access(NC_INDEPENDENT);
            nc_plot_var.put_all(tmp.dataPtr(), {local_start_nt,local_start_y,local_start_x},
                                           {local_nt, local_ny, local_nx});
            }
            {
            FArrayBox tmp;
            tmp.resize(tmp_bx_2d,1,amrex::The_Pinned_Arena());
            tmp.template copy<RunOn::Device>((*vec_sustr[lev])[mfi.index()],0,0,1);
            Gpu::streamSynchronize();

            auto nc_plot_var = ncf.var("sustr");
            //nc_plot_var.par_access(NC_INDEPENDENT);
            nc_plot_var.put_all(tmp.dataPtr(), {local_start_nt,local_start_y,local_start_x},
                                           {local_nt, local_ny, local_nx});
            }
        } // in subdomain
    } // mfi

    // Writing v (we loop over cons to get cell-centered box)
    for (MFIter mfi(*cons_new[lev],false); mfi.isValid(); ++mfi)
    {
        Box bx = mfi.validbox();

        if (subdomain.contains(bx))
        {
            //
            // We only include one grow cell at subdomain boundaries, not internal grid boundaries
            //
            Box tmp_bx(bx); tmp_bx.surroundingNodes(1);
            if (tmp_bx.smallEnd()[0] == subdomain.smallEnd()[0]) tmp_bx.growLo(0,1);
            if (tmp_bx.bigEnd()[0]   == subdomain.bigEnd()[0])   tmp_bx.growHi(0,1);
            // amrex::Print() << "    BX " << bx << std::endl;
            // amrex::Print() << "TMP_BX " << tmp_bx << std::endl;

            Box tmp_bx_2d(tmp_bx);
            tmp_bx_2d.makeSlab(2,0);

            //
            // These are the dimensions of the data we write for only this box
            //
            long unsigned local_nx = tmp_bx.length()[0];
            long unsigned local_ny = tmp_bx.length()[1];
            long unsigned local_nz = tmp_bx.length()[2];

            // We do the "+1" because the offset needs to start at 0
            long unsigned local_start_x  = static_cast<long unsigned>(tmp_bx.smallEnd()[0]+1);
            long unsigned local_start_y  = static_cast<long unsigned>(tmp_bx.smallEnd()[1]);
            long unsigned local_start_z  = static_cast<long unsigned>(tmp_bx.smallEnd()[2]);

            {
            FArrayBox tmp;
            tmp.resize(tmp_bx,1,amrex::The_Pinned_Arena());
            tmp.template copy<RunOn::Device>((*yvel_new[lev])[mfi.index()],0,0,1);
            Gpu::streamSynchronize();

            auto nc_plot_var = ncf.var("v");
            //nc_plot_var.par_access(NC_INDEPENDENT);
             nc_plot_var.put_all(tmp.dataPtr(), {local_start_nt,local_start_z,local_start_y,local_start_x},
                                            {local_nt, local_nz, local_ny, local_nx});
            }

            {
            FArrayBox tmp;
            tmp.resize(tmp_bx_2d,1,amrex::The_Pinned_Arena());
            tmp.template copy<RunOn::Device>((*vec_vbar[lev])[mfi.index()],0,0,1);
            Gpu::streamSynchronize();

            auto nc_plot_var = ncf.var("vbar");
            //nc_plot_var.par_access(NC_INDEPENDENT);
            nc_plot_var.put_all(tmp.dataPtr(), {local_start_nt,local_start_y,local_start_x},
                                           {local_nt, local_ny, local_nx});
            }
            {
            FArrayBox tmp;
            tmp.resize(tmp_bx_2d,1,amrex::The_Pinned_Arena());
            tmp.template copy<RunOn::Device>((*vec_svstr[lev])[mfi.index()],0,0,1);
            Gpu::streamSynchronize();

            auto nc_plot_var = ncf.var("svstr");
            //nc_plot_var.par_access(NC_INDEPENDENT);
            nc_plot_var.put_all(tmp.dataPtr(), {local_start_nt,local_start_y,local_start_x},
                                           {local_nt, local_ny, local_nx});
            }

        } // in subdomain
    } // mfi

    mask_arrays_for_write(lev, 0.0_rt);

    ncf.close();

    REMORA::total_nc_plot_file_step += 1;
}


