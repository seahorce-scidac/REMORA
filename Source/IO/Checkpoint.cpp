#include <ROMSX.H>
#include "AMReX_PlotFileUtil.H"

using namespace amrex;

// utility to skip to next line in Header
void
ROMSX::GotoNextLine (std::istream& is)
{
    constexpr std::streamsize bl_ignore_max { 100000 };
    is.ignore(bl_ignore_max, '\n');
}

void
ROMSX::WriteCheckpointFile () const
{
    // chk00010            write a checkpoint file with this root directory
    // chk00010/Header     this contains information you need to save (e.g., finest_level, t_new, etc.) and also
    //                     the BoxArrays at each level
    // chk00010/Level_0/
    // chk00010/Level_1/
    // etc.                these subdirectories will hold the MultiFab data at each level of refinement

    // checkpoint file name, e.g., chk00010
    const std::string& checkpointname = amrex::Concatenate(check_file,istep[0],5);

    amrex::Print() << "Writing checkpoint " << checkpointname << "\n";

    const int nlevels = finest_level+1;

    // ---- prebuild a hierarchy of directories
    // ---- dirName is built first.  if dirName exists, it is renamed.  then build
    // ---- dirName/subDirPrefix_0 .. dirName/subDirPrefix_nlevels-1
    // ---- if callBarrier is true, call ParallelDescriptor::Barrier()
    // ---- after all directories are built
    // ---- ParallelDescriptor::IOProcessor() creates the directories
    amrex::PreBuildDirectorHierarchy(checkpointname, "Level_", nlevels, true);

    // write Header file
   if (ParallelDescriptor::IOProcessor()) {

       std::string HeaderFileName(checkpointname + "/Header");
       VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
       std::ofstream HeaderFile;
       HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
       HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |
                                               std::ofstream::trunc |
                                               std::ofstream::binary);
       if( ! HeaderFile.good()) {
           amrex::FileOpenFailed(HeaderFileName);
       }

       HeaderFile.precision(17);

       // write out title line
       HeaderFile << "Checkpoint file for ROMSX\n";

       // write out finest_level
       HeaderFile << finest_level << "\n";

       // write the number of components
       // for each variable we store

       // conservative, cell-centered vars
       HeaderFile << NCONS << "\n";

       // x-velocity on faces
       HeaderFile << 1 << "\n";

       // y-velocity on faces
       HeaderFile << 1 << "\n";

       // z-velocity on faces
       HeaderFile << 1 << "\n";

       HeaderFile << 2 << "\n";

       HeaderFile << 2 << "\n";

       HeaderFile << 3 << "\n";

       HeaderFile << 3 << "\n";

       // write out array of istep
       for (int i = 0; i < istep.size(); ++i) {
           HeaderFile << istep[i] << " ";
       }
       HeaderFile << "\n";

       // write out array of dt
       for (int i = 0; i < dt.size(); ++i) {
           HeaderFile << dt[i] << " ";
       }
       HeaderFile << "\n";

       // write out array of t_new
       for (int i = 0; i < t_new.size(); ++i) {
           HeaderFile << t_new[i] << " ";
       }
       HeaderFile << "\n";

       // write the BoxArray at each level
       for (int lev = 0; lev <= finest_level; ++lev) {
           boxArray(lev).writeOn(HeaderFile);
           HeaderFile << '\n';
       }
   }

   // write the MultiFab data to, e.g., chk00010/Level_0/
   // Here we make copies of the MultiFab with no ghost cells
   for (int lev = 0; lev <= finest_level; ++lev)
   {
       BoxList bl2d = grids[lev].boxList();
       for (auto& b : bl2d) {
           b.setRange(2,0);
       }
       BoxArray ba2d(std::move(bl2d));

       MultiFab cons(grids[lev],dmap[lev],NCONS,0);
       MultiFab::Copy(cons,*cons_new[lev],0,0,NCONS,0);
       VisMF::Write(cons, amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "Cell"));

       MultiFab xvel(convert(grids[lev],IntVect(1,0,0)),dmap[lev],1,0);
       MultiFab::Copy(xvel,*xvel_new[lev],0,0,1,0);
       VisMF::Write(xvel, amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "XFace"));

       MultiFab yvel(convert(grids[lev],IntVect(0,1,0)),dmap[lev],1,0);
       MultiFab::Copy(yvel,*yvel_new[lev],0,0,1,0);
       VisMF::Write(yvel, amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "YFace"));

       MultiFab zvel(convert(grids[lev],IntVect(0,0,1)),dmap[lev],1,0);
       MultiFab::Copy(zvel,*zvel_new[lev],0,0,1,0);
       VisMF::Write(zvel, amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "ZFace"));

       MultiFab mf_ru(grids[lev],dmap[lev],2,NGROW);
       MultiFab::Copy(mf_ru,*vec_ru[lev],0,0,2,NGROW);
       VisMF::Write(mf_ru, amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "XRHS"));

       MultiFab mf_rv(grids[lev],dmap[lev],2,NGROW);
       MultiFab::Copy(mf_rv,*vec_rv[lev],0,0,2,NGROW);
       VisMF::Write(mf_rv, amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "YRHS"));

       MultiFab mf_ubar(ba2d,dmap[lev],3,NGROW);
       MultiFab::Copy(mf_ubar,*(vec_ubar[lev]),0,0,3,NGROW);
       VisMF::Write(mf_ubar, amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "XBar"));

       MultiFab mf_vbar(ba2d,dmap[lev],3,NGROW);
       MultiFab::Copy(mf_vbar,*(vec_vbar[lev]),0,0,3,NGROW);
       VisMF::Write(mf_vbar, amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "YBar"));

       VisMF::Write(*(vec_rufrc[lev]), amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "rufrc"));
       VisMF::Write(*(vec_rvfrc[lev]), amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "rvfrc"));

       VisMF::Write(*(vec_sustr[lev]), amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "sustr"));
       VisMF::Write(*(vec_svstr[lev]), amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "svstr"));

       VisMF::Write(*(vec_rdrag[lev]), amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "rdrag"));

       VisMF::Write(*(vec_bustr[lev]), amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "bustr"));
       VisMF::Write(*(vec_bvstr[lev]), amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "bvstr"));

       VisMF::Write(*(vec_DU_avg1[lev]), amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "DU_avg1"));
       VisMF::Write(*(vec_DU_avg2[lev]), amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "DU_avg2"));
       VisMF::Write(*(vec_DV_avg1[lev]), amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "DV_avg1"));
       VisMF::Write(*(vec_DV_avg2[lev]), amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "DV_avg2"));

       VisMF::Write(*(vec_zeta[lev]), amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "zeta"));
       VisMF::Write(*(vec_Zt_avg1[lev]), amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "Zt_avg1"));
   }

#ifdef ROMSX_USE_PARTICLES
   particleData.Checkpoint(checkpointname);
#endif

#ifdef ROMSX_USE_NETCDF
   // Write bdy_data files
   if ( ParallelDescriptor::IOProcessor() && (solverChoice.ic_bc_type == IC_BC_Type::Real) )
   {

     // Vector dimensions
     int num_time = bdy_data_xlo.size();
     int num_var  = bdy_data_xlo[0].size();

     // Open header file and write to it
     std::ofstream bdy_h_file(amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "bdy_H"));
     bdy_h_file << std::setprecision(1) << std::fixed;
     bdy_h_file << num_time << "\n";
     bdy_h_file << num_var  << "\n";
     bdy_h_file << start_bdy_time << "\n";
     bdy_h_file << bdy_time_interval << "\n";
     bdy_h_file << bdy_width << "\n";
     for (int ivar(0); ivar<num_var; ++ivar) {
       bdy_h_file << bdy_data_xlo[0][ivar].box() << "\n";
       bdy_h_file << bdy_data_xhi[0][ivar].box() << "\n";
       bdy_h_file << bdy_data_ylo[0][ivar].box() << "\n";
       bdy_h_file << bdy_data_yhi[0][ivar].box() << "\n";
     }

     // Open data file and write to it
     std::ofstream bdy_d_file(amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "bdy_D"));
     for (int itime(0); itime<num_time; ++itime) {
       for (int ivar(0); ivar<num_var; ++ivar) {
         bdy_data_xlo[itime][ivar].writeOn(bdy_d_file,0,1);
         bdy_data_xhi[itime][ivar].writeOn(bdy_d_file,0,1);
         bdy_data_ylo[itime][ivar].writeOn(bdy_d_file,0,1);
         bdy_data_yhi[itime][ivar].writeOn(bdy_d_file,0,1);
       }
     }
   }
#endif
}

void
ROMSX::ReadCheckpointFile ()
{
    amrex::Print() << "Restart from checkpoint " << restart_chkfile << "\n";

    // Header
    std::string File(restart_chkfile + "/Header");

    VisMF::IO_Buffer io_buffer(VisMF::GetIOBufferSize());

    Vector<char> fileCharPtr;
    ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
    std::string fileCharPtrString(fileCharPtr.dataPtr());
    std::istringstream is(fileCharPtrString, std::istringstream::in);

    std::string line, word;

    int chk_ncomp;

    // read in title line
    std::getline(is, line);

    // read in finest_level
    is >> finest_level;
    GotoNextLine(is);

    // read the number of components
    // for each variable we store

    // conservative, cell-centered vars
    is >> chk_ncomp;
    GotoNextLine(is);
    AMREX_ASSERT(chk_ncomp == NCONS);

    // x-velocity on faces
    is >> chk_ncomp;
    GotoNextLine(is);
    AMREX_ASSERT(chk_ncomp == 1);

    // y-velocity on faces
    is >> chk_ncomp;
    GotoNextLine(is);
    AMREX_ASSERT(chk_ncomp == 1);

    // z-velocity on faces
    is >> chk_ncomp;
    GotoNextLine(is);
    AMREX_ASSERT(chk_ncomp == 1);

    is >> chk_ncomp;
    GotoNextLine(is);
    AMREX_ASSERT(chk_ncomp == 2);

    is >> chk_ncomp;
    GotoNextLine(is);
    AMREX_ASSERT(chk_ncomp == 2);

    is >> chk_ncomp;
    GotoNextLine(is);
    AMREX_ASSERT(chk_ncomp == 3);

    is >> chk_ncomp;
    GotoNextLine(is);
    AMREX_ASSERT(chk_ncomp == 3);

    // read in array of istep
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
            istep[i++] = std::stoi(word);
        }
    }

    // read in array of dt
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
            dt[i++] = std::stod(word);
        }
    }

    // read in array of t_new
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
            t_new[i++] = std::stod(word);
        }
    }

    for (int lev = 0; lev <= finest_level; ++lev) {

        // read in level 'lev' BoxArray from Header
        BoxArray ba;
        ba.readFrom(is);
        GotoNextLine(is);

        // create a distribution mapping
        DistributionMapping dm { ba, ParallelDescriptor::NProcs() };

        MakeNewLevelFromScratch (lev, t_new[lev], ba, dm);
    }

    // read in the MultiFab data
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        BoxList bl2d = grids[lev].boxList();
        for (auto& b : bl2d) {
            b.setRange(2,0);
        }
        BoxArray ba2d(std::move(bl2d));

        MultiFab cons(grids[lev],dmap[lev],NCONS,0);
        VisMF::Read(cons, amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "Cell"));
        MultiFab::Copy(*cons_new[lev],cons,0,0,NCONS,0);

        MultiFab xvel(convert(grids[lev],IntVect(1,0,0)),dmap[lev],1,0);
        VisMF::Read(xvel, amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "XFace"));
        MultiFab::Copy(*xvel_new[lev],xvel,0,0,1,0);

        MultiFab yvel(convert(grids[lev],IntVect(0,1,0)),dmap[lev],1,0);
        VisMF::Read(yvel, amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "YFace"));
        MultiFab::Copy(*yvel_new[lev],yvel,0,0,1,0);

        MultiFab zvel(convert(grids[lev],IntVect(0,0,1)),dmap[lev],1,0);
        VisMF::Read(zvel, amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "ZFace"));
        MultiFab::Copy(*zvel_new[lev],zvel,0,0,1,0);

       MultiFab mf_ru(grids[lev],dmap[lev],2,NGROW);
       VisMF::Read(mf_ru, amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "XRHS"));
       MultiFab::Copy(*(vec_ru[lev]),mf_ru,0,0,2,(vec_ru[lev])->nGrowVect());

       MultiFab mf_rv(grids[lev],dmap[lev],2,NGROW);
       VisMF::Read(mf_rv, amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "YRHS"));
       MultiFab::Copy(*(vec_rv[lev]),mf_rv,0,0,2,(vec_rv[lev])->nGrowVect());

       MultiFab mf_ubar(ba2d,dmap[lev],3,NGROW);
       VisMF::Read(mf_ubar, amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "XBar"));
       MultiFab::Copy(*(vec_ubar[lev]),mf_ubar,0,0,3,(vec_ubar[lev])->nGrowVect());

       MultiFab mf_vbar(ba2d,dmap[lev],3,NGROW);
       VisMF::Read(mf_vbar, amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "YBar"));
       MultiFab::Copy(*(vec_vbar[lev]),mf_vbar,0,0,3,(vec_vbar[lev])->nGrowVect());

       VisMF::Read(*(vec_rufrc[lev]), amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "rufrc"));
       VisMF::Read(*(vec_rvfrc[lev]), amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "rvfrc"));

       VisMF::Read(*(vec_sustr[lev]), amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "sustr"));
       VisMF::Read(*(vec_svstr[lev]), amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "svstr"));

       VisMF::Read(*(vec_rdrag[lev]), amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "rdrag"));

       VisMF::Read(*(vec_bustr[lev]), amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "bustr"));
       VisMF::Read(*(vec_bvstr[lev]), amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "bvstr"));

       VisMF::Read(*(vec_DU_avg1[lev]), amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "DU_avg1"));
       VisMF::Read(*(vec_DU_avg2[lev]), amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "DU_avg2"));
       VisMF::Read(*(vec_DV_avg1[lev]), amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "DV_avg1"));
       VisMF::Read(*(vec_DV_avg2[lev]), amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "DV_avg2"));

       VisMF::Read(*(vec_zeta[lev]), amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "zeta"));
       VisMF::Read(*(vec_Zt_avg1[lev]), amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "Zt_avg1"));
    }

#ifdef ROMSX_USE_PARTICLES
   particleData.Restart((amrex::ParGDBBase*)GetParGDB(),restart_chkfile);
#endif

#ifdef ROMSX_USE_NETCDF
    // Read bdy_data files
    if ( solverChoice.ic_bc_type == IC_BC_Type::Real)
    {
        int ioproc = ParallelDescriptor::IOProcessorNumber();  // I/O rank
        int num_time;
        int num_var;
        Vector<Box> bx_v;
        if (ParallelDescriptor::IOProcessor()) {
            // Open header file and read from it
            std::ifstream bdy_h_file(amrex::MultiFabFileFullPrefix(0, restart_chkfile, "Level_", "bdy_H"));
            bdy_h_file >> num_time;
            bdy_h_file >> num_var;
            bdy_h_file >> start_bdy_time;
            bdy_h_file >> bdy_time_interval;
            bdy_h_file >> bdy_width;
            bx_v.resize(4*num_var);
            for (int ivar(0); ivar<num_var; ++ivar) {
                bdy_h_file >> bx_v[4*ivar  ];
                bdy_h_file >> bx_v[4*ivar+1];
                bdy_h_file >> bx_v[4*ivar+2];
                bdy_h_file >> bx_v[4*ivar+3];
            }

            // IO size the FABs
            bdy_data_xlo.resize(num_time);
            bdy_data_xhi.resize(num_time);
            bdy_data_ylo.resize(num_time);
            bdy_data_yhi.resize(num_time);
            for (int itime(0); itime<num_time; ++itime) {
                bdy_data_xlo[itime].resize(num_var);
                bdy_data_xhi[itime].resize(num_var);
                bdy_data_ylo[itime].resize(num_var);
                bdy_data_yhi[itime].resize(num_var);
                for (int ivar(0); ivar<num_var; ++ivar) {
                    bdy_data_xlo[itime][ivar].resize(bx_v[4*ivar  ]);
                    bdy_data_xhi[itime][ivar].resize(bx_v[4*ivar+1]);
                    bdy_data_ylo[itime][ivar].resize(bx_v[4*ivar+2]);
                    bdy_data_yhi[itime][ivar].resize(bx_v[4*ivar+3]);
                }
            }

            // Open data file and read from it
            std::ifstream bdy_d_file(amrex::MultiFabFileFullPrefix(0, restart_chkfile, "Level_", "bdy_D"));
            for (int itime(0); itime<num_time; ++itime) {
                for (int ivar(0); ivar<num_var; ++ivar) {
                    bdy_data_xlo[itime][ivar].readFrom(bdy_d_file);
                    bdy_data_xhi[itime][ivar].readFrom(bdy_d_file);
                    bdy_data_ylo[itime][ivar].readFrom(bdy_d_file);
                    bdy_data_yhi[itime][ivar].readFrom(bdy_d_file);
                }
            }
        } // IO

        // Broadcast the data
        ParallelDescriptor::Barrier();
        ParallelDescriptor::Bcast(&start_bdy_time,1,ioproc);
        ParallelDescriptor::Bcast(&bdy_time_interval,1,ioproc);
        ParallelDescriptor::Bcast(&bdy_width,1,ioproc);
        ParallelDescriptor::Bcast(&num_time,1,ioproc);
        ParallelDescriptor::Bcast(&num_var,1,ioproc);

        // Everyone size their boxes
        bx_v.resize(4*num_var);

        ParallelDescriptor::Bcast(bx_v.dataPtr(),bx_v.size(),ioproc);

        // Everyone but IO size their FABs
        if (!ParallelDescriptor::IOProcessor()) {
          bdy_data_xlo.resize(num_time);
          bdy_data_xhi.resize(num_time);
          bdy_data_ylo.resize(num_time);
          bdy_data_yhi.resize(num_time);
          for (int itime(0); itime<num_time; ++itime) {
            bdy_data_xlo[itime].resize(num_var);
            bdy_data_xhi[itime].resize(num_var);
            bdy_data_ylo[itime].resize(num_var);
            bdy_data_yhi[itime].resize(num_var);
            for (int ivar(0); ivar<num_var; ++ivar) {
              bdy_data_xlo[itime][ivar].resize(bx_v[4*ivar  ]);
              bdy_data_xhi[itime][ivar].resize(bx_v[4*ivar+1]);
              bdy_data_ylo[itime][ivar].resize(bx_v[4*ivar+2]);
              bdy_data_yhi[itime][ivar].resize(bx_v[4*ivar+3]);
            }
          }
        }

        for (int itime(0); itime<num_time; ++itime) {
            for (int ivar(0); ivar<num_var; ++ivar) {
                ParallelDescriptor::Bcast(bdy_data_xlo[itime][ivar].dataPtr(),bdy_data_xlo[itime][ivar].box().numPts(),ioproc);
                ParallelDescriptor::Bcast(bdy_data_xhi[itime][ivar].dataPtr(),bdy_data_xhi[itime][ivar].box().numPts(),ioproc);
                ParallelDescriptor::Bcast(bdy_data_ylo[itime][ivar].dataPtr(),bdy_data_ylo[itime][ivar].box().numPts(),ioproc);
                ParallelDescriptor::Bcast(bdy_data_yhi[itime][ivar].dataPtr(),bdy_data_yhi[itime][ivar].box().numPts(),ioproc);
            }
        }
    } // init real
#endif
}
