#include <REMORA.H>
#include "AMReX_PlotFileUtil.H"

using namespace amrex;

// utility to skip to next line in Header
void
REMORA::GotoNextLine (std::istream& is)
{
    constexpr std::streamsize bl_ignore_max { 100000 };
    is.ignore(bl_ignore_max, '\n');
}

void
REMORA::WriteCheckpointFile ()
{

    // chk00010            write a checkpoint file with this root directory
    // chk00010/Header     this contains information you need to save (e.g., finest_level, t_new, etc.) and also
    //                     the BoxArrays at each level
    // chk00010/Level_0/
    // chk00010/Level_1/
    // etc.                these subdirectories will hold the MultiFab data at each level of refinement

    // checkpoint file name, e.g., chk00010
    const std::string& checkpointname = amrex::Concatenate(check_file,istep[0],file_min_digits);

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
       HeaderFile << "Checkpoint file for REMORA\n";

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

       MultiFab cons(grids[lev],dmap[lev],NCONS,cons_new[lev]->nGrowVect());
       MultiFab::Copy(cons,*cons_new[lev],0,0,NCONS,cons_new[lev]->nGrowVect());
       VisMF::Write(cons, amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "Cell"));

       MultiFab::Copy(cons,*cons_old[lev],0,0,NCONS,cons_old[lev]->nGrowVect());
       VisMF::Write(cons, amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "Cell_old"));

       MultiFab xvel(convert(grids[lev],IntVect(1,0,0)),dmap[lev],1,xvel_new[lev]->nGrowVect());
       MultiFab::Copy(xvel,*xvel_new[lev],0,0,1,xvel_new[lev]->nGrowVect());
       VisMF::Write(xvel, amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "XFace"));

       MultiFab::Copy(xvel,*xvel_old[lev],0,0,1,xvel_old[lev]->nGrowVect());
       VisMF::Write(xvel, amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "XFace_old"));

       MultiFab yvel(convert(grids[lev],IntVect(0,1,0)),dmap[lev],1,yvel_new[lev]->nGrowVect());
       MultiFab::Copy(yvel,*yvel_new[lev],0,0,1,yvel_new[lev]->nGrowVect());
       VisMF::Write(yvel, amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "YFace"));

       MultiFab::Copy(yvel,*yvel_old[lev],0,0,1,yvel_old[lev]->nGrowVect());
       VisMF::Write(yvel, amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "YFace_old"));

       MultiFab zvel(convert(grids[lev],IntVect(0,0,1)),dmap[lev],1,zvel_new[lev]->nGrowVect());
       MultiFab::Copy(zvel,*zvel_new[lev],0,0,1,zvel_new[lev]->nGrowVect());
       VisMF::Write(zvel, amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "ZFace"));

       MultiFab::Copy(zvel,*zvel_old[lev],0,0,1,zvel_old[lev]->nGrowVect());
       VisMF::Write(zvel, amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "ZFace_old"));

       MultiFab mf_ru(convert(grids[lev],IntVect(1,0,0)),dmap[lev],2,(vec_ru[lev])->nGrowVect());
       MultiFab::Copy(mf_ru,*vec_ru[lev],0,0,2,(vec_ru[lev])->nGrowVect());
       VisMF::Write(mf_ru, amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "XRHS"));

       MultiFab mf_rv(convert(grids[lev],IntVect(0,1,0)),dmap[lev],2,(vec_rv[lev])->nGrowVect());
       MultiFab::Copy(mf_rv,*vec_rv[lev],0,0,2,(vec_rv[lev])->nGrowVect());
       VisMF::Write(mf_rv, amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "YRHS"));

       MultiFab mf_ubar(convert(ba2d,IntVect(1,0,0)),dmap[lev],3,(vec_ubar[lev])->nGrowVect());
       MultiFab::Copy(mf_ubar,*(vec_ubar[lev]),0,0,3,(vec_ubar[lev])->nGrowVect());
       VisMF::Write(mf_ubar, amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "XBar"));

       MultiFab mf_vbar(convert(ba2d,IntVect(0,1,0)),dmap[lev],3,(vec_vbar[lev])->nGrowVect());
       MultiFab::Copy(mf_vbar,*(vec_vbar[lev]),0,0,3,(vec_vbar[lev])->nGrowVect());
       VisMF::Write(mf_vbar, amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "YBar"));

       MultiFab mf_ru2d(convert(ba2d,IntVect(1,0,0)),dmap[lev],2,(vec_ru2d[lev])->nGrowVect());
       MultiFab::Copy(mf_ru2d,*(vec_ru2d[lev]),0,0,2,(vec_ru2d[lev])->nGrowVect());
       VisMF::Write(mf_ru2d, amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "XRHS2d"));

       MultiFab mf_rv2d(convert(ba2d,IntVect(0,1,0)),dmap[lev],2,(vec_rv2d[lev])->nGrowVect());
       MultiFab::Copy(mf_rv2d,*(vec_rv2d[lev]),0,0,2,(vec_rv2d[lev])->nGrowVect());
       VisMF::Write(mf_rv2d, amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "YRHS2d"));

       MultiFab mf_mskr(ba2d,dmap[lev],1,(vec_mskr[lev])->nGrowVect());
       MultiFab::Copy(mf_mskr,*(vec_mskr[lev]),0,0,1,(vec_mskr[lev])->nGrowVect());
       VisMF::Write(mf_mskr, amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "Mskr"));

       MultiFab mf_msku(convert(ba2d,IntVect(1,0,0)),dmap[lev],1,(vec_msku[lev])->nGrowVect());
       MultiFab::Copy(mf_msku,*(vec_msku[lev]),0,0,1,(vec_msku[lev])->nGrowVect());
       VisMF::Write(mf_msku, amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "Msku"));

       MultiFab mf_mskv(convert(ba2d,IntVect(0,1,0)),dmap[lev],1,(vec_mskv[lev])->nGrowVect());
       MultiFab::Copy(mf_mskv,*(vec_mskv[lev]),0,0,1,(vec_mskv[lev])->nGrowVect());
       VisMF::Write(mf_mskv, amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "Mskv"));

       VisMF::Write(*(vec_rdrag[lev]), amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "rdrag"));

       VisMF::Write(*(vec_DU_avg1[lev]), amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "DU_avg1"));
       VisMF::Write(*(vec_DU_avg2[lev]), amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "DU_avg2"));
       VisMF::Write(*(vec_DV_avg1[lev]), amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "DV_avg1"));
       VisMF::Write(*(vec_DV_avg2[lev]), amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "DV_avg2"));

       VisMF::Write(*(vec_zeta[lev]), amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "zeta"));
       VisMF::Write(*(vec_Zt_avg1[lev]), amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "Zt_avg1"));

       VisMF::Write(*(vec_hOfTheConfusingName[lev]), amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "h"));

       VisMF::Write(*(vec_tke[lev]), amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "tke"));
       VisMF::Write(*(vec_gls[lev]), amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "gls"));
       VisMF::Write(*(vec_Lscale[lev]), amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "Lscale"));
       VisMF::Write(*(vec_Akk[lev]), amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "Akk"));
       VisMF::Write(*(vec_Akp[lev]), amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "Akp"));
       VisMF::Write(*(vec_Akv[lev]), amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "Akv"));
       VisMF::Write(*(vec_Akt[lev]), amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "Akt"));
   }

#ifdef REMORA_USE_PARTICLES
   particleData.Checkpoint(checkpointname);
#endif

#ifdef REMORA_USE_NETCDF
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
REMORA::ReadCheckpointFile ()
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
#ifdef AMREX_USE_FLOAT
            dt[i++] = std::stof(word);
#else
            dt[i++] = std::stod(word);
#endif
        }
    }

    // read in array of t_new
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
#ifdef AMREX_USE_FLOAT
            t_new[i++] = std::stof(word);
#else
            t_new[i++] = std::stod(word);
#endif
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

        MultiFab cons(grids[lev],dmap[lev],NCONS,cons_new[lev]->nGrowVect());
        VisMF::Read(cons, amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "Cell"));
        MultiFab::Copy(*cons_new[lev],cons,0,0,NCONS,cons_new[lev]->nGrowVect());

        VisMF::Read(cons, amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "Cell_old"));
        MultiFab::Copy(*cons_old[lev],cons,0,0,NCONS,cons_old[lev]->nGrowVect());

        MultiFab xvel(convert(grids[lev],IntVect(1,0,0)),dmap[lev],1,xvel_new[lev]->nGrowVect());
        VisMF::Read(xvel, amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "XFace"));
        MultiFab::Copy(*xvel_new[lev],xvel,0,0,1,xvel_new[lev]->nGrowVect());

        VisMF::Read(xvel, amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "XFace_old"));
        MultiFab::Copy(*xvel_old[lev],xvel,0,0,1,xvel_old[lev]->nGrowVect());

        MultiFab yvel(convert(grids[lev],IntVect(0,1,0)),dmap[lev],1,yvel_new[lev]->nGrowVect());
        VisMF::Read(yvel, amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "YFace"));
        MultiFab::Copy(*yvel_new[lev],yvel,0,0,1,yvel_new[lev]->nGrowVect());

        VisMF::Read(yvel, amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "YFace_old"));
        MultiFab::Copy(*yvel_old[lev],yvel,0,0,1,yvel_old[lev]->nGrowVect());

        MultiFab zvel(convert(grids[lev],IntVect(0,0,1)),dmap[lev],1,zvel_new[lev]->nGrowVect());
        VisMF::Read(zvel, amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "ZFace"));
        MultiFab::Copy(*zvel_new[lev],zvel,0,0,1,zvel_new[lev]->nGrowVect());

        VisMF::Read(zvel, amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "ZFace_old"));
        MultiFab::Copy(*zvel_old[lev],zvel,0,0,1,zvel_old[lev]->nGrowVect());

       MultiFab mf_ru(convert(grids[lev],IntVect(1,0,0)),dmap[lev],2,(vec_ru[lev])->nGrowVect());
       VisMF::Read(mf_ru, amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "XRHS"));
       MultiFab::Copy(*(vec_ru[lev]),mf_ru,0,0,2,(vec_ru[lev])->nGrowVect());

       MultiFab mf_rv(convert(grids[lev],IntVect(0,1,0)),dmap[lev],2,(vec_rv[lev])->nGrowVect());
       VisMF::Read(mf_rv, amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "YRHS"));
       MultiFab::Copy(*(vec_rv[lev]),mf_rv,0,0,2,(vec_rv[lev])->nGrowVect());

       MultiFab mf_ubar(convert(ba2d,IntVect(1,0,0)),dmap[lev],3,(vec_ubar[lev])->nGrowVect());
       VisMF::Read(mf_ubar, amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "XBar"));
       MultiFab::Copy(*(vec_ubar[lev]),mf_ubar,0,0,3,(vec_ubar[lev])->nGrowVect());

       MultiFab mf_vbar(convert(ba2d,IntVect(0,1,0)),dmap[lev],3,(vec_vbar[lev])->nGrowVect());
       VisMF::Read(mf_vbar, amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "YBar"));
       MultiFab::Copy(*(vec_vbar[lev]),mf_vbar,0,0,3,(vec_vbar[lev])->nGrowVect());

       MultiFab mf_ru2d(convert(ba2d,IntVect(1,0,0)),dmap[lev],2,(vec_ru2d[lev])->nGrowVect());
       VisMF::Read(mf_ru2d, amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "XRHS2d"));
       MultiFab::Copy(*(vec_ru2d[lev]),mf_ru2d,0,0,2,(vec_ru2d[lev])->nGrowVect());

       MultiFab mf_rv2d(convert(ba2d,IntVect(0,1,0)),dmap[lev],2,(vec_rv2d[lev])->nGrowVect());
       VisMF::Read(mf_rv2d, amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "YRHS2d"));
       MultiFab::Copy(*(vec_rv2d[lev]),mf_rv2d,0,0,2,(vec_rv2d[lev])->nGrowVect());

       MultiFab mf_mskr(ba2d,dmap[lev],1,(vec_mskr[lev])->nGrowVect());
       VisMF::Read(mf_mskr, amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "Mskr"));
       MultiFab::Copy(*(vec_mskr[lev]),mf_mskr,0,0,1,(vec_mskr[lev])->nGrowVect());

       MultiFab mf_msku(convert(ba2d,IntVect(1,0,0)),dmap[lev],1,(vec_msku[lev])->nGrowVect());
       VisMF::Read(mf_msku, amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "Msku"));
       MultiFab::Copy(*(vec_msku[lev]),mf_msku,0,0,1,(vec_msku[lev])->nGrowVect());

       MultiFab mf_mskv(convert(ba2d,IntVect(0,1,0)),dmap[lev],1,(vec_mskv[lev])->nGrowVect());
       VisMF::Read(mf_mskv, amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "Mskv"));
       MultiFab::Copy(*(vec_mskv[lev]),mf_mskv,0,0,1,(vec_mskv[lev])->nGrowVect());

       VisMF::Read(*(vec_rdrag[lev]), amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "rdrag"));

       VisMF::Read(*(vec_DU_avg1[lev]), amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "DU_avg1"));
       VisMF::Read(*(vec_DU_avg2[lev]), amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "DU_avg2"));
       VisMF::Read(*(vec_DV_avg1[lev]), amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "DV_avg1"));
       VisMF::Read(*(vec_DV_avg2[lev]), amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "DV_avg2"));

       VisMF::Read(*(vec_zeta[lev]), amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "zeta"));
       VisMF::Read(*(vec_Zt_avg1[lev]), amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "Zt_avg1"));

       VisMF::Read(*(vec_hOfTheConfusingName[lev]), amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "h"));

       VisMF::Read(*(vec_tke[lev]), amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "tke"));
       VisMF::Read(*(vec_gls[lev]), amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "gls"));
       VisMF::Read(*(vec_Lscale[lev]), amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "Lscale"));
       VisMF::Read(*(vec_Akk[lev]), amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "Akk"));
       VisMF::Read(*(vec_Akp[lev]), amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "Akp"));
       VisMF::Read(*(vec_Akt[lev]), amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "Akt"));
       VisMF::Read(*(vec_Akv[lev]), amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "Akv"));

       stretch_transform(lev);

    }

#ifdef REMORA_USE_PARTICLES
   particleData.Restart((amrex::ParGDBBase*)GetParGDB(),restart_chkfile);
#endif

#ifdef REMORA_USE_NETCDF
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
