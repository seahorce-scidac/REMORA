
    int rrx = remora.refRatio(0)[0];
    int rry = remora.refRatio(0)[1];

    // Must not be doubly periodic, and must have terrain
    //AMREX_ALWAYS_ASSERT( !NSPeriodic || !EWPeriodic);
    AMREX_ALWAYS_ASSERT( !m_solverChoice.flat_bathymetry);

    mf_h.setVal(geomdata.ProbHi(2));

    if (lev==0 && m_solverChoice.init_l1ad_T) {
        const int Lm = geom.Domain().size()[0] * rrx;
        const int Mm = geom.Domain().size()[1] * rry;
        BoxArray old_ba(mf_h.boxArray());
        BoxArray tmp_ba = old_ba.refine(IntVect(rrx,rry,1));
        DistributionMapping tmp_dm(tmp_ba);
        MultiFab mf_h_hires(tmp_ba, tmp_dm, 2, mf_h.nGrow());
        mf_h_hires.setVal(geomdata.ProbHi(2));

        for ( MFIter mfi(mf_h_hires, TilingIfNotGPU()); mfi.isValid(); ++mfi )
        {
            Array4<Real> const& h  = (mf_h_hires).array(mfi);

            Box bx = mfi.tilebox();
            Box gbx3 = bx;
            gbx3.grow(IntVect(NGROW+1,NGROW+1,0));

            Box gbx3D = gbx3;
            gbx3D.makeSlab(2,0);

            Gpu::streamSynchronize();

            if (EWPeriodic) {

                ParallelFor(gbx3D, [=] AMREX_GPU_DEVICE (int i, int j, int )
                {
                    int jFort = j+1; // (+1 is to match the Fortran indexing in ROMS)

                    Real val1 = (jFort<=Mm/2.0) ? jFort : Mm+1-jFort;
                    val1 -= 0.5;
                    Real adj = geomdata.CellSize()[0]/(1000.0_rt * rrx);

                    h(i,j,0) = std::min(-geomdata.ProbLo(2),(84.5_rt+66.526_rt*std::tanh((val1*adj-10.0_rt)/7.0_rt)));
                });
            } else if (NSPeriodic) {

                ParallelFor(gbx3D, [=] AMREX_GPU_DEVICE (int i, int j, int )
                {
                    int iFort = i+1; // (+1 is to match the Fortran indexing in ROMS)

                    Real val1 = (iFort <= Lm/2.0) ? iFort : Lm+1-iFort;
                    val1 -= 0.5;
                    Real adj = geomdata.CellSize()[1]/(1000.0_rt * rry);

                    h(i,j,0) = std::min(-geomdata.ProbLo(2),(84.5_rt+66.526_rt*std::tanh((val1*adj-10.0_rt)/7.0_rt)));
                });
            }
        } // mfi
        average_down(mf_h_hires, mf_h, 0, mf_h.nComp(), IntVect(rrx,rry,1));
    } else {
        const int Lm = geom.Domain().size()[0];
        const int Mm = geom.Domain().size()[1];
        for ( MFIter mfi(mf_h, TilingIfNotGPU()); mfi.isValid(); ++mfi )
        {
            Array4<Real> const& h  = (mf_h).array(mfi);

            Box bx = mfi.tilebox();
            Box gbx3 = bx;
            gbx3.grow(IntVect(NGROW+1,NGROW+1,0));

            Box gbx3D = gbx3;
            gbx3D.makeSlab(2,0);

            Gpu::streamSynchronize();

            if (EWPeriodic) {

                ParallelFor(gbx3D, [=] AMREX_GPU_DEVICE (int i, int j, int )
                {
                    int jFort = j+1; // (+1 is to match the Fortran indexing in ROMS)

                    Real val1 = (jFort<=Mm/2.0) ? jFort : Mm+1-jFort;
                    val1 -= 0.5;
                    Real adj = geomdata.CellSize()[0]/1000.0_rt;

                    h(i,j,0) = std::min(-geomdata.ProbLo(2),(84.5_rt+66.526_rt*std::tanh((val1*adj-10.0_rt)/7.0_rt)));
                });
            } else if (NSPeriodic) {

                ParallelFor(gbx3D, [=] AMREX_GPU_DEVICE (int i, int j, int )
                {
                    int iFort = i+1; // (+1 is to match the Fortran indexing in ROMS)

                    Real val1 = (iFort <= Lm/2.0) ? iFort : Lm+1-iFort;
                    val1 -= 0.5;
                    Real adj = geomdata.CellSize()[1]/1000.0_rt;

                    h(i,j,0) = std::min(-geomdata.ProbLo(2),(84.5_rt+66.526_rt*std::tanh((val1*adj-10.0_rt)/7.0_rt)));
                });
            }
        } // mfi
    }
