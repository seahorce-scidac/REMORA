
    // Must not be doubly periodic, and must have terrain
    AMREX_ALWAYS_ASSERT( !NSPeriodic || !EWPeriodic);
    AMREX_ALWAYS_ASSERT( !m_solverChoice.flat_bathymetry);

    mf_h.setVal(geomdata.ProbHi(2));

    const int Lm = geom.Domain().size()[0];
    const int Mm = geom.Domain().size()[1];

    for ( MFIter mfi(mf_h, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
      Array4<Real> const& h  = (mf_h).array(mfi);

      Box bx = mfi.tilebox();
      Box gbx = bx;
      gbx.grow(mf_h.nGrowVect());

      Box gbxD = gbx;
      gbxD.makeSlab(2,0);

      Gpu::streamSynchronize();

      if (NSPeriodic) {

          ParallelFor(gbxD, [=] AMREX_GPU_DEVICE (int i, int j, int )
          {
              int iFort = i+1; // (+1 is to match the Fortran indexing in ROMS)

              Real val1 = (iFort <= Lm/2.0) ? iFort : Lm+1-iFort;
              val1 -= 0.5;
              Real adj = geomdata.CellSize()[1]/1000.0_rt;

              h(i,j,0) = std::min(-geomdata.ProbLo(2),(84.5_rt+66.526_rt*std::tanh((val1*adj-10.0_rt)/7.0_rt)));
          });

      } else if (EWPeriodic) {

          ParallelFor(gbxD, [=] AMREX_GPU_DEVICE (int i, int j, int )
          {
              int jFort = j+1; // (+1 is to match the Fortran indexing in ROMS)

              Real val1 = (jFort<=Mm/2.0) ? jFort : Mm+1-jFort;
              val1 -= 0.5;
              Real adj = geomdata.CellSize()[0]/1000.0_rt;

              h(i,j,0) = std::min(-geomdata.ProbLo(2),(84.5_rt+66.526_rt*std::tanh((val1*adj-10.0_rt)/7.0_rt)));
          });
      }
    } // mfi
