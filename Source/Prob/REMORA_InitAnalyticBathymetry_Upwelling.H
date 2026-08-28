
    // Must not be doubly periodic, and must have terrain
    AMREX_ALWAYS_ASSERT( !NSPeriodic || !EWPeriodic);

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

              Real val1 = (iFort <= Lm/two) ? iFort : Lm+1-iFort;
              val1 -= Real(0.5);
              Real adj = geomdata.CellSize()[0]/Real(1000.0);

              h(i,j,0) = std::min(-geomdata.ProbLo(2),(Real(84.5)+Real(66.526)*std::tanh((val1*adj-Real(10.0))/Real(7.0))));
          });

      } else if (EWPeriodic) {

          ParallelFor(gbxD, [=] AMREX_GPU_DEVICE (int i, int j, int )
          {
              int jFort = j+1; // (+1 is to match the Fortran indexing in ROMS)

              Real val1 = (jFort<=Mm/two) ? jFort : Mm+1-jFort;
              val1 -= Real(0.5);
              Real adj = geomdata.CellSize()[1]/Real(1000.0);

              h(i,j,0) = std::min(-geomdata.ProbLo(2),(Real(84.5)+Real(66.526)*std::tanh((val1*adj-Real(10.0))/Real(7.0))));
          });
      }
    } // mfi
