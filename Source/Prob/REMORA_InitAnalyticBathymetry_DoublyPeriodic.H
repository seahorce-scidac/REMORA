
    mf_h.setVal(geomdata.ProbHi(2));

    const int Lm = geom.Domain().size()[0];
    const int Mm = geom.Domain().size()[1];

    for ( MFIter mfi(mf_h, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
      Array4<Real> const& h  = (mf_h).array(mfi);

      Box bx = mfi.tilebox();
      Box gbx2 = bx;
      // See the note in REMORA_InitAnalyticBathymetry_Seamount.H
      gbx2.grow(mf_h.nGrowVect());

      // auto N = geom.Domain().length(2); // Number of vertical "levels" aka, NZ

      Box gbx2D = gbx2;
      gbx2D.makeSlab(2,0);

      Gpu::streamSynchronize();
      ParallelFor(gbx2D, [=] AMREX_GPU_DEVICE (int i, int j, int )
      {
          Real val1;
          [[maybe_unused]] Real val2;
          int iFort = i+1;
          int jFort = j+1;
          if (NSPeriodic) {
              if (iFort<=Lm/two) {
                  val1=iFort;
              } else {
                  val1=Lm+1-iFort;
              }
          h(i,j,0) = std::min(-geomdata.ProbLo(2),(Real(84.5)+Real(66.526)*std::tanh((val1-Real(10.0))/Real(7.0))));
          }
          else if(EWPeriodic) {
              if (jFort<=Mm/two) {
                  val1=jFort;
              } else {
                  val1=Mm+1-jFort;
              }
          h(i,j,0) = std::min(-geomdata.ProbLo(2),(Real(84.5)+Real(66.526)*std::tanh((val1-Real(10.0))/Real(7.0))));
          }
      });
    } // mfi
