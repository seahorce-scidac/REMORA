
    mf_h.setVal(geomdata.ProbHi(2));

    const int Lm = geom.Domain().size()[0];
    const int Mm = geom.Domain().size()[1];

    for ( MFIter mfi(mf_h, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
      Array4<Real> const& h  = (mf_h).array(mfi);

      Box bx = mfi.tilebox();
      Box gbx2 = bx;
      gbx2.grow(IntVect(NGROW,NGROW,0));

      // auto N = geom.Domain().length(2); // Number of vertical "levels" aka, NZ

      Box gbx2D = gbx2;
      gbx2D.makeSlab(2,0);

      if(!m_solverChoice.flat_bathymetry) {
          Gpu::streamSynchronize();
          ParallelFor(gbx2D, [=] AMREX_GPU_DEVICE (int i, int j, int )
          {
              Real val1;
              [[maybe_unused]] Real val2;
              int iFort = i+1;
              int jFort = j+1;
              if (NSPeriodic) {
                  if (iFort<=Lm/2.0) {
                      val1=iFort;
                  } else {
                      val1=Lm+1-iFort;
                  }
              h(i,j,0) = std::min(-geomdata.ProbLo(2),(84.5_rt+66.526_rt*std::tanh((val1-10.0_rt)/7.0_rt)));
              }
              else if(EWPeriodic) {
                  if (jFort<=Mm/2.0) {
                      val1=jFort;
                  } else {
                      val1=Mm+1-jFort;
                  }
              h(i,j,0) = std::min(-geomdata.ProbLo(2),(84.5_rt+66.526_rt*std::tanh((val1-10.0_rt)/7.0_rt)));
              }
          });
      } else { // Flat
          Gpu::streamSynchronize();
          ParallelFor(gbx2, [=] AMREX_GPU_DEVICE (int i, int j, int k)
          {
              h(i,j,0,0) = -geomdata.ProbLo(2);
              if (k==0) {
                  h(i,j,0,1) = h(i,j,0,0);
              }
          });
      }
    } // mfi
