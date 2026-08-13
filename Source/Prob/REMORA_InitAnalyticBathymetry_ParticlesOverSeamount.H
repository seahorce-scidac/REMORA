
    mf_h.setVal(geom.ProbHi(2));
    const int Lm = geom.Domain().size()[0];
    const int Mm = geom.Domain().size()[1];

    for ( MFIter mfi(mf_h, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {

      Array4<Real> const& h  = (mf_h).array(mfi);

      Box bx = mfi.tilebox();
      Box gbx2 = bx;
      gbx2.grow(IntVect(NGROW+1,NGROW+1,0));

      Gpu::streamSynchronize();

      if (!m_solverChoice.flat_bathymetry) {
          ParallelFor(makeSlab(gbx2,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int )
          {
              Real val1;
              int iFort = i+1;
              int jFort = j+1;
              if(NSPeriodic) {
                  if (iFort<=Lm/two)
                      val1=iFort;
                  else
                      val1=Lm+1-iFort;
              h(i,j,0) = std::min(-geomdata.ProbLo(2),(Real(84.5)+Real(66.526)*std::tanh((val1-Real(10.0))/Real(7.0))));

              } else if(EWPeriodic) {

                  if (jFort<=Mm/2.0)
                      val1=jFort;
                  else
                      val1=Mm+1-jFort;
              h(i,j,0) = std::min(-geomdata.ProbLo(2),(Real(84.5)+Real(66.526)*std::tanh((val1-Real(10.0))/Real(7.0))));
              }
          });
      } else {

          ParallelFor(gbx2, [=] AMREX_GPU_DEVICE (int i, int j, int k)
          {
              h(i,j,0) = -geomdata.ProbLo(2);
              if (k==0) {
                  h(i,j,0,1) = h(i,j,0);
              }
          });
      }
    } // mfi
