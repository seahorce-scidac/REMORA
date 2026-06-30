    mf_h.setVal(geom.ProbHi(2));

    for ( MFIter mfi(mf_h, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {

      Array4<Real> const& h  = (mf_h).array(mfi);

      Box bx = mfi.tilebox();
      Box gbx2 = grow(bx,IntVect(NGROW+1,NGROW+1,0));

      amrex::Real Xsize = 320000.0_rt;
      amrex::Real Esize = 320000.0_rt;
      amrex::Real depth = 5000.0_rt;

      if(!m_solverChoice.flat_bathymetry)
      {
          ParallelFor(makeSlab(gbx2,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int)
          {
              const auto prob_lo         = geomdata.ProbLo();
              const auto dx              = geomdata.CellSize();

              const Real x = prob_lo[0] + (i + 0.5_rt) * dx[0];
              const Real y = prob_lo[1] + (j + 0.5_rt) * dx[1];

              Real val1, val2;
              val1 = (x-0.5_rt*Xsize)/40000.0_rt;
              val2 = (y-0.5_rt*Esize)/40000.0_rt;
              h(i,j,0) = depth - 4500.0_rt * std::exp(-(val1*val1+val2*val2));
          });

      } else {

          ParallelFor(makeSlab(gbx2,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int )
          {
              h(i,j,0) = -geomdata.ProbLo(2);
              h(i,j,0,1) = h(i,j,0);
          });
      }
    }
