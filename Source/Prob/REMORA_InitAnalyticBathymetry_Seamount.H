    mf_h.setVal(geom.ProbHi(2));

    for ( MFIter mfi(mf_h, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {

      Array4<Real> const& h  = (mf_h).array(mfi);

      Box bx = mfi.tilebox();
      // Grow to whatever the target actually has: with hires_grid_level this fab is
      // vec_h_full_domain, whose grow cells are max(cum_ref_ratios[lev], NGROW+1).
      Box gbx2 = grow(bx,mf_h.nGrowVect());

      amrex::Real Xsize = Real(320000.0);
      amrex::Real Esize = Real(320000.0);
      amrex::Real depth = Real(5000.0);

      if(!m_solverChoice.flat_bathymetry)
      {
          ParallelFor(makeSlab(gbx2,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int)
          {
              const auto prob_lo         = geomdata.ProbLo();
              const auto dx              = geomdata.CellSize();

              const Real x = prob_lo[0] + (i + Real(0.5)) * dx[0];
              const Real y = prob_lo[1] + (j + Real(0.5)) * dx[1];

              Real val1, val2;
              val1 = (x-Real(0.5)*Xsize)/Real(40000.0);
              val2 = (y-Real(0.5)*Esize)/Real(40000.0);
              h(i,j,0) = depth - Real(4500.0) * std::exp(-(val1*val1+val2*val2));
          });

      } else {

          ParallelFor(makeSlab(gbx2,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int )
          {
              h(i,j,0) = -geomdata.ProbLo(2);
              h(i,j,0,1) = h(i,j,0);
          });
      }
    }
