    const auto prob_lo         = geomdata.ProbLo();
    const auto prob_hi         = geomdata.ProbHi();

    //If we had wind stress and bottom stress we would need to set these:
    Real pi = 3.14159265359_rt;

    const Real yextent = prob_hi[1] - prob_lo[1];
    const Real windamp = -0.05_rt / m_solverChoice.rho0;
    const Real val1  = 2.0_rt * pi / yextent;
    for ( MFIter mfi((mf_sustr), TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        const Box& bx = mfi.tilebox();
        [[maybe_unused]] const Box& xbx = surroundingNodes(bx,0);
        const Box& xbx2 = mfi.grownnodaltilebox(0, IntVect(NGROW,NGROW,0));

        Array4<Real> const& sustr = mf_sustr.array(mfi);
        ParallelFor(xbx2, [=] AMREX_GPU_DEVICE(int i, int j, int /*k*/) noexcept
        {
            const auto dx = geomdata.CellSize();

            const Real y  = prob_lo[1] + (j + 0.5) * dx[1];// - ycent;

            sustr(i,j,0) = windamp * std::cos(val1 * y);
        });
    }
    mf_svstr.setVal(0.0_rt);
