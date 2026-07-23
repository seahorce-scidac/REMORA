    const auto prob_lo         = geomdata.ProbLo();
    const auto prob_hi         = geomdata.ProbHi();
    const auto dx              = geomdata.CellSize();

    //If we had wind stress and bottom stress we would need to set these:
    Real pi = Real(3.14159265359);

    const Real yextent = prob_hi[1] - prob_lo[1];
    const Real windamp = Real(-0.05) / m_solverChoice.rho0;
    const Real val1  = two * pi / yextent;
    const Real prob_lo_y = prob_lo[1];
    const Real dy = dx[1];
    for ( MFIter mfi((mf_sustr), TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        const Box& bx = mfi.tilebox();
        [[maybe_unused]] const Box& xbx = surroundingNodes(bx,0);
        const Box& xbx2 = mfi.grownnodaltilebox(0, IntVect(NGROW,NGROW,0));

        Array4<Real> const& sustr = mf_sustr.array(mfi);
        ParallelFor(xbx2, [=] AMREX_GPU_DEVICE(int i, int j, int /*k*/) noexcept
        {
            const Real y  = prob_lo_y + (j + half) * dy;// - ycent;

            sustr(i,j,0) = windamp * std::cos(val1 * y);
        });
    }
    mf_svstr.setVal(zero);
