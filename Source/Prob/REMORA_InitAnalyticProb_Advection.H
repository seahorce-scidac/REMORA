
    ParmParse pp_prob("remora.prob");
    Real u_0 = 0.0; pp_prob.query("u_0", u_0);
    Real v_0 = 0.0; pp_prob.query("v_0", v_0);

    bool l_use_salt = m_solverChoice.use_salt;

    auto geomdata = geom.data();
    const int khi = geomdata.Domain().bigEnd()[2];

    auto T0 = m_solverChoice.T0;
    auto S0 = m_solverChoice.S0;

    for (MFIter mfi(mf_cons, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box &bx = mfi.tilebox();
        AMREX_ALWAYS_ASSERT(bx.length()[2] == khi+1);

        Array4<      Real> const& state = mf_cons.array(mfi);
        Array4<      Real> const& x_vel = mf_xvel.array(mfi);
        Array4<      Real> const& y_vel = mf_yvel.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            const auto prob_lo         = geomdata.ProbLo();
            const auto prob_hi         = geomdata.ProbHi();
            const auto dx              = geomdata.CellSize();

            state(i, j, k, Temp_comp) = 1.;

            state(i,j,k,Temp_comp)=T0; //+8.0*std::exp(z/50.0_rt);

            // Set tracer = 0 everywhere
            const Real xcent = 0.5*(prob_lo[0] + prob_hi[0]);
            const Real ycent = 0.5*(prob_lo[1] + prob_hi[1]);

            const Real x  = prob_lo[0] + (i + 0.5) * dx[0] - xcent;
            const Real y  = prob_lo[1] + (j + 0.5) * dx[1] - ycent;
            const Real r2 = x*x + y*y;
            const Real rad = 0.1 * (prob_hi[0]-prob_lo[0]);
            const Real radsq = rad*rad;
            const Real rad_inner = 0.05 * (prob_hi[0]-prob_lo[0]);
            [[maybe_unused]] const Real rad_inner_sq = rad_inner*rad_inner;

            if (l_use_salt) {
                state(i,j,k,Salt_comp)= S0;
            }

            // Single circle of tracer (default)
            state(i, j, k, Tracer_comp) = std::exp(-r2/(2.*radsq));

            // Donut of tracer
            //state(i, j, k, Tracer_comp) = 1.25 * (std::exp(-r2/(2.*radsq)) - std::exp(-r2/(2*rad_inner_sq)));
        });

        // Construct a box that is on x-faces
        const Box& xbx = surroundingNodes(bx,0);
        // Set the x-velocity
        ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
              x_vel(i, j, k) = u_0;
        });

        // Construct a box that is on y-faces
        const Box& ybx = surroundingNodes(bx,1);

        // Set the y-velocity
        ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
              y_vel(i, j, k) = v_0;
        });
    }
    Gpu::streamSynchronize();
