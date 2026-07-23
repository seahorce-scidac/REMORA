
    ParmParse pp_prob("remora.prob");
    Real u_0     = zero; pp_prob.query("u_0", u_0);
    Real v_0     = zero; pp_prob.query("v_0", v_0);
    Real z0      = Real(0.1); pp_prob.query("z0",  z0);        // Surface Roughness
    Real zRef    = Real(80.); pp_prob.query("zRef", zRef);     // Reference Height
    Real velRef  = zero; pp_prob.query("velRef", velRef); // Reference Wind Speed

    bool l_use_salt = m_solverChoice.use_salt;

    auto T0 = m_solverChoice.T0;
    auto S0 = m_solverChoice.S0;

    for (MFIter mfi(mf_cons, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box &bx = mfi.tilebox();

        auto geomdata = geom.data();

        Array4<      Real> const& state = mf_cons.array(mfi);
        Array4<      Real> const& x_vel = mf_xvel.array(mfi);
        Array4<      Real> const& y_vel = mf_yvel.array(mfi);

        Array4<const Real> const& z_r = remora.vec_z_r[lev]->const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            const auto prob_lo         = geomdata.ProbLo();
            const auto prob_hi         = geomdata.ProbHi();
            const auto dx              = geomdata.CellSize();
            const Real z = z_r(i,j,k);

            state(i, j, k, Temp_comp) = one;

            state(i,j,k,Temp_comp)=T0+Real(8.0)*std::exp(z/Real(50.0));
            if (l_use_salt) {
                state(i,j,k,Salt_comp)=S0;
            }

            const Real xcent = half*(prob_lo[0] + prob_hi[0]);
            const Real ycent = half*(prob_lo[1] + prob_hi[1]);

            const Real x  = prob_lo[0] + (i + half) * dx[0] - xcent;
            const Real y  = prob_lo[1] + (j + half) * dx[1] - ycent*half;
            const Real r2 = x*x + y*y;
            const Real rad = Real(0.1) * (prob_hi[0]-prob_lo[0]);
            const Real radsq = rad*rad;

            // Set tracer = 0 everywhere
            state(i, j, k, Tracer_comp) = std::exp(-r2/(two*radsq));
        });

        const Box& xbx = surroundingNodes(bx,0);
        const Box& ybx = surroundingNodes(bx,1);
        ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            const Real z = -z_r(i,j,k);
            x_vel(i, j, k) = u_0 + velRef * std::log((z + z0)/z0) / std::log((zRef +z0)/z0);
        });
        ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            y_vel(i, j, k) = zero;
        });
    }
