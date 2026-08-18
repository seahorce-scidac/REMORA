
    // Whether this run carries a passive dye tracer. nComp() cannot answer this: with
    // biology active and no dye, Tracer_comp is the first biology component, so testing
    // nComp() would write the dye initial value into NO3.
    const bool l_do_dye = (remora.num_passive_scalars() > 0);
    ParmParse pp_prob("remora.prob");
    Real u_0     = Real(0.0); pp_prob.query("u_0", u_0);
    Real v_0     = Real(0.0); pp_prob.query("v_0", v_0);
    Real z0      = Real(0.1); pp_prob.query("z0", z0);       // Surface Roughness
    Real zRef    = Real(80.); pp_prob.query("zRef", zRef);  // Reference Height
    Real uRef    = Real(0.0); pp_prob.query("uRef", uRef);   // Reference Wind Speed

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

            // Set tracer = 0 everywhere
            const Real xcent = Real(0.5)*(prob_lo[0] + prob_hi[0]);
            const Real ycent = Real(0.5)*(prob_lo[1] + prob_hi[1]);

            const Real x  = prob_lo[0] + (i + Real(0.5)) * dx[0] - xcent;
            const Real y  = prob_lo[1] + (j + Real(0.5)) * dx[1] - ycent;
            const Real r2 = x*x + y*y;
            const Real rad = Real(0.1) * (prob_hi[0]-prob_lo[0]);
            const Real radsq = rad*rad;

            if (l_do_dye) { state(i, j, k,Tracer_comp) = (r2 < radsq) ? one : zero; }
        });

        // Construct a box that is on x-faces
        const Box& xbx = surroundingNodes(bx,0);
        // Set the x-velocity
        ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
              // const auto prob_lo         = geomdata.ProbLo();
              // const auto dx              = geomdata.CellSize();

              // const Real x = prob_lo[0] + (i + 0.5) * dx[0];
              // const Real y = prob_lo[1] + (j + 0.5) * dx[1];
              const Real z = -z_r(i,j,k);

              // Set the x-velocity
              x_vel(i, j, k) = u_0 + uRef * std::log((z + z0)/z0) / std::log((zRef +z0)/z0);
        });

        // Construct a box that is on y-faces
        const Box& ybx = surroundingNodes(bx,1);

        // Set the y-velocity
        ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
              //const auto prob_lo         = geomdata.ProbLo();
              //const auto dx              = geomdata.CellSize();

              // const Real x = prob_lo[0] + (i + 0.5) * dx[0];
              // const Real y = prob_lo[1] + (j + 0.5) * dx[1];
              y_vel(i, j, k) = v_0; });
    }
    Gpu::streamSynchronize();
