    ParmParse pp_prob("remora.prob");

    bool traditional = true;
    pp_prob.query("traditional", traditional);

    // Flat-bottom depth for the traditional dogbone. Settable so a run can be matched to a
    // ROMS grid file with a different h; set remora.prob_lo's z to -depth to match.
    Real depth = Real(10.0);
    pp_prob.query("depth", depth);

    if (traditional) {
        mf_h.setVal(depth);
    } else {
        // Coordinates from the Geometry rather than remora.vec_xr/vec_yr, as in the mask hook:
        // the hires_grid_level path evaluates this on a full-domain MultiFab that is not built
        // on grids[lev], so the per-level coordinate arrays cannot be indexed with its MFIter.
        // Matches set_grid_coords_from_grid_scale's (i + 1/2) * dx convention.
        const Real dx = geom.CellSize(0);
        const Real dy = geom.CellSize(1);

        mf_h.setVal(zero);
        for (MFIter mfi(mf_h, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box &bx = mfi.growntilebox();

            Array4<      Real> const& h = mf_h.array(mfi);

            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int ) noexcept
            {
                const Real x_r = (i + Real(0.5)) * dx;
                const Real y_r = (j + Real(0.5)) * dy;
                h(i,j,0) = Real(10.0) + (x_r - Real(4200.0)) * (x_r - Real(4200.0)) * Real(5.0) / (Real(4200.0) * Real(4200.0)) -
                            (y_r - Real(375.0)) * (y_r - Real(375.0)) * Real(5.0) / (Real(375.0) * Real(375.0));
            });
        }
    }
