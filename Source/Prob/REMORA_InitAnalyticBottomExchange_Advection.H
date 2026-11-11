    Real q = Real(0);
    pp.get("bottom_q", q);
    ParmParse pp_prob("remora.prob");
    Real q_amplitude = Real(0);
    Real tracer_amplitude = Real(0);
    Real ramp_time = Real(0);
    pp_prob.query("bottom_q_amplitude", q_amplitude);
    pp_prob.query("bottom_tracer_amplitude", tracer_amplitude);
    pp_prob.query("bottom_ramp_time", ramp_time);
    const Real temporal = ramp_time > Real(0) ? Real(1) + time / ramp_time : Real(1);
    auto& raw = q_and_concentrations;
    raw.setVal(q * temporal, 0, 1, raw.nGrowVect());
    const auto& names = remora.GetTracerNames();
    for (int n = 0; n < names.size(); ++n) {
        Real incoming = std::numeric_limits<Real>::quiet_NaN();
        pp.query(("bottom_" + names[n]).c_str(), incoming);
        raw.setVal(incoming, n + 1, 1, raw.nGrowVect());
    }
    const int dye_comp = remora.num_passive_scalars() > 0
                       ? remora.component_index("tracer") + 1 : -1;
    const auto lo = geom.ProbLoArray();
    const auto hi = geom.ProbHiArray();
    const auto dx = geom.CellSizeArray();
    const Real two_pi = Real(6.2831853071795864769);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(raw, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto a = raw.array(mfi);
        ParallelFor(mfi.growntilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            const Real x = (Real(i) + Real(0.5)) * dx[0] / (hi[0] - lo[0]);
            const Real y = (Real(j) + Real(0.5)) * dx[1] / (hi[1] - lo[1]);
            const Real mode = std::cos(two_pi * x) * std::cos(two_pi * y);
            a(i,j,k,0) *= Real(1) + q_amplitude * mode;
            if (dye_comp >= 0) {
                a(i,j,k,dye_comp) *= Real(1) + tracer_amplitude * mode;
            }
        });
    }
