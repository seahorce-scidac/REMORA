    auto T0 = m_solverChoice.T0;
    auto S0 = m_solverChoice.S0;

    mf_cons.setVal(T0,Temp_comp,1);
    mf_cons.setVal(S0,Salt_comp,1);
    // Only a passive dye lives at Tracer_comp; with biology and no dye that slot is NO3,
    // which init_biology_ic owns.
    if (remora.num_passive_scalars() > 0) { mf_cons.setVal(zero,Tracer_comp,1); }
    mf_xvel.setVal(zero);
    mf_yvel.setVal(zero);
