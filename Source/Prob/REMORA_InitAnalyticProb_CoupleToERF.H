    auto T0 = m_solverChoice.T0;
    auto S0 = m_solverChoice.S0;

    mf_cons.setVal(T0,Temp_comp,1);
    mf_cons.setVal(S0,Salt_comp,1);
    mf_cons.setVal(0.0_rt,Tracer_comp,1);
    mf_xvel.setVal(0.0_rt);
    mf_yvel.setVal(0.0_rt);
