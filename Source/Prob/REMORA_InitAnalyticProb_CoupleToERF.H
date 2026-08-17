    auto T0 = m_solverChoice.T0;
    auto S0 = m_solverChoice.S0;

    mf_cons.setVal(T0,Temp_comp,1);
    mf_cons.setVal(S0,Salt_comp,1);
    if (mf_cons.nComp() > Tracer_comp) { mf_cons.setVal(zero,Tracer_comp,1); }
    mf_xvel.setVal(zero);
    mf_yvel.setVal(zero);
