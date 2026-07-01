    Real time = remora.get_t_old(lev);

    //If we had wind stress and bottom stress we would need to set these:
    Real tdays=time/Real(24.0*60.0*60.0);
    Real dstart=0.0_rt;
    Real userval=0.01_rt;
    Real mf_fcor=1.0e-4_rt;
    //It's possible these should be set to be nonzero only at the boundaries they affect

    // Don't allow doubly periodic in this case

    Real windamp_u = userval*std::sin(0.9_rt*mf_fcor*86400.0_rt*(tdays-dstart))/m_solverChoice.rho0;
    mf_sustr.setVal(windamp_u);
    mf_svstr.setVal(0.0_rt);
