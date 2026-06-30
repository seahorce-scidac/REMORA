    //If we had wind stress and bottom stress we would need to set these:
    Real pi = 3.14159265359_rt;
    Real tdays=remora.get_t_old(lev)/Real(24.0*60.0*60.0);
    Real dstart=0.0_rt;
    Real windamp;
    //It's possible these should be set to be nonzero only at the boundaries they affect

    // Don't allow doubly periodic in this case
    AMREX_ALWAYS_ASSERT( !NSPeriodic || !EWPeriodic);

    // Flow in x-direction (EW):
    if (NSPeriodic) {
        mf_sustr.setVal(0.0_rt);
    }
    else if (EWPeriodic) {
        if ((tdays-dstart)<=2.0)
            windamp=-0.1_rt*Real(sin(pi*(tdays-dstart)/4.0_rt))/Real(m_solverChoice.rho0);
        else
            windamp=-0.1_rt/m_solverChoice.rho0;
        mf_sustr.setVal(windamp);
    }

    // Flow in y-direction (NS):
    if (NSPeriodic) {
        if ((tdays-dstart)<=2.0)
            windamp=-0.1_rt*Real(sin(pi*(tdays-dstart)/4.0_rt))/Real(m_solverChoice.rho0);
        else
            windamp=-0.1_rt/m_solverChoice.rho0;
        mf_svstr.setVal(windamp);
    }
    else if(EWPeriodic) {
        mf_svstr.setVal(0.0_rt);
    }
