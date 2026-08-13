    //If we had wind stress and bottom stress we would need to set these:
    Real pi = Real(3.14159265359);
    Real tdays=remora.get_t_old(lev)/(Real(24.0)*Real(60.0)*Real(60.0));
    Real dstart=zero;
    Real windamp;
    //It's possible these should be set to be nonzero only at the boundaries they affect

    // Don't allow doubly periodic in this case
    //AMREX_ALWAYS_ASSERT( !NSPeriodic || !EWPeriodic);

    // Flow in x-direction (EW):
    if (EWPeriodic) {
        if ((tdays-dstart)<=two)
            windamp=Real(-0.1)*Real(sin(pi*(tdays-dstart)/Real(4.0)))/Real(m_solverChoice.rho0);
        else
            windamp=Real(-0.1)/m_solverChoice.rho0;
        mf_sustr.setVal(windamp);
    } else if (NSPeriodic) {
        mf_sustr.setVal(zero);
    }

    // Flow in y-direction (NS):
    if(EWPeriodic) {
        mf_svstr.setVal(zero);
    } else if (NSPeriodic) {
        if ((tdays-dstart)<=two)
            windamp=Real(-0.1)*Real(sin(pi*(tdays-dstart)/Real(4.0)))/Real(m_solverChoice.rho0);
        else
            windamp=Real(-0.1)/m_solverChoice.rho0;
        mf_svstr.setVal(windamp);
    }
    mf_sustr.setVal(zero);
    mf_svstr.setVal(zero);
