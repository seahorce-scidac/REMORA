    Real time = remora.get_t_old(lev);

    //If we had wind stress and bottom stress we would need to set these:
    Real tdays=time/(Real(24.0)*Real(60.0)*Real(60.0));
    Real dstart=zero;
    // Zero leaves the wind off, which a test needs if its exact solution is rest.
    amrex::ParmParse pp_prob("remora.prob");
    bool quiescent = false;
    pp_prob.query("quiescent", quiescent);

    Real userval = quiescent ? Real(0.0) : Real(0.01);
    Real mf_fcor=Real(1.0e-4);
    //It's possible these should be set to be nonzero only at the boundaries they affect

    // Don't allow doubly periodic in this case

    Real windamp_u = userval*std::sin(Real(0.9)*mf_fcor*Real(86400.0)*(tdays-dstart))/m_solverChoice.rho0;
    mf_sustr.setVal(windamp_u);
    mf_svstr.setVal(zero);
