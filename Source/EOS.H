#ifndef _EOS_H_
#define _EOS_H_

#include <ROMSX_Constants.H>
#include <AMReX.H>
#include <cmath>

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real getTgivenRandTh(const amrex::Real rho, const amrex::Real temp)
{
    amrex::Real p_loc = p_0 * std::pow(R_d * temp * ip_0, Gamma);
    return p_loc / (R_d * rho);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real getThgivenRandT(const amrex::Real rho, const amrex::Real T)
{
    amrex::Real p_loc = rho * R_d * T;
    return T * std::pow((p_0/p_loc),rdOcp);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real getRhoGivenThetaPress (const amrex::Real theta, const amrex::Real p) {
    return std::pow(p_0, rdOcp) * std::pow(p, iGamma) / (R_d * theta);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real getdPdRgivenConstantTheta(const amrex::Real rho, const amrex::Real theta)
{
    return Gamma * p_0 * std::pow( (R_d * theta * ip_0), Gamma) * std::pow(rho, Gamma-1.0) ;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real getExnergivenP(const amrex::Real P)
{
    // Exner function pi in terms of P
    return std::pow(P * ip_0, rdOcp);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real getTempgivenP(const amrex::Real p)
{
    // diagnostic relation for the full pressure
    return std::pow(p*std::pow(p_0, Gamma-1), iGamma) * iR_d;
}
#endif

