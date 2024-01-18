#ifndef _EOS_H_
#define _EOS_H_

#include <REMORA_Constants.H>
#include <AMReX.H>
#include <cmath>

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real getTgivenRandTh(const amrex::Real rho, const amrex::Real temp)
{
    amrex::Real p_loc = p_0 * std::pow(R_d * temp * ip_0, Gamma);
    return p_loc / (R_d * rho);
}
#endif

