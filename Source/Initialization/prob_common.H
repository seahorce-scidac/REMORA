#ifndef _PROBCOMMON_H_
#define _PROBCOMMON_H_

#include <DataStruct.H>

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>

void amrex_probinit (
  const amrex_real* problo,
  const amrex_real* probhi);

void
init_custom_prob (
  const amrex::Box& bx,
  amrex::Array4<amrex::Real      > const& state,
  amrex::Array4<amrex::Real      > const& x_vel,
  amrex::Array4<amrex::Real      > const& y_vel,
  amrex::Array4<amrex::Real      > const& z_vel,
  amrex::Array4<amrex::Real const> const& z_w,
  amrex::Array4<amrex::Real const> const& z_r,
  amrex::Array4<amrex::Real const> const& Hz,
  amrex::Array4<amrex::Real const> const& h,
  amrex::Array4<amrex::Real const> const& Zt_avg1,
  amrex::GeometryData const& geomdata,
  SolverChoice const& m_solverChoice);

void init_custom_bathymetry (const amrex::Geometry& geom,
                        amrex::MultiFab& mf_h,
                        amrex::MultiFab& mf_Zt_avg1,
                        SolverChoice const& m_solverChoice);

void init_custom_vmix (const amrex::Geometry& geom,
                        amrex::MultiFab& mf_Akv,
                        amrex::MultiFab& mf_Akt,
                        amrex::MultiFab& mf_z_w,
                        SolverChoice const& m_solverChoice);

void init_custom_hmix (const amrex::Geometry& geom,
                        amrex::MultiFab& mf_visc2_p,
                        amrex::MultiFab& mf_visc2_r,
                        amrex::MultiFab& mf_diff2,
                        SolverChoice const& m_solverChoice);

void init_custom_smflux (const amrex::Geometry& geom,
                        amrex::Real time,
                        amrex::MultiFab& mf_sustr,
                        amrex::MultiFab& mf_svstr,
                        SolverChoice const& m_solverChoice);

#endif
