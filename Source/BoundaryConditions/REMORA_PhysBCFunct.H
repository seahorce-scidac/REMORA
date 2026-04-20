#ifndef REMORA_PhysBCFunct_H_
#define REMORA_PhysBCFunct_H_
#include <AMReX_Config.H>

#include <AMReX_BCRec.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_FilCC_C.H>
#include <AMReX_FilND_C.H>
#include <AMReX_FilFC_C.H>
#include "AMReX_TypeTraits.H"
#include "AMReX_Orientation.H"

#include <REMORA_TimeInterpolatedData.H>
#include <REMORA_IndexDefines.H>
#include <REMORA_DataStruct.H>

using PlaneVector = amrex::Vector<amrex::FArrayBox>;

class REMORAPhysBCFunct
{
public:
    /** \brief Constructor for physical boundary condition functor */
    REMORAPhysBCFunct (const int lev,
                      const amrex::Geometry& geom, const amrex::Vector<amrex::BCRec>& domain_bcs_type,
                      const amrex::Gpu::DeviceVector<amrex::BCRec>& domain_bcs_type_d,
                      const amrex::Vector<amrex::GpuArray<amrex::Real,AMREX_SPACEDIM*2>>& bc_extdir_vals,
                      int ncons
                    )
        : m_lev(lev), m_geom(geom), m_domain_bcs_type(domain_bcs_type),
          m_domain_bcs_type_d(domain_bcs_type_d),
          m_bc_extdir_vals(bc_extdir_vals),
          m_ncons(ncons)
    {
        m_bc_extdir_vals_d.resize(m_bc_extdir_vals.size());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                         m_bc_extdir_vals.begin(), m_bc_extdir_vals.end(),
                         m_bc_extdir_vals_d.begin());
    }

    ~REMORAPhysBCFunct () {}

    /** \brief apply boundary condition to mf */
    void operator() (amrex::MultiFab& mf, const amrex::MultiFab& mask, int icomp, int ncomp, amrex::IntVect const& nghost,
                     amrex::Real time, int bccomp, int n_not_fill=0,
                     const amrex::MultiFab& mf_calc = amrex::MultiFab(),
                     const amrex::MultiFab& mf_msku = amrex::MultiFab(),
                     const amrex::MultiFab& mf_mskv = amrex::MultiFab());

    /** \brief apply x-velocity type boundary conditions */
    void impose_xvel_bcs (const amrex::Array4<amrex::Real>& dest_arr, const amrex::Box& bx, const amrex::Box& domain,
                          const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> dxInv, const amrex::Array4<const amrex::Real>& msku,
                          const amrex::Array4<const amrex::Real>& calc_arr, amrex::Real time, int bccomp);

    /** \brief apply y-velocity type boundary conditions */
    void impose_yvel_bcs (const amrex::Array4<amrex::Real>& dest_arr, const amrex::Box& bx, const amrex::Box& domain,
                          const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> dxInv, const amrex::Array4<const amrex::Real>& mskv,
                          const amrex::Array4<const amrex::Real>& calc_arr, amrex::Real time, int bccomp);


    /** \brief apply z-velocity type boundary conditions */
    void impose_zvel_bcs (const amrex::Array4<amrex::Real>& dest_arr,
                          const amrex::Box& bx, const amrex::Box& domain,
                          const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> dxInv, const amrex::Array4<const amrex::Real>& mskr,
                          amrex::Real time, int bccomp);

    /** \brief apply scalar type boundary conditions */
    void impose_cons_bcs (const amrex::Array4<amrex::Real>& mf, const amrex::Box& bx, const amrex::Box& valid_bx,
                          const amrex::Box& domain,
                          const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> dxInv, const amrex::Array4<const amrex::Real>& mskr,
                          const amrex::Array4<const amrex::Real>& msku, const amrex::Array4<const amrex::Real>& mskv,
                          const amrex::Array4<const amrex::Real>& calc_arr,int icomp, int ncomp, amrex::Real time, int bccomp, int n_not_fill);


private:
    int                  m_lev;
    amrex::Geometry      m_geom;
    amrex::Vector<amrex::BCRec>            m_domain_bcs_type;
    amrex::Gpu::DeviceVector<amrex::BCRec> m_domain_bcs_type_d;
    amrex::Vector<amrex::GpuArray<amrex::Real, AMREX_SPACEDIM*2>> m_bc_extdir_vals;
    amrex::Gpu::DeviceVector<amrex::GpuArray<amrex::Real, AMREX_SPACEDIM*2>> m_bc_extdir_vals_d;
    int m_ncons;
};

#endif
