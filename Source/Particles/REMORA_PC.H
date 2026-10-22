#ifndef REMORA_PC_H_
#define REMORA_PC_H_

#ifdef REMORA_USE_PARTICLES

#include <string>
#include <REMORA_Constants.H>
#include <AMReX_Particles.H>

struct REMORAParticlesIntIdxAoS
{
    enum {
        k = 0,
        ncomps
    };
};

struct REMORAParticlesRealIdxAoS
{
    enum {
        ncomps = 0
    };
};

struct REMORAParticlesIntIdxSoA
{
    enum {
        ncomps = 0
    };
};

struct REMORAParticlesRealIdxSoA
{
    enum {
        vx = 0,
        vy,
        vz,
        mass,
        ncomps
    };
};

namespace REMORAParticleInitializations
{
    /* list of particle initializations */
    const std::string init_box_uniform  = "box";
}

namespace REMORAParticleNames
{
    const std::string tracers = "tracer_particles";
    const std::string hydro = "hydro_particles";
}

struct REMORAParticlesAssignor
{
    template <typename P>
    AMREX_GPU_HOST_DEVICE
    amrex::IntVect operator() (  P const& p,
                                 amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                                 amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi,
                                 const amrex::Box& domain ) const noexcept
    {
        amrex::IntVect iv(
            AMREX_D_DECL(   int(amrex::Math::floor((p.pos(0)-plo[0])*dxi[0])),
                            int(amrex::Math::floor((p.pos(1)-plo[1])*dxi[1])),
                            p.idata(REMORAParticlesIntIdxAoS::k) ) );
        iv[0] += domain.smallEnd()[0];
        iv[1] += domain.smallEnd()[1];
        return iv;
    }
};

// Particles are tracked vertically by a fractional w-level index s in [0,nz]: s = k + f puts
// a particle the fraction f of the way from w-level k to w-level k+1 of its column, so s = 0
// is the bed and s = nz the free surface. Moving at fixed s follows the terrain-following
// levels, which is what keeps a particle in the water column over sloping bathymetry.

/**
 * Height of w-level a_k at horizontal position (a_x, a_y), bilinear in the nodal heights
 * a_z_nd (z_phys_nd), whose (i,j,k) is the low corner of cell (i,j,k)
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real level_height (  amrex::ParticleReal a_x, amrex::ParticleReal a_y, int a_k,
                            amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& a_plo,
                            amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& a_dxi,
                            amrex::Array4<amrex::Real const> const& a_z_nd )
{
    const auto lx = static_cast<amrex::Real>((a_x-a_plo[0])*a_dxi[0]);
    const auto ly = static_cast<amrex::Real>((a_y-a_plo[1])*a_dxi[1]);
    const int i = static_cast<int>(amrex::Math::floor(lx));
    const int j = static_cast<int>(amrex::Math::floor(ly));
    const amrex::Real sx = lx - static_cast<amrex::Real>(i);
    const amrex::Real sy = ly - static_cast<amrex::Real>(j);
    return a_z_nd(i  ,j  ,a_k) * (one-sx) * (one-sy) +
           a_z_nd(i+1,j  ,a_k) *      sx  * (one-sy) +
           a_z_nd(i  ,j+1,a_k) * (one-sx) *      sy  +
           a_z_nd(i+1,j+1,a_k) *      sx  *      sy;
}

/**
 * Fractional w-level index of height a_z in the column at (a_x, a_y). A height below the
 * bed or above the free surface is clamped onto it. a_k, the cell the particle was last
 * in, is where the search starts.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real height_to_level_index ( amrex::ParticleReal a_x, amrex::ParticleReal a_y,
                                    amrex::ParticleReal a_z, int a_k, int a_nz,
                                    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& a_plo,
                                    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& a_dxi,
                                    amrex::Array4<amrex::Real const> const& a_z_nd )
{
    const auto z = static_cast<amrex::Real>(a_z);
    int k = amrex::max(0, amrex::min(a_k, a_nz-1));
    amrex::Real zlo = level_height(a_x, a_y, k  , a_plo, a_dxi, a_z_nd);
    amrex::Real zhi = level_height(a_x, a_y, k+1, a_plo, a_dxi, a_z_nd);
    while (z > zhi && k < a_nz-1) {
        ++k;
        zlo = zhi;
        zhi = level_height(a_x, a_y, k+1, a_plo, a_dxi, a_z_nd);
    }
    while (z < zlo && k > 0) {
        --k;
        zhi = zlo;
        zlo = level_height(a_x, a_y, k, a_plo, a_dxi, a_z_nd);
    }
    const amrex::Real f = amrex::max(zero, amrex::min(one, (z - zlo) / (zhi - zlo)));
    return static_cast<amrex::Real>(k) + f;
}

/** Cell holding fractional w-level index a_s; the free surface s = nz belongs to the top cell */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int level_index_cell (amrex::Real a_s, int a_nz)
{
    return amrex::max(0, amrex::min(static_cast<int>(amrex::Math::floor(a_s)), a_nz-1));
}

/** Height of fractional w-level index a_s at (a_x, a_y); the inverse of height_to_level_index */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real level_index_to_height ( amrex::ParticleReal a_x, amrex::ParticleReal a_y,
                                    amrex::Real a_s, int a_nz,
                                    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& a_plo,
                                    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& a_dxi,
                                    amrex::Array4<amrex::Real const> const& a_z_nd )
{
    const int k = level_index_cell(a_s, a_nz);
    const amrex::Real f = a_s - static_cast<amrex::Real>(k);
    const amrex::Real zlo = level_height(a_x, a_y, k  , a_plo, a_dxi, a_z_nd);
    const amrex::Real zhi = level_height(a_x, a_y, k+1, a_plo, a_dxi, a_z_nd);
    return zlo + f * (zhi - zlo);
}

/**
 * Trilinear interpolation of a_arr to horizontal position (a_x, a_y) and fractional w-level
 * index a_s. a_nodal gives a_arr's staggering. Data that is cell-centered in z is held
 * constant over the half cells at the bed and the surface, so the vertical stencil stays
 * within the valid cells [0,nz) (or w-levels [0,nz] for z-nodal data).
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real interpolate_to_level_index ( amrex::ParticleReal a_x, amrex::ParticleReal a_y,
                                         amrex::Real a_s, int a_nz, amrex::IntVect const& a_nodal,
                                         amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& a_plo,
                                         amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& a_dxi,
                                         amrex::Array4<amrex::Real const> const& a_arr )
{
    const amrex::Real lx = static_cast<amrex::Real>((a_x-a_plo[0])*a_dxi[0])
                         - static_cast<amrex::Real>(!a_nodal[0]) * myhalf;
    const amrex::Real ly = static_cast<amrex::Real>((a_y-a_plo[1])*a_dxi[1])
                         - static_cast<amrex::Real>(!a_nodal[1]) * myhalf;
    const int i0 = static_cast<int>(amrex::Math::floor(lx));
    const int j0 = static_cast<int>(amrex::Math::floor(ly));
    const amrex::Real sx[] = {one - (lx - static_cast<amrex::Real>(i0)), lx - static_cast<amrex::Real>(i0)};
    const amrex::Real sy[] = {one - (ly - static_cast<amrex::Real>(j0)), ly - static_cast<amrex::Real>(j0)};

    const int kmax = a_nodal[2] ? a_nz : a_nz-1;
    const amrex::Real lz = amrex::max(zero, amrex::min(static_cast<amrex::Real>(kmax),
                                      a_s - static_cast<amrex::Real>(!a_nodal[2]) * myhalf));
    const int k0 = amrex::min(static_cast<int>(amrex::Math::floor(lz)), amrex::max(kmax-1, 0));
    const int k1 = amrex::min(k0+1, kmax);
    const amrex::Real sz[] = {one - (lz - static_cast<amrex::Real>(k0)), lz - static_cast<amrex::Real>(k0)};

    amrex::Real val = zero;
    for (int jj = 0; jj < 2; ++jj) {
        for (int ii = 0; ii < 2; ++ii) {
            val += sx[ii] * sy[jj] * ( sz[0] * a_arr(i0+ii, j0+jj, k0)
                                     + sz[1] * a_arr(i0+ii, j0+jj, k1) );
        }
    }
    return val;
}

/**
 * Bilinear interpolation of the 2D cell-centered field a_arr (stored at k = 0, like pm and pn)
 * to horizontal position (a_x, a_y)
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real interpolate_horizontal ( amrex::ParticleReal a_x, amrex::ParticleReal a_y,
                                     amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& a_plo,
                                     amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& a_dxi,
                                     amrex::Array4<amrex::Real const> const& a_arr )
{
    const amrex::Real lx = static_cast<amrex::Real>((a_x-a_plo[0])*a_dxi[0]) - myhalf;
    const amrex::Real ly = static_cast<amrex::Real>((a_y-a_plo[1])*a_dxi[1]) - myhalf;
    const int i = static_cast<int>(amrex::Math::floor(lx));
    const int j = static_cast<int>(amrex::Math::floor(ly));
    const amrex::Real sx = lx - static_cast<amrex::Real>(i);
    const amrex::Real sy = ly - static_cast<amrex::Real>(j);
    return a_arr(i  ,j  ,0) * (one-sx) * (one-sy) +
           a_arr(i+1,j  ,0) *      sx  * (one-sy) +
           a_arr(i  ,j+1,0) * (one-sx) *      sy  +
           a_arr(i+1,j+1,0) *      sx  *      sy;
}

class REMORAPC : public amrex::ParticleContainer<  REMORAParticlesRealIdxAoS::ncomps,  // AoS real attributes
                                                REMORAParticlesIntIdxAoS::ncomps,   // AoS integer attributes
                                                REMORAParticlesRealIdxSoA::ncomps,  // SoA real attributes
                                                REMORAParticlesIntIdxSoA::ncomps,   // SoA integer attributes
                                                amrex::DefaultAllocator,
                                                REMORAParticlesAssignor >
{
    public:

        /*! Constructor */
        REMORAPC ( amrex::ParGDBBase* a_gdb,
                const std::string& a_name = "particles" )
            : amrex::ParticleContainer< REMORAParticlesRealIdxAoS::ncomps,   // AoS real attributes
                                        REMORAParticlesIntIdxAoS::ncomps,    // AoS integer attributes
                                        REMORAParticlesRealIdxSoA::ncomps,   // SoA real attributes
                                        REMORAParticlesIntIdxSoA::ncomps,    // SoA integer attributes
                                        amrex::DefaultAllocator,
                                        REMORAParticlesAssignor> (a_gdb)
        {
            BL_PROFILE("REMORAPCPC::REMORAPC()");
            m_name = a_name;
            readInputs();
        }

        /*! Constructor */
        REMORAPC ( const amrex::Geometry&            a_geom,
                const amrex::DistributionMapping& a_dmap,
                const amrex::BoxArray&            a_ba,
                const std::string&                a_name = "particles" )
            : amrex::ParticleContainer< REMORAParticlesRealIdxAoS::ncomps,   // AoS real attributes
                                        REMORAParticlesIntIdxAoS::ncomps,    // AoS real attributes
                                        REMORAParticlesRealIdxSoA::ncomps,   // SoA real attributes
                                        REMORAParticlesIntIdxSoA::ncomps,    // SoA integer attributes
                                        amrex::DefaultAllocator,
                                        REMORAParticlesAssignor> ( a_geom, a_dmap, a_ba )
        {
            BL_PROFILE("REMORAPCPC::REMORAPC()");
            m_name = a_name;
            readInputs();
        }

        /*! Initialize particles in domain */
        virtual void InitializeParticles (const std::unique_ptr<amrex::MultiFab>& a_ptr = nullptr);

        /*! Evolve particles for one time step */
        virtual void EvolveParticles (   int a_lev,
                                         amrex::Real a_dt_lev,
                                         amrex::Vector<amrex::MultiFab const*>& a_flow_vel,
                                         const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& a_z_phys_nd,
                                         amrex::MultiFab const& a_pm,
                                         amrex::MultiFab const& a_pn );

        /*! Get real-type particle attribute names */
        virtual amrex::Vector<std::string> varNames () const
        {
            BL_PROFILE("REMORAPCPC::varNames()");
            return {AMREX_D_DECL("xvel","yvel","zvel"),"mass"};
        }

        /*! Get real-type particle attribute names */
        virtual amrex::Vector<std::string> meshPlotVarNames () const
        {
            BL_PROFILE("REMORAPCPC::varNames()");
            return {"mass_density"};
        }

        /*! Uses midpoint method to advance particles using flow velocity. */
        virtual void AdvectWithFlow (  int a_lev,
                                       amrex::Real a_dt,
                                       amrex::Vector<amrex::MultiFab const*>& a_umac,
                                       const std::unique_ptr<amrex::MultiFab>& a_z_height,
                                       amrex::MultiFab const& a_pm,
                                       amrex::MultiFab const& a_pn );

        /*! Compute mass density */
        virtual void massDensity ( amrex::MultiFab&, const int&, const int& a_comp = 0) const;

        /*! Compute mesh variable from particles */
        virtual void computeMeshVar(    const std::string&  a_var_name,
                                        amrex::MultiFab&    a_mf,
                                        const int           a_lev) const
        {
            if (a_var_name == "mass_density") {
                massDensity( a_mf, a_lev );
            } else {
                a_mf.setVal(zero);
            }
        }

        /*! Specify if particles should advect with flow */
        inline void setAdvectWithFlow (bool a_flag)
        {
            BL_PROFILE("REMORAPCPC::setAdvectWithFlow()");
            m_advect_w_flow = a_flag;
        }

        // the following functions should ideally be private or protected, but need to be
        // public due to CUDA extended lambda capture rules

        /*! Default particle initialization */
        void initializeParticlesUniformDistributionInBox (const std::unique_ptr<amrex::MultiFab>& a_ptr,
                                                          const amrex::RealBox& particle_box);

    protected:

        bool m_advect_w_flow;               /*!< advect with flow velocity */

        amrex::RealBox m_particle_box;      /*!< box within which to place particles */

        std::string m_name;                 /*!< name of this particle species */

        std::string m_initialization_type;  /*!< initial particle distribution type */
        int m_ppc_init;                     /*!< initial number of particles per cell */

        /*! read inputs from file */
        virtual void readInputs ();

    private:

        bool place_randomly_in_cells; /*!< place particles at random positions? */
};

#endif
#endif
