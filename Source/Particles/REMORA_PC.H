#ifndef REMORA_PC_H_
#define REMORA_PC_H_

#ifdef REMORA_USE_PARTICLES

#include <string>
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

template <typename P>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void update_location_idata (  P& a_p,
                              amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& a_plo,
                              amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& a_dxi,
                              const amrex::Array4<amrex::Real const>&  a_height_arr )
{
    amrex::IntVect iv( int(amrex::Math::floor((a_p.pos(0)-a_plo[0])*a_dxi[0])),
                       int(amrex::Math::floor((a_p.pos(1)-a_plo[1])*a_dxi[1])),
                       a_p.idata(REMORAParticlesIntIdxAoS::k) );

    if (a_height_arr) {
        amrex::Real lx = (a_p.pos(0)-a_plo[0])*a_dxi[0] - static_cast<amrex::Real>(iv[0]);
        amrex::Real ly = (a_p.pos(1)-a_plo[1])*a_dxi[1] - static_cast<amrex::Real>(iv[1]);
        auto zlo = a_height_arr(iv[0]  ,iv[1]  ,iv[2]  ) * (1.0-lx) * (1.0-ly) +
                   a_height_arr(iv[0]+1,iv[1]  ,iv[2]  ) *      lx  * (1.0-ly) +
                   a_height_arr(iv[0]  ,iv[1]+1,iv[2]  ) * (1.0-lx) * ly +
                   a_height_arr(iv[0]+1,iv[1]+1,iv[2]  ) *      lx  * ly;
        auto zhi = a_height_arr(iv[0]  ,iv[1]  ,iv[2]+1) * (1.0-lx) * (1.0-ly) +
                   a_height_arr(iv[0]+1,iv[1]  ,iv[2]+1) *      lx  * (1.0-ly) +
                   a_height_arr(iv[0]  ,iv[1]+1,iv[2]+1) * (1.0-lx) * ly +
                   a_height_arr(iv[0]+1,iv[1]+1,iv[2]+1) *      lx  * ly;

        if (a_p.pos(2) > zhi) {
            a_p.idata(REMORAParticlesIntIdxAoS::k) += 1;
        } else if (a_p.pos(2) <= zlo) {
            a_p.idata(REMORAParticlesIntIdxAoS::k) -= 1;
        }
    }
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
        virtual void EvolveParticles (   int lev,
                                         amrex::Real dt,
                                         amrex::Vector<amrex::MultiFab const*>& a_flow_vel,
                                         const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& );

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
        virtual void AdvectWithFlow (  int lev,
                                       amrex::Real dt,
                                       amrex::Vector<amrex::MultiFab const*>&,
                                       const std::unique_ptr<amrex::MultiFab>& );

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
                a_mf.setVal(0.0);
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
