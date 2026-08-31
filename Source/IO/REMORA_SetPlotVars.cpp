#include <REMORA.H>
#include "AMReX_Interp_3D_C.H"
#include "AMReX_PlotFileUtil.H"

using namespace amrex;

template<typename V, typename T>
bool containerHasElement(const V& iterable, const T& query) {
    return std::find(iterable.begin(), iterable.end(), query) != iterable.end();
}

/**
 * @param pp_plot_var_names_3d   list of variable names to plot read in from parameter file
 */
void
REMORA::set3DPlotVariables (const std::string& pp_plot_var_names_3d)
{
    ParmParse pp(pp_prefix);

    if (pp.contains("plot_vars")) {
        amrex::Abort("You must use plot_vars_3d rather than plot_vars");
    }

    if (pp.contains(pp_plot_var_names_3d.c_str()))
    {
        std::string nm;

        int nPltVars = pp.countval(pp_plot_var_names_3d.c_str());

        for (int i = 0; i < nPltVars; i++)
        {
            pp.get(pp_plot_var_names_3d.c_str(), nm, i);

            // Add the named variable to our list of plot variables
            // if it is not already in the list
            if (!containerHasElement(plot_var_names_3d, nm)) {
                plot_var_names_3d.push_back(nm);
            }
        }
    } else {
        //
        // The default is to add none of the variables to the list
        //
        plot_var_names_3d.clear();
    }

    const bool expand_fennel = containerHasElement(plot_var_names_3d, "fennel") &&
                               REMORABiology::has_biology(biology_model);
    if (expand_fennel) {
        // "fennel" means the biology tracers only, not any passive scalars ahead of them
        for (int icomp = Bio_comp; icomp < ncons; ++icomp) {
            if (!containerHasElement(plot_var_names_3d, cons_names[icomp])) {
                plot_var_names_3d.push_back(cons_names[icomp]);
            }
        }
    }

    // Get state variables in the same order as we define them,
    // since they may be in any order in the input list
    Vector<std::string> tmp_plot_names;

    for (int i = 0; i < ncons; ++i) {
        if ( containerHasElement(plot_var_names_3d, cons_names[i]) ) {
            tmp_plot_names.push_back(cons_names[i]);
        }
    }
    // Check for velocity since it's not in cons_names
    // If we are asked for any velocity component, we will need them all
    if (containerHasElement(plot_var_names_3d, "x_velocity") ||
        containerHasElement(plot_var_names_3d, "y_velocity") ||
        containerHasElement(plot_var_names_3d, "z_velocity")) {
        tmp_plot_names.push_back("x_velocity");
        tmp_plot_names.push_back("y_velocity");
        tmp_plot_names.push_back("z_velocity");
    }

    // If we are asked for any location component, we will provide them all
    if (containerHasElement(plot_var_names_3d, "x_cc") ||
        containerHasElement(plot_var_names_3d, "y_cc") ||
        containerHasElement(plot_var_names_3d, "z_cc")) {
        tmp_plot_names.push_back("x_cc");
        tmp_plot_names.push_back("y_cc");
        tmp_plot_names.push_back("z_cc");
    }

    for (int i = 0; i < derived_names_3d.size(); ++i) {
        if ( containerHasElement(plot_var_names_3d, derived_names_3d[i]) ) {
               tmp_plot_names.push_back(derived_names_3d[i]);
        } // if
    } // i

#ifdef REMORA_USE_PARTICLES
    const auto& particles_namelist( particleData.getNamesUnalloc() );
    for (auto it = particles_namelist.cbegin(); it != particles_namelist.cend(); ++it) {
        std::string tmp( (*it)+"_count" );
        if (containerHasElement(plot_var_names_3d, tmp) ) {
            tmp_plot_names.push_back(tmp);
        }
    }
#endif

    // Check to see if we found all the requested variables
    for (auto plot_name : plot_var_names_3d) {
      if (plot_name == "fennel" && expand_fennel) {
          continue;
      }
      if (!containerHasElement(tmp_plot_names, plot_name)) {
           Warning("\nWARNING: Requested to plot variable '" + plot_name + "' in 3D list but it is not available");
      }
    }
    plot_var_names_3d = tmp_plot_names;
}

/**
 * @param pp_plot_var_names_2d   list of variable names to plot read in from parameter file
 */
void
REMORA::set2DPlotVariables (const std::string& pp_plot_var_names_2d)
{
    ParmParse pp(pp_prefix);

    if (pp.contains(pp_plot_var_names_2d.c_str()))
    {
        std::string nm;

        int nPltVars = pp.countval(pp_plot_var_names_2d.c_str());

        for (int i = 0; i < nPltVars; i++)
        {
            pp.get(pp_plot_var_names_2d.c_str(), nm, i);

            // Add the named variable to our list of plot variables
            // if it is not already in the list
            if (!containerHasElement(plot_var_names_2d, nm)) {
                plot_var_names_2d.push_back(nm);
            }
        }
    } else {
        //
        // The default is to add none of the variables to the list
        //
        plot_var_names_2d.clear();
    }

    // If horizontal mixing is scaled_to_grid, automatically output the spatially
    // varying coefficients used by the run as 2D fields.
    if (solverChoice.horiz_mixing_type == HorizMixingType::scaled_to_grid) {
        if (!containerHasElement(plot_var_names_2d, "visc2")) {
            plot_var_names_2d.push_back("visc2");
        }
        for (int n = 0; n < ncons; ++n) {
            const std::string nm = std::string("diff2_") + cons_names[n];
            if (!containerHasElement(plot_var_names_2d, nm)) {
                plot_var_names_2d.push_back(nm);
            }
        }
    }

    Vector<std::string> tmp_plot_names;

    if (containerHasElement(plot_var_names_2d, "zeta")) tmp_plot_names.push_back("zeta");
    if (containerHasElement(plot_var_names_2d, "h")) tmp_plot_names.push_back("h");
    if (containerHasElement(plot_var_names_2d, "f")) tmp_plot_names.push_back("f");
    if (containerHasElement(plot_var_names_2d, "ubar")) tmp_plot_names.push_back("ubar");
    if (containerHasElement(plot_var_names_2d, "vbar")) tmp_plot_names.push_back("vbar");
    if (containerHasElement(plot_var_names_2d, "sustr")) tmp_plot_names.push_back("sustr");
    if (containerHasElement(plot_var_names_2d, "bustr")) tmp_plot_names.push_back("bustr");
    if (containerHasElement(plot_var_names_2d, "svstr")) tmp_plot_names.push_back("svstr");
    if (containerHasElement(plot_var_names_2d, "bvstr")) tmp_plot_names.push_back("bvstr");

    for (int i = 0; i < derived_names_2d.size(); ++i) {
        if ( containerHasElement(plot_var_names_2d, derived_names_2d[i]) ) {
               tmp_plot_names.push_back(derived_names_2d[i]);
        } // if
    } // i

    // Horizontal mixing coefficients (2D rho points)
    if (containerHasElement(plot_var_names_2d, "visc2")) {
        tmp_plot_names.push_back("visc2");
    }
    for (int n = 0; n < ncons; ++n) {
        const std::string nm = std::string("diff2_") + cons_names[n];
        if (containerHasElement(plot_var_names_2d, nm)) {
            tmp_plot_names.push_back(nm);
        }
    }
    for (int n = 0; n < ncons; ++n) {
        const std::string nm = std::string("stflux_") + cons_names[n];
        if (containerHasElement(plot_var_names_2d, nm)) {
            tmp_plot_names.push_back(nm);
        }
    }

    if (containerHasElement(plot_var_names_2d, "lrflux")) tmp_plot_names.push_back("lrflux");
    if (containerHasElement(plot_var_names_2d, "lhflux")) tmp_plot_names.push_back("lhflux");
    if (containerHasElement(plot_var_names_2d, "srflux")) tmp_plot_names.push_back("srflux");
    if (containerHasElement(plot_var_names_2d, "shflux")) tmp_plot_names.push_back("shflux");

    if (containerHasElement(plot_var_names_2d, "mask_rho")) tmp_plot_names.push_back("mask_rho");
    if (containerHasElement(plot_var_names_2d, "mask_u")) tmp_plot_names.push_back("mask_u");
    if (containerHasElement(plot_var_names_2d, "mask_v")) tmp_plot_names.push_back("mask_v");

    // Check to see if we found all the requested variables
    for (auto plot_name : plot_var_names_2d) {
      if (!containerHasElement(tmp_plot_names, plot_name)) {
           Warning("\nWARNING: Requested to plot variable '" + plot_name + "' in 2D list but it is not available");
      }
    }
    plot_var_names_2d = tmp_plot_names;
}

/**
 * @param pp_plot_var_names     variables to add to plot list
 */
void
REMORA::append3DPlotVariables (const std::string& pp_plot_var_names_3d)
{
    ParmParse pp(pp_prefix);

    if (pp.contains("plot_vars")) {
        amrex::Abort("You must use plot_vars_3d rather than plot_vars");
    }

    // This runs after the particle containers are set up, so the only thing it may add is
    // a particle mesh variable, which set3DPlotVariables could not know about yet. Adding
    // anything else back would resurrect a name set3DPlotVariables deliberately dropped as
    // unavailable, and nothing downstream would ever fill it: the plotfile would carry the
    // uninitialized sentinel from WritePlotFile under a real-looking variable name. The
    // "fennel" keyword needs no handling here, since set3DPlotVariables already expanded it.
#ifdef REMORA_USE_PARTICLES
    Vector<std::string> particle_mesh_plot_names;
    particleData.GetMeshPlotVarNames( particle_mesh_plot_names );

    if (pp.contains(pp_plot_var_names_3d.c_str())) {
        std::string nm;
        int nPltVars = pp.countval(pp_plot_var_names_3d.c_str());
        for (int i = 0; i < nPltVars; i++) {
            pp.get(pp_plot_var_names_3d.c_str(), nm, i);
            if (containerHasElement(particle_mesh_plot_names, nm) &&
                !containerHasElement(plot_var_names_3d, nm)) {
                plot_var_names_3d.push_back(nm);
            }
        }
    }
#else
    amrex::ignore_unused(pp_plot_var_names_3d);
#endif
}

/**
 * @param pp_plot_var_names_2d    variables to add to plot list
 */
void
REMORA::append2DPlotVariables (const std::string& pp_plot_var_names_2d)
{
    // Same rule as append3DPlotVariables: only a particle mesh variable may be added back.
#ifdef REMORA_USE_PARTICLES
    ParmParse pp(pp_prefix);

    Vector<std::string> particle_mesh_plot_names;
    particleData.GetMeshPlotVarNames( particle_mesh_plot_names );

    if (pp.contains(pp_plot_var_names_2d.c_str())) {
        std::string nm;
        int nPltVars = pp.countval(pp_plot_var_names_2d.c_str());
        for (int i = 0; i < nPltVars; i++) {
            pp.get(pp_plot_var_names_2d.c_str(), nm, i);
            if (containerHasElement(particle_mesh_plot_names, nm) &&
                !containerHasElement(plot_var_names_2d, nm)) {
                plot_var_names_2d.push_back(nm);
            }
        }
    }
#else
    amrex::ignore_unused(pp_plot_var_names_2d);
#endif
}
