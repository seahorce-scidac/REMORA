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

    // If horizontal mixing is scaled to grid size, automatically output the
    // spatially-varying coefficients used by the run.
    if (solverChoice.horiz_mixing_type == HorizMixingType::scaled_to_grid) {
        if (!containerHasElement(plot_var_names_3d, "visc2")) {
            plot_var_names_3d.push_back("visc2");
        }
        for (int n = 0; n < NCONS; ++n) {
            const std::string nm = std::string("diff2_") + cons_names[n];
            if (!containerHasElement(plot_var_names_3d, nm)) {
                plot_var_names_3d.push_back(nm);
            }
        }
    }

    // Get state variables in the same order as we define them,
    // since they may be in any order in the input list
    Vector<std::string> tmp_plot_names;

    for (int i = 0; i < NCONS; ++i) {
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

    for (int i = 0; i < derived_names.size(); ++i) {
        if ( containerHasElement(plot_var_names_3d, derived_names[i]) ) {
               tmp_plot_names.push_back(derived_names[i]);
        } // if
    } // i

    // Horizontal mixing coefficients (cell-centered rho points)
    if (containerHasElement(plot_var_names_3d, "visc2")) {
        tmp_plot_names.push_back("visc2");
    }
    for (int n = 0; n < NCONS; ++n) {
        const std::string nm = std::string("diff2_") + cons_names[n];
        if (containerHasElement(plot_var_names_3d, nm)) {
            tmp_plot_names.push_back(nm);
        }
    }

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
      if (!containerHasElement(tmp_plot_names, plot_name)) {
           Warning("\nWARNING: Requested to plot variable '" + plot_name + "' but it is not available");
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

    // Get state variables in the same order as we define them,
    // since they may be in any order in the input list
    Vector<std::string> tmp_plot_names;

    for (int i = 0; i < NCONS; ++i) {
        if ( containerHasElement(plot_var_names_2d, cons_names[i]) ) {
            tmp_plot_names.push_back(cons_names[i]);
        }
    }
    // Check for velocity since it's not in cons_names
    // If we are asked for any velocity component, we will need them all
    if (containerHasElement(plot_var_names_2d, "x_velocity") ||
        containerHasElement(plot_var_names_2d, "y_velocity") ||
        containerHasElement(plot_var_names_2d, "z_velocity")) {
        tmp_plot_names.push_back("x_velocity");
        tmp_plot_names.push_back("y_velocity");
        tmp_plot_names.push_back("z_velocity");
    }

    // If we are asked for any location component, we will provide them all
    if (containerHasElement(plot_var_names_2d, "x_cc") ||
        containerHasElement(plot_var_names_2d, "y_cc") ||
        containerHasElement(plot_var_names_2d, "z_cc")) {
        tmp_plot_names.push_back("x_cc");
        tmp_plot_names.push_back("y_cc");
        tmp_plot_names.push_back("z_cc");
    }

    for (int i = 0; i < derived_names.size(); ++i) {
        if ( containerHasElement(plot_var_names_2d, derived_names[i]) ) {
               tmp_plot_names.push_back(derived_names[i]);
        } // if
    } // i

#ifdef REMORA_USE_PARTICLES
    const auto& particles_namelist( particleData.getNamesUnalloc() );
    for (auto it = particles_namelist.cbegin(); it != particles_namelist.cend(); ++it) {
        std::string tmp( (*it)+"_count" );
        if (containerHasElement(plot_var_names_2d, tmp) ) {
            tmp_plot_names.push_back(tmp);
        }
    }
#endif

    // Check to see if we found all the requested variables
    for (auto plot_name : plot_var_names_2d) {
      if (!containerHasElement(tmp_plot_names, plot_name)) {
           Warning("\nWARNING: Requested to plot variable '" + plot_name + "' but it is not available");
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

    if (pp.contains(pp_plot_var_names_3d.c_str())) {
        std::string nm;
        int nPltVars = pp.countval(pp_plot_var_names_3d.c_str());
        for (int i = 0; i < nPltVars; i++) {
            pp.get(pp_plot_var_names_3d.c_str(), nm, i);
            // Add the named variable to our list of plot variables
            // if it is not already in the list
            if (!containerHasElement(plot_var_names_3d, nm)) {
                plot_var_names_3d.push_back(nm);
            }
        }
    }

    // Same auto-append for scaled_to_grid as in set3DPlotVariables.
    if (solverChoice.horiz_mixing_type == HorizMixingType::scaled_to_grid) {
        if (!containerHasElement(plot_var_names_3d, "visc2")) {
            plot_var_names_3d.push_back("visc2");
        }
        for (int n = 0; n < NCONS; ++n) {
            const std::string nm = std::string("diff2_") + cons_names[n];
            if (!containerHasElement(plot_var_names_3d, nm)) {
                plot_var_names_3d.push_back(nm);
            }
        }
    }

    Vector<std::string> tmp_plot_names(0);
#ifdef REMORA_USE_PARTICLES
    Vector<std::string> particle_mesh_plot_names;
    particleData.GetMeshPlotVarNames( particle_mesh_plot_names );
    for (int i = 0; i < particle_mesh_plot_names.size(); i++) {
        std::string tmp(particle_mesh_plot_names[i]);
        if (containerHasElement(plot_var_names_3d, tmp) ) {
            tmp_plot_names.push_back(tmp);
        }
    }
#endif

    for (int i = 0; i < tmp_plot_names.size(); i++) {
        plot_var_names_3d.push_back( tmp_plot_names[i] );
    }
}

/**
 * @param pp_plot_var_names_2d    variables to add to plot list
 */
void
REMORA::append2DPlotVariables (const std::string& pp_plot_var_names_2d)
{
    ParmParse pp(pp_prefix);

    if (pp.contains(pp_plot_var_names_2d.c_str())) {
        std::string nm;
        int nPltVars = pp.countval(pp_plot_var_names_2d.c_str());
        for (int i = 0; i < nPltVars; i++) {
            pp.get(pp_plot_var_names_2d.c_str(), nm, i);
            // Add the named variable to our list of plot variables
            // if it is not already in the list
            if (!containerHasElement(plot_var_names_2d, nm)) {
                plot_var_names_2d.push_back(nm);
            }
        }
    }

    Vector<std::string> tmp_plot_names(0);
#ifdef REMORA_USE_PARTICLES
    Vector<std::string> particle_mesh_plot_names;
    particleData.GetMeshPlotVarNames( particle_mesh_plot_names );
    for (int i = 0; i < particle_mesh_plot_names.size(); i++) {
        std::string tmp(particle_mesh_plot_names[i]);
        if (containerHasElement(plot_var_names_2d, tmp) ) {
            tmp_plot_names.push_back(tmp);
        }
    }
#endif

    for (int i = 0; i < tmp_plot_names.size(); i++) {
        plot_var_names_2d.push_back( tmp_plot_names[i] );
    }
}
