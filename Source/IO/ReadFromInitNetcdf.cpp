#include "NCFile.H"
#include "AMReX_FArrayBox.H"
#include "DataStruct.H"

using namespace amrex;

#ifdef REMORA_USE_NETCDF
void
read_data_from_netcdf (int /*lev*/,
                       const Box& domain,
                       const std::string& fname,
                       FArrayBox& NC_temp_fab, FArrayBox& NC_salt_fab,
                       FArrayBox& NC_xvel_fab, FArrayBox& NC_yvel_fab,
                       FArrayBox& NC_ubar_fab, FArrayBox& NC_vbar_fab,
                       FArrayBox& NC_zeta_fab,
                       FArrayBox& NC_mskr_fab, FArrayBox& NC_msku_fab,
                       FArrayBox& NC_mskv_fab)
{
    amrex::Print() << "Loading initial solution data from NetCDF file " << fname << std::endl;

    Vector<FArrayBox*> NC_fabs;
    Vector<std::string> NC_names;
    Vector<enum NC_Data_Dims_Type> NC_dim_types;

    NC_fabs.push_back(&NC_temp_fab); NC_names.push_back("temp");     NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE); // 0
    NC_fabs.push_back(&NC_salt_fab); NC_names.push_back("salt");     NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE); // 1
    NC_fabs.push_back(&NC_xvel_fab); NC_names.push_back("u");        NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE); // 2
    NC_fabs.push_back(&NC_yvel_fab); NC_names.push_back("v");        NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE); // 3
    NC_fabs.push_back(&NC_ubar_fab), NC_names.push_back("ubar");     NC_dim_types.push_back(NC_Data_Dims_Type::Time_SN_WE); // 4
    NC_fabs.push_back(&NC_vbar_fab); NC_names.push_back("vbar");     NC_dim_types.push_back(NC_Data_Dims_Type::Time_SN_WE); // 5
    NC_fabs.push_back(&NC_zeta_fab); NC_names.push_back("zeta");     NC_dim_types.push_back(NC_Data_Dims_Type::Time_SN_WE); // 6
    NC_fabs.push_back(&NC_mskr_fab); NC_names.push_back("mask_rho"); NC_dim_types.push_back(NC_Data_Dims_Type::Time_SN_WE); // 7
    NC_fabs.push_back(&NC_msku_fab); NC_names.push_back("mask_u");   NC_dim_types.push_back(NC_Data_Dims_Type::Time_SN_WE); // 8
    NC_fabs.push_back(&NC_mskv_fab); NC_names.push_back("mask_v");   NC_dim_types.push_back(NC_Data_Dims_Type::Time_SN_WE); // 9

    // Read the netcdf file and fill these FABs
    BuildFABsFromNetCDFFile<FArrayBox,Real>(domain, fname, NC_names, NC_dim_types, NC_fabs);
}

void
read_bathymetry_from_netcdf (int /*lev*/,
                             const Box& domain,
                             const std::string& fname,
                             FArrayBox& NC_h_fab,
                             FArrayBox& NC_pm_fab, FArrayBox& NC_pn_fab)
{
    amrex::Print() << "Loading initial bathymetry from NetCDF file " << fname << std::endl;

    Vector<FArrayBox*> NC_fabs;
    Vector<std::string> NC_names;
    Vector<enum NC_Data_Dims_Type> NC_dim_types;

    NC_fabs.push_back(&NC_h_fab )   ; NC_names.push_back("h")    ; NC_dim_types.push_back(NC_Data_Dims_Type::SN_WE); // 0
    NC_fabs.push_back(&NC_pm_fab)   ; NC_names.push_back("pm")   ; NC_dim_types.push_back(NC_Data_Dims_Type::SN_WE); // 1
    NC_fabs.push_back(&NC_pn_fab)   ; NC_names.push_back("pn")   ; NC_dim_types.push_back(NC_Data_Dims_Type::SN_WE); // 2

    // Read the netcdf file and fill these FABs
    BuildFABsFromNetCDFFile<FArrayBox,Real>(domain, fname, NC_names, NC_dim_types, NC_fabs);
}

void
read_coriolis_from_netcdf (const Box& domain,
                           const std::string& fname,
                           FArrayBox& NC_fcor_fab)
{
    amrex::Print() << "Loading initial coriolis from NetCDF file " << fname << std::endl;

    Vector<FArrayBox*> NC_fabs;
    Vector<std::string> NC_names;
    Vector<enum NC_Data_Dims_Type> NC_dim_types;

    NC_fabs.push_back(&NC_fcor_fab )   ; NC_names.push_back("f")    ; NC_dim_types.push_back(NC_Data_Dims_Type::SN_WE); // 0

    // Read the netcdf file and fill these FABs
    BuildFABsFromNetCDFFile<FArrayBox,Real>(domain, fname, NC_names, NC_dim_types, NC_fabs);
}
#endif // ROMSX_USE_NETCDF
