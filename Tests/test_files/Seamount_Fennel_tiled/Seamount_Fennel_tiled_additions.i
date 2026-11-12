# Added to Tests/test_files/Seamount/Seamount.i by compose_test_input in
# CTestList.cmake, so the Seamount case itself is defined once. Settings here
# override the base's, since ParmParse reads the last occurrence.
#
# Seamount's terrain-following layers sit at different depths over the bump, so
# temperature and layer thickness at each level vary in both x and y, and so does
# the biology after one step. The test runs this input untiled and with small CPU
# tiles and requires the two runs to agree bitwise, so a kernel that mixes up
# columns across tile edges in either direction fails.
remora.max_step = 3
remora.plot_int = 3
# Fill new arrays with signaling NaN and trap invalid operations, so a biology
# kernel that reads scratch it has not written fails instead of using garbage.
amrex.init_snan = 1
amrex.fpe_trap_invalid = 1
remora.bulk_fluxes = true
remora.surface_radiation_flux = 150.0
remora.uwind = 5.0
remora.vwind = 2.0
remora.biology_ic_type = constant
remora.plot_vars_3d = salt temp x_velocity y_velocity z_velocity biology
remora.biology_model = fennel
# The options Upwelling_Fennel exercises, so their branches are checked under tiling.
remora.fennel.carbon = true
remora.fennel.oxygen = true
remora.fennel.river_don = true
remora.fennel.talk_nonconserv = true
remora.fennel.co2_schmidt = rw14
remora.fennel.oxygen_schmidt = rw14
# The secular atmospheric pCO2 needs a real date for the model clock.
remora.fennel.pco2air_type = secular
remora.time_ref = 20200301.0
remora.biology_initial_NO3 = 8.0
remora.biology_initial_NH4 = 0.1
remora.biology_initial_chlorophyll = 0.02
remora.biology_initial_phytoplankton = 0.08
remora.biology_initial_zooplankton = 0.06
remora.biology_initial_LdetritusN = 0.02
remora.biology_initial_SdetritusN = 0.04
remora.biology_initial_RdetritusN = 0.05
remora.biology_initial_LdetritusC = 0.002
remora.biology_initial_SdetritusC = 0.06
remora.biology_initial_TIC = 2100.0
remora.biology_initial_alkalinity = 2350.0
remora.biology_initial_RdetritusC = 0.3
remora.biology_initial_oxygen = 446.0
