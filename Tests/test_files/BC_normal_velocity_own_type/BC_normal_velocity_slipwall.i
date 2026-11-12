# Both velocity components are slip walls on the East side; see BC_normal_velocity_own_type.i.
FILE = BC_per_variable.i
#                        West      South       East      North
remora.bc.u.type      =  inflow    noslipwall  slipwall  slipwall
remora.bc.v.type      =  inflow    noslipwall  slipwall  slipwall
remora.bc.ubar.type   =  inflow    noslipwall  slipwall  slipwall
remora.bc.vbar.type   =  inflow    noslipwall  slipwall  slipwall
remora.plot_file      =  plt_slipwall
