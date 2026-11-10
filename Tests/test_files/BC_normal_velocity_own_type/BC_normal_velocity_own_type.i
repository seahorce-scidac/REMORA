# The East side is a slip wall for the normal velocity u and ubar. The tangential v and vbar
# use outflow, which gives them the same foextrap condition a slip wall would. The result
# must therefore match BC_normal_velocity_slipwall.i, where v and vbar are slip walls too.
FILE = BC_per_variable.i
#                        West      South       East      North
remora.bc.u.type      =  inflow    noslipwall  slipwall  slipwall
remora.bc.v.type      =  inflow    noslipwall  outflow   slipwall
remora.bc.ubar.type   =  inflow    noslipwall  slipwall  slipwall
remora.bc.vbar.type   =  inflow    noslipwall  outflow   slipwall
remora.plot_file      =  plt
