
# Flather on every side, but zeta on the East radiates without reading the boundary file.
# Flather there would relax ubar toward an elevation that is never filled, so setup must stop.
#                        West      South       East          North
remora.bc.ubar.type   =  flather   flather     flather       flather
remora.bc.vbar.type   =  flather   flather     flather       flather
remora.bc.zeta.type   =  chapman   chapman     orlanski_rad  chapman
