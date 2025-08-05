The DogboneAnalytic problem is an analytic-initialized version of the
standard dogbone test in `ROMS <https://www.myroms.org/wiki/Test_Cases>`_. The
``inputs`` file matches the unrefined ``dogbone_whole`` case. ``inputs_ml`` does
the same but with refinement where the absolute x-velocities are in excess
of 0.05. ``inputs_ml_quad`` refines on the lower left quadrant of the grid.
