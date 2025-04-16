
 .. role:: cpp(code)
    :language: c++

.. _Climatology:


Nudging to climatology
======================

REMORA supports nudging to climatology both on the boundaries and in an interior sponge region. Nudging at the boundaries is discussed in :ref:`Boundary Conditions`<sec:domainBCs>`.

The climatology history and nudging coefficients are read from NetCDF files.

Nudging for Tracers
-------------------

Tracers :math:`S`  are nudged towards solution by:

.. math::
  S = S + C_{\mathrm{tracer}} \Delta t \left(S_{\mathrm{clim}} - S\right)

where :math:`S_{mathrm{clim}}` is the value given in the climatology and :math:`C_{\mathrm{tracer}}` is the nudging coefficient for the relevant tracer.

Nudging for 2D velocity
-----------------------

The 2D velocity components are udpated in the right-hand side term ``rhs_ubar`` according to:

.. math::
   \mathrm{RHS}_{\overline{u}} = \mathrm{RHS}_{\overline{u}} + C_{\mathrm{M2}} \frac{D}{mn} \left(\overline{u}_{\mathrm{clim}} - \overline{u}_{\mathrm{krhs}}\right)

or similar for v. :math:`C_{\mathrm{M2}}` is the nudging coefficient, :math:`D` is the water column depth, and :math:`\overline{u}_{\mathrm{clim}}` is the climatology value read from file. :math:`\mathrm{RHS}_{\overline{u}}` is the variable ``rhs_ubar`` and :math:`\overline{u}_{\mathrm{krhs}}` is the value of :math:`\overline{u}` at the time index ``krhs``.

Nudging for 3D velocity
-----------------------

The 3D velocity components are udpated in the right-hand side term ``ru`` according to:

.. math::
   \mathrm{RHS}_{u} = \mathrm{RHS}_{u} + C_{\mathrm{M3}} \frac{H_{z}}{mn} \left(u_{\mathrm{clim}} - u_{\mathrm{old}}\right)

or similar for v. :math:`C_{\mathrm{M3}}` is the nudging coefficient and :math:`u_{\mathrm{clim}}` is the climatology value read from file. :math:`\mathrm{RHS}_{u}` is the variable ``ru`` and :math:`u_{\mathrm{old}}` is the value of :math:`u` at the prior time step.

