
 .. role:: cpp(code)
    :language: c++


.. _Equations:

Prognostic Equations
====================

.. _`ROMS documentation`: https://www.myroms.org/wiki/Equations_of_Motion

These equations are as given in the `ROMS documentation`_

The momentum balance in the :math:`x`- and :math:`y`-directions are:

.. math::
   {\frac {\partial u}{\partial t}}+{\vec {v}}\cdot \nabla u-fv=-{\frac {\partial \phi }{\partial x}}-{\partial \over \partial z}\left({\overline{u'w'}}-\nu {\partial u \over \partial z}\right)+{\cal {F}}_{u}+{\cal {D}}_{u}

   {\frac {\partial v}{\partial t}}+{\vec {v}}\cdot \nabla v+fu=-{\frac {\partial \phi }{\partial y}}-{\partial \over \partial z}\left({\overline{v'w'}}-\nu {\partial v \over \partial z}\right)+{\cal {F}}_{v}+{\cal {D}}_{v}

The time evolution of a scalar concentration field, :math:`C(x,y,z,t)`, e.g., salinity, temperature, or nutrient species,
is governed by the advective-diffusive equation:

.. math::
   {\frac {\partial C}{\partial t}}+{\vec {v}}\cdot \nabla C=-{\partial \over \partial z}\left({\overline{C'w'}}-\nu _{\theta }{\partial C \over \partial z}\right)+{\cal {F}}_{C}+{\cal {D}}_{C}

The equation of state is given by:

.. math::
   \rho =\rho (T,S,P)

In the Boussinesq approximation, density variations are neglected in the momentum equations except in their contribution to the buoyancy force in the vertical momentum equation. Under the hydrostatic approximation, it is further assumed that the vertical pressure gradient balances the buoyancy force:

.. math::
   {\frac {\partial \phi }{\partial z}}=-{\frac {\rho g}{\rho _{o}}}

The final equation expresses the continuity equation for an incompressible fluid:

.. math::
   {\frac {\partial u}{\partial x}}+{\frac {\partial v}{\partial y}}+{\frac {\partial w}{\partial z}}    = 0

In the ocean, vertical mixing due to molecular viscosity is extremely weak compared to the turbulent mixing,
so the terms involving :math:`\nu` and :math:`\nu_\theta` can be neglected.

To close these questions, parametrizations of the Reynolds stresses and turbulent tracer fluxes are introduced as functions of the other
variables and vertical turbulent eddy viscosity and eddy diffusivity coefficients :math:`K_m` and :math:`K_C`, respectively.

.. math::
    \overline{u^\prime w^\prime} = -K_M \frac{\partial u}{\partial z}; \hspace{0.5in} \overline{v^\prime w^\prime} = -K_M \frac{\partial v}{\partial z}; \hspace{0.5in} \overline{C^\prime w^\prime} = -K_C \frac{\partial C}{\partial z};

An overbar represents a time average and a prime represents a fluctuation about the mean.

See :ref:`sec:VerticalMixing` for how :math:`K_m` (:math:`K_v`) and :math:`K_C` (:math:`K_t`) are computed in REMORA.

Vertical Boundary Conditions
----------------------------

The vertical boundary conditions can be presribed as follows.

At the top:

.. math::
   K_m \frac{\partial{u}}{\partial{z}} &= \tau_s^x(x,y,t)

   K_m \frac{\partial{v}}{\partial{z}} &= \tau_s^y(x,y,t)

   K_C \frac{\partial{C}}{\partial{z}} &= \frac{Q_C}{\rho_{o} c_{P}}

   w &= \frac{\partial{\zeta}}{\partial{t}}

At the bottom:

.. math::
   K_m \frac{\partial{u}}{\partial{z}} &= \tau_b^x(x,y,t)

   K_m \frac{\partial{v}}{\partial{z}} &= \tau_b^y(x,y,t)

   K_C \frac{\partial{C}}{\partial{z}} &= 0

   -w + \vec{v} \cdot \nabla h &= 0

where :math:`Q_C` is the surface concentration flux of temperature or salinity, :math:`\tau_s^x` and :math:`\tau_s^y` are the surface wind stresses, and :math:`\tau_b^x` and :math:`\tau_b^y` are the bottom stresses. The bottom stress can have a linear, quadratic, or log layer form. The surface stresses are set either set explicitly or calculated from prescribed winds based on the bulk flux parameterization (optional). The temperature flux is also set with the parametrization, which accounts for atmospheric conditions and solar irradiation. Calculated net evaporation/precipitation determines the salt flux.
