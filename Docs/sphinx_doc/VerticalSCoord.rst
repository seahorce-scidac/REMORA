
 .. role:: cpp(code)
    :language: c++

Vertical S-coordinate
=====================

.. _'ROMS documentation': https://www.myroms.org/wiki/Vertical_S-coordinate

These equations are as given in the 'ROMS documentation'_

ROMSX has a generalized vertical, terrain-following, coordinate system.  Currently, the vertical transformation equation :math:'z=z\left( x,y,\sigma ,t\right) ', is available which can support vertical stretching 1D-functions when several constraints are satisfied.

Transformation Equation
-----------------------

.. math::
   z\left( x,y,\sigma ,t\right) &=\zeta \left( x,y,t\right) +\left[ \zeta \left( x,y,t\right) +h\left( x,y\right) \right] S\left( x,y,\sigma \right) ,\\
   S\left( x,y,\sigma \right) &=\frac{h_c \sigma +h \left( x,y\right) C\left( \sigma \right) }{h_c +h\left( x,y\right)}

where :math:'S\left( x,y,\sigma \right)' is a nonlinear vertical transformation functional, :math:'\zeta \left( x,y,t\right)' is the time-varying free-surface, :math:'h\left( x,y\right) ' is the unperturbed water column thickness and :math:'z=-h\left( x,y\right) ' corresponds to the ocean bottom, :math:'\sigma ' is a fractional vertical stretching coordinate ranging from :math:'-1\leq \sigma \leq 0,C\left( \sigma \right) ' is a nondimensional, monotonic, vertical stretching function ranging from :math:'-1\leq C\left( \sigma \right) \leq 0', and :math:'h_c' is a positive thickness controlling the stretching.  <!-- In sediment applications, :math:'h=h\left( x,y,t\right) ' is changed at every time-step since it is affected by erosion and deposition processes. -->

We find it convenient to define:

.. math::
   H_z \equiv \frac{\partial z}{\partial \sigma }

where :math:'H_z = H_z \left( x,y,\sigma ,t\right) ' are the vertical grid thicknesses. In ROMSX, :math:'H_z' is computed discretely as :math:'\delta z/\delta \sigma ' since this leads to the vertical sum of :math:'H_z' being exactly the total water column thickness :math:'D'.

Notice that,

.. math::
   S\left( x,y,\sigma \right) = \begin{cases}{0,& \text{if } \sigma = 0, & C\left( \sigma \right) = 0, & \text{at the free-surface;}}\\{-1 , & \text{if } \sigma = -1, & C\left( \sigma \right) = -1, & \text{at the ocean bottom}} \end{cases}

This transformation offers several advantages:

* Regardless of the design of :math:'C\left( \sigma \right) ', it behaves like equally-spaced sigma-coordinates in shallow regions, where :math:'h\left( x,y\right) \ll h_c'. This is advantageous because it avoids excessive resolution and associated CFL limitation in such areas.

* Near-surface refinement behaves more or less like geopotential coordinates in deep regions (level thicknesses, :math:'H_z', do not depend or weakly depend on bathymetry), while near-bottom like sigma coordinates (:math:'H_z' is roughly proportional to depth). This reduces the extreme r-factors near the bottom and reduces pressure gradient errors.

* The **true** sigma-coordinate system is recovered as :math:'h_c \rightarrow \infty '. This is useful when configuring applications with **flat** bathymetry and **uniform** level thickness. Practically, you can achieve this by setting tcline to **1.0d+16** in DataStruct.h. This will set :math:'h_c =1.0\times 10^{16}'. Although not necessary, we also recommend that you set :math:'\theta _S = 0' and :math:'\theta _B =0'.

In an undisturbed ocean state, :math:'\zeta \equiv 0', this transformation yields the following unperturbed depths, :math:'\hat{z}',

.. math::
   \hat{z} \left( x,y,\sigma \right) \equiv h\left( x,y\right) S\left( x,y,\sigma \right) =h\left( x,y\right) \left[ \frac{h_c \sigma +h\left( x,y\right) C\left( \sigma \right)}{h_c +h\left( x,y\right) } \right]

and

.. math::
   d\hat{z} =d\sigma h\left( x,y\right) \left[ \frac{h_c}{h_c +h\left( x,y\right) } \right]

As a consequence, the uppermost grid box retains very little dependency from bathymetry in deep areas, where :math:'h_c \ll h\left( x,y\right) '. For example if :math:'h_c =250m' and :math:'h\left( x,y\right) ' changes from :math:'2000' to :math:'6000m', the uppermost grid box changes only by a factor of 1.08 (less than 10%).



Vertical Stretching
-------------------

The above generic vertical transformation design facilitates numerous vertical stretching functions :math:'C\left( \sigma \right) '. This function is defined in terms of several parameters which are specified in the standard input file, inputs.

Stretching Function Properties
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* :math:'C\left( \sigma \right) ' is a dimensionless, nonlinear, monotonic function;
* :math:'C\left( \sigma \right) ' is a continuous differentiable function, or a differentiable piecewise function with smooth transition;
* :math:'C\left( \sigma \right) ' must be discretized in terms of the fractional stretched vertical coordinate :math:'\sigma ',

  .. math::
     \sigma \left( k \right) = \begin{cases}{\frac{k-N}{N}, & \text{at vertical }W\text{-points}, & k=0,\ldots ,N,}\\{\frac{k-N-0.5}{N}, & \text{at vertical}\rho \text{-points}, & k=1,\ldots ,N} \end{cases}

* :math:'C\left( \sigma \right) ' must be constrained by :math:'-1 \leq C\left( \sigma \right) \leq 0', that is,

  .. math::
     C\left( \sigma \right) = \begin{cases}{0, & \text{if } \sigma = 0, & C\left( \sigma \right) = 0, & \text{at the free-surface};}\\{-1, & \text{if } \sigma = -1, & C\left( \sigma \right) = -1, & \text{at the ocean bottom}.} \end{cases}

Available Stretching Functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
4. A. Shchepetkin (2010) UCLA-ROMS current function, Vstretching=4. :math:'C\left( \sigma \right) ' is defined as a continuous, double stretching function:

   Surface refinement function:
   .. math::
      C\left( \sigma \right) = \frac{1-\cosh \left( \theta _S \sigma \right) }{\cosh \left( \theta _S \right) -1}, & \text{for } \theta _S > 0, & C\left( \sigma \right) = - \sigma ^2, & \text{for }\theta _S \leq 0

   Bottom refinement function:
   .. math::
      C\left( \sigma \right) = \frac{\exp \left( \theta _B C\left( \sigma \right) \right) -1}{1-\exp \left( -\theta _B \right) }, & \text{for }\theta _B >0

   Notice that the bottom function is the second stretching of an already stretched transform. The resulting stretching function is continuous with respect to :math:'\theta _S' and :math:'\theta _B' as their values approach zero. The range of meaningful values for :math:'\theta _S' and :math:'\theta _B' are:
   .. math::
      0\leq \theta _S \leq 10 & \text{and} & 0\leq \theta _B \leq 4

   However, users need to pay attention to extreme r-factor (rx1) values near the bottom.

