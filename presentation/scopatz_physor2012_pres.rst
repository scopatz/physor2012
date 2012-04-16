Physor 2012
==============================

.. container:: main-title

    INTERPOLATIONS OF NUCLIDE-SPECIFIC SCATTERING KERNELS GENERATED WITH SERPENT

.. container:: main-names


    April 17th, 2012, Physor, Knoxville, TN

    Anthony Scopatz - The University of Chicago

    Erich Schneider - The University of Texas at Austin

    scopatz@flash.uchicago.edu


Goals
==============================
**Goal 1:** 

    Find the neutron scattering kernel

.. container:: align-center

    .. math::

        \sigma_{s,g\to h,i}

    for all incident and exiting energy groups and species:

.. container:: align-center

    .. math::

        g,h \in G, \, i \in I


Goals
==============================
**Goal 1:** 

    Find the neutron scattering kernel

.. container:: align-center

    .. math::

        \sigma_{s,g\to h,i}

    for all incident and exiting energy groups and species:

.. container:: align-center

    .. math::

        g,h \in G, \, i \in I

**Goal 2:** 

    Show that we may interpolate between related kernels.


Motivation
==============================
- Multigroup solvers that employ a perturbation-based approach to generate group 
  constants from a small number of pre-computed cases require nuclide specific 
  scattering cross sections.

Motivation
==============================
- Multigroup solvers that employ a perturbation-based approach to generate group 
  constants from a small number of pre-computed cases require nuclide specific 
  scattering cross sections.

* If we may interpolate between kernels, then the method by which 
  we obtain the scattering kernel (typically Monte Carlo) may be performed 
  less often. 

Motivation
==============================
- Multigroup solvers that employ a perturbation-based approach to generate group 
  constants from a small number of pre-computed cases require nuclide specific 
  scattering cross sections.

* If we may interpolate between kernels, then the method by which 
  we obtain the scattering kernel (typically Monte Carlo) may be performed 
  less often. 

- The major application for this is *the Fuel Cycle*.
