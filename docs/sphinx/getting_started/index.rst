.. SPDX-FileCopyrightText: 2026 INRIA
..
.. SPDX-License-Identifier: Apache-2.0

.. _siconos_getting_started:

Quickstart
==========

This manual gives a short overview of the features of the software and some guidelines for your first steps in the simulation of nonsmooth systems
with Siconos. You will learn the basics needed to describe a nonsmooth problem and run a simulation with Siconos, with both APIs, C++ and Python.
For a more detailed description of Siconos and its features, see the :ref:`full_documentation` table of contents.


What is Siconos?
----------------

Siconos is open-source scientific software primarily aimed at modeling and simulating nonsmooth dynamical systems :cite:`Acary.Brogliato2008`:

* **Mechanical systems** (rigid or deformable) with unilateral contact, Coulomb friction and impacts
  (nonsmooth mechanics, contact dynamics, multibody system dynamics, granular materials).
* **Switched electrical circuits**, such as circuits with ideal and piecewise linear components: power converters, rectifiers, Phase-Locked Loops (PLL) or Analog-to-Digital converters.
* Sliding mode control systems.
* Gene regulatory networks in biology.

Other applications are found in systems and control (hybrid systems, differential inclusions,
optimal control with state constraints), optimization (complementarity systems and variational inequalities),
fluid mechanics, computer graphics, etc.

See the :ref:`siconos_examples` for an overview of the many problems handled with Siconos.


Try it
------

You can use one of the end-user Docker images available in the `siconos-tutorials registry <https://gricad-gitlab.univ-grenoble-alpes.fr/nonsmooth/siconos-tutorials/container_registry>`_.

For example, start a Jupyter container with Siconos, all the working examples and some tutorial notebooks:

.. code-block:: bash


   docker run -p 8888:8888 -ti gricad-registry.univ-grenoble-alpes.fr/nonsmooth/siconos-tutorials/siconoslab-main

   # Then open http://localhost:8888/... in your browser.
   # The exact address is printed by the command above.


Another (less stable) way to get started with Siconos is to try the tutorial notebooks available in the siconos-notebooks directory:

.. image:: https://plmbinder.math.cnrs.fr/binder/badge_logo.svg
   :target: https://plmbinder.math.cnrs.fr/binder/v2/git/https%3A%2F%2Fgricad-gitlab.univ-grenoble-alpes.fr%2Fnonsmooth%2Fsiconos-tutorials.git/HEAD

This page provides an interactive Python interface (notebooks) to Siconos that you can use online, in your web browser.


Siconos usage in a few steps
----------------------------

#.  :ref:`download` the software
#.  :ref:`siconos_install_guide`
#.  Write and run your first simulation.

There are two ways to use Siconos:

* As a library (C++): you write a C++ file (driver) that describes and solves your problem,
  and then run siconos to build and execute your program.
* As a Python package.

The Python API is generated (with pybind11) from the C++ code, so both APIs are nearly equivalent, although the C++ one may be more complete.

For new users, we recommend the Python API, which is easier to understand.

Below are two examples (Python and C++) of a Siconos program. They simply build and print a first-order dynamical
system (see :cpp:class:`siconos::modeling::FirstOrderLinearDS`).

.. admonition:: tip

   Have a look at the many `examples in the siconos-tutorials repository <https://gricad-gitlab.univ-grenoble-alpes.fr/nonsmooth/siconos-tutorials/>`_ for inspiration, and use the search bar to find details on all (Python, C or C++) objects.

**Example of Siconos Python API usage**

.. code-block:: python

   # File run.py
   # import siconos package
   import siconos.modeling as sm
   # import numpy package
   import numpy as np

   # Create a dynamical system
   size = 10
   x0 = np.random.random(size)
   ds = sm.FirstOrderLinearDS(x0, sm.copy_t)
   print(ds)

And:

.. code-block:: bash

   siconos run.py
   # --> calls python and anything else required

.. admonition:: tip

   For more on the *siconos* script options, try:

   .. code-block:: shell

      siconos --help
      siconos --info

**Example of Siconos C++ API usage**

Write a C++ file, e.g. run.cpp:

.. code-block:: c++

   // File run.cpp
   #include "SiconosKernel.hpp"
   int main()
   {
    unsigned int size = 10;
    auto x0 = std::make_shared<siconos::algebra::SiconosVector>(size);
    auto A = std::make_shared<siconos::algebra::SiconosDenseMatrix>(size, size);
    A->randomize();

    auto ds = std::make_shared<siconos::modeling::FirstOrderLinearDS>(x0, A);
    ds->display();
   }

Then compile, link and execute it (in one go, with the siconos script):

.. code-block:: shell

   siconos run.cpp

For a new simulation, start from the :ref:`template for a C++ driver file <template_siconos_driver>`
or adapt one of the examples available at https://gricad-gitlab.univ-grenoble-alpes.fr/nonsmooth/siconos-tutorials.git.


More
----

* :doc:`users_guide/index`
* :doc:`reference/index`
* :doc:`tutorial_python/index`
* :doc:`tutorial_cpp/siconos_tutorial`
* :ref:`running_siconos`
