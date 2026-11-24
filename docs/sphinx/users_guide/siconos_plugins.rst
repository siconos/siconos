.. SPDX-FileCopyrightText: 2026 INRIA
..
.. SPDX-License-Identifier: Apache-2.0

.. _siconos_plugins:


User-defined plugins
--------------------


Siconos provides a 'plugin' system that lets users give their own functions to describe specific behaviours of
some class components.

For example, consider a linear Lagrangian dynamical system, :math:`M\ddot q + C \dot q + K q =  F_{Ext}(t) + p`.

Suppose you want to set :math:`F_{Ext}(t) = cos(t)`. You can define a function that computes this cosine and 'plug' it into
the dynamical system: each time the system computes its external forces, your function is called.


Plugins are currently available for :cpp:class:`DynamicalSystems`, :cpp:class:`Relation` and their derived classes.

**How to provide user-defined functions?**

When the plugin mechanism is available for an operator, the class has two methods:

* setCompute<VARNAME>Function: sets the function to use (any C++ function, lambda function, etc.)
* compute<VARNAME>

Then:

1. Look up the class documentation, find setCompute<VARNAME>Function and identify the arguments of the required user-defined function.

   All types are defined in :ref:`FunctionTypes.hpp <_file_kernel_src_modeling_FunctionTypes.hpp>`.

   Example for :cpp:class:`siconos::modeling::LagrangianDS`: :cpp:func:`siconos::modeling::LagrangianDS::setComputeFextFunction` and :cpp:func:`siconos::modeling::LagrangianDS::computeFext`.

2. Define your function as you would any other function.

3. Call setCompute<VARNAME>Function.

That's it!

Example:

.. code-block:: c++

   siconos::modeling::LagrangianDS lds{q0, velocity0, siconos::algebra::alias_t};
   lds.setComputeFextFunction(
      [](double time, Eigen::Ref<siconos::algebra::MapVectorType> result) {
        auto i = 0;
        for (auto& v : result) v = time;
      });
   ...
   // ...
   ds.computeFext(2.)
   // --> calls the user function with time == 2.

.. tip::

   Some good places to find complete examples of user-defined functions:

   - Test files in the Siconos repository, e.g. `LagrangianDSTest.cpp <https://gricad-gitlab.univ-grenoble-alpes.fr/nonsmooth/siconos/-/blob/main/kernel/src/modelingTools/test/LagrangianDSTest.cpp?ref_type=heads>`_
   - Search for setCompute... in the `siconos-tutorials repository <https://gricad-gitlab.univ-grenoble-alpes.fr/nonsmooth/siconos-tutorials>`_
