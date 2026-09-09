.. _siconos_plugins:


User-defined plugins
--------------------


Siconos proposes a 'plugin' system that allows users to provide their own function(s) to describe some specific behavior for
some classes components.

For example, consider a lagrangian linear dynamical system, where :math:`M\ddot q + C \dot q + K q =  F_{Ext}(t) + p`.

Suppose you want to set :math:`F_{Ext}(t) = cos(t)`, then you can define a C function to compute this cosine and 'plug' it to
the dynamical system so that each time the system needs to compute its external forces, your cosine function will be called.


At the time, plug-in are available for :cpp:class:`DynamicalSystems`, :cpp:class:`Relation` and derived classes.

**How to provide user-defined functions ?**

When pluggin mechanism is available for an operator, two methods exist in the class:

* setCompute<VARNAME>Function: to define the function to be used (any other c++ function, lambda function ...)
* compute<VARNAME>

So 

1. Check the documentation class, find de the setCompute<VARNAME>Function and identify the arguments of the required user-defined function

   All types are defined in :ref:`FunctionTypes.hpp <_file_kernel_src_modeling_FunctionTypes.hpp>`.

   Example for :cpp:class:`siconos::modeling::LagrangianDS`: :cpp:func:`siconos::modeling::LagrangianDS::setComputeFextFunction` and :cpp:func:`siconos::modeling::LagrangianDS::computeFext`.

2. Define your function as you would for any other function

3. call setCompute<VARNAME>Function

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
   // --> call external_forces with time == 2.
   
.. tip::
   
   Some good places to find complete examples of user-defined functions usage:

   - Tests files in siconos repository; e.g. `LagrangianDSTest.cpp <https://gricad-gitlab.univ-grenoble-alpes.fr/nonsmooth/siconos/-/blob/main/kernel/src/modelingTools/test/LagrangianDSTest.cpp?ref_type=heads>`_
   - Try to find setCompute... inside `siconos-tutorials repository <https://gricad-gitlab.univ-grenoble-alpes.fr/nonsmooth/siconos-tutorials>`_ 

