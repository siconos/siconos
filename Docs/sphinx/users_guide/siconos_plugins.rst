.. _siconos_plugins:

User-defined plugins
====================

*Work in progress ...*


Siconos proposes a 'plugin' system that allows users to provide their own function(s) to describe some specific behavior for
some classes components.

For example, consider a lagrangian linear dynamical system, where :math:`M\ddot q + C \dot q + K q =  F_{Ext}(t,z) + p`.

Suppose you want to set :math:`F_{Ext}(t,z) = cos(t)`, then you can define a C function to compute this cosine and 'plug' it to
the dynamical system so that each time the system needs to compute its external forces, your cosine function will be called.


At the time, plug-in are available for :doxysiconos:`DynamicalSystems` and :doxysiconos:`Relation`. For both of them and for their derived classes, a list
of the variables that can be plugged is given in :ref:`ds_plugins` and :ref:`relation_plugins`.

* find the variable you want to plug and check what is the expected list of arguments for a function plugged to this variable
(check :ref:`ds_plugins` or :ref:`relation_plugins`).

* write a C function::

    extern "C" external_forces(double time, int size, double* fext, int zsize, double *z)
    {
    for(int i=0;i<size;++i)
       (*fext)(i) = cos(time);
    }

* connect your function to the variable. For each 'plugable' variable, a setComputeVARNAMEFunction exists ::

    ds->setComputeFExtFunction('myPlugin', 'external_forces');
    // ...
    ds->computeFExt(2.)
    // --> call external_forces with time == 2.
    

Plugins overview
----------------

.. csv-table:: plugins in siconos classes
   :header: "Class Name", "operator", "plugin name", "signature"
   :widths: 10 5 5 40

   :doxysiconos:`DynamicalSystem`, ":math:`g(\dot x, x, t, z)`", g, "``(double time, int size, double* fext, int zsize, double *z)``"
   :doxysiconos:`LagrangianLinearTIDS`, ":math:`F_{Ext}(t,z)`", FExt, "``(double time, int size, double* fext, int zsize, double *z)``"
   :doxysiconos:`FirstOrderR`, ":math:`h(x,t,\lambda,z)`", h, "``(double time, int x.size, double * x, int lambda.size, double * lambda, double * y, int z.size, double * z)``"
