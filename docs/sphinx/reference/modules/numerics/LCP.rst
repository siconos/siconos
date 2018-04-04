LCP
---

Usage
^^^^^

In python the LinearComplementarityProblem from C API is renamed :index:`LCP`.

The solution of the problem exposed in :ref:`intro-lcp`::

  import siconos.numerics as Numerics
  lcp = Numerics.LCP([[2., 1.], [1., 2.]], [-5., -6.])

can be reached with the ``Numerics`` module by first providing a guess. This guess is made of the two vectors ``z`` and ``w``. As it will be an input as well
as an output parameter for the solver, we must use a numpy array that can be modified in place. Standard python sequences cannot be modified
in this interface and are not suitable for this kind of parameter::

  from numpy import array
  z = array([0., 0.])
  w = array([0., 0.])

We must also provide a solver options object::

  SO = Numerics.SolverOptions(lcp, Numerics.SICONOS_LCP_LEMKE)

That object brings the decision to use a Lemke method. It also allows
the manipulation of the different solver parameters such as the wanted
precision. The ``integer`` option parameters are in an ``iparam`` array and
the ``double`` option parameters are in a ``dparam`` array.

For the Lemke method the default precision is the first element of ``SO.dparam`` array::

  >>> SO.dparam
  array([  1.00000000e-06,   0.00000000e+00,   0.00000000e+00,
           0.00000000e+00,   0.00000000e+00])


The solution may be now computed with a Lemke solver::

  info = Numerics.lcp_lexicolemke(lcp, z, w, SO)

The ``info`` output is an integer. For all ``Numerics`` solvers, ``0`` means a successfull resolution. ``z`` and ``w`` now contain
a correct solution for the asked precision::

  >>> z
  array([ 1.33333333,  2.33333333])
  >>> w
  array([ 0.,  0.])


We may then compute the error::

  d = Numerics.lcp_compute_error(lcp, z, w, 1e-6)
  [ ... ]
  


LCP API
^^^^^^^

.. automodule:: siconos.numerics
  :members: :eval:`starting_with(['LCP', 'lcp'])`
