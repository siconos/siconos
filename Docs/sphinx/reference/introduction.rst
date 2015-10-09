.. _siconos_python_introduction:


Introduction
============
	      
The different ``siconos`` libraries may be imported as python
modules. To import all Siconos modules, one may write for example:
 
.. testcode::

  import siconos.numerics
  import siconos.kernel
  import siconos.io
  import siconos.mechanics
  import siconos.control

The python bindings for siconos try to be as close as possible to the ``C`` and
``C++`` API of the different libraries, but some discrepancies
exist. These are explained in the sequel of this introduction.

.. _intro-lcp:

Dense matrices and vectors
--------------------------

Standard `python <https://www.python.org/>`_ sequences as well as
`numpy <http://www.numpy.org/>`_ arrays and matrices may be used as
input for matrices and vectors in the ``C`` and ``C++`` API. For
example, with ``numerics`` module a linear complementarity problem :

.. math::

      \begin{cases}
        w=Mz+q \\
        0 \leq w \perp z\geq 0
      \end{cases}

with 

.. math::

   M = \begin{pmatrix}
   2 & 1 \\
   1 & 2 \\
   \end{pmatrix}

and

.. math::

   q = \begin{pmatrix}
        -5 \\
        -6
   \end{pmatrix}

may be declared like this:

.. testcode::

  import siconos.numerics as numerics
  lcp = numerics.LCP([[2., 1.], [1., 2.]], [-5., -6.])

In this ``python`` declaration, the ``M`` matrix of the linear
complementary problem is:

.. testcode::

  M = [[2., 1.],[1., 2.]]

It is a simple ``python`` list of list of float numbers where the rows
of the matrix are the inner lists.

The ``q`` vector of the linear complementary problem is:

.. testcode::

  q = [-5., -6.]

It is a ``python`` list of float numbers.

numpy arrays may also be used in inputs:

.. testcode::

  import siconos.numerics as numerics
  import numpy as np
  
  M = np.array([[2., 1.],[1., 2.]])
  q = np.array([-5., -6.])
  
  lcp = numerics.LCP(M, q)

This is important: where the parameter is an input as well as an
output parameter, numpy array *must* be used!


With the ``kernel`` and ``mechanics`` modules every ``SimpleMatrix``
and ``SiconosVector`` may be replaced by ``python`` standard sequences
or ``numpy`` arrays. For example, we can build a Lagrangian dynamical
system with 3 degrees of freedom like this:

.. testcode::

  from siconos.kernel import LagrangianDS

  position = [0, 0, 0]
  velocity = [0, 0, 0]
  mass = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]

  lds = LagrangianDS(position, velocity, mass)

It is possible to use ``SimpleMatrix`` and ``SiconosVector`` arguments:

.. testcode::

  from siconos.kernel import LagrangianDS, SimpleMatrix, SiconosVector

  position = SiconosVector(3)
  position.zero()

  velocity = SiconosVector(3)
  velocity.zero()

  mass = SimpleMatrix(3,3)
  mass.eye()

  lds =  LagrangianDS(position, velocity, mass)
  
Please note that ``kernel.SimpleMatrix`` and ``kernel.SiconosVector``
objects cannot be used as arguments to ``numerics`` module functions.
The instantiation of previous ``numerics.LCP`` can only be done with
standard ``python`` sequences or ``numpy`` arrays.



Sparse matrices
---------------

`Scipy <http://www.scipy.org/>`_ sparse matrices must be used in input
where cs_sparse is needed in the ``C`` API.

Here is for example the conversion from a sparse compressed column
matrix into a sparse block matrix with blocks of 3 rows and 3 columns:

.. testcode::

  import scipy.sparse
  import siconos.numerics as numerics

  # create a sparse compressed column matrix
  m = scipy.sparse.csc_matrix([[1,2,3,4,5,6],[1,2,3,4,5,6],[1,2,3,4,5,6],[1,2,3,4,5,6],[1,2,3,4,5,6],[1,2,3,4,5,6]])

  result = numerics.sparseToSBM(3, m)
  
  # result[0] is the info result and should be 0
  # result[1] is the sparse block matrix
  sbm = result[1]

  # print the matrix
  numerics.printSBM(sbm)

Omitted parameters
------------------

This concerns the ``C`` API of the ``numerics`` library:

 1. Where the size of an input vector may be inferred, the size must not be given in the arguments list. 

 2. Output only parameters given in the argument list in the ``C`` API are ``python`` return parameters

Here is an example that shows both cases:

.. testcode::

  # the C signature:
  #   void frictionContact3D_AlartCurnierFunction(
  #     unsigned int problemSize,
  #     double *reaction,
  #     double *velocity,
  #     double *mu,
  #     double *rho,
  #     double *result,
  #     double *A,
  #     double *B)  # 
 
  from numpy import array
  from siconos.numerics import frictionContact3D_AlartCurnierFunction
 
  mu = array([0.1])
  reactions = array([1., 1., 1.])
  velocities = array([1., 1., 1.])
  rho = array([1., 1., 1.])
 
  # problemSize is omitted in python call as it can be infered from the
  # size of given vectors (reactions, velocities, rho)
  # result A and B are only given for output so are Python return parameters 
 
  result,A,B = frictionContact3D_AlartCurnierFunction(reactions, velocities, mu, rho)




C++ Visitors
------------

``Siconos`` C++ visitors are not binded. The class of a returned object on the ``python`` side is always the ``true`` class and never a more general class, so the visitor pattern is not relevant here:

.. testcode::

  import siconos.kernel as K
  dsA = K.LagrangianDS([0],[0],[[1]])
  dsB = K.FirstOrderLinearDS([0],[[1]])
  model = K.Model(0, 0)
  model.nonSmoothDynamicalSystem().insertDynamicalSystem(dsA)
  model.nonSmoothDynamicalSystem().insertDynamicalSystem(dsB)

  assert(type(model.nonSmoothDynamicalSystem().dynamicalSystem(dsA.number())) == K.LagrangianDS)
  assert(type(model.nonSmoothDynamicalSystem().dynamicalSystem(dsB.number())) == K.FirstOrderLinearDS)

Shared pointers
---------------

For ``Siconos`` C++ libraries (``kernel``, ``io``, ``mechanics``, ``control``) the
shared pointer mechanisms is totally hidden and the namespaces SP,
SPC, SPA are not present in the ``Python`` modules.

Other differences specific to ``Siconos Libraries`` are documented in relevant sections.




	      
