Introduction
============

This is the documentation of ``Siconos Front-End`` which is a `python <https://www.python.org/>`_ interface to the `Siconos <http://siconos.gforge.inria.fr/>`_ libraries.

The different ``Siconos`` libraries may be imported as python modules. To import all Siconos modules, one may write for example::
  
  import Siconos.Numerics
  import Siconos.Kernel
  import Siconos.IO
  import Siconos.Mechanics

``Siconos Front-End`` try to be as close as possible to the ``C`` and ``C++``
API of the different libraries, but some discrepancies exist. These are explained in the sequel of this introduction.

.. _intro-lcp:

Dense matrices and vectors
--------------------------

Standard `python <https://www.python.org/>`_ sequences as well as `numpy <http://www.numpy.org/>`_ arrays and matrices may be used as input for matrices and vectors in the ``C`` and ``C++`` API. For example, with ``Numerics`` module a linear complementary problem :

.. math::
   \begin{eqnarray*}
      \begin{cases}
        w=Mz+q \\
        0 \leq w \perp z\geq 0
      \end{cases}
    \end{eqnarray*}

with 

.. math::

   \begin{eqnarray*}
        M = \begin{pmatrix}
        2 & 1 \\
        1 & 2 \\
        \end{pmatrix}
    \end{eqnarray*}

and

.. math::
   
   \begin{eqnarray*}
        q = \begin{pmatrix}
        -5 \\
        -6
        \end{pmatrix}
    \end{eqnarray*}

may be declared like this::
  
  import Siconos.Numerics as Numerics
  lcp = Numerics.LCP([[2., 1.], [1., 2.]], [-5., -6.])

In this ``python`` declaration, the ``M`` matrix of the linear
complementary problem is::

  M = [[2., 1.],[1., 2.]]

It is a simple ``python`` list of list of float numbers where the rows
of the matrix are the inner lists.

The ``q`` vector of the linear complementary problem is::

  q = [-5., -6.]

It is a ``python`` list of float numbers.

numpy arrays may also be used in inputs::

  import Siconos.Numerics as Numerics
  import numpy as np
  
  M = np.array([[2., 1.],[1., 2.]])
  q = np.array([-5., -6.])
  
  lcp = Numerics.LCP(M, q)

This is important: where the parameter is an input as well as an
output parameter, numpy array *must* be used!


With the ``Kernel`` and ``Mechanics`` modules every ``SimpleMatrix`` and
``SiconosVector`` may be replaced by ``python`` standard sequences or
``numpy`` arrays. For example, we can build a lagrangian dynamical system with 3 degrees of freedom like this::

  from Siconos.Kernel import LagrangianDS

  position = [0, 0, 0]
  velocity = [0, 0, 0]
  mass = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]

  lds = LagrangianDS(position, velocity, mass)

It is possible to use ``SimpleMatrix`` and ``SiconosVector`` arguments::

  from Siconos.Kernel import LagrangianDS, SimpleMatrix, SiconosVector

  position = SiconosVector(3)
  position.zero()

  velocity = SiconosVector(3)
  velocity.zero()

  mass = SimpleMatrix(3,3)
  mass.eye()

  lds =  LagrangianDS(position, velocity, mass)
  
Please note that ``Kernel.SimpleMatrix`` and ``Kernel.SiconosVector`` objects cannot be used as arguments to ``Numerics`` module functions.
The instanciation of previous ``Numerics.LCP`` can only be done with standard ``python`` sequences or ``numpy`` arrays.



Sparse matrices
---------------

`Scipy <http://www.scipy.org/>`_ sparse matrices must be used in input where cs_sparse is needed in the ``C`` API.

Here is for example the conversion from a sparse compressed column matrix into a sparse block matrix with blocks of 3 rows and 3 columns::

  import scipy.sparse
  import Siconos.Numerics as Numerics

  # create a sparse compressed column matrix
  m = scipy.sparse.csc_matrix([[1,2,3,4,5,6],[1,2,3,4,5,6],[1,2,3,4,5,6],[1,2,3,4,5,6],[1,2,3,4,5,6],[1,2,3,4,5,6]])

  result = Numerics.sparseToSBM(3, m)
  
  # result[0] is the info result and should be 0
  # result[1] is the sparse block matrix
  sbm = result[1]

  # print the matrix
  Numerics.printSBM(sbm)

Size of matrices and vectors in arguments
-----------------------------------------

This concerns the ``C`` API of the ``Numerics`` library. Where the size of a an input vector may be infered, it must not be given in the arguments list::
     
  example...

Output parameters
-----------------

This concerns the ``C`` API of the ``Numerics`` library. In the
``Numerics`` module, output parameters given in the argument list in
the ``C`` API are ``python`` return parameters::

  example...


C++ Visitors
------------

``Siconos`` C++ visitors are not binded. The class of a returned object on the ``python`` side is always the ``true`` class and never a more general class, so the visitor pattern is not relevant here::

     example...


Shared pointers
---------------

For ``Siconos`` C++ libraries (``Kernel``, ``IO``, ``Mechanics``)
the shared pointer mechanisms is totally hidden.

Other differences specific to ``Siconos Libraries`` are documented in relevant sections.



