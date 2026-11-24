.. SPDX-FileCopyrightText: 2026 INRIA
..
.. SPDX-License-Identifier: Apache-2.0

.. _siconos_dependencies:

Siconos required and optional dependencies
==========================================

.. warning::

   This page might be outdated. For an up-to-date list of the Siconos dependencies, look at the configurations used for continuous integration.
   The Dockerfiles and requirements.txt files in the `ci_gitlab/dockerfiles directory <https://gricad-gitlab.univ-grenoble-alpes.fr/nonsmooth/siconos/-/tree/main/ci_gitlab/dockerfiles?ref_type=heads>`_
   of the Siconos sources show the list of required dependencies and how to install them.



* a compiler suite, with C++, C and gfortran compilers.

  C++20 support is required.

* cmake, version > 3.19 - https://cmake.org

  cmake is easy to install or update, whatever your system, e.g.:

  .. code-block:: bash

     python3 -m pip install cmake

* Eigen https://libeigen.gitlab.io/

* boost, version >= 1.75 (http://www.boost.org)

* BLAS and LAPACK (see :ref:`about_blas_lapack`)


To generate the documentation, you need:

* doxygen
* sphinx

For the python bindings:

* python (>= 3.8)

  We strongly recommend using Python virtual environments (e.g. `Micromamba <https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html>`_, `venv <https://docs.python.org/3/library/venv.html>`_...)

  See the details in :ref:`about_python`.


To run the tests:

* cppunit


.. _about_python:

About Python
------------

We strongly recommend using Python virtual environments (e.g. `Micromamba <https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html>`_, `venv <https://docs.python.org/3/library/venv.html>`_...)


Siconos venv example
""""""""""""""""""""

Installation:

.. code-block:: bash

   python3 -m venv $HOME/siconosenv
   source $HOME/siconosenv/bin/activate
   pip install -U -r requirements.txt


A requirements.txt file for Siconos is available in the source directory, in `ci_gitlab/dockerfiles <https://gricad-gitlab.univ-grenoble-alpes.fr/nonsmooth/siconos/-/tree/main/ci_gitlab/dockerfiles?ref_type=heads>`_

Activation:

.. code-block:: bash

   source $HOME/siconosenv/bin/activate


Siconos micromamba env example
""""""""""""""""""""""""""""""

Download the file `siconoslabenv.yml <https://gricad-gitlab.univ-grenoble-alpes.fr/nonsmooth/siconos/-/blob/main/ci_gitlab/sicolabenv.yml?ref_type=heads>`_ from the Siconos sources.

Change the name "base" to sicoenv in this file.


Installation:

.. code-block:: bash

   micromamba install -y -f sicolabenv.yml

Activation:

.. code-block:: bash

   micromamba activate sicoenv



.. warning::
   To get a proper mamba configuration and the right channel, you may need to run the following commands first:

   .. code-block::

      micromamba config prepend channels conda-forge
      micromamba config set channel_priority strict
      micromamba self-update



.. _about_blas_lapack:

About BLAS and LAPACK
---------------------

The BLAS (Basic Linear Algebra Subprograms, http://www.netlib.org/blas/) are routines that provide standard building blocks for performing basic vector and matrix operations, while LAPACK (http://www.netlib.org/lapack/#_presentation) provides routines for solving systems of simultaneous linear equations, least-squares solutions of linear systems of equations, eigenvalue problems, and singular value problems.
Different implementations are available, such as:

* openblas (http://www.openblas.net),
* the one from MKL (https://software.intel.com/en-us/intel-mkl),
* the Accelerate framework on macOS (https://developer.apple.com/library/prerelease/mac/documentation/Accelerate/Reference/BLAS_Ref/index.html)...

For Siconos we recommend:

* Accelerate on macOS,
* OpenBLAS + LAPACKE on Linux.

Warning: we do not support ATLAS.

Power users can still choose their favourite BLAS and LAPACK(E) vendor.

The BLAS/LAPACK setup of your system is checked when cmake runs.

If this fails, or if you need a specific implementation, the following variables can be given to cmake to help the search (see :ref:`siconos_detailed_install`):

* BLA_VENDOR: BLAS implementation type (also used as a hint for LAPACK),
  one of:
  * OpenBLAS
  * Matlab
  * Intel10_32 (intel mkl v10 32 bit)
  * Intel10_64lp (intel mkl v10+ 64 bit, threaded code, lp64 model)
  * Intel10_64lp_seq (intel mkl v10+ 64 bit, sequential code, lp64 model)
  * Intel10_64ilp (intel mkl v10+ 64 bit, threaded code, ilp64 model)
  * Intel10_64ilp_seq (intel mkl v10+ 64 bit, sequential code, ilp64 model)
  * Apple
  * Generic

* BLAS_ROOT: location of the BLAS implementation.
* LAPACK_ROOT: location of the LAPACK implementation.


About Boost
-----------

Boost provides many useful C++ libraries, in particular uBLAS, a C++ template class library that provides BLAS level 1, 2 and 3 functionality
for dense, packed and sparse matrices.

uBLAS is used in Siconos to define and implement matrices and vectors.

About Boost: http://www.boost.org/

About Ublas: http://www.boost.org/libs/numeric/ublas/doc/index.htm

Installation: a suitable Boost version comes with most Linux distributions, so no extra installation is usually required.

To get and install Boost, see
Boost Getting Started.

Note that we also use boost-bindings:
"Boost Bindings is a bindings library (not just) for Boost.Ublas. It offers an easy way of calling BLAS, LAPACK, UMFPACK, MUMPS and many other mature legacy numerical codes from within C++."

They are distributed and installed with Siconos, but you can also get the latest version here:
http://mathema.tician.de/software/boost-bindings

GMP
---

"GMP is a free library for arbitrary precision arithmetic, operating on signed integers, rational numbers, and floating point numbers ... "

This library usually comes with gcc. If not, see http://gmplib.org/ for download and installation instructions.
