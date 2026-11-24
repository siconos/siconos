.. SPDX-FileCopyrightText: 2026 INRIA
..
.. SPDX-License-Identifier: Apache-2.0

.. _siconos_algebra:

Linear Algebra in Siconos
=========================

.. _kernel_matrix_storage:

Matrices and vectors in the kernel component
--------------------------------------------

In the Siconos kernel, all matrices and vectors are interfaces to `Eigen <https://libeigen.gitlab.io/eigen/docs-nightly/group__TutorialMatrixClass.html>`_ vectors and matrices
and are defined in the files :file:`SiconosMatrix.hpp` and :file:`SiconosVector.hpp`.

Please refer to the Eigen documentation for details.

Siconos wrappers are:

* :cpp:class:`SiconosVector`: dense vector of real numbers (double precision)
* :cpp:class:`BlockVector`: vector of :cpp:class:`SiconosVector`
* :cpp:class:`SiconosDenseMatrix`: dense matrix of doubles
* :cpp:class:`SiconosSparseMatrix`: sparse matrix of doubles, column-major storage

and some fixed-size objects, like :cpp:class:`SiconosMatrix33`, :cpp:class:`SiconosMatrix66`, etc.

Note that a BlockVector is no more than a collection of pointers to SiconosVector objects,
so in most cases you build one by inserting existing vectors.


.. _numerics_matrix_storage:

Matrix storage in the numerics component
----------------------------------------

The numerics component provides different ways to store 'matrix-like' objects, all handled through a C structure, :cpp:struct:`NumericsMatrix`.

A number (:cpp:member:`NumericsMatrix::storageType`) identifies the type of storage, and only one of the pointers NumericsMatrix.matrixX, X = storageType = 0, 1 or 2, is not NULL and holds the values of the matrix.

The following storage types are currently available:



* "classical" (i.e. dense) column-major storage in a double*, :cpp:member:`NumericsMatrix::matrix0`

* sparse block storage, in a structure of type :cpp:struct:`SparseBlockStructuredMatrix` (warning: only for square matrices!), :cpp:member:`NumericsMatrix::matrix1`

* sparse storage (csc, csr or triplet), based on CSparse (by T. Davis), :cpp:member:`NumericsMatrix::matrix2`

As an example, consider the following 8x8 matrix A:



.. math::

    \begin{equation*} A=\left[\begin{array}{cccc|cc|cc} 1 & 2 & 0 & 4 & 3 &-1 & 0 & 0\\ 2 & 1 & 0 & 0 & 4 & 1 & 0 & 0\\ 0 & 0 & 1 &-1 & 0 & 0 & 0 & 0\\ 5 & 0 &-1 & 6 & 0 & 6 & 0 & 0\\ \hline 0 & 0 & 0 & 0 & 1 & 0 & 0 & 5\\ 0 & 0 & 0 & 0 & 0 & 2 & 0 & 2\\ \hline 0 & 0 & 2 & 1 & 0 & 0 & 2 & 2\\ 0 & 0 & 2 & 2 & 0 & 0 & -1& 2\\ \end{array}\right] \end{equation*}

How can we store this matrix?

The classical (dense) storage gives:



* M.storageType = 0

* M.size0 = 8, M.size1 = 8

* M.dense = [1 2 0 5 0 0 0 0 2 1 0 0 ...]

  matrix0 being a double* of size 64.

With the second storage type, :cpp:struct:`SparseBlockStructuredMatrix`, we have:

* M.storageType = 1

* M.size0 = 8, M.size1 = 8

* M.sparse_block a :cpp:struct:`SparseBlockStructuredMatrix` in which we save:

  * the number of nonzero blocks, 6 (matrix1->nbblocks), and the number of diagonal blocks, 3 (matrix1->size);

  * the vector matrix1->blocksize, the cumulative sum of the sizes of the diagonal blocks, equal to [4,6,8]:

    blocksize[i] = blocksize[i-1] + ni, ni being the size of the diagonal block at (block) row i.

    Note that the last element of blocksize is the actual size of the matrix;

  * the positions of the nonzero blocks in the vectors matrix1->ColumnIndex and matrix1->RowIndex, equal to [0,1,1,2,0,2] and [0,0,1,1,2,2];

  * the list of nonzero blocks, in matrix1->block, stored in Fortran order (column-major) as

    matrix1->block[0] = [1,2,0,5,2,1,0,0,0,0,1,-1,4,0,-1,6]

    matrix1->block[1] = [3,4,0,0,-1,1,0,6]

    ...

    matrix1->block[5] = [2,-1,2,2]

TODO: write proper documentation for the CSparse storage and complete the example above.

.. _numerics_matrix_storage_1NumericsMatrixTools:
.. rubric:: Functions on NumericsMatrix:

.. _numerics_matrix_storage_1NMAlloc:
.. rubric:: Create, fill and delete NumericsMatrix functions:

* :cpp:func:`NM_create`: create an empty shell (no matrix data allocated)

* :cpp:func:`NM_create_with_data`: create a NumericsMatrix and attach caller-owned data

* :cpp:func:`NM_free()`: free a NumericsMatrix

:cpp:func:`NM_create_with_data` accepts a *data* parameter which, if non-NULL, contains the matrix data. The allocation of the underlying storage (``matrix0``, ``matrix1`` or ``matrix2``) remains the caller's responsibility.

.. _numerics_matrix_storage_1NM_LA:
.. rubric:: Linear Algebra:

The following linear algebra operations are supported:

* BLAS-like functions:

  * matrix-vector product: :cpp:func:`NM_gemv` and :cpp:func:`NM_tgemv` (transpose)

  * matrix-matrix product: :cpp:func:`NM_gemm`

  * partial matrix-vector product: :cpp:func:`NM_row_prod`

* LAPACK-like functions: NM_gesv() solves a linear system Ax = b

.. _numerics_matrix_storage_1NM_IO:
.. rubric:: Input / Output:

* :cpp:func:`NM_display`: display a :cpp:struct:`NumericsMatrix`

* :cpp:func:`NM_display_row_by_row`: display a :cpp:struct:`NumericsMatrix` row by row

* :cpp:func:`NM_write_in_filename`, :cpp:func:`NM_write_in_file`: save to the file system

* :cpp:func:`NM_read_in_filename`, :cpp:func:`NM_read_in_file`: fill a :cpp:struct:`NumericsMatrix` from a file

* :cpp:func:`NM_new_from_file`: create a new :cpp:struct:`NumericsMatrix` from a file
