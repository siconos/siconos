.. _siconos_algebra:

Linear Algebra in Siconos
=========================

.. _kernel_matrix_storage:

Matrices and vectors handling in kernel component
-------------------------------------------------

In Siconos kernel, all matrices and vectors are just interfaces to `Eigen <https://libeigen.gitlab.io/eigen/docs-nightly/group__TutorialMatrixClass.html>`_ vectors and matrices
and defined in the files :file:`SiconosMatrix.hpp` and :file:`SiconosVector.hpp`.

Please refer to Eigen documentation for details.

Siconos wrappers are:

* :cpp:class:`SiconosVector` : dense vector of real numbers (type : double precision)
* :cpp:class:`BlockVector` : vector of :cpp:class:`SiconosVector`
* :cpp:class:`SiconosDenseMatrix` : dense matrix of double
* :cpp:class:`SiconosSparseMatrix` : sparse matrix of double, column-major storage.

and some fixed-sized objects, like :cpp:class:`SiconosMatrix33`, :cpp:class:`SiconosMatrix66` ....

Notice that BlockVector are no more that a collection of pointers to SiconosVector.
Then in most cases, to build such an object, you just need to insert some existing objects.
The usual ways of construction are described below.


.. _numerics_matrix_storage:

Matrix Storage in numerics component
------------------------------------

Numerics component proposes different ways to store 'matrix-like' objects, all handled through a C structure, :cpp:struct:`NumericsMatrix` .

Numerics component proposes different ways to store 'matrix-like' objects, all handled through a C structure, :cpp:struct:`NumericsMatrix` .

A number (:cpp:member:`NumericsMatrix::storageType` ) identify the type of storage while only one pointer among NumericsMatrix.matrixX, X = storageType = 0, 1 or 2, is not NULL and hold the values of the matrix.

At the time, the following storages are available:



* "classical" (i.e. dense) column-major storage in a double*, :cpp:member:`NumericsMatrix::matrix0`

* sparse block storage, in a structure of type :cpp:struct:`SparseBlockStructuredMatrix` (warning: only for square matrices!!), :cpp:member:`NumericsMatrix::matrix1`

* sparse storage (csc, csr or triplet), based on CSparse (from T.Davis), :cpp:member:`NumericsMatrix::matrix2`

As an example, consider the following matrix A of size 8X8:



.. math::

    \begin{equation*} A=\left[\begin{array}{cccc|cc|cc} 1 & 2 & 0 & 4 & 3 &-1 & 0 & 0\\ 2 & 1 & 0 & 0 & 4 & 1 & 0 & 0\\ 0 & 0 & 1 &-1 & 0 & 0 & 0 & 0\\ 5 & 0 &-1 & 6 & 0 & 6 & 0 & 0\\ \hline 0 & 0 & 0 & 0 & 1 & 0 & 0 & 5\\ 0 & 0 & 0 & 0 & 0 & 2 & 0 & 2\\ \hline 0 & 0 & 2 & 1 & 0 & 0 & 2 & 2\\ 0 & 0 & 2 & 2 & 0 & 0 & -1& 2\\ \end{array}\right] \end{equation*}

How can we store this matrix ?

The first classical storage results in:



* M.storageType = 0

* M.size0 = 8, M.size1 = 8

* M.matrix0 = [1 2 0 5 0 0 0 0 2 1 0 0 ...]
  
  matrix0 being a double* of size 64.

For the second way of storage, :cpp:struct:`SparseBlockStructuredMatrix` we have:

* M.storageType = 1

* M.size0 = 8, M.size1 = 8

* M.matrix1 a :cpp:struct:`SparseBlockStructuredMatrix` in which we save:
  
  * the number of non null blocks, 6 (matrix1->nbblocks) and the number of diagonal blocks, 3 (matrix1->size).
  
  * the vector matrix1->blocksize which collects the sum of diagonal blocks sizes until the present one, is equal to [4,6,8],
    
    blocksize[i] = blocksize[i-1] + ni, ni being the size of the diagonal block at row(block) i.
    
    Note that the last element of blocksize corresponds to the real size of the matrix.
  
  * the list of positions of non null blocks in vectors matrix1->ColumnIndex and matrix1->RowIndex, equal to [0,1,1,2,0,2] and [0,0,1,1,2,2]
  
  * the list of non null blocks, in matrix1->block, stored in Fortran order (column-major) as
    
    matrix1->block[0] = [1,2,0,5,2,1,0,0,0,0,1,-1,4,0,-1,6]
    
    matrix1->block[1] = [3,4,0,0,-1,1,0,6]
    
    ...
    
    matrix1->block[5] = [2,-1,2,2]

Todo write proper doc for CSparse storage and complete the example above.

.. _numerics_matrix_storage_1NumericsMatrixTools:
.. rubric:: Functions on NumericsMatrix:

.. _numerics_matrix_storage_1NMAlloc:
.. rubric:: Create, fill and delete NumericsMatrix functions:

* :cpp:func:`NM_create` : allocation without initial values

* :cpp:func:`NM_create_from_data` : allocation and set default values from external data

* :cpp:func:`NM_fill` : needs a pre-defined :cpp:struct:`NumericsMatrix` , set default values from external data

* :cpp:func:`NM_free()` : free a NumericsMatrix

These last two functions accept a *data* parameter, which if non-NULL contains the matrix data.

.. _numerics_matrix_storage_1NM_LA:
.. rubric:: Linear Algebra:

The following linear algebra operation are supported:

* BLAS-like functions:
  
  * product matrix - vector: :cpp:func:`NM_gemv` and :cpp:func:`NM_tgemv` (transpose)
  
  * product matrix - matrix: :cpp:func:`NM_gemm`
  
  * partial product matrix - vector: :cpp:func:`NM_row_prod`

-LAPACK-like functions -NM_gesv(): solve a linear system Ax = b

.. _numerics_matrix_storage_1NM_IO:
.. rubric:: Input / Output:

* :cpp:func:`NM_display` : display a :cpp:struct:`NumericsMatrix`

* :cpp:func:`NM_display_row_by_row` : display a :cpp:struct:`NumericsMatrix` row by row

* :cpp:func:`NM_write_in_filename` , :cpp:func:`NM_write_in_file` : save to filesystem

* :cpp:func:`NM_read_in_filename` , :cpp:func:`NM_read_in_file` : fill a :cpp:struct:`NumericsMatrix` from a file

* :cpp:func:`NM_new_from_file` : create new :cpp:struct:`NumericsMatrix` from a file

