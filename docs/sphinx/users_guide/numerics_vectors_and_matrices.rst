.. index:: single: Matrix Storage in numerics component
.. _doxid-_numerics_matrix_page:

Matrix Storage in numerics component
====================================

Numerics component proposes different ways to store 'matrix-like' objects, all handled through a C structure, :ref:`NumericsMatrix <doxid-struct_numerics_matrix>` .

Numerics component proposes different ways to store 'matrix-like' objects, all handled through a C structure, :ref:`NumericsMatrix <doxid-struct_numerics_matrix>` .

A number ( :ref:`NumericsMatrix.storageType <doxid-struct_numerics_matrix_1a1967b9a134e02375a0a5ecbe804a2b49>` ) identify the type of storage while only one pointer among NumericsMatrix.matrixX, X = storageType = 0, 1 or 2, is not NULL and hold the values of the matrix.

At the time, the following storages are available:



* "classical" (i.e. dense) column-major storage in a double*, :ref:`NumericsMatrix.matrix0 <doxid-struct_numerics_matrix_1aed33596859b4c613f48bce4fd9fd707c>`

* sparse block storage, in a structure of type :ref:`SparseBlockStructuredMatrix <doxid-struct_sparse_block_structured_matrix>` (warning: only for square matrices!!), :ref:`NumericsMatrix.matrix1 <doxid-struct_numerics_matrix_1a3ebb848fffd1648ed0609abd17eb9441>`

* sparse storage (csc, csr or triplet), based on CSparse (from T.Davis), :ref:`NumericsMatrix.matrix2 <doxid-struct_numerics_matrix_1ad0498bfc6c8c84c46a860bd77ef1abb5>`

As an example, consider the following matrix A of size 8X8:



.. math::

    \begin{equation*} A=\left[\begin{array}{cccc|cc|cc} 1 & 2 & 0 & 4 & 3 &-1 & 0 & 0\\ 2 & 1 & 0 & 0 & 4 & 1 & 0 & 0\\ 0 & 0 & 1 &-1 & 0 & 0 & 0 & 0\\ 5 & 0 &-1 & 6 & 0 & 6 & 0 & 0\\ \hline 0 & 0 & 0 & 0 & 1 & 0 & 0 & 5\\ 0 & 0 & 0 & 0 & 0 & 2 & 0 & 2\\ \hline 0 & 0 & 2 & 1 & 0 & 0 & 2 & 2\\ 0 & 0 & 2 & 2 & 0 & 0 & -1& 2\\ \end{array}\right] \end{equation*}

How can we store this matrix ?

The first classical storage results in:



* M.storageType = 0

* M.size0 = 8, M.size1 = 8

* M.matrix0 = [1 2 0 5 0 0 0 0 2 1 0 0 ...]
  
  matrix0 being a double* of size 64.

For the second way of storage, :ref:`SparseBlockStructuredMatrix <doxid-struct_sparse_block_structured_matrix>` we have:

* M.storageType = 1

* M.size0 = 8, M.size1 = 8

* M.matrix1 a :ref:`SparseBlockStructuredMatrix <doxid-struct_sparse_block_structured_matrix>` in which we save:
  
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

.. _doxid-_numerics_matrix_page_1NumericsMatrixTools:
.. rubric:: Functions on NumericsMatrix:

.. _doxid-_numerics_matrix_page_1NMAlloc:
.. rubric:: Create, fill and delete NumericsMatrix functions:

* :ref:`NM_create() <doxid-_numerics_matrix_8h_1a7bf697d892ef962778d10c62696735a7>` : allocation without initial values

* :ref:`NM_create_from_data() <doxid-_numerics_matrix_8h_1acf0c380f241e12e90046341ffe283f6f>` : allocation and set default values from external data

* :ref:`NM_fill() <doxid-_numerics_matrix_8h_1ab980dffc6a809393994c57c4c2dfc748>` : needs a pre-defined :ref:`NumericsMatrix <doxid-struct_numerics_matrix>` , set default values from external data

* :ref:`NM_free() <doxid-_numerics_matrix_8h_1a21829f090afbac3cb9b1d4dc3e8e3312>` : free a :ref:`NumericsMatrix <doxid-struct_numerics_matrix>`

These last two functions accept a *data* parameter, which if non-NULL contains the matrix data.

.. _doxid-_numerics_matrix_page_1NM_LA:
.. rubric:: Linear Algebra:

The following linear algebra operation are supported:

* BLAS-like functions:
  
  * product matrix - vector: :ref:`NM_gemv() <doxid-_numerics_matrix_8h_1a7dc236c34ee1ed9dd3bc680ff332241a>` and :ref:`NM_tgemv() <doxid-_numerics_matrix_8h_1a61088ce617e69fb5e25f0246758f2ff8>` (transpose)
  
  * product matrix - matrix: :ref:`NM_gemm() <doxid-_numerics_matrix_8h_1a1d00d2d368f5eea0c1dce711033fecf3>`
  
  * partial product matrix - vector: :ref:`NM_row_prod() <doxid-_numerics_matrix_8h_1a4e37dc94ecee8a398f44481683c91b4d>`

-LAPACK-like functions -NM_gesv(): solve a linear system Ax = b

.. _doxid-_numerics_matrix_page_1NM_IO:
.. rubric:: Input / Output:

* :ref:`NM_display() <doxid-_numerics_matrix_8h_1ab5b41fe722c5aedbb2a7e80fad32a3c9>` : display a :ref:`NumericsMatrix <doxid-struct_numerics_matrix>`

* :ref:`NM_display_row_by_row() <doxid-_numerics_matrix_8h_1a08c490d545a730ae230b9e1b7e56e42f>` : display a :ref:`NumericsMatrix <doxid-struct_numerics_matrix>` row by row

* :ref:`NM_write_in_filename() <doxid-_numerics_matrix_8h_1a38d9c8dfa9ba1ac0013e9c99506e6629>` , :ref:`NM_write_in_file() <doxid-_numerics_matrix_8h_1af5dcc5a62ff0f9035d07cb30addedd14>` : save to filesystem

* :ref:`NM_read_in_filename() <doxid-_numerics_matrix_8h_1a72b1385f5ff6e4e158e900984ab8647e>` , :ref:`NM_read_in_file() <doxid-_numerics_matrix_8h_1a79412e3c65f0775299ac4690abbcb63d>` : fill a :ref:`NumericsMatrix <doxid-struct_numerics_matrix>` from a file

* :ref:`NM_new_from_file() <doxid-_numerics_matrix_8h_1ac5c99652db6920701a65489f857f97ee>` : create new :ref:`NumericsMatrix <doxid-struct_numerics_matrix>` from a file

