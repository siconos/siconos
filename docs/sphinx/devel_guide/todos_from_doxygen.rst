.. index:: single: Todo List
.. _doxid-todo:

Todo List
=========

*  :cpp:struct:`GlobalFrictionContactProblem`

  * Implement ProdTransSBM
  * Improve the splitting Algorithm with a smaller granularity
  * Use a global projection perhaps

* :cpp:struct:`LinearComplementarityProblem`

  * Optimize the memory allocation (try to avoid copying JacH into A)
  * Add rules for the computation of the penalization rho
  * Add a globalization strategy based on the decrease of a merit function (nonmonotone LCP). Reference: Ferris, Kanzow 2002

* :cpp:struct:`LinearComplementarityProblem`

  * Use the relaxation parameter
  * Add a test
  * Add a vector of relaxation parameters with automatic sizing (see the SOR algorithm for linear solvers)

* :cpp:struct:`LinearComplementarityProblem`

  * Size the regularization parameter and apply it only to zero diagonal terms

*  :ref:`Matrix Storage in numerics component <numerics_matrix_storage>`

   Write proper documentation for the CSparse storage and complete the example above.

* :cpp:struct:`MixedLinearComplementarityProblem`

  * Size the regularization parameter and apply it only to zero diagonal terms

*  ``NonSmoothDrivers.h``

   * solve_qp does not exist.
