.. index:: single: Todo List
.. _doxid-todo:

Todo List
=========
  
*  :cpp:struct:`GlobalFrictionContactProblem`

  * Implement ProdTransSBM
  * Improve the splitting Algorithm with a smaller granularity
  * Use a global projection perhaps

* :cpp:struct:`LinearComplementarityProblem`

  * Optimizing the memory allocation (Try to avoid the copy of JacH into A)
  * Add rules for the computation of the penalization rho
  * Add a globalization strategy based on a decrease of a merit function. (Nonmonotone LCP) Reference in Ferris Kanzow 2002

* :cpp:struct:`LinearComplementarityProblem`

  * Optimizing the memory allocation (Try to avoid the copy of JacH into A)
  * Add rules for the computation of the penalization rho
  * Add a globalization strategy based on a decrease of a merit function. (Nonmonotone LCP) Reference in Ferris Kanzow 2002

* :cpp:struct:`LinearComplementarityProblem`

  * use the relax parameter
  * add test
  * add a vector of relaxation parameter wtith an auto sizing (see SOR algorithm for linear solver.)

* :cpp:struct:`LinearComplementarityProblem`

  * Sizing the regularization paramter and apply it only on null diagnal term

*  :ref:`Matrix Storage in numerics component <numerics_matrix_storage>`

   write proper doc for CSparse storage and complete the example above.

* :cpp:struct:`MixedLinearComplementarityProblem`

  * Sizing the regularization parameter and apply it only on null diagnal term

*  ``NonSmoothDrivers.h``

   * solve_qp does not exist

