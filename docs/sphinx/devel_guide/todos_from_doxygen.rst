.. index:: single: Todo List
.. _doxid-todo:

Todo List
=========
  
*  :ref:`GlobalFrictionContactProblem <doxid-struct_global_friction_contact_problem`

  * Implement ProdTransSBM
  * Improve the splitting Algorithm with a smaller granularity
  * Use a global projection perhaps

* :ref:`LinearComplementarityProblem <doxid-struct_linear_complementarity_problem>`

  * Optimizing the memory allocation (Try to avoid the copy of JacH into A)
  * Add rules for the computation of the penalization rho
  * Add a globalization strategy based on a decrease of a merit function. (Nonmonotone LCP) Reference in Ferris Kanzow 2002

* :ref:`LinearComplementarityProblem <doxid-struct_linear_complementarity_problem>`

  * Optimizing the memory allocation (Try to avoid the copy of JacH into A)
  * Add rules for the computation of the penalization rho
  * Add a globalization strategy based on a decrease of a merit function. (Nonmonotone LCP) Reference in Ferris Kanzow 2002

* :ref:`LinearComplementarityProblem <doxid-struct_linear_complementarity_problem>`

  * use the relax parameter
  * add test
  * add a vector of relaxation parameter wtith an auto sizing (see SOR algorithm for linear solver.)

* :ref:`LinearComplementarityProblem <doxid-struct_linear_complementarity_problem>`

  * Sizing the regularization paramter and apply it only on null diagnal term

*  :ref:`Matrix Storage in numerics component <doxid-_numerics_matrix_page>`

   write proper doc for CSparse storage and complete the example above.

* :ref:`MixedLinearComplementarityProblem <doxid-struct_mixed_linear_complementarity_problem>`

  * Sizing the regularization parameter and apply it only on null diagnal term

*  ``NonSmoothDrivers.h``

   * solve_qp does not exist

