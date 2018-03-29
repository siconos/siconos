.. index:: pair: page; Todo List
.. _doxid-todo:

Todo List
=========

.. list-table::
    :widths: 20 80

    *
        - :target:`doxid-todo_1_todo000008` Member :ref:`gfc3d_nsgs <doxid-gfc3d___solvers_8h_1a8744cb3657964afe57ebe1f68e15305a>` ( :ref:`GlobalFrictionContactProblem <doxid-struct_global_friction_contact_problem>` *problem, double *reaction, double *velocity, double *globalVelocity, int *info, SolverOptions *options)

        - Implement ProdTransSBM
          
          Improve the splitting Algorithm with a smaller granularity
          
          Use a global projection perhaps

    *
        - :target:`doxid-todo_1_todo000007` Member :ref:`lcp_newton_FB <doxid-_l_c_p___solvers_8h_1ab44e0b227d2e6f1f96a906f2f023c238>` ( :ref:`LinearComplementarityProblem <doxid-struct_linear_complementarity_problem>` *problem, double *z, double *w, int *info, SolverOptions *options)

        - Optimizing the memory allocation (Try to avoid the copy of JacH into A)
          
          Add rules for the computation of the penalization rho
          
          Add a globalization strategy based on a decrease of a merit function. (Nonmonotone LCP) Reference in Ferris Kanzow 2002

    *
        - :target:`doxid-todo_1_todo000006` Member :ref:`lcp_newton_min <doxid-_l_c_p___solvers_8h_1a759c437f42dfae616b69cf80eb595fe5>` ( :ref:`LinearComplementarityProblem <doxid-struct_linear_complementarity_problem>` *problem, double *z, double *w, int *info, SolverOptions *options)

        - Optimizing the memory allocation (Try to avoid the copy of JacH into A)
          
          Add rules for the computation of the penalization rho
          
          Add a globalization strategy based on a decrease of a merit function. (Nonmonotone LCP) Reference in Ferris Kanzow 2002

    *
        - :target:`doxid-todo_1_todo000005` Member :ref:`lcp_psor <doxid-_l_c_p___solvers_8h_1aa26a3cd3c19ad5869fc5f93c31f3c6ec>` ( :ref:`LinearComplementarityProblem <doxid-struct_linear_complementarity_problem>` *problem, double *z, double *w, int *info, SolverOptions *options)

        - use the relax parameter
          
          add test
          
          add a vector of relaxation parameter wtith an auto sizing (see SOR algorithm for linear solver.)

    *
        - :target:`doxid-todo_1_todo000004` Member :ref:`lcp_rpgs <doxid-_l_c_p___solvers_8h_1acd882d2802211d19d8dc9b324cda02c7>` ( :ref:`LinearComplementarityProblem <doxid-struct_linear_complementarity_problem>` *problem, double *z, double *w, int *info, SolverOptions *options)

        - Sizing the regularization paramter and apply it only on null diagnal term

    *
        - :target:`doxid-todo_1_todo000003` Page :ref:`Matrix Storage in numerics component <doxid-_numerics_matrix_page>`

        - write proper doc for CSparse storage and complete the example above.

    *
        - :target:`doxid-todo_1_todo000002` Member :ref:`mlcp_driver <doxid-_non_smooth_drivers_8h_1a514d2fda3f57ebad039803d5c7aba26f>` ( :ref:`MixedLinearComplementarityProblem <doxid-struct_mixed_linear_complementarity_problem>` *problem, double *z, double *w, SolverOptions *options)

        - Sizing the regularization parameter and apply it only on null diagnal term

    *
        - :target:`doxid-todo_1_todo000001` File ``NonSmoothDrivers.h``

        - solve_qp does not exist

