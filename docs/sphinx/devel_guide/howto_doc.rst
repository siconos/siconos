.. _howto:


Test doxygen/sphinx links
-------------------------

Below are some examples of internal links to Siconos objects in Sphinx documents (rst).

* Link to a Siconos class or struct:

  .. code-block:: rst

     Try to link to :cpp:class:`siconos::modeling::DynamicalSystem`

  **Result**:

  Try to link to :cpp:class:`siconos::modeling::DynamicalSystem`


.. * Link to a file (programs listing):

..   .. code-block:: rst

..      Try to link to :ref:`pgm_kernel_src_modeling_DynamicalSystem.hpp`

..   **Result** :

..   Try to link to :ref:`pgm_kernel_src_modeling_DynamicalSystem.hpp`

* Link to a file (documentation):

  .. code-block:: rst

     Try to link to :ref:`file_kernel_src_modeling_DynamicalSystem.hpp`

  **Result**:

  Try to link to :ref:`file_kernel_src_modeling_DynamicalSystem.hpp`


* Link to a class method:

  .. code-block:: rst

     Try to link to :cpp:func:`Simulation::nextStep`

  **Result**:

  Try to link to :cpp:func:`Simulation::nextStep`


* Link to a function:

  .. code-block:: rst

     Try to link to :cpp:func:`cs_dl_norm`

  **Result**:

  Try to link to :cpp:func:`cs_dl_norm`
