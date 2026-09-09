.. _howto:


Test doxygen/sphinx links
-------------------------

Below some examples on how to write internal links to siconos objects in sphinx documents (rst).

* Link to a siconos class or a struct:

  .. code-block:: rst
		  
     Try to link to :cpp:class:`siconos::modeling::DynamicalSystem`

  **Result** : 
      
  Try link to :cpp:class:`siconos::modeling::DynamicalSystem`


.. * Link to a file (programs listing): 

..   .. code-block:: rst

..      Try to link to :ref:`pgm_kernel_src_modeling_DynamicalSystem.hpp`
      
..   **Result** : 

..   Try to link to :ref:`pgm_kernel_src_modeling_DynamicalSystem.hpp`

* Link to a file (documentation): 

  .. code-block:: rst

     Try to link to :ref:`file_kernel_src_modeling_DynamicalSystem.hpp`
      
  **Result** : 

  Try to link to :ref:`file_kernel_src_modeling_DynamicalSystem.hpp`


* Link to a class method : 

  .. code-block:: rst

     Try to link to :cpp:func:`Simulation::nextStep`

  **Result** :

  Try to link to :cpp:func:`Simulation::nextStep`


* Link to a function : 

  .. code-block:: rst

     Try to link to :cpp:func:`cs_dl_norm`

  **Result** :

  Try to link to :cpp:func:`cs_dl_norm`
