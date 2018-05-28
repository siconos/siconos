.. _howto:


Test doxygen/sphinx links
-------------------------



To link to a siconos class or a struct:

.. code-block:: rst

   Try to link to :class:`DynamicalSystem`

Result : 
      
Try link to :class:`DynamicalSystem`


To a file (programs listing): 

.. code-block:: rst

   Try to link to :ref:`pgm_kernel_src_modelingTools_DynamicalSystem.hpp`
      
Result : 

Try to link to :ref:`pgm_kernel_src_modelingTools_DynamicalSystem.hpp`
Try to link to :ref:`pgm_externals_SuiteSparse_CXSparse_cs.h`

To a file (documentation): 

.. code-block:: rst

   Try to link to :ref:`file_kernel_src_modelingTools_DynamicalSystem.hpp`
      
Result : 

Try to link to :ref:`file_kernel_src_modelingTools_DynamicalSystem.hpp`
Try to link to :ref:`file_externals_SuiteSparse_CXSparse_cs.h`


To link to a method : 

.. code-block:: rst

   Try to link to :func:`Simulation::nextStep`

Result :

Try to link to :func:`Simulation::nextStep`


To link to a function : 

.. code-block:: rst

   Try to link to :func:`cs_dl_norm`

Result :

Try to link to :func:`cs_dl_norm`
