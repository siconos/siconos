.. _template_siconos_driver:


C++ template for siconos driver
===============================

.. highlight:: c++
	       
.. code::

   // Header file
   #include <SiconosKernel.h>
   using namespace std;
   // main program
   int main(int argc, char* argv[])
   {
   // Exception handling
   try
   {
    // == User-defined parameters ==
    
    // ================= Creation of the model =======================
    
    // == Creation of the NonSmoothDynamicalSystem ==
    // -- DynamicalSystem(s) --
    // -- Interaction --
    // - Relations - 
    // - NonSmoothLaw -
    // -- NonSmoothDynamicalSystem --	
    // == Creation of the Simulation ==
    // -- TimeDiscretisation --
    // -- Simulation (time stepping or event-driven)
    // -- OneStepIntegrator --
    // -- OneStepNSProblem --
    // ================================= Computation =================================

    // --- Time loop ---
  }
  
  // == Catch exceptions ==
  catch(...)
    {
      Siconos::exception::process();
    }
  // == get elapsed time ==
  cout << "Computation Time " << time.elapsed()  << endl;  
  }
