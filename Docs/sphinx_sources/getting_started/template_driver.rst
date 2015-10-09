.. _template_siconos_driver:


C++ template for siconos driver
===============================

.. highlight:: c++
	       
.. code::

   // Header file
   #include "SiconosKernel.h"
   using namespace std;
   // main program
   int main(int argc, char* argv[])
   {
   // == Start timer ==
   boost::timer time;
   time.restart();
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
    // -- Model --
    // == Creation of the Simulation ==
    // -- TimeDiscretisation --
    // -- Simulation (time stepping or event-driven)
    // -- OneStepIntegrator --
    // -- OneStepNSProblem --
    // ================================= Computation =================================

    // --- Initialisation of the simulation ---
	
    // --- Time loop ---
  }
  
  // == Catch exceptions ==
  catch(SiconosException e)
    {cout << e.report() << endl;}
  catch(...)
    {cout << "Exception caught in mySample.cpp" << endl;}
  
  // == get elapsed time ==
  cout << "Computation Time " << time.elapsed()  << endl;  
  }
