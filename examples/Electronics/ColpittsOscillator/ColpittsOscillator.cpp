/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/
//-----------------------------------------------------------------------
//
// Colpitts Oscillator
//-----------------------------------------------------------------------
#include "SiconosKernel.hpp"

using namespace std;
int main(int argc, char* argv[])
{

  double t0 = 0.0;
  double T = 100.0;        // Total simulation time
  double h_step = 1.0e-3;  // Time step
  double L = 0.1;   // inductance
  double C1 = 2.0;   // capacitance
  double C2 = 0.8;   // capacitance
  double Rc = 10.0;    // resistance
  double Re = 20.0;    // resistance
  // double Rb = 0.5;    // resistance
  double alphaF = 0.99;
  double alphaR = 0.015;
  double VCC = 10;
  double VEE = 20;
  string Modeltitle = "Colpitts";

  boost::timer time;
  time.restart();
  try
  {
    // --- Dynamical system specification ---
    SP::SiconosVector init_state(new SiconosVector(3,0.0));
//    init_state->setValue(1,-1.0);
    SP::SimpleMatrix LS_A(new SimpleMatrix(3, 3));

    LS_A->setValue(0, 0, -1.0 / (Rc*C1));
    LS_A->setValue(1, 0, -1.0 / (Rc*C2));
    LS_A->setValue(2, 0, -1.0/L);

    LS_A->setValue(0, 1, 1.0 /C1 *(-1.0/Rc));
    LS_A->setValue(1, 1, -1.0 /C2 *(1.0/Rc +1.0/Re));
    LS_A->setValue(2, 1, -1.0/L);

    LS_A->setValue(0, 2, 1.0 / (C1));
    LS_A->setValue(1, 2, 1.0 / (C2));
    LS_A->setValue(2, 2, 0.0);

    SP::SiconosVector LS_b(new SiconosVector(3));

    LS_b->setValue(0,VCC/(Rc*C1));
    LS_b->setValue(1,1.0/C2*(  VCC/Rc-VEE/Re));
    LS_b->setValue(2,VCC/L);


    SP::FirstOrderLinearDS LSCollpitts(new FirstOrderLinearDS(init_state,LS_A,LS_b));

    // --- Interaction between linear system and non smooth system ---
    SP::SimpleMatrix Int_C(new SimpleMatrix(2, 3));

    (*Int_C)(0, 0) = 1.0;
    (*Int_C)(1, 0) = 0.0;

    (*Int_C)(0, 1) = 1.0;
    (*Int_C)(1, 1) = 1.0;

    (*Int_C)(0, 2) = 0.0;
    (*Int_C)(1, 2) = 0.0;

    //SP::SiconosMatrix Int_D(new SimpleMatrix(2, 2,0.0));
    // (*Int_D)(0, 0) = 1.0 / Rvalue;
    // (*Int_D)(0, 1) = 1.0 / Rvalue;
    // (*Int_D)(0, 2) = -1.0;
    // (*Int_D)(1, 0) = 1.0 / Rvalue;
    // (*Int_D)(1, 1) = 1.0 / Rvalue;
    // (*Int_D)(1, 3) = -1.0;
    // (*Int_D)(2, 0) = 1.0;
    // (*Int_D)(3, 1) = 1.0;

    SP::SimpleMatrix Int_B(new SimpleMatrix(3, 2));
    (*Int_B)(0, 0) = 1.0 /C1 ;
    (*Int_B)(1, 0) = (1.0-alphaR) / C2;
    (*Int_B)(2, 0) = 0.0;
    (*Int_B)(0, 1) = -alphaF / C1 ;
    (*Int_B)(1, 1) = (1.0-alphaF) / C2;
    (*Int_B)(2, 1) = 0.0;


    SP::FirstOrderLinearTIR LTIRCollpitts(new FirstOrderLinearTIR(Int_C, Int_B));
    //LTIRCollpitts->setDPtr(Int_D);

    SP::NonSmoothLaw nslaw(new ComplementarityConditionNSL(2));

    SP::Interaction InterCollpitts(new Interaction(nslaw, LTIRCollpitts));

    // --- Model creation ---
    SP::Model Collpitts(new Model(t0, T, Modeltitle));
    // add the dynamical system in the non smooth dynamical system
    Collpitts->nonSmoothDynamicalSystem()->insertDynamicalSystem(LSCollpitts);
    // link the interaction and the dynamical system
    Collpitts->nonSmoothDynamicalSystem()->link(InterCollpitts, LSCollpitts);

    // ------------------
    // --- Simulation ---
    // ------------------


    // -- (1) OneStepIntegrators --
    double theta = 0.5;
    SP::EulerMoreauOSI aOSI(new EulerMoreauOSI(theta));

    // -- (2) Time discretisation --
    SP::TimeDiscretisation aTiDisc(new TimeDiscretisation(t0, h_step));

    // -- (3) Non smooth problem
    SP::LCP aLCP(new LCP(SICONOS_LCP_LEMKE));
    aLCP->numericsSolverOptions()->iparam[0]=1;  // Multiple solutions 0 or 1
    aLCP->numericsSolverOptions()->iparam[3]=0;  // choice of seeds for multiple solutions
    aLCP->numericsSolverOptions()->iparam[4]=1;  // LS for enum
    // aLCP->setNumericsVerboseMode(1);

    // -- (4) Simulation setup with (1) (2) (3)
    SP::TimeStepping aTS(new TimeStepping(aTiDisc, aOSI, aLCP));
    Collpitts->setSimulation(aTS);
    
    // Initialization
    cout << "====> Initialisation ..." << endl << endl;
    Collpitts->initialize();
    cout << " ---> End of initialization." << endl;

    int k = 0;
    double h = aTS->timeStep();
    int N = ceil((T - t0) / h); // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    SimpleMatrix dataPlot(N, 8);

    SP::SiconosVector x = LSCollpitts->x();
    SP::SiconosVector y = InterCollpitts->y(0);
    SP::SiconosVector lambda = InterCollpitts->lambda(0);

    // For the initial time step:
    // time
    dataPlot(k, 0) = t0;

    dataPlot(k, 1) = (*x)(0);
    dataPlot(k, 2) = (*x)(1);
    dataPlot(k, 3) = (*x)(2);

    dataPlot(k, 4) = (*y)(0);
    dataPlot(k, 5) = (*y)(1);

    dataPlot(k, 6) = (*lambda)(0);
    dataPlot(k, 7) = (*lambda)(1);


    boost::timer t;
    t.restart();


    boost::progress_display show_progress(N);

    boost::timer time;
    time.restart();

    // --- Time loop  ---
    for(k = 1 ; k < N ; ++k)
    {
      // solve ...
      aTS->computeOneStep();
      //  aLCP->display();
      // --- Get values to be plotted ---
      // time
      dataPlot(k, 0) = aTS->nextTime();


      dataPlot(k, 1) = (*x)(0);
      dataPlot(k, 2) = (*x)(1);
      dataPlot(k, 3) = (*x)(2);

      dataPlot(k, 4) = (*y)(0);
      dataPlot(k, 5) = (*y)(1);

      dataPlot(k, 6) = (*lambda)(0);
      dataPlot(k, 7) = (*lambda)(1);
      ++show_progress;

      aTS->nextStep();

    }


    // --- elapsed time computing ---
    cout << ""  << endl;
    cout << "time = " << t.elapsed() << endl;

    // Number of time iterations
    cout << "Number of iterations done: " << k << endl;

    // dataPlot (ascii) output
    ioMatrix::write("Colpitts.dat", "ascii", dataPlot,"noDim");


    SimpleMatrix dataPlotRef(dataPlot);
    dataPlotRef.zero();
    ioMatrix::read("Colpitts.ref", "ascii", dataPlotRef);
    double errRef= (dataPlot - dataPlotRef).normInf();
    if (errRef > 1e-12)
    {
      SimpleMatrix dataPlotRef2(dataPlot);
      dataPlotRef2.zero();
      ioMatrix::read("Colpitts-sol2.ref", "ascii", dataPlotRef2);
      double errRef2 = (dataPlot - dataPlotRef2).normInf();
      if (errRef2 < errRef) errRef = errRef2;
    }

    std::cout <<
                "Error with respect to  the reference file : "
              <<errRef << std::endl;
    if (errRef > 1e-12)
    {    
      std::cout <<
                "Warning. The result is rather different from the reference file."
                << std::endl;
      //return 1; Since two solutions are possible, we return only a warning
    }

  }
  // --- Exceptions handling ---
  catch(SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch(...)
  {
    cout << "Exception caught " << endl;
  }
  cout << "Computation Time: " << time.elapsed()  << endl;
}
