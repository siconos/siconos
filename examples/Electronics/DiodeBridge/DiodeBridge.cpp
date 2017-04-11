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
//  DiodeBridge  : sample of an electrical circuit involving :
//  - a linear dynamical system consisting of an LC oscillator (1 µF , 10 mH)
//  - a non smooth system (a 1000 Ohm resistor supplied through a 4
//    diodes bridge) in parallel with the oscillator
//
//  Expected behavior :

//  The initial state (Vc = 10 V , IL = 0) of the oscillator provides
//  an initial energy.
//  The period is 2 Pi sqrt(LC) ~ 0,628 ms.
//      The non smooth system is a full wave rectifier :
//  each phase (positive and negative) of the oscillation allows current to flow
//  through the resistor in a constant direction, resulting in an energy loss :
//  the oscillation damps.
//
//  State variables :
//  - the voltage across the capacitor (or inductor)
//  - the current through the inductor
//
//  Since there is only one dynamical system, the interaction is defined by :
//  - complementarity laws between diodes current and
//    voltage. Depending on the diode position in the bridge, y stands
//    for the reverse voltage across the diode or for the diode
//    current (see figure in the template file)
//  - a linear time invariant relation between the state variables and
//    y and lambda (derived from Kirchhoff laws)
//
//-----------------------------------------------------------------------
#include "SiconosKernel.hpp"

using namespace std;
int main(int argc, char* argv[])
{

  double t0 = 0.0;
  double T = 5.0e-3;        // Total simulation time
  std::string h_step = std::string("1.0e-6");  // Time step
  double Lvalue = 1e-2;   // inductance
  double Cvalue = 1e-6;   // capacitance
  double Rvalue = 1e3;    // resistance
  double Vinit = 10.0;    // initial voltage
  string Modeltitle = "DiodeBridge";

  boost::timer time;
  time.restart();
  try
  {
    // --- Dynamical system specification ---
    SP::SiconosVector init_state(new SiconosVector(2));
    init_state->setValue(0, Vinit);

    SP::SimpleMatrix LS_A(new SimpleMatrix(2, 2));
    LS_A->setValue(0, 1, -1.0 / Cvalue);
    LS_A->setValue(1, 0, 1.0 / Lvalue);

    SP::FirstOrderLinearDS LSDiodeBridge(new FirstOrderLinearDS(init_state,
                                         LS_A));

    // --- Interaction between linear system and non smooth system ---
    SP::SimpleMatrix Int_C(new SimpleMatrix(4, 2));
    (*Int_C)(2, 0) = -1.0;
    (*Int_C)(3, 0) = 1.0;

    SP::SimpleMatrix Int_D(new SimpleMatrix(4, 4));
    (*Int_D)(0, 0) = 1.0 / Rvalue;
    (*Int_D)(0, 1) = 1.0 / Rvalue;
    (*Int_D)(0, 2) = -1.0;
    (*Int_D)(1, 0) = 1.0 / Rvalue;
    (*Int_D)(1, 1) = 1.0 / Rvalue;
    (*Int_D)(1, 3) = -1.0;
    (*Int_D)(2, 0) = 1.0;
    (*Int_D)(3, 1) = 1.0;

    SP::SimpleMatrix Int_B(new SimpleMatrix(2, 4));
    (*Int_B)(0, 2) = -1.0 / Cvalue ;
    (*Int_B)(0, 3) = 1.0 / Cvalue;

    SP::FirstOrderLinearTIR LTIRDiodeBridge(new FirstOrderLinearTIR(Int_C, Int_B));
    LTIRDiodeBridge->setDPtr(Int_D);

    SP::NonSmoothLaw nslaw(new ComplementarityConditionNSL(4));

    SP::Interaction InterDiodeBridge(new Interaction(nslaw, LTIRDiodeBridge));

    // --- Model creation ---
    SP::Model DiodeBridge(new Model(t0, T, Modeltitle));
    // add the dynamical system in the non smooth dynamical system
    DiodeBridge->nonSmoothDynamicalSystem()->insertDynamicalSystem(LSDiodeBridge);
    // link the interaction and the dynamical system
    DiodeBridge->nonSmoothDynamicalSystem()->link(InterDiodeBridge, LSDiodeBridge);

    // ------------------
    // --- Simulation ---
    // ------------------


    // -- (1) OneStepIntegrators --
    double theta = 0.5;
    SP::EulerMoreauOSI aOSI(new EulerMoreauOSI(theta));

    // -- (2) Time discretisation --
    SP::TimeDiscretisation aTiDisc(new TimeDiscretisation(t0, h_step));

    // -- (3) Non smooth problem
    SP::LCP aLCP(new LCP());

    // -- (4) Simulation setup with (1) (2) (3)
    SP::TimeStepping aTS(new TimeStepping(aTiDisc, aOSI, aLCP));
    DiodeBridge->setSimulation(aTS);
    
    // Initialization
    cout << "====> Initialisation ..." << endl << endl;
    DiodeBridge->initialize();
    cout << " ---> End of initialization." << endl;

    int k = 0;
    double h = aTS->timeStep();
    int N = ceil((T - t0) / h); // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    SimpleMatrix dataPlot(N, 8);

    SP::SiconosVector x = LSDiodeBridge->x();
    SP::SiconosVector y = InterDiodeBridge->y(0);
    SP::SiconosVector lambda = InterDiodeBridge->lambda(0);

    // For the initial time step:
    // time
    dataPlot(k, 0) = t0;

    // inductor voltage
    dataPlot(k, 1) = (*x)(0);

    // inductor current
    dataPlot(k, 2) = (*x)(1);

    // diode R1 current
    dataPlot(k, 3) = (*y)(0);

    // diode R1 voltage
    dataPlot(k, 4) = -(*lambda)(0);

    // diode F2 voltage
    dataPlot(k, 5) = -(*lambda)(1);

    // diode F1 current
    dataPlot(k, 6) = (*lambda)(2);

    // resistor current
    dataPlot(k, 7) = (*y)(0) + (*lambda)(2);

    boost::timer t;
    t.restart();

    // --- Time loop  ---
    for (k = 1 ; k < N ; ++k)
    {
      // solve ...
      aTS->computeOneStep();
      //  aLCP->display();
      // --- Get values to be plotted ---
      // time
      dataPlot(k, 0) = aTS->nextTime();

      // inductor voltage
      dataPlot(k, 1) = (*x)(0);

      // inductor current
      dataPlot(k, 2) = (*x)(1);

      // diode R1 current
      dataPlot(k, 3) = (*y)(0);

      // diode R1 voltage
      dataPlot(k, 4) = -(*lambda)(0);

      // diode F2 voltage
      dataPlot(k, 5) = -(*lambda)(1);

      // diode F1 current
      dataPlot(k, 6) = (*lambda)(2);

      // resistor current
      dataPlot(k, 7) = (*y)(0) + (*lambda)(2);

      aTS->nextStep();

    }


    // --- elapsed time computing ---
    cout << "time = " << t.elapsed() << endl;

    // Number of time iterations
    cout << "Number of iterations done: " << k << endl;

    // dataPlot (ascii) output
    ioMatrix::write("DiodeBridge.dat", "ascii", dataPlot, "noDim");


    SimpleMatrix dataPlotRef(dataPlot);
    dataPlotRef.zero();
    ioMatrix::read("result.ref", "ascii", dataPlotRef);
    std::cout << (dataPlot - dataPlotRef).normInf() << std::endl;
    if ((dataPlot - dataPlotRef).normInf() > 1e-12)
    {
      std::cout <<
                "Warning. The result is rather different from the reference file."
                << std::endl;
      return 1;
    }

  }
  // --- Exceptions handling ---
  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught " << endl;
  }
  cout << "Computation Time: " << time.elapsed()  << endl;
}
