
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
//  CircuitRLCD  : sample of an electrical circuit involving :
//  - a linear dynamical system consisting of an LC oscillator (1 µF , 10 mH)
//  - a non smooth system (a 1000 Ohm resistor in series with a diode) in parallel
//    with the oscillator
//
//  Expected behavior :
//  The initial state of the oscillator provides an initial energy.
//  The period is 2 Pi sqrt(LC) ~ 0,628 ms.
//  A positive voltage across the capacitor allows current to flow
//  through the resistor-diode branch , resulting in an energy loss :
//  the oscillation damps.
//
//  State variables :
//  - the voltage across the capacitor (or inductor)
//  - the current through the inductor
//
//  Since there is only one dynamical system, the interaction is defined by :
//  - a complementarity law between diode current and voltage where y stands
//    for the reverse voltage across the diode and lambda stands for the
//    the diode current
//  - a linear time invariant relation between the state variables and
//    y and lambda (derived from Kirchhoff laws)
//
//-----------------------------------------------------------------------

#include "SiconosKernel.hpp"

using namespace std;

int main(int argc, char* argv[])
{
  double t0 = 0.0;
  double T = 5e-3;        // Total simulation time
  double h_step = 10.0e-6;// Time step
  double Lvalue = 1e-2;   // inductance
  double Cvalue = 1e-6;   // capacitance
  double Rvalue = 1e3;    // resistance
  double Vinit = 10.0;    // initial voltage
  string Modeltitle = "CircuitRLCD";

  try
  {
    // --- Dynamical system specification ---
    SP::SiconosVector init_state(new SiconosVector(2));
    init_state->setValue(0, Vinit);
    init_state->setValue(1, 0.0);

    SP::SimpleMatrix LS_A(new SimpleMatrix(2, 2));
    LS_A->setValue(0 , 1, -1.0 / Cvalue);
    LS_A->setValue(1 , 0, 1.0 / Lvalue);

    SP::FirstOrderLinearTIDS LSCircuitRLCD(new FirstOrderLinearTIDS(init_state, LS_A));

    // --- Interaction between linear system and non smooth system ---
    SP::SimpleMatrix Int_C(new SimpleMatrix(1, 2));
    Int_C->setValue(0 , 0 , -1.0);

    SP::SimpleMatrix Int_D(new SimpleMatrix(1, 1));
    Int_D->setValue(0 , 0, Rvalue);

    SP::SimpleMatrix Int_B(new SimpleMatrix(2, 1));
    Int_B->setValue(0 , 0, -1.0 / Cvalue);

    SP::FirstOrderLinearTIR LTIRCircuitRLCD(new FirstOrderLinearTIR(Int_C, Int_B));
    SP::NonSmoothLaw NSLaw(new ComplementarityConditionNSL(1));

    LTIRCircuitRLCD->setDPtr(Int_D);

    SP::Interaction InterCircuitRLCD(new Interaction(NSLaw, LTIRCircuitRLCD));

    // --- Model creation ---
    SP::Model CircuitRLCD(new Model(t0, T, Modeltitle));

    // add the dynamical system in the non smooth dynamical system
    CircuitRLCD->nonSmoothDynamicalSystem()->insertDynamicalSystem(LSCircuitRLCD);

    // link the interaction and the dynamical system
    CircuitRLCD->nonSmoothDynamicalSystem()->link(InterCircuitRLCD, LSCircuitRLCD);
    CircuitRLCD->nonSmoothDynamicalSystem()->display();

    // ------------------
    // --- Simulation ---
    // ------------------
    double theta = 0.5000000000001;

    // -- (1) OneStepIntegrators --
    SP::EulerMoreauOSI OSI_RLCD(new EulerMoreauOSI(theta));

    // -- (2) Time discretisation --
    SP::TimeDiscretisation TiDiscRLCD(new TimeDiscretisation(t0, h_step));
    // --- (3) one step non smooth problem
    SP::LCP LCP_RLCD(new LCP());

    // -- (4) Simulation setup with (1) (2) (3)
    SP::TimeStepping StratCircuitRLCD(new TimeStepping(TiDiscRLCD, OSI_RLCD, LCP_RLCD));
    CircuitRLCD->setSimulation(StratCircuitRLCD);
    
    cout << "====> Initialisation ..." << endl << endl;
    CircuitRLCD->initialize();
    cout << " -----> End of initialization." << endl;

    double h = StratCircuitRLCD->timeStep();
    int N = ceil((T - t0) / h); // Number of time steps
    int k = 0;

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    SimpleMatrix dataPlot(N, 6);

    // For the initial time step:

    // time
    dataPlot(k, 0) = CircuitRLCD->t0();

    // inductor voltage
    dataPlot(k, 1) = (*LSCircuitRLCD->x())(0);

    // inductor current
    dataPlot(k, 2) = (*LSCircuitRLCD->x())(1);

    // diode voltage
    dataPlot(k, 3) = - (*InterCircuitRLCD->y(0))(0);

    // diode current
    dataPlot(k, 4) = (InterCircuitRLCD->getLambda(0))(0);

    dataPlot(k, 5) = (LSCircuitRLCD->getR())(0);

    boost::timer t;
    t.restart();

    // --- Time loop  ---
    for (k = 1 ; k < N ; ++k)
    {
      // solve ...
      StratCircuitRLCD->computeOneStep();

      // --- Get values to be plotted ---
      // time
      dataPlot(k, 0) = StratCircuitRLCD->nextTime();

      // inductor voltage
      dataPlot(k, 1) = (*LSCircuitRLCD->x())(0);

      // inductor current
      dataPlot(k, 2) = (*LSCircuitRLCD->x())(1);

      // diode voltage
      dataPlot(k, 3) = - (*InterCircuitRLCD->y(0))(0);

      // diode current
      dataPlot(k, 4) = (InterCircuitRLCD->getLambda(0))(0);

      //dataPlot(k,5) = (LSCircuitRLCD->getR())(0);
      //    dataPlot(k,5) = OSI_RLCD->computeResidu();
      dataPlot(k, 5) = 0;
      // transfer of state i+1 into state i and time incrementation
      StratCircuitRLCD->nextStep();

    }
    // Number of time iterations
    cout << "Number of iterations done: " << k - 1 << endl;
    cout << "Computation Time " << t.elapsed()  << endl;

    // dataPlot (ascii) output
    ioMatrix::write("CircuitRLCD.dat", "ascii", dataPlot, "noDim");

    SimpleMatrix dataPlotRef(dataPlot);
    dataPlotRef.zero();
    ioMatrix::read("CircuitRLCD.ref", "ascii", dataPlotRef);
    std::cout << "Error w.r.t reference file   " <<(dataPlot - dataPlotRef).normInf() << std::endl;
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
}
