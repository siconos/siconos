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
//  This circuit corresponds to the circuit (b) p 20, equation 1.39 of
//  "Nonsmooth modeling and simulation for switched circuits"
//  Acary, V., Brogliato, B. and Bonnefon, B.
//  LNEE 69, Springer Verlag 2010
//-----------------------------------------------------------------------

#include "SiconosKernel.hpp"

using namespace std;

int main(int argc, char* argv[])
{
  double t0 = 0.0;
  double T = 5.0;        // Total simulation time
  double h_step = 5e-2;// Time step
  double Lvalue = 1.0;   // inductance
  double Cvalue = 1.0 / (4 * M_PI * M_PI); // capacitance
  double Rvalue = 10.0 ;    // resistance
  string Modeltitle = "CircuitRLCD";

  try
  {
    // --- Dynamical system specification ---
    SP::SiconosVector init_state(new SiconosVector(2));
    init_state->setValue(0, 1.0);
    init_state->setValue(1, 1.0);

    SP::SimpleMatrix LS_A(new SimpleMatrix(2, 2));
    LS_A->setValue(0 , 1, -1.0);
    LS_A->setValue(1 , 0, 1.0 / (Lvalue * Cvalue));

    SP::FirstOrderLinearTIDS LSCircuitRLCD(new FirstOrderLinearTIDS(init_state, LS_A));

    // --- Interaction between linear system and non smooth system ---

    SP::SimpleMatrix Int_C(new SimpleMatrix(1, 2));
    Int_C->setValue(0 , 0 , 1.0 / Cvalue);

    SP::SimpleMatrix Int_D(new SimpleMatrix(1, 1));
    Int_D->setValue(0 , 0, Rvalue);

    SP::SimpleMatrix Int_B(new SimpleMatrix(2, 1));
    Int_B->setValue(0 , 0, 1.0);

    SP::FirstOrderLinearTIR LTIRCircuitRLCD(new FirstOrderLinearTIR(Int_C, Int_B));
    SP::NonSmoothLaw NSLaw(new ComplementarityConditionNSL(1));

    LTIRCircuitRLCD->setDPtr(Int_D);

    SP::Interaction InterCircuitRLCD(new Interaction(NSLaw, LTIRCircuitRLCD));

    // --- Model creation ---
    SP::Model CircuitRLCD(new Model(t0, T, Modeltitle));
    CircuitRLCD->nonSmoothDynamicalSystem()->insertDynamicalSystem(LSCircuitRLCD);
    CircuitRLCD->nonSmoothDynamicalSystem()->link(InterCircuitRLCD, LSCircuitRLCD);
    // ------------------
    // --- Simulation ---
    // ------------------
    double theta = 0.500000000000;
    double gamma = 0.500000000000;

    // -- (1) OneStepIntegrators --
    SP::EulerMoreauOSI OSI_RLCD(new EulerMoreauOSI(theta, gamma));
    OSI_RLCD->setUseGammaForRelation(true);
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
    SimpleMatrix dataPlot(N, 9);

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

    // Storage function
    dataPlot(k, 6) = 1.0 / 2.0 * (
                       1.0 / (Cvalue) * (*LSCircuitRLCD->x())(0) * (*LSCircuitRLCD->x())(0)
                       + Lvalue * (*LSCircuitRLCD->x())(1) * (*LSCircuitRLCD->x())(1)
                     ) ;
    // Dissipation function
    //dataPlot(k,7) =Rvalue*(InterCircuitRLCD->getLambda(0))(0) ;
    dataPlot(k, 7) = 0.0;

    // Total energy
    dataPlot(k, 8) = dataPlot(k, 6) + dataPlot(k, 7)  ;


    boost::timer t;
    t.restart();
    boost::progress_display show_progress(N);
    // --- Time loop  ---
    for (k = 1 ; k < N ; ++k)
    {
      // solve ...
      StratCircuitRLCD->computeOneStep();
      //LCP_RLCD->display();
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

      // Storage function
      dataPlot(k, 6) = 1.0 / 2.0 * (1.0 / Cvalue * (*LSCircuitRLCD->x())(0) * (*LSCircuitRLCD->x())(0)
                                    + Lvalue * (*LSCircuitRLCD->x())(1) * (*LSCircuitRLCD->x())(1)) ;

      // Dissipation function
      dataPlot(k, 7) = h * Rvalue * (InterCircuitRLCD->getLambda(0))(0) * (InterCircuitRLCD->getLambda(0))(0)
                       + dataPlot(k - 1, 7) ;

      // Total energy function
      // Total energy
      dataPlot(k, 8) = dataPlot(k, 6) + dataPlot(k, 7)  ;



      // transfer of state i+1 into state i and time incrementation
      StratCircuitRLCD->nextStep();
      ++show_progress;

    }
    // Number of time iterations
    cout << endl << "Number of iterations done: " << k - 1 << endl;
    cout << "Computation Time " << t.elapsed()  << endl;

    // dataPlot (ascii) output
    dataPlot.resize(k, 9);
    ioMatrix::write("CircuitRLCD.dat", "ascii", dataPlot, "noDim");
    std::cout << "Comparison with a reference file" << std::endl;
    SimpleMatrix dataPlotRef(dataPlot);
    dataPlotRef.zero();
    ioMatrix::read("CircuitRLCD.ref", "ascii", dataPlotRef);
    
    double error = (dataPlot - dataPlotRef).normInf();
    std::cout << "error = " << error << std::endl;
    if ( error > 1e-10)
    {
      std::cout << "Warning. The results is rather different from the reference file." << std::endl;
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
