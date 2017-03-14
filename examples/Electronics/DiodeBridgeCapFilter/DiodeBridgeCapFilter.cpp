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
//  DiodeBridgeCapFilter  : sample of an electrical circuit involving :
//  - a 1st linear dynamical system LSDiodeBridge1 consisting of
//        an LC oscillator (1 µF , 10 mH)
//  - a non smooth system : a 4 diodes bridge used as a full wave rectifier
//        of the supplied voltage across the LC oscillator, providing power
//    to the resistor load of the 2nd dynamical system
//      - a 2nd linear dynamical system LSDiodeBridge2 consisting of
//        a filtering capacitor in parallel with a load resistor
//
//  Expected behavior :
//  The initial state (Vc = 10 V , IL = 0) of the oscillator provides
//      an initial energy.
//  The oscillator period is 2 Pi sqrt(LC) ~ 0,628 ms.
//      The non smooth system is a full wave rectifier :
//  each phase (positive and negative) of the oscillation allows current
//      to flow in a constant direction through the load.
//      The capacitor filter acts as a tank providing energy to the load resistor
//      when the voltage across the oscillator weakens.
//      The load resistor consumes energy : the oscillation damps.
//
//  State variables LSDiodeBridge1:
//  - the voltage across the capacitor (or inductor)
//  - the current through the inductor
//
//  State variable LSDiodeBridge2:
//  - the voltage across the filtering capacitor
//
//  The interaction between the two dynamical systems is defined by :
//  - complementarity laws between diodes current and voltage. Depending on
//        the diode position in the bridge, y stands for the reverse voltage across
//    the diode or for the diode current (see figure in the template file)
//  - a linear time invariant relation between the state variables and y and
//    lambda (derived from the Kirchhoff laws)
//
//-----------------------------------------------------------------------

#include "SiconosKernel.hpp"

using namespace std;

int main(int argc, char* argv[])
{

  double t0 = 0.0;
  double T = 5e-3;           // Total simulation time
  double h_step = 1e-6;    // Time step
  double Lvalue = 1e-2;      // inductance
  double Cvalue = 1e-6;      // capacitance LC oscillator
  double Rvalue = 1e3;       // load resistance
  double Cfilt  = 300.0e-9;  // filtering capacitor
  double VinitLS1 = 10.0;    // initial voltage LC oscillator
  double VinitLS2 = 0.0;     // initial voltage Cfilt
  string Modeltitle = "DiodeBridgeCapFilter";

  try
  {

    // --- Linear system 1 (LC oscillator) specification ---
    SP::SiconosVector init_stateLS1(new SiconosVector(2));
    (*init_stateLS1)(0) = VinitLS1;

    SP::SimpleMatrix LS1_A(new SimpleMatrix(2, 2));
    (*LS1_A)(0 , 1) = -1.0 / Cvalue;
    (*LS1_A)(1 , 0) = 1.0 / Lvalue;

    // cout << " LS1 matrice A = " << endl;
    // LS1_A->display();
    SP::FirstOrderLinearDS LS1DiodeBridgeCapFilter(new FirstOrderLinearDS(init_stateLS1, LS1_A));

    // --- Linear system 2 (load and filter) specification ---
    SP::SiconosVector init_stateLS2(new SiconosVector(1));
    (*init_stateLS2)(0) = VinitLS2;

    SP::SimpleMatrix LS2_A(new SimpleMatrix(1, 1));
    (*LS2_A)(0 , 0) = -1.0 / (Rvalue * Cfilt);

    // cout << " LS2 matrice A = " << endl;
    // LS2_A->display();
    SP::FirstOrderLinearDS LS2DiodeBridgeCapFilter(new FirstOrderLinearDS(init_stateLS2, LS2_A));

    // --- Interaction between linear systems and non smooth system ---

    SP::SimpleMatrix Int_C(new SimpleMatrix(4, 3));
    (*Int_C)(0 , 2) = 1.0;
    (*Int_C)(2 , 0) = -1.0;
    (*Int_C)(2 , 2) = 1.0;
    (*Int_C)(3 , 0) = 1.0;

    SP::SimpleMatrix Int_D(new SimpleMatrix(4, 4));
    (*Int_D)(0 , 1) = -1.0;
    (*Int_D)(1 , 0) = 1.0;
    (*Int_D)(1 , 2) = 1.0;
    (*Int_D)(1 , 3) = -1.0;
    (*Int_D)(2 , 1) = -1.0;
    (*Int_D)(3 , 1) = 1.0;

    SP::SimpleMatrix Int_B(new SimpleMatrix(3, 4));
    (*Int_B)(0 , 2) = -1.0 / Cvalue;
    (*Int_B)(0 , 3) = 1.0 / Cvalue;
    (*Int_B)(2 , 0) = 1.0 / Cfilt;
    (*Int_B)(2 , 2) = 1.0 / Cfilt;

    SP::FirstOrderLinearTIR LTIRDiodeBridgeCapFilter(new FirstOrderLinearTIR(Int_C, Int_B));
    LTIRDiodeBridgeCapFilter->setDPtr(Int_D);
    SP::NonSmoothLaw nslaw(new ComplementarityConditionNSL(4));

    SP::Interaction InterDiodeBridgeCapFilter(new Interaction(nslaw, LTIRDiodeBridgeCapFilter));

    // --- Model creation ---
    SP::Model DiodeBridgeCapFilter(new Model(t0, T, Modeltitle));
    // add the dynamical system in the non smooth dynamical system
    DiodeBridgeCapFilter->nonSmoothDynamicalSystem()->insertDynamicalSystem(LS1DiodeBridgeCapFilter);
    DiodeBridgeCapFilter->nonSmoothDynamicalSystem()->insertDynamicalSystem(LS2DiodeBridgeCapFilter);
    // link the interaction and the dynamical system
    DiodeBridgeCapFilter->nonSmoothDynamicalSystem()->link(InterDiodeBridgeCapFilter, LS1DiodeBridgeCapFilter, LS2DiodeBridgeCapFilter);


    // ------------------
    // --- Simulation ---
    // ------------------

    // -- (1) OneStepIntegrators --
    double theta = 1.0;
    SP::EulerMoreauOSI aOSI(new EulerMoreauOSI(theta));
    // -- (2) Time discretisation --
    SP::TimeDiscretisation aTiDisc(new TimeDiscretisation(t0, h_step));
    // -- (3) Non smooth problem
    SP::LCP aLCP(new LCP());
    // -- (4) Simulation setup with (1) (2) (3)
    SP::TimeStepping aTS(new TimeStepping(aTiDisc, aOSI, aLCP));
    DiodeBridgeCapFilter->setSimulation(aTS);
    
    // Initialization
    cout << "====> Initialisation ..." << endl << endl;
    DiodeBridgeCapFilter->initialize();
    cout << " ---> End of initialization." << endl;


    int k = 0;
    double h = aTS->timeStep();
    int N = ceil((T - t0) / h); // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    SimpleMatrix dataPlot(N, 7);

    // For the initial time step:

    // time
    dataPlot(k, 0) = DiodeBridgeCapFilter->t0();

    // inductor voltage
    dataPlot(k, 1) = (*LS1DiodeBridgeCapFilter->x())(0);

    // inductor current
    dataPlot(k, 2) = (*LS1DiodeBridgeCapFilter->x())(1);

    // diode R1 current
    dataPlot(k, 3) = (InterDiodeBridgeCapFilter->getLambda(0))(0);

    // diode R1 voltage
    dataPlot(k, 4) = -(*InterDiodeBridgeCapFilter->y(0))(0);

    // diode F2 voltage
    dataPlot(k, 5) = -(InterDiodeBridgeCapFilter->getLambda(0))(1);

    // diode F1 current
    dataPlot(k, 6) = (InterDiodeBridgeCapFilter->getLambda(0))(2);

    // --- Compute elapsed time ---
    boost::timer t;
    t.restart();
    // --- Time loop  ---
    while (k < N - 1)
    {
      // get current time step
      k++;

      // solve ...
      aTS->computeOneStep();

      // --- Get values to be plotted ---
      // time
      dataPlot(k, 0) = aTS->nextTime();

      // inductor voltage
      dataPlot(k, 1) = (*LS1DiodeBridgeCapFilter->x())(0);

      // inductor current
      dataPlot(k, 2) = (*LS1DiodeBridgeCapFilter->x())(1);

      // diode R1 current
      dataPlot(k, 3) = (InterDiodeBridgeCapFilter->getLambda(0))(0);

      // diode R1 voltage
      dataPlot(k, 4) = -(*InterDiodeBridgeCapFilter->y(0))(0);

      // diode F2 voltage
      dataPlot(k, 5) = -(InterDiodeBridgeCapFilter->getLambda(0))(1);

      // diode F1 current
      dataPlot(k, 6) = (InterDiodeBridgeCapFilter->getLambda(0))(2);

      aTS->nextStep();

    }


    // --- elapsed time computing ---
    cout << "time = " << t.elapsed() << endl;

    // Number of time iterations
    cout << "Number of iterations done: " << k << endl;

    // dataPlot (ascii) output
    ioMatrix::write("DiodeBridgeCapFilter.dat", "ascii", dataPlot, "noDim");

    // Comparison with a reference file
    std::cout << "Comparison with a reference file" << std::endl;
    SimpleMatrix dataPlotRef(dataPlot);
    dataPlotRef.zero();
    ioMatrix::read("DiodeBridgeCapFilter.ref", "ascii", dataPlotRef);
    SP::SiconosVector err(new SiconosVector(dataPlot.size(1)));
    (dataPlot - dataPlotRef).normInfByColumn(err);
    err->display();
    double error = 0.0;
    for (unsigned int i = 0; i < 3; ++i)
    {
      if (error < (*err)(i))
        error = (*err)(i);
    }

    std::cout << "Error = "<< error << std::endl;
    if (error > 1e-12)
    {
      (dataPlot - dataPlotRef).display();
      std::cout << "Warning. The results is rather different from the reference file." << std::endl;
      std::cout << "Error = "<< error << std::endl;
      return 1;
    }




  }

  // --- Exceptions handling ---
  catch (SiconosException e)
  {
    cout << "SiconosException" << endl;
    cout << e.report() << endl;
  }
  catch (std::exception& e)
  {
    cout << "Exception: " << e.what() << endl;
    exit(-1);
  }
  catch (...)
  {
    cout << "Exception caught " << endl;
  }
}
