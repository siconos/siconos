/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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

/* !\file SMCExampleImplicit.cpp
  \brief Two independent systems of dimension one controlled to slide
  on \f$x = 0\f$. An implicit scheme is used
  O. Huber
  */

#include "SiconosKernel.hpp"
#include "SiconosControl.hpp"
using namespace std;
class MyDS : public FirstOrderLinearDS
{
public:
  MyDS(SP::SiconosVector x0, SP::SiconosMatrix A) : FirstOrderLinearDS(x0, A)
  {
    _b.reset(new SiconosVector(x0->size()));
  };
  void computeb(double time)
  {
    printf("computeB\n");
    double t = sin(50 * time);
    _b->setValue(0,t);
    _b->setValue(1,-t);
    printf("b[0] = %g, b[1] = %g\n", _b->getValue(0), _b->getValue(1));
  };
};
// main program
int main(int argc, char* argv[])
{
  // User-defined parameters
  unsigned int ndof = 2;          // Number of degrees of freedom of your system
  double t0 = 0.0;                // Starting time
  double T = 1;                   // Total simulation time
  double h = 1.0e-4;              // Time step for simulation
  double hControl = 1.0e-2;       // Time step for control
  double Xinit = 1.0;

  if (h > hControl)
  {
    RuntimeException::selfThrow("hControl must be bigger than h");
  }

  // ================= Creation of the model =======================
  // Steps:
  // - create a Dynamical System
  // - add a Simulation to the model

  // -------------------------
  // --- Dynamical systems ---
  // -------------------------

  // First System:
  // dx/dt = Ax + u(t) + r
  // x(0) = x0
  // Note: r = Blambda, B defines in relation below.

  // Matrix declaration
  SP::SiconosMatrix A(new SimpleMatrix(ndof, ndof, 0));
  SP::SiconosVector x0(new SiconosVector(ndof));
  (*x0)(0) = Xinit;
  (*x0)(1) = -Xinit;
  SP::SimpleMatrix sensorC(new SimpleMatrix(2, 2));
  sensorC->eye();
  SP::SimpleMatrix sensorD(new SimpleMatrix(2, 2, 0));
  SP::SimpleMatrix Csurface(new SimpleMatrix(1, 2, 0));
  (*Csurface)(0, 1) = 1;
  SP::SimpleMatrix Brel(new SimpleMatrix(2, 1, 0));
  (*Brel)(1, 0) = 2;

  // Dynamical Systems
  SP::FirstOrderLinearDS processDS(new MyDS(x0, A));
  // -------------
  // --- Model process ---
  // -------------
  SP::ControlSimulation sim(new ControlZOHSimulation(t0, T, h));
  //sim->setSaveOnlyMainSimulation(true);
  sim->addDynamicalSystem(processDS);

  // ------------------
  // --- Simulation ---
  // ------------------
  // Control stuff
  // use a controlSensor
  SP::LinearSensor sens(new LinearSensor(processDS, sensorC, sensorD));
  sim->addSensor(sens, hControl);
  // add the sliding mode controller
  SP::LinearSMC act(new LinearSMC(sens));
  act->setCsurface(Csurface);
  act->setB(Brel);
  sim->addActuator(act, hControl);
  // =========================== End of model definition ===========================

  // ================================= Computation =================================

  // --- Simulation initialization ---

  cout << "====> Simulation initialisation ..." << endl << endl;
  // initialise the process and the ControlManager
  sim->initialize();

  // ==== Simulation loop =====
  cout << "====> Start computation ... " << endl << endl;
  sim->run();
  // --- Output files ---
  cout << "====> Output file writing ..." << endl;
  SimpleMatrix& dataPlot = *sim->data();
  ioMatrix::write("SMCExampleImplicit.dat", "ascii", dataPlot, "noDim");

  // // Comparison with a reference file
  // SimpleMatrix dataPlotRef(dataPlot);
  // dataPlotRef.zero();
  // ioMatrix::read("SMCExampleImplicit.ref", "ascii", dataPlotRef);
  // std::cout << (dataPlot - dataPlotRef).normInf() << std::endl;

  // if ((dataPlot - dataPlotRef).normInf() > 1e-12)
  // {
  //   std::cout << "Warning. The results is rather different from the reference file." << std::endl;
  //   return 1;
  // }

}
