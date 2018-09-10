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

#include "SiconosKernel.hpp"
#include "SiconosControl.hpp"
#include "PID.hpp"
#include "LinearSensor.hpp"

using namespace std;

#define WITH_CONTROL

int main(int argc, char* argv[])
{

  // ================= Creation of the model =======================

  // User-defined main parameters
  unsigned int nDof = 2;           // degrees of freedom for the system
  double t0 = 0;                   // initial computation time
  double T = 100;                  // final computation time
  double h = 0.05;                // time step
  double hControl = 2*h;
  double position_init = 10;      // initial position for lowest bead.
  double velocity_init = 0.0;      // initial velocity for lowest bead.
  double xFinal = 0;              // final value
  // -------------------------
  // --- Dynamical systems ---
  // -------------------------

  cout << "====> Model loading ..." << endl << endl;

  SP::SiconosMatrix A(new SimpleMatrix(nDof, nDof));
  A->zero();
  (*A)(0, 1) = 1;
  (*A)(1, 0) = -1;

  SP::SimpleMatrix B(new SimpleMatrix(nDof, 1));
  (*B)(nDof-1, 0) = 1.0;

  // -- Initial positions and velocities --
  SP::SiconosVector x0(new SiconosVector(nDof));
  (*x0)(0) = position_init;
  (*x0)(1) = velocity_init;

  // -- The dynamical system --
  SP::FirstOrderLinearTIDS doubleIntegrator(new FirstOrderLinearTIDS(x0, A));

  // -------------
  // --- Model ---
  // -------------
  SP::ControlSimulation sim(new ControlZOHSimulation(t0, T, h));

  // add the dynamical system in the non smooth dynamical system
  sim->addDynamicalSystem(doubleIntegrator);

  // use a controlSensor
  SP::SimpleMatrix C(new SimpleMatrix(1, 2, 0));
  (*C)(0, 0) = 1.0;
  SP::LinearSensor sens(new LinearSensor(doubleIntegrator, C));
  sim->addSensor(sens, hControl);

  // add the Observer
  SP::SimpleMatrix L(new SimpleMatrix(2, 1));
  (*L)(0, 0) = -7.5125146;
  (*L)(1, 0) = -50.04168751;
  SP::SiconosVector xHat0(new SiconosVector(2));
  (*xHat0)(0) = position_init;
  (*xHat0)(1) = -5;
  SP::Observer obs(new SlidingReducedOrderObserver(sens, *xHat0, C, L));
  sim->addObserver(obs, hControl);
  // add the PID controller
  SP::SiconosVector K(new SiconosVector(3, 0));
  (*K)(0) = .25;
  (*K)(1) = .125;
  (*K)(2) = 2;
#ifdef WITH_CONTROL
  SP::PID act(new PID(sens));
  act->setRef(xFinal);
  act->setK(K);
  act->setB(B);
  sim->addActuator(act, hControl);
#endif

  cout << "=== End of model loading === " << endl;
  // =========================== End of model definition ===========================

  // ================================= Computation =================================

  // --- Simulation initialization ---

  cout << "====> Initialisation ..." << endl << endl;
  // Initialize the model and the controlManager
  sim->initialize();
 // --- Time loop ---
  cout << "====> Start computation ... " << endl << endl;
  // ==== Simulation loop - Writing without explicit event handling =====
  sim->run();
  // --- Output files ---
  cout << "====> Output file writing ..." << endl;
  SimpleMatrix& dataPlot = *sim->data();
  ioMatrix::write("SlidingReducedOrderObserver.dat", "ascii", dataPlot, "noDim");

  double error=0.0, eps=1e-12;
  if ((error=ioMatrix::compareRefFile(dataPlot, "SlidingReducedOrderObserver.ref",
                                      eps)) >= 0.0
      && error > eps)
    return 1;
  else
    return 0;

}
