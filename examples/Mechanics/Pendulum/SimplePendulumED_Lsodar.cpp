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


// =============================== Double Pendulum Example ===============================
//
// Author: Vincent Acary
//
// Keywords: LagrangianDS, LagrangianLinear relation, MoreauJeanOSI TimeStepping, LCP.
//
// =============================================================================================

#include "SiconosKernel.hpp"
#include <stdlib.h>
using namespace std;
#define PI 3.141592653589793
#define gravity  9.8100

// User-defined main parameters
unsigned int nDof = 2;           // degrees of freedom for robot arm
double L = 1.0;   // Length of the pendulum
double InitAngle = PI / 3.0; // Initial inclination angle
double m = 1.0;             // Mass of the pendulum
double t0 = 0;                   // initial computation time
double T = 10.0;                   // final computation time
double h = 0.01;                // time step
unsigned int N =  ceil(T / h) + 1;   // Number of points to be saved
double e = 0.9;                  // nslaw
double _rho = 0.99;

int main(int argc, char* argv[])
{
  //---------------------------- calculate the computation time --------------------------------------------------
  boost::timer time;
  time.restart();
  try
  {
    // ================= Creation of the model =======================



    // -> mind to set the initial conditions below.

    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    // --- DS: Simple Pendulum ---

    // Initial position (angles in radian)
    SP::SiconosVector q0(new SiconosVector(nDof));
    SP::SiconosVector v0(new SiconosVector(nDof));
    (*q0).zero();
    (*v0).zero();
    (*q0)(0) = L * sin(InitAngle);
    (*q0)(1) = L * cos(InitAngle);

    SP::SimpleMatrix Mass(new SimpleMatrix(nDof, nDof));
    (*Mass)(0, 0) = m;
    (*Mass)(1, 1) = m;
    SP::LagrangianDS simplependulum(new LagrangianLinearTIDS(q0, v0, Mass));


    SP::SiconosVector ForceExtern(new SiconosVector(nDof));
    (*ForceExtern)(0) = 0.0;
    (*ForceExtern)(1) = m * gravity;
    simplependulum->setFExtPtr(ForceExtern);

    // -------------------
    // --- Interactions---
    // -------------------


    SP::NonSmoothLaw nslaw(new NewtonImpactNSL(e));
    SP::Relation relation(new LagrangianScleronomousR("SimplePendulumBilateralConstraintPlugin:h0", "SimplePendulumBilateralConstraintPlugin:G0", "SimplePendulumBilateralConstraintPlugin:G0dot"));
    SP::Interaction inter(new Interaction(nslaw, relation));

    // -------------
    // --- Model ---
    // -------------

    SP::Model Pendulum(new Model(t0, T));
    Pendulum->nonSmoothDynamicalSystem()->insertDynamicalSystem(simplependulum);
    Pendulum->nonSmoothDynamicalSystem()->link(inter, simplependulum);

    // ----------------
    // --- Simulation ---
    // ----------------

    //1. Time discretization
    SP::TimeDiscretisation TimeDiscret(new TimeDiscretisation(t0, h));
    //2. Integration solver for one step
    SP::OneStepIntegrator OSI(new LsodarOSI());
    //3. Nonsmooth problem
    SP::OneStepNSProblem impact(new LCP());
    SP::OneStepNSProblem acceleration(new LCP());
    //4. Simulation with (1), (2), (3)
    SP::Simulation EDscheme(new EventDriven(TimeDiscret));
    EDscheme->insertIntegrator(OSI);
    EDscheme->insertNonSmoothProblem(impact, SICONOS_OSNSP_ED_IMPACT);
    EDscheme->insertNonSmoothProblem(acceleration, SICONOS_OSNSP_ED_SMOOTH_ACC);
    Pendulum->setSimulation(EDscheme); // initialize the model

    // =========================== End of model definition ===========================
    // --- Simulation Initialization ---
    cout << "====> Simulation initialisation ..." << endl << endl;
    EDscheme->setPrintStat(true);
    Pendulum->initialize(); // initialize the model
    cout << "End of simulation initialisation" << endl;

    // SP::LsodarOSI lsodar = std11::static_pointer_cast<LsodarOSI>(OSI);
    // lsodar->setMinMaxStepSizes(9.5e-4,1.0e-3);
    // lsodar->setTol(1,1.0e-3,1.0e-6);
    // lsodar->setMaxOrder(2, 2);


    // ================================= Computation =================================

    SP::EventsManager eventsManager = EDscheme->eventsManager(); // ponters point to the "eventsManager" object
    SP::SiconosVector _q = simplependulum->q();              // pointer points to the position vector of the rocking block
    SP::SiconosVector _qdot = simplependulum->velocity();       // pointer points to the velocity of the rocking block
    SP::SiconosVector _qddot = simplependulum->acceleration();       // pointer points to the velocity of the rocking block
    SP::SiconosVector _g = inter->y(0);
    SP::SiconosVector _gdot = inter->y(1);
    SP::SiconosVector _lambda = inter->lambda(2);
    SP::InteractionsGraph indexSet0 = Pendulum->nonSmoothDynamicalSystem()->topology()->indexSet(0);
    SP::InteractionsGraph indexSet1 = Pendulum->nonSmoothDynamicalSystem()->topology()->indexSet(1);
    SP::InteractionsGraph indexSet2 = Pendulum->nonSmoothDynamicalSystem()->topology()->indexSet(2);
    cout << "Size of IndexSet0: " << indexSet0->size() << endl;
    cout << "Size of IndexSet1: " << indexSet1->size() << endl;
    cout << "Size of IndexSet2: " << indexSet2->size() << endl;
    //-------------------- Save the output during simulation ---------------------------------------------------------
    SimpleMatrix DataPlot(N, 10);
    //------------- At the initial time -----------------------------------------------------------------------------
    DataPlot(0, 0) = Pendulum->t0();
    DataPlot(0, 1) = (*_q)(0); // Position X
    DataPlot(0, 2) = (*_q)(1); // Position Y
    DataPlot(0, 3) = (*_qdot)(0); // Velocity Vx
    DataPlot(0, 4) = (*_qdot)(1); // Velocity Vy
    DataPlot(0, 5) = (*_qddot)(0); // Acceleration ax
    DataPlot(0, 6) = (*_qddot)(1); // Acceleration ay
    DataPlot(0, 7) = (*_g)(0);     // Contraint in position
    DataPlot(0, 8) = (*_gdot)(0);  // Constraint in velocity
    DataPlot(0, 9) = (*_lambda)(0); // Reaction force

    //----------------------------------- Simulation starts ----------------------------------------------------------
    cout << "====> Start computation ... " << endl << endl;
    bool NSEvent = false;
    unsigned int NumberNSEvent = 0;
    unsigned int k = 1;
    boost::progress_display show_progress(N);
    while ((EDscheme->hasNextEvent()) && (k < N))
    {
      EDscheme->advanceToEvent(); // lead the simulation run from one event to the next
      //---------- detect the statue of the current event ------------------------------------
      if (eventsManager->nextEvent()->getType() == 2) // the current event is non-smooth
      {
        NSEvent = true;
      };
      EDscheme->processEvents();  // process the current event
      //------------------- get data at the beginning of non-smooth events ---------------------------
      if (NSEvent)
      {
        DataPlot(k, 0) = EDscheme->startingTime(); // instant at non-smooth event
        SP::SiconosVector _qMemory = simplependulum->qMemory()->getSiconosVector(1);
        SP::SiconosVector _qdotMemory = simplependulum->velocityMemory()->getSiconosVector(1);
        DataPlot(k, 1) = (*_qMemory)(0);
        DataPlot(k, 2) = (*_qMemory)(1);
        DataPlot(k, 3) = (*_qdotMemory)(0);
        DataPlot(k, 4) = (*_qdotMemory)(1);
        k++;
        ++NumberNSEvent;
        ++show_progress;
        NSEvent = false;                        // The next event is maybe smooth
      };
      //-------------------- get data at smooth events or at the end of non-smooth events ---------------
      DataPlot(k, 0) = EDscheme->startingTime();
      DataPlot(k, 1) = (*_q)(0); // Position X
      DataPlot(k, 2) = (*_q)(1); // Position Y
      DataPlot(k, 3) = (*_qdot)(0); // Velocity Vx
      DataPlot(k, 4) = (*_qdot)(1); // Velocity Vy
      DataPlot(k, 5) = (*_qddot)(0); // Acceleration ax
      DataPlot(k, 6) = (*_qddot)(1); // Acceleration ay
      DataPlot(k, 7) = (*_g)(0);     // Contraint in position
      DataPlot(k, 8) = (*_gdot)(0);  // Constraint in velocity
      DataPlot(k, 9) = (*_lambda)(0); // Reaction force
      // go to the next time step
      k++;
      ++show_progress;
    }
    //----------------------- At the end of the simulation --------------------------
    cout << " " << endl;
 
    cout << "End of the simulation" << endl;
    cout << "Number of events processed during simulation: " << (k + 1) << endl;
    cout << "Number of non-smooth events: " << NumberNSEvent << endl;
    cout << "====> Output file writing ..." << endl << endl;
    ioMatrix::write("SimplependulumResult.dat", "ascii", DataPlot, "noDim");
    // Comparison with a reference file
    SimpleMatrix dataPlotRef(DataPlot);
    dataPlotRef.zero();
    ioMatrix::read("result_LsodarOSI.ref", "ascii", dataPlotRef);
    std::cout << (DataPlot - dataPlotRef).normInf() << std::endl;
    if ((DataPlot - dataPlotRef).normInf() > 1e-12)
    {
      std::cout << "Warning. The results is rather different from the reference file." << std::endl;
      return 1;
    }
  }

  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught" << endl;
  }
  cout << "Computation Time: " << time.elapsed()  << endl;
}
