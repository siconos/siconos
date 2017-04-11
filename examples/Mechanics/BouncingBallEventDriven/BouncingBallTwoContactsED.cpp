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

/*!\file BouncingBallED.cpp
  \brief \ref EMBouncingBall - C++ input file, Event-Driven version - V. Acary, F. Perignon.

  A Ball bouncing on the ground.
  Direct description of the model.
  Simulation with an Event-Driven scheme.
*/

#include "SiconosKernel.hpp"
#include <boost/numeric/ublas/matrix.hpp>

using namespace std;

int main(int argc, char* argv[])
{
  boost::timer time;
  time.restart();
  try
  {
    // ================= Creation of the model =======================

    // User-defined main parameters
    unsigned int nDof = 3;           // degrees of freedom for the ball
    double t0 = 0;                   // initial computation time
    double T = 10.0;                   // final computation time
    double h = 0.01;                // time step
    double position_init = 1.0;      // initial position for lowest bead.
    double velocity_init = 10.0;      // initial velocity for lowest bead.
    double Heightbox = 1.5;
    double R = 0.1; // Ball radius
    double m = 1; // Ball mass
    double g = 10.0; // Gravity

    // -------------------------
    // --- Dynamical systems ---
    // -------------------------
    cout << "====> Model loading ..." << endl << endl;
    SP::SiconosMatrix Mass(new SimpleMatrix(nDof, nDof));
    (*Mass)(0, 0) = m;
    (*Mass)(1, 1) = m;
    (*Mass)(2, 2) = 3. / 5 * m * R * R;

    // -- Initial positions and velocities --
    SP::SiconosVector q0(new SiconosVector(nDof));
    SP::SiconosVector v0(new SiconosVector(nDof));
    (*q0)(0) = position_init;
    (*v0)(0) = velocity_init;
    // -- The dynamical system --
    SP::LagrangianLinearTIDS ball(new LagrangianLinearTIDS(q0, v0, Mass));
    // -- Set external forces (weight) --
    SP::SiconosVector weight(new SiconosVector(nDof));
    (*weight)(0) = -m * g;
    ball->setFExtPtr(weight);
    //
    // --------------------
    // --- Interactions ---
    // --------------------

    // -- nslaw --
    double e = 0.8; // Warning this example does not work with e=0.0

    // Interaction ball-floor-ceiling
    //
    SP::SimpleMatrix H1(new SimpleMatrix(1, nDof));
    (*H1)(0, 0) = 1.0;
    SP::SiconosVector E1(new SiconosVector(1));
    (*E1)(0) = 0.0;//-1.0*R;
    //
    SP::SimpleMatrix H2(new SimpleMatrix(1, nDof));
    (*H2)(0, 0) = -1.0;
    SP::SiconosVector E2(new SiconosVector(1));
    (*E2)(0) = Heightbox ;//- R;
    // impact law
    SP::NonSmoothLaw  nslaw(new NewtonImpactNSL(e));
    // Interaction at contact 1 (ball-floor)
    SP::Relation relation1(new LagrangianLinearTIR(H1, E1));
    SP::Interaction inter1(new Interaction(nslaw, relation1));
    // Interaction at contact 2 (ball-ceiling)
    SP::Relation relation2(new LagrangianLinearTIR(H2, E2));
    SP::Interaction inter2(new Interaction(nslaw, relation2));
    // --------------------------------
    // --- NonSmoothDynamicalSystem ---
    // --------------------------------

    // -------------
    // --- Model ---
    // -------------

    SP::Model bouncingBall(new Model(t0, T));
    // add the dynamical system in the non smooth dynamical system
    bouncingBall->nonSmoothDynamicalSystem()->insertDynamicalSystem(ball);

    // link the interaction and the dynamical system
    bouncingBall->nonSmoothDynamicalSystem()->link(inter1, ball);
    bouncingBall->nonSmoothDynamicalSystem()->link(inter2, ball);

    // ----------------
    // --- Simulation ---
    // ----------------

    // -- (1) OneStepIntegrators --
    SP::OneStepIntegrator OSI(new LsodarOSI());

    // -- (2) Time discretisation --
    SP::TimeDiscretisation t(new TimeDiscretisation(t0, h));

    // -- (3) Non smooth problem --
    SP::OneStepNSProblem impact(new LCP());
    SP::OneStepNSProblem acceleration(new LCP());

    // -- (4) Simulation setup with (1) (2) (3)
    SP::EventDriven s(new EventDriven(t));
    s->insertIntegrator(OSI);
    s->insertNonSmoothProblem(impact, SICONOS_OSNSP_ED_IMPACT);
    s->insertNonSmoothProblem(acceleration, SICONOS_OSNSP_ED_SMOOTH_ACC);
    cout << "SICONOS_OSNSP_ED_IMPACT: " << SICONOS_OSNSP_ED_IMPACT << endl;
    cout << "SICONOS_OSNSP_ED_ACCELERATION :" << SICONOS_OSNSP_ED_SMOOTH_ACC << endl;
    bouncingBall->setSimulation(s);
    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    // --- Simulation initialization ---
    cout << "====> Simulation initialisation ..." << endl << endl;
    s->setPrintStat(true);
    bouncingBall->initialize();
    OSI->display();
    int N = 1850; // Number of saved points: depends on the number of events ...
    int ll = 0;

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 9;
    SimpleMatrix dataPlot(N, outputSize);
    SP::SiconosVector q = ball->q();        // ball position
    SP::SiconosVector v = ball->velocity(); // ball velocity
    SP::SiconosVector gamma = ball->acceleration(); // ball velocity
    SP::SiconosVector f = ball->p(2);       // resultant force deduced from the LCP at acceleration level
    SP::SiconosVector p = ball->p(1);       // resultant force deduced from the LCP at velocity level


    SP::SiconosVector y1 = inter1->y(0);
    SP::SiconosVector y2 = inter2->y(0);
    //   SiconosVector * y = bouncingBall->nonSmoothDynamicalSystem()->interaction(0)->y(0);

    SP::EventsManager eventsManager = s->eventsManager();
    OSI->display();
    // For the initial time step:
    // time

    dataPlot(0, 0) = bouncingBall->t0();
    dataPlot(0, 1) = (*q)(0);
    dataPlot(0, 2) = (*v)(0);
    dataPlot(0, 3) = (*p)(0);
    dataPlot(0, 4) = 0;
    dataPlot(0, 5) = (*y1)(0);
    dataPlot(0, 6) = (*y2)(0);
    dataPlot(0, 7) = (*gamma)(0);
    dataPlot(0, 8) = (*f)(0);

    // --- Time loop ---
    cout << "====> Start computation ... " << endl << endl;
    bool nonSmooth = false;
    unsigned int numberOfEvent = 0 ;
    double k = 1;
    boost::progress_display show_progress(N);
    s->setPrintStat(true);
    //    s->setTolerance(1e-10);
    while (s->hasNextEvent())
    {
      s->advanceToEvent(); // run simulation from one event to the next
      if (eventsManager->nextEvent()->getType() == 2)
        nonSmooth = true;

      s->processEvents();  // process events
      // If the treated event is non smooth, the pre-impact state has been solved in memory vectors during process.
      if (nonSmooth) // if the event is nonsmooth
      {
        dataPlot(k,0) = s->startingTime(); // get the time at nonsmooth event
        dataPlot(k,1) = (*ball->qMemory()->getSiconosVector(1))(0);
        dataPlot(k,2) = (*ball->velocityMemory()->getSiconosVector(1))(0);
        k++;
        nonSmooth = false;
        ++show_progress;
        dataPlot(k,4) = 1;
        ++ll;
        //         cout << "========================================" << endl;
        //         cout << "Nonsmooth event" << endl;
      }
      dataPlot(k, 0) = s->startingTime();
      dataPlot(k, 1) = (*q)(0);
      dataPlot(k, 2) = (*v)(0);
      dataPlot(k, 3) = (*p)(0);
      dataPlot(k, 5) = (*y1)(0);
      dataPlot(k, 6) = (*y2)(0);
      dataPlot(k, 4) = 0;
      dataPlot(k, 7) = (*gamma)(0);
      dataPlot(k, 8) = (*f)(0);




      cout << "========================================" << endl;
      cout << " time: " << s->startingTime() << endl;
      cout << "ball position: " << (*q)(0) << endl;
      cout << "ball velocity: " << (*v)(0) << endl;
      cout << "gap at contact 1: " << (*y1)(0) << endl;
      cout << "gap at contact 2: " << (*y2)(0) << endl;
      //
      k++;
      ++numberOfEvent;
      ++show_progress;
    }

    // --- Output files ---
    cout << endl;
    cout << "===== End of Event Driven simulation. " << numberOfEvent << " events have been processed. ==== " << endl << endl;
    cout << "Number of nonsmooth events = " << ll << endl;
    cout << "====> Output file writing ..." << endl << endl;
    dataPlot.resize(k, outputSize);
    ioMatrix::write("BouncingBallTwoContactsED.dat", "ascii", dataPlot, "noDim");

    // Comparison with a reference file
    SimpleMatrix dataPlotRef(dataPlot);
    dataPlotRef.zero();
    ioMatrix::read("BouncingBallTwoContactsED.ref", "ascii", dataPlotRef);

    std:: cout << " Error ="<< (dataPlot - dataPlotRef).normInf() << std::endl;

    if ((dataPlot - dataPlotRef).normInf() > 1e-12)
    {
      std::cout << "Warning. The results is rather different from the reference file." << std::endl;
      return 1;
    }

  }

  catch (SiconosException e)
  {
    cout << e.report() << endl;
    return 1;
  }
  catch (...)
  {
    cout << "Exception caught." << endl;
    return 1;
  }
  cout << "Computation Time: " << time.elapsed()  << endl;
}
