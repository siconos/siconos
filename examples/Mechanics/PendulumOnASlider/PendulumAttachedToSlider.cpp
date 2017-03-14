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


/*!\file
 *C++ input file, MoreauJeanOSI-Time-Stepping version
 */

#include "SiconosKernel.hpp"
#define WITH_FRICTION
#define DISPLAY_INTER
using namespace std;

int main(int argc, char* argv[])
{
  try
  {
    // ================= Creation of the model =======================

    // parameters according to Table 1
    unsigned int nDof = 3; // degrees of freedom for robot arm
    double t0 = 0;         // initial computation time
    double T = 3.0;       // final computation time
    double h = 1e-3;       // time step : do not decrease, because of strong penetrations

    // contact parameters
    double eN = 0.5;
    double eT = 0.;
    double mu = 0.8;

    // initial conditions
    SP::SiconosVector q0(new SiconosVector(nDof));
    SP::SiconosVector v0(new SiconosVector(nDof));
    q0->zero();
    v0->zero();
    (*q0)(0) = 0.1;
    (*q0)(2) = 0.1;
    (*v0)(0) = 2.;

    // -------------------------
    // --- Dynamical systems ---
    // -------------------------
    cout << "====> Model loading ..." << endl << endl;

    SP::LagrangianDS pendulum(new LagrangianDS(q0, v0, "PendulumPlugin:mass"));
    pendulum->setComputeFGyrFunction("PendulumPlugin", "FGyr");
    pendulum->setComputeJacobianFGyrqFunction("PendulumPlugin", "jacobianFGyrq");
    pendulum->setComputeJacobianFGyrqDotFunction("PendulumPlugin", "jacobianFGyrqDot");
    pendulum->setComputeFIntFunction("PendulumPlugin", "FInt");
    pendulum->setComputeJacobianFIntqFunction("PendulumPlugin", "jacobianFIntq");
    pendulum->setComputeJacobianFIntqDotFunction("PendulumPlugin", "jacobianFIntqDot");

    // -------------------
    // --- Interactions---
    // -------------------
    SP::NonSmoothLaw nslaw1(new NewtonImpactFrictionNSL(eN, eT, mu, 2));
    SP::Relation relation1(new LagrangianScleronomousR("PendulumPlugin:g1", "PendulumPlugin:W1"));
    SP::Interaction inter1(new Interaction(nslaw1, relation1));

    SP::NonSmoothLaw nslaw2(new NewtonImpactFrictionNSL(eN, 0.0, 0.0, 2));
    //SP::NonSmoothLaw nslaw2(new NewtonImpactNSL(eN));
    SP::Relation relation2(new LagrangianScleronomousR("PendulumPlugin:g2", "PendulumPlugin:W2"));
    SP::Interaction inter2(new Interaction(nslaw2, relation2));

    SP::NonSmoothLaw nslaw3(new NewtonImpactFrictionNSL(eN, 0.0, 0.0, 2));
    SP::Relation relation3(new LagrangianScleronomousR("PendulumPlugin:g3", "PendulumPlugin:W3"));
    SP::Interaction inter3(new Interaction(nslaw3, relation3));

    // -------------
    // --- Model ---
    // -------------
    SP::Model pendulumWithSlider(new Model(t0, T));
    pendulumWithSlider->nonSmoothDynamicalSystem()->insertDynamicalSystem(pendulum);
    pendulumWithSlider->nonSmoothDynamicalSystem()->link(inter1, pendulum);
    pendulumWithSlider->nonSmoothDynamicalSystem()->link(inter2, pendulum);
    pendulumWithSlider->nonSmoothDynamicalSystem()->link(inter3, pendulum);

    // ----------------
    // --- Simulation ---
    // ----------------
    SP::MoreauJeanOSI OSI(new MoreauJeanOSI(0.5, 1.));

    SP::TimeDiscretisation t(new TimeDiscretisation(t0, h));
    SP::OneStepNSProblem impact(new FrictionContact(2, SICONOS_FRICTION_2D_ENUM));
    //SP::OneStepNSProblem impact(new LCP(SICONOS_LCP_ENUM));

    impact->numericsSolverOptions()->dparam[0] = 1e-08;
    impact->numericsSolverOptions()->iparam[0] = 100;
    impact->numericsSolverOptions()->iparam[2] = 1; // random
    SP::TimeStepping s(new TimeStepping(t));
    s->insertIntegrator(OSI);
    s->insertNonSmoothProblem(impact, SICONOS_OSNSP_TS_VELOCITY);
    s->setNewtonTolerance(1e-10);
    s->setNewtonMaxIteration(200);

    SP::Topology topo = pendulumWithSlider->nonSmoothDynamicalSystem()->topology();

    // ================================= Computation =================================

    // --- Simulation initialization ---
    cout << "====> Initialisation ..." << endl << endl;
    pendulumWithSlider->setSimulation(s);
    pendulumWithSlider->initialize();
    int N = ceil((T - t0) / h) + 1; // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 14;
    SimpleMatrix dataPlot(N + 1, outputSize);

    SP::SiconosVector q = pendulum->q();
    SP::SiconosVector v = pendulum->velocity();

    dataPlot(0, 0) = pendulumWithSlider->t0();
    dataPlot(0, 1) = (*q)(0);
    dataPlot(0, 2) = (*q)(1);
    dataPlot(0, 3) = (*q)(2);
    dataPlot(0, 4) = (*v)(0);
    dataPlot(0, 5) = (*v)(1);
    dataPlot(0, 6) = (*v)(2);
    dataPlot(0, 7) = (*inter1->y(0))(0) ; // g1
    dataPlot(0, 8) = (*inter2->y(0))(0) ; // g2
    dataPlot(0, 9) = (*inter3->y(0))(0) ; // g3

    // --- Time loop ---
    cout << "====> Start computation ... " << endl << endl;

    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 1;
    boost::progress_display show_progress(N);

    boost::timer time;
    time.restart();
    SP::InteractionsGraph indexSet1 = topo->indexSet(1);

    while ((s->hasNextEvent()) && (k<= 3000))
//    while ((s->hasNextEvent()))
    {
      std::cout <<"t = " <<s->nextTime()-h  <<std::endl;
      //std::cout <<"=====================================================" <<std::endl;
      cout << "q[0] = " << (*q)(0)  << endl;
      cout << "q[1] = " << (*q)(1)  << endl;
      cout << "q[2] = " << (*q)(2)  << endl;
      cout << "v[0] = " << (*v)(0)  << endl;
      cout << "v[1] = " << (*v)(1)  << endl;
      cout << "v[2] = " << (*v)(2)  << endl;
      impact->display();



      //std::cout << "=============== Step k ="<< k<< std::endl;
      s->advanceToEvent();
      impact->setNumericsVerboseMode(0);
      // --- Get values to be plotted ---
      dataPlot(k, 0) = s->nextTime();
      dataPlot(k, 1) = (*q)(0) ;
      dataPlot(k, 2) = (*q)(1);
      dataPlot(k, 3) = (*q)(2);
      dataPlot(k, 4) = (*v)(0);
      dataPlot(k, 5) = (*v)(1);
      dataPlot(k, 6) = (*v)(2);
      dataPlot(k, 7) = (*inter1->y(0))(0) ; // g1
      dataPlot(k, 8) = (*inter2->y(0))(0) ; // g2
      dataPlot(k, 9) = (*inter3->y(0))(0) ; // g3
      dataPlot(k, 10) = s->getNewtonNbIterations();
      dataPlot(k, 11) = indexSet1->size();

      if (indexSet1->size() > 5)
      {
        impact->display();
      }
      //      if (s->nextTime() > 0.035 and (*inter1->lambda(1))(0) >0.0)
#ifdef DISPLAY_INTER
        std::cout << "=============== Step k =" << k << std::endl;
        std::cout << "Time " << s->nextTime() << std::endl;

        impact->display();
        //std::cout << " (*inter1->lambda(1))(0) " << (*inter1->lambda(1))(0) << std:: endl;
        //std::cout << " (*inter2->lambda(1))(0) " << (*inter2->lambda(1))(0) << std:: endl;
#endif

      s->processEvents();
      ++show_progress;
      k++;
    }

    cout << endl << "End of computation - Number of iterations done: " << k - 1 << endl;
    cout << "Computation Time " << time.elapsed()  << endl;

    // --- Output files ---
    cout << "====> Output file writing ..." << endl;
    dataPlot.resize(k, outputSize);
    ioMatrix::write("result.dat", "ascii", dataPlot, "noDim");
  }

  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught in SliderCrankD1MinusLinearOSI.cpp" << endl;
  }
}
