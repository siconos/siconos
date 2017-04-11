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


double gravity = 10.0;
double m1 = 1.0;
double l1 = 1.0 ;





int main(int argc, char* argv[])
{
  try
  {

    // ================= Creation of the model =======================

    // User-defined main parameters
    unsigned int nDof = 1;           // degrees of freedom for robot arm
    double t0 = 0;                   // initial computation time
    double T = 50.0;                   // final computation time
    double h = 0.0005;                // time step
    double criterion = 0.00005;
    unsigned int maxIter = 2000;
    double e = 0.9;                  // nslaw

    // -> mind to set the initial conditions below.

    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    // unsigned int i;

    // --- DS: Double Pendulum ---

    // Initial position (angles in radian)
    SP::SiconosVector q0(new SiconosVector(nDof));
    SP::SiconosVector v0(new SiconosVector(nDof));
    (*q0).zero();
    (*v0).zero();
    (*q0)(0) = 1;

    SP::LagrangianDS simplependulum(new LagrangianDS(q0, v0));

    SP::SimpleMatrix Mass(new SimpleMatrix(nDof, nDof));
    (*Mass)(0, 0) = m1 * l1;
    simplependulum->setMassPtr(Mass);


    // external plug-in
    //simplependulum->setComputeMassFunction("SimplePendulumPlugin","mass");


    simplependulum->setComputeFIntFunction("SimplePendulumPlugin", "FInt");
    simplependulum->setComputeJacobianFIntqDotFunction("SimplePendulumPlugin", "jacobianVFInt");
    simplependulum->setComputeJacobianFIntqFunction("SimplePendulumPlugin", "jacobianFIntq");

    // -------------------
    // --- Interactions---
    // -------------------


    // -- relations --


    //     SimpleMatrix H(1,2);
    //     SiconosVector b(1);
    //     H.zero();
    //     H(0,0) =1.0;
    //     H(0,1) =0.0;

    //     b(0) = 0.0;


    //     NonSmoothLaw nslaw(new NewtonImpactNSL(e));
    //     Relation relation(new LagrangianLinearTIR(H,b));
    //     Interaction inter =  new Interaction("floor-mass1", allDS,1,1, nslaw, relation);)


    string G = "SimplePendulumPlugin:G0";
    SP::NonSmoothLaw nslaw(new NewtonImpactNSL(e));
    SP::Relation relation(new LagrangianScleronomousR("SimplePendulumPlugin:h0", G));
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

    // -- Time discretisation --
    SP::TimeDiscretisation t(new TimeDiscretisation(t0, h));

    SP::TimeStepping s(new TimeStepping(t));
    //    s->setUseRelativeConvergenceCriteron(true);
    //  s->setRelativeConvergenceTol(0.5);


    // -- OneStepIntegrators --

    //double theta=0.500001;
    double theta = 0.500001;

    SP::OneStepIntegrator OSI(new MoreauJeanOSI(theta));
    s->insertIntegrator(OSI);

    SP::OneStepNSProblem osnspb(new LCP());
    s->insertNonSmoothProblem(osnspb);

    Pendulum->setSimulation(s);

    cout << "=== End of model loading === " << endl;

    // =========================== End of model definition ===========================  dataPlot(k,7) = (*inter->y(0))(0);


    // ================================= Computation =================================


    // --- Simulation initialization ---
    Pendulum->initialize();
    cout << "End of simulation initialisation" << endl;

    int k = 0;
    int N = ceil((T - t0) / h);
    cout << "Number of time step" << N << endl;
    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 11;
    SimpleMatrix dataPlot(N + 1, outputSize);
    // For the initial time step:
    // time
    dataPlot(k, 0) =  Pendulum->t0();
    dataPlot(k, 1) = (*simplependulum->q())(0);
    dataPlot(k, 2) = (*simplependulum->velocity())(0);
    dataPlot(k, 3) =  l1 * sin((*simplependulum->q())(0));
    dataPlot(k, 4) = -l1 * cos((*simplependulum->q())(0));
    dataPlot(k, 5) =  l1 * cos((*simplependulum->q())(0)) * ((*simplependulum->velocity())(0));
    // --- Compute elapsed time ---
    boost::timer tt;
    //    EventsManager eventsManager = s->eventsManager();
    tt.restart();
    // --- Time loop ---
    cout << "Start computation ... " << endl;
    cout << "Number of time step" << N << "\n";
    while (s->hasNextEvent())
    {
      k++;
      if (!(div(k, 1000).rem))  cout << "Step number " << k << "\n";

      // Solve problem
      s->newtonSolve(criterion, maxIter);
      // Data Output
      dataPlot(k, 0) =  s->nextTime();
      dataPlot(k, 1) = (*simplependulum->q())(0);
      dataPlot(k, 2) = (*simplependulum->velocity())(0);
      dataPlot(k, 3) =  l1 * sin((*simplependulum->q())(0));
      dataPlot(k, 4) = -l1 * cos((*simplependulum->q())(0));
      dataPlot(k, 5) =  l1 * cos((*simplependulum->q())(0)) * ((*simplependulum->velocity())(0));
      s->nextStep();
    }

    cout << "time = " << tt.elapsed() << endl;
    cout << "End of computation - Number of iterations done: " << k << endl;

    // --- Output files ---
    ioMatrix::write("SimplePendulumResult.dat", "ascii", dataPlot, "noDim");

  }

  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught in \'sample/MultiBeadsColumn\'" << endl;
  }
}
