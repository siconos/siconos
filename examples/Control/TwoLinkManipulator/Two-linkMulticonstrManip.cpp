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


// =============================== Robot arm sample (HuMAnsPa10) ===============================
//
// see modelRobot1.jpg for complete system view.
//
// Keywords: LagrangianDS, LagrangianLinear relation, MoreauJeanOSI TimeStepping, LCP.
//
// =============================================================================================

#define _USE_MATH_DEFINES
#include "SiconosKernel.hpp"
#include <math.h>

#define PI 3.14159265

using namespace std;


int main(int argc, char* argv[])
{
  boost::timer time;
  time.restart();
  try
  {

    // ================= Creation of the model =======================

    // User-defined main parameters
    unsigned int nDof = 2;           // degrees of freedom for robot arm
    double t0 = 0;                   // initial computation time
    double T = 30;                   // final computation time
    double h = 1e-3;                // time step
    double criterion = 1e-8;
    unsigned int maxIter = 200000;
    double e = 0.7;                  // nslaw
    double e2 = 0.0;
    //double L = 0.0;
    int test = 0;
    int nimpact = 0;

    // -> mind to set the initial conditions below.

    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    // --- DS: manipulator arm ---

    // The dof are angles between ground and arm and between differents parts of the arm. (See corresponding .pdf for more details)

    // Initial position (angles in radian)
    SiconosVector q0(nDof), v0(nDof);
    q0.zero();
    v0.zero();
    q0(0) = 1.5;
    q0(1) = -0.9;
    SP::SiconosVector z(new SiconosVector(nDof * 12 + 1));
    (*z)(0) = q0(0);
    (*z)(1) = q0(1);
    (*z)(2) = v0(0);
    (*z)(3) = v0(1);
    (*z)(4) = 0;
    (*z)(5) = 0;
    (*z)(6) = 0;
    (*z)(7) = 0;
    (*z)(8) = 0;
    (*z)(9) = 0;
    (*z)(10) = 0;
    (*z)(11) = PI;
    (*z)(12) = 0;
    (*z)(13) = 0;
    (*z)(14) = 0;
    (*z)(15) = 0;
    (*z)(16) = 0;
    (*z)(17) = 0;
    (*z)(22) = 0;
    (*z)(23) = 0;
    (*z)(24) = 0;

    SP::LagrangianDS  arm(new LagrangianDS(createSPtrSiconosVector(q0), createSPtrSiconosVector(v0)));

    // external plug-in
    arm->setComputeMassFunction("Two-linkMultiPlugin", "mass");
    arm->setComputeFGyrFunction("Two-linkMultiPlugin", "FGyr");
    arm->setComputeJacobianFGyrqDotFunction("Two-linkMultiPlugin", "jacobianVFGyr");
    arm->setComputeJacobianFGyrqFunction("Two-linkMultiPlugin", "jacobianFGyrq");
    arm->setComputeFIntFunction("Two-linkMultiPlugin", "U");
    arm->setComputeJacobianFIntqDotFunction("Two-linkMultiPlugin", "jacobFintV");
    arm->setComputeJacobianFIntqFunction("Two-linkMultiPlugin", "jacobFintQ");
    arm->setzPtr(z);

    // -------------------
    // --- Interactions---
    // -------------------

    //  - one with Lagrangian non linear relation to define contact with ground
    //  Both with newton impact nslaw.

    // -- relations --

    SP::NonSmoothLaw nslaw(new NewtonImpactNSL(e));
    //SP::Relation relation(new LagrangianScleronomousR("Two-linkMultiPlugin:h0", "Two-linkMultiPlugin:G0"));
    //SP::Interaction inter(new Interaction(nslaw, relation));
    SP::Relation relation01(new LagrangianScleronomousR("Two-linkMultiPlugin:h01", "Two-linkMultiPlugin:G01"));
    SP::Interaction inter01(new Interaction(nslaw, relation01));
    SP::Relation relation02(new LagrangianScleronomousR("Two-linkMultiPlugin:h02", "Two-linkMultiPlugin:G02"));
    SP::Interaction inter02(new Interaction(nslaw, relation02));
    
    SP::Relation relation31(new LagrangianScleronomousR("Two-linkMultiPlugin:h31", "Two-linkMultiPlugin:G31"));
    SP::Interaction inter31(new Interaction(nslaw, relation31));
    SP::Relation relation32(new LagrangianScleronomousR("Two-linkMultiPlugin:h32", "Two-linkMultiPlugin:G32"));
    SP::Interaction inter32(new Interaction(nslaw, relation32));

    SP::SimpleMatrix H10(new SimpleMatrix(1, 2));
    SP::SiconosVector b10(new SiconosVector(1));
    H10->zero();
    (*H10)(0, 0) = -1;
    (*b10)(0) = PI;

    SP::NonSmoothLaw nslaw2(new NewtonImpactNSL(e2));
    SP::Relation relation10(new LagrangianLinearTIR(H10, b10));
    SP::Interaction inter10(new Interaction(nslaw2, relation10));

    SP::SimpleMatrix H11(new SimpleMatrix(1, 2));
    SP::SiconosVector b11(new SiconosVector(1));
    H11->zero();
    (*H11)(0, 0) = 1;
    (*b11)(0) = 0;

    SP::Relation relation11(new LagrangianLinearTIR(H11, b11));
    SP::Interaction inter11(new Interaction(nslaw2, relation11));

    SP::SimpleMatrix H20(new SimpleMatrix(1, 2));
    SP::SiconosVector b20(new SiconosVector(1));
    H20->zero();
    (*H20)(0, 1) = -1;
    (*b20)(0) = 0.0001;

    SP::Relation relation20(new LagrangianLinearTIR(H20, b20));
    SP::Interaction inter20(new Interaction(nslaw2, relation20));

    SP::SimpleMatrix H21(new SimpleMatrix(1, 2));
    SP::SiconosVector b21(new SiconosVector(1));
    H21->zero();
    (*H21)(0, 1) = 1;
    (*b21)(0) = PI - 0.0001;


    SP::Relation relation21(new LagrangianLinearTIR(H21, b21));
    SP::Interaction inter21(new Interaction(nslaw2, relation21));

    // -------------
    // --- Model ---
    // -------------

    SP::Model Manipulator(new Model(t0, T));
    // add the dynamical system in the non smooth dynamical system
    Manipulator->nonSmoothDynamicalSystem()->insertDynamicalSystem(arm);

    // link the interaction and the dynamical system
    Manipulator->nonSmoothDynamicalSystem()->link(inter01, arm);
    Manipulator->nonSmoothDynamicalSystem()->link(inter02, arm);
    Manipulator->nonSmoothDynamicalSystem()->link(inter31, arm);
    Manipulator->nonSmoothDynamicalSystem()->link(inter32, arm);
    Manipulator->nonSmoothDynamicalSystem()->link(inter10, arm);
    Manipulator->nonSmoothDynamicalSystem()->link(inter20, arm);
    Manipulator->nonSmoothDynamicalSystem()->link(inter11, arm);
    Manipulator->nonSmoothDynamicalSystem()->link(inter21, arm);

    // ----------------
    // --- Simulation ---
    // ----------------

    // -- Time discretisation --
    SP::TimeDiscretisation t(new TimeDiscretisation(t0, h));

    SP::TimeStepping s(new TimeStepping(t));

    // -- OneStepIntegrators --
    SP::OneStepIntegrator OSI(new MoreauJeanOSI(0.500001));
    s->insertIntegrator(OSI);
    ;
    // -- OneStepNsProblem --
    SP::OneStepNSProblem osnspb(new LCP());
    s->insertNonSmoothProblem(osnspb);
    Manipulator->setSimulation(s);
    cout << "=== End of model loading === " << endl;

    // =========================== End of model definition ===========================


    // ================================= Computation
    // --- Simulation initialization ---



    Manipulator->initialize();
    cout << "End of model initialisation" << endl;

    int k = 0;
    unsigned int N = ceil((T - t0) / h); // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 15;
    SimpleMatrix dataPlot(N + 1, outputSize);
    // For the initial time step:
    // time

    SP::SiconosVector q = arm->q();
    SP::SiconosVector v = arm->velocity();
    SP::SiconosVector p = arm->p(1);
    // EventsManager * eventsManager = s->eventsManager();

    dataPlot(k, 0) =  Manipulator->t0();
    dataPlot(k, 1) = (*q)(0);
    dataPlot(k, 2) = (*q)(1);
    dataPlot(k, 3) = (*inter02->y(0))(0);
    dataPlot(k, 4) = (*v)(0);
    dataPlot(k, 5) = (*v)(1);
    dataPlot(k, 6) = (*inter01->y(0))(0) - 2;
    dataPlot(k, 7) = nimpact; //(*inter->y(1))(1);
    dataPlot(k, 8) = (*z)(6);
    dataPlot(k, 9) = (*z)(4); //L
    dataPlot(k, 10) = (*z)(24);
    dataPlot(k, 11) = test;
    dataPlot(k, 12) = (*p)(1);
    dataPlot(k, 13) = (*z)(22);
    dataPlot(k, 14) = (*z)(23);

    boost::progress_display show_progress(N);

    boost::timer time;
    time.restart();
    while (k < N)
    {
      (*z)(0) = (*q)(0);
      (*z)(1) = (*q)(1);
      (*z)(2) = (*v)(0);
      (*z)(3) = (*v)(1);
      (*z)(16) = (*z)(14);
      (*z)(17) = (*z)(15);
      (*z)(20) = (*z)(18);
      (*z)(21) = (*z)(19);

      // get current time step
      k++;

      dataPlot(k, 0) =  s->nextTime();
      dataPlot(k, 1) = (*q)(0);
      dataPlot(k, 2) = (*q)(1);
      dataPlot(k, 3) = (*inter02->y(0))(0);
      dataPlot(k, 4) = (*v)(0);
      dataPlot(k, 5) = (*v)(1);
      dataPlot(k, 6) = (*inter01->y(0))(0) - 2;
      dataPlot(k, 7) = nimpact; //(*inter->y(1))(1);
      dataPlot(k, 8) = (*z)(6);
      if (test == 3) dataPlot(k, 9) = (*z)(4) / h;
      else dataPlot(k, 9) = (*z)(4);
      if (test == 5) dataPlot(k, 10) = (*z)(24) / h;
      else dataPlot(k, 10) = (*z)(24);
      dataPlot(k, 11) = test;
      dataPlot(k, 13) = (*z)(22);
      dataPlot(k, 14) = (*z)(23);

      s->newtonSolve(criterion, maxIter);
      dataPlot(k, 12) = (*p)(1);
      (*z)(4) = (inter02->getLambda(1))(0);
      (*z)(24) = (inter31->getLambda(1))(0);
      s->nextStep();

      //    controller during impacts accumulation phase before the first impact
      if ((dataPlot(k, 3) <= 0.01) && (test == 0) && (dataPlot(k, 6) < 0.3))
      {
        (*z)(8) = dataPlot(k, 0);
        (*z)(5) =  0.7 + 0.5 * cos(PI / 2 + 2 * PI * ((*z)(8)) / (*z)(11));
        (*z)(7) = (*z)(9);
        arm->setComputeFIntFunction("Two-linkMultiPlugin", "U10");
        test = 1;
      }

      //  controller during impacts accumulation phase after the first impact
      if ((dataPlot(k, 12) > 0) && (test == 1))
      {
        (*z)(8) = dataPlot(k, 0);
        arm->setComputeFIntFunction("Two-linkMultiPlugin", "U11");
        test = 2;
      }
      if ((dataPlot(k, 12) > 0) && (test == 2))
        nimpact = nimpact + 1;

      // controller during constraint-motion phase.
      if ((dataPlot(k, 12) > 0) && (test == 2) && (dataPlot(k, 7) - dataPlot(k - 1, 7) == 1) && (fabs((*inter02->y(0))(0)) < 1e-6))
      {
        // L= dataPlot(k,0)-(*z)(8);
        (*z)(8) = dataPlot(k, 0);
        arm->setComputeFIntFunction("Two-linkMultiPlugin", "U20");
        test = 3;
        nimpact = 0;
      }
      //  controller during impacts accumulation phase after the first impact
      if (((*z)(24) > 0) && (test == 3))
      {
        arm->setComputeFIntFunction("Two-linkMultiPlugin", "U21");
        test = 4;
      }
      if (((*z)(24) > 0) && (test == 4))
        nimpact = nimpact + 1;
      // controller during constraint-motion phase.
      if (((*z)(24) > 0) && (test == 4) && (dataPlot(k, 7) - dataPlot(k - 1, 7) == 1)) // && (fabs((*inter0->y(1))(0))<1e-6))
      {
        (*z)(8) = dataPlot(k, 0);
        arm->setComputeFIntFunction("Two-linkMultiPlugin", "U22");
        test = 5;
        nimpact = 0;
      }
      // change of control law with a particular design of the desired trajectory that guarantee the take-off
      if ((trunc((dataPlot(k, 0) + h) / (*z)(11)) > trunc((dataPlot(k, 0)) / (*z)(11))) && (test == 5))
      {
        (*z)(8) = dataPlot(k, 0) + h;
        (*z)(10) = (*z)(12);
        arm->setComputeFIntFunction("Two-linkMultiPlugin", "U3");
        test = 6;
        // L = 0;
      }

      //  controller during free-motion phase
      if (((*z)(13) - 0.1 >= 0) && (test == 6))
      {
        arm->setComputeFIntFunction("Two-linkMultiPlugin", "U");
        test = 0;
        (*z)(13) = 0;
      }
      ++show_progress;
    }
    cout << endl << "End of computation - Number of iterations done: " << k << endl;
    cout << "Computation Time " << time.elapsed()  << endl;
    // --- Output files ---
    ioMatrix::write("result.dat", "ascii", dataPlot, "noDim");

  }

  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught in TwolinkMulticonstrManip" << endl;
  }

}
