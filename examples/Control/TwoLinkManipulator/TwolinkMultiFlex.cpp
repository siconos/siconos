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
    unsigned int nDof = 4;           // degrees of freedom for robot arm
    double t0 = 0;                   // initial computation time
    double T = 30;                   // final computation time
    double h = 1e-3;                // time step
    double criterion = 1e-8;
    unsigned int maxIter = 20000;
    double e = 0.0;
    double L = 0.0;
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
    q0(2) = 1.5;
    q0(3) = -0.9;
    SP::SiconosVector z(new SiconosVector(nDof * 6 + 1));
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
    (*z)(13) = -1;
    (*z)(14) = 0;//q0(2);
    (*z)(15) = 0;//q0(3);
    (*z)(16) = 0;//v0(2);
    (*z)(17) = 0;//v0(3);
    (*z)(18) = 0;
    (*z)(19) = 0;
    (*z)(20) = 0;
    (*z)(21) = 0;
    (*z)(22) = 0;
    (*z)(23) = 0;
    (*z)(24) = 0;


    SP::LagrangianDS  arm(new LagrangianDS(createSPtrSiconosVector(q0), createSPtrSiconosVector(v0)));

    // external plug-in
    arm->setComputeMassFunction("TwolinkMultiFlexPlugin", "mass");
    arm->setComputeFGyrFunction("TwolinkMultiFlexPlugin", "FGyr");
    arm->setComputeJacobianFGyrqDotFunction("TwolinkMultiFlexPlugin", "jacobianVFGyr");
    arm->setComputeJacobianFGyrqFunction("TwolinkMultiFlexPlugin", "jacobianFGyrq");
    arm->setComputeFIntFunction("TwolinkMultiFlexPlugin", "U");
    arm->setComputeJacobianFIntqDotFunction("TwolinkMultiFlexPlugin", "jacobFintV");
    arm->setComputeJacobianFIntqFunction("TwolinkMultiFlexPlugin", "jacobFintQ");
    arm->setzPtr(z);

    // -------------------
    // --- Interactions---
    // -------------------

    //  - one with Lagrangian non linear relation to define contact with ground
    //  Both with newton impact nslaw.
    // -- relations --

    SP::NonSmoothLaw nslaw(new NewtonImpactNSL(e));
    SP::Relation relation01(new LagrangianScleronomousR("TwolinkMultiFlexPlugin:h01", "TwolinkMultiFlexPlugin:G01"));
    SP::Interaction inter01(new Interaction(nslaw, relation01));
    SP::Relation relation02(new LagrangianScleronomousR("TwolinkMultiFlexPlugin:h02", "TwolinkMultiFlexPlugin:G02"));
    SP::Interaction inter02(new Interaction(nslaw, relation02));
    
    SP::Relation relation31(new LagrangianScleronomousR("TwolinkMultiFlexPlugin:h31", "TwolinkMultiFlexPlugin:G31"));
    SP::Interaction inter31(new Interaction(nslaw, relation31));
    SP::Relation relation32(new LagrangianScleronomousR("TwolinkMultiFlexPlugin:h32", "TwolinkMultiFlexPlugin:G32"));
    SP::Interaction inter32(new Interaction(nslaw, relation32));
    // -------------
    // --- Model ---
    // -------------

    SP::Model Manipulator(new Model(t0, T));
   // add the dynamical system in the non smooth dynamical system
    Manipulator->nonSmoothDynamicalSystem()->insertDynamicalSystem(arm);

    // link the interaction and the dynamical system
    Manipulator->nonSmoothDynamicalSystem()->link(inter01, arm);
    Manipulator->nonSmoothDynamicalSystem()->link(inter02, arm);
    // link the interaction and the dynamical system
    Manipulator->nonSmoothDynamicalSystem()->link(inter31, arm);
    Manipulator->nonSmoothDynamicalSystem()->link(inter31, arm);

    // ----------------
    // --- Simulation ---
    // ----------------

    // -- Time discretisation --
    SP::TimeDiscretisation t(new TimeDiscretisation(t0, h));

    SP::TimeStepping s(new TimeStepping(t));

    // -- OneStepIntegrators --
    SP::OneStepIntegrator OSI(new MoreauJeanOSI(0.500001));
    s->insertIntegrator(OSI);

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
    unsigned int outputSize = 16;
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
    dataPlot(k, 10) = (*z)(18);
    dataPlot(k, 11) = test;
    dataPlot(k, 12) = (*p)(1);
    dataPlot(k, 13) = (*z)(14);
    dataPlot(k, 14) = (*z)(15);
    dataPlot(k, 15) = (*z)(16);

    bool stop = 0;

    boost::progress_display show_progress(N);

    boost::timer time;
    time.restart();
    while (k < N)
    {
      (*z)(0) = (*q)(0);
      (*z)(1) = (*q)(1);
      (*z)(2) = (*v)(0);
      (*z)(3) = (*v)(1);
      // (*z)(14) = (*q)(2);
      //         (*z)(15) = (*q)(3);
      //  (*z)(16) = (*v)(2);
      //         (*z)(17) = (*v)(3);

      // get current time step
      k++;


      //      relation->computeOutput(s->nextTime());
      //  if(k==1106) stop = 1;

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
      if (test == 5) dataPlot(k, 10) = (*z)(18) / h;
      else dataPlot(k, 10) = (*z)(18);
      dataPlot(k, 11) = test;
      dataPlot(k, 13) = (*z)(14);
      dataPlot(k, 14) = (*z)(15);
      dataPlot(k, 15) = (*z)(16); //(19)*(*z)(19);

      s->newtonSolve(criterion, maxIter);
      dataPlot(k, 12) = (*p)(1);
      (*z)(4) = (inter02->getLambda(1))(0);
      (*z)(18) = (inter31->getLambda(1))(0);
      //  if(k==41000)
      //    {
      //      (*v)(0) = 0.0;(*v)(1) = 0.0;(*v)(2) = 0.0;(*v)(3) = 0.0;
      //    }
      s->nextStep();


      //    controller during impacts accumulation phase before the first impact
      if ((dataPlot(k - 1, 14) <= 0.1) && (test == 0) && (dataPlot(k, 13) < 0.6))
      {
        (*z)(8) = dataPlot(k, 0);
        (*z)(5) = (*z)(14);
        (*z)(10) = dataPlot(k, 3);
        (*z)(7) = (*z)(9);
        arm->setComputeFIntFunction("TwolinkMultiFlexPlugin", "U1");
        test = 1;
      }

      //controller during impacts accumulation phase after the first impact
      if (((*z)(4) > 0) && (test == 1))
      {
        (*z)(8) = dataPlot(k, 0);
        arm->setComputeFIntFunction("TwolinkMultiFlexPlugin", "U2");
        test = 2;
      }
      if (((*z)(4) > 0) && (test == 2))
        nimpact = nimpact + 1;

      // controller during constraint-motion phase.
      if (((*z)(4) > 0) && (test == 2) && (dataPlot(k, 7) - dataPlot(k - 3, 7) == 3)) //  && (fabs((*inter0->y(1))(1))<1e-6))
      {
        // L= dataPlot(k,0)-(*z)(8);
        (*z)(8) = dataPlot(k, 0);
        arm->setComputeFIntFunction("TwolinkMultiFlexPlugin", "U3");
        test = 3;
        nimpact = 0;
      }
      //  controller during impacts accumulation phase after the first impact
      if ((dataPlot(k, 10) > 0) && (test == 3))
      {
        arm->setComputeFIntFunction("TwolinkMultiFlexPlugin", "U4");
        test = 4;
      }

      if (((*z)(18) > 0) && (test == 4))
        nimpact = nimpact + 1;
      // controller during constraint-motion phase.
      if (((*z)(18) > 0) && (test == 4) && (dataPlot(k, 7) - dataPlot(k - 3, 7) == 3)) // && (fabs((*inter0->y(1))(0))<1e-6))
      {
        (*z)(8) = dataPlot(k, 0);
        arm->setComputeFIntFunction("TwolinkMultiFlexPlugin", "U5");
        test = 5;
        nimpact = 0;
      }
      // change of control law with a particular design of the desired trajectory that guarantee the take-off
      if ((trunc((dataPlot(k, 0) + h) / (*z)(11)) > trunc((dataPlot(k, 0)) / (*z)(11))) && (test == 5))
      {
        (*z)(8) = dataPlot(k, 0) + h;
        (*z)(10) = (*z)(12);
        arm->setComputeFIntFunction("TwolinkMultiFlexPlugin", "U6");
        test = 6;
        // L = 0;
      }

      //  controller during free-motion phase
      if (((*z)(13) >= 0) && (test == 6))
      {
        arm->setComputeFIntFunction("TwolinkMultiFlexPlugin", "U");
        test = 0;
        (*z)(13) = 0;
      }

      if (stop) break;
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
