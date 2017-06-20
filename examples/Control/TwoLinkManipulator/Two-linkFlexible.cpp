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
#include "SiconosKernel.hpp"
#include <math.h>

#define PI 3.14159265

#ifdef _MSC_VER
double trunc(double d){ return (d>0) ? floor(d) : ceil(d) ; }
#endif

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
    q0(0) = 0.9;
    q0(1) = -1.5;
    q0(2) = 0.9;
    q0(3) = -1.5;
    SP::SiconosVector z(new SiconosVector(nDof * 5));
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
    (*z)(14) = 0;//q0(2);
    (*z)(15) = 0;//q0(3);
    (*z)(16) = 0;//v0(2);
    (*z)(17) = 0;//v0(3);
    (*z)(18) = 0;
    (*z)(19) = 0;


    SP::LagrangianDS  arm(new LagrangianDS(createSPtrSiconosVector(q0), createSPtrSiconosVector(v0)));

    // external plug-in
    arm->setComputeMassFunction("Two-linkFlexiblePlugin", "mass");
    arm->setComputeFGyrFunction("Two-linkFlexiblePlugin", "FGyr");
    arm->setComputeJacobianFGyrqDotFunction("Two-linkFlexiblePlugin", "jacobianVFGyr");
    arm->setComputeJacobianFGyrqFunction("Two-linkFlexiblePlugin", "jacobianFGyrq");
    arm->setComputeFIntFunction("Two-linkFlexiblePlugin", "U");
    arm->setComputeJacobianFIntqDotFunction("Two-linkFlexiblePlugin", "jacobFintV");
    arm->setComputeJacobianFIntqFunction("Two-linkFlexiblePlugin", "jacobFintQ");
    arm->setzPtr(z);

    // -------------------
    // --- Interactions---
    // -------------------

    //  - one with Lagrangian non linear relation to define contact with ground
    //  Both with newton impact nslaw.

    // -- relations --

    SP::NonSmoothLaw nslaw(new NewtonImpactNSL(e));
    SP::Relation relation1(new LagrangianScleronomousR("Two-linkFlexiblePlugin:h01", "Two-linkFlexiblePlugin:G01"));
    SP::Interaction inter1(new Interaction(nslaw, relation1));
    SP::Relation relation2(new LagrangianScleronomousR("Two-linkFlexiblePlugin:h02", "Two-linkFlexiblePlugin:G02"));
    SP::Interaction inter2(new Interaction(nslaw, relation2));
    //  Relation * relation0 = new LagrangianScleronomousR("Two-linkFlexiblePlugin:h3","Two-linkFlexiblePlugin:G3");
    //  Interaction * inter0 = new Interaction("wall-arm", allDS,1,2, nslaw, relation0);
  
    // -------------
    // --- Model ---
    // -------------

    SP::Model Manipulator(new Model(t0, T));

    // add the dynamical system in the non smooth dynamical system
    Manipulator->nonSmoothDynamicalSystem()->insertDynamicalSystem(arm);

    // link the interaction and the dynamical system
    Manipulator->nonSmoothDynamicalSystem()->link(inter1, arm);
    Manipulator->nonSmoothDynamicalSystem()->link(inter2, arm);

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
    osnspb->numericsSolverOptions()->dparam[0] = 1e-8;



    s->insertNonSmoothProblem(osnspb);
    Manipulator->setSimulation(s);
    cout << "=== End of model loading === " << endl;

    // =========================== End of model definition ===========================


    // ================================= Computation
    // --- Simulation initialization ---



    Manipulator->initialize();
    cout << "End of model initialisation" << endl;

    unsigned int k = 0;
    unsigned int N = ceil((T - t0) / h); // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 14;
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
    dataPlot(k, 3) = (*inter2->y(0))(0);
    dataPlot(k, 4) = (*v)(0);
    dataPlot(k, 5) = (*v)(1);
    dataPlot(k, 6) = (*inter1->y(0))(0) - 2;
    dataPlot(k, 7) = nimpact; //(*inter->y(1))(1);
    dataPlot(k, 8) = (*z)(6);
    dataPlot(k, 9) = (*z)(4); //L
    dataPlot(k, 10) = test;
    dataPlot(k, 11) = (*p)(1);
    dataPlot(k, 12) = (*z)(14);
    dataPlot(k, 13) = (*z)(15);

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

      dataPlot(k, 0) =  s->nextTime();
      dataPlot(k, 1) = (*q)(0);
      dataPlot(k, 2) = (*q)(1);
      dataPlot(k, 3) = (*inter2->y(0))(0);
      dataPlot(k, 4) = (*v)(0);
      dataPlot(k, 5) = (*v)(1);
      dataPlot(k, 6) = (*inter1->y(0))(0) - 2;
      dataPlot(k, 7) = nimpact; //(*inter->y(1))(1);
      dataPlot(k, 8) = (*z)(6);
      if (test == 3) dataPlot(k, 9) = (*z)(4) / h;
      else dataPlot(k, 9) = (*z)(4);
      dataPlot(k, 10) = test;
      dataPlot(k, 12) = (*z)(14);
      dataPlot(k, 13) = (*z)(15);


      s->newtonSolve(criterion, maxIter);
      dataPlot(k, 11) = (*p)(1);
      (*z)(4) = (inter2->getLambda(1))(0);
      s->nextStep();

      //    controller during impacts accumulation phase before the first impact
      if ((dataPlot(k, 13) <= 0.05) && (test == 0) && (dataPlot(k, 6) < 0.3))
      {
        (*z)(8) = dataPlot(k, 0);
        (*z)(5) = (*z)(14);
        // (*z)(10)= dataPlot(k,3);
        (*z)(7) = (*z)(9);
        arm->setComputeFIntFunction("Two-linkFlexiblePlugin", "U1");
        test = 1;
      }

      //controller during impacts accumulation phase after the first impact
      if ((dataPlot(k - 1, 11) > 0) && (test == 1))
      {
        //(*z)(8) = dataPlot(k-1,0);
        arm->setComputeFIntFunction("Two-linkFlexiblePlugin", "U2");
        test = 2;
      }
      if ((dataPlot(k, 11) > 0) && (test == 2))
        nimpact = nimpact + 1;

      // controller during constraint-motion phase.
      if ((dataPlot(k, 11) > 0) && (test == 2) && (dataPlot(k, 7) - dataPlot(k - 1, 7) == 1)) //  && (fabs((*inter0->y(1))(0))<1e-6))
      {
        // L= dataPlot(k,0)-(*z)(8);
        (*z)(8) = dataPlot(k, 0);
        arm->setComputeFIntFunction("Two-linkFlexiblePlugin", "U3");
        test = 3;
        nimpact = 0;
      }

      // change of control law with a particular design of the desired trajectory that guarantee the take-off
      if ((trunc((dataPlot(k, 0) + h) / (*z)(11)) > trunc((dataPlot(k, 0)) / (*z)(11))) && (test == 3))
      {
        (*z)(8) = dataPlot(k, 0) + h;
        (*z)(10) = (*z)(12);
        arm->setComputeFIntFunction("Two-linkFlexiblePlugin", "U4");
        test = 4;
        // L = 0;
      }

      //  controller during free-motion phase
      if (((*z)(13) - 0.01 >= 0) && (test == 4))
      {
        arm->setComputeFIntFunction("Two-linkFlexiblePlugin", "U");
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
    cout << "Exception caught in TwolinkFlexManipulator" << endl;
  }
}





















// SiconosVector * z = new SiconosVector(nDof*4);
//     (*z)(0) = q0(0);
//     (*z)(1) = q0(1);
//     (*z)(2) = v0(0);
//     (*z)(3) = v0(1);
//     (*z)(4) = 0;
//     (*z)(5) = 0;
//     (*z)(6) = 0;
//     (*z)(7) = 0;
//     (*z)(8) = 0;
//     (*z)(9) = 0;
//     (*z)(10) = 0;
//     (*z)(11) = 0;
//     (*z)(12) = 0;
//     (*z)(13) = 0;
//     // (*z)(14) = q0(2);
// //     (*z)(15) = q0(3);


//     LagrangianDS * arm = new LagrangianDS(1, q0, v0);

//     // external plug-in
//     arm->setComputeMassFunction("Two-linkFlexiblePlugin","mass");
//     arm->setComputeFGyrFunction("Two-linkFlexiblePlugin","FGyr");
//     arm->setComputeJacobianFGyrFunction(1,"Two-linkFlexiblePlugin","jacobianVFGyr");
//     arm->setComputeJacobianFGyrFunction(0,"Two-linkFlexiblePlugin","jacobianFGyrq");
//     arm->setComputeFExtFunction("Two-linkFlexiblePlugin","U");
//     arm->setzPtr(z);

//     allDS.insert(arm);

//     // -------------------
//     // --- Interactions---
//     // -------------------

//     //  - one with Lagrangian non linear relation to define contact with ground
//     //  Both with newton impact nslaw.

//     InteractionsSet allInteractions;

//     // -- relations --

//     NonSmoothLaw * nslaw = new NewtonImpactNSL(e);
//     Relation * relation = new LagrangianScleronomousR("Two-linkFlexiblePlugin:h0","Two-linkFlexiblePlugin:G0");
//     Interaction * inter = new Interaction("floor-arm", allDS,0,2, nslaw, relation);


//    //  SimpleMatrix H1(2,4);
// //     SiconosVector b1(2);
// //     H1.zero();
// //     H1(0,0) =-1;
// //     H1(1,0) =1;

// //     b1(0) = PI;
// //     b1(1) = 0;


// //     NonSmoothLaw * nslaw2 = new NewtonImpactNSL(e2);
// //     Relation * relation1 = new LagrangianLinearTIR(H1,b1);
// //     Interaction * inter1 =  new Interaction("floor-arm2", allDS,1,2, nslaw2, relation1);

// //     SimpleMatrix H2(2,4);
// //     SiconosVector b2(2);
// //     H2.zero();
// //     H2(0,1) =-1;
// //     H2(1,1) =1;

// //     b2(0) = 0.0001;
// //     b2(1) = PI-0.0001;

// //     Relation * relation2 = new LagrangianLinearTIR(H2,b2);
// //     Interaction * inter2 =  new Interaction("singular-points", allDS,2,2, nslaw2, relation2);



//     allInteractions.insert(inter);
//    //  allInteractions.insert(inter1);
// //     allInteractions.insert(inter2);


// // Interaction * inter = new Interaction("floor-arm", allDS,0,4, nslaw, relation);




// //     //  NonSmoothLaw * nslaw = new NewtonImpactNSL(e);
// // //      Relation * relation = new LagrangianScleronomousR("Two-linkPlugin:h0","Two-linkPlugin:G0");
// // //      Interaction * inter = new Interaction("floor-arm", allDS,0,1, nslaw, relation);

// //     allInteractions.insert(inter);
// //     allInteractions.insert(inter1);
// //     allInteractions.insert(inter2);
//     // --------------------------------
//     // --- NonSmoothDynamicalSystem ---
//     // -------------------------------

//     NonSmoothDynamicalSystem * nsds = new NonSmoothDynamicalSystem(allDS, allInteractions);

//     // -------------
//     // --- Model ---
//     // -------------

//     Model * Manipulator = new Model(t0,T);
//     Manipulator->setNonSmoothDynamicalSystemPtr(nsds); // set NonSmoothDynamicalSystem of this model

//     // ----------------
//     // --- Simulation ---
//     // ----------------

//     // -- Time discretisation --
//     TimeDiscretisation * t = new TimeDiscretisation(h,Manipulator);
















