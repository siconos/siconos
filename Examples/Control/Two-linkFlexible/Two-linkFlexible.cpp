/* Siconos-sample version 2.1.1, Copyright INRIA 2005-2007.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
 */


// =============================== Robot arm sample (HuMAnsPa10) ===============================
//
// see modelRobot1.jpg for complete system view.
//
// Keywords: LagrangianDS, LagrangianLinear relation, Moreau TimeStepping, LCP.
//
// =============================================================================================
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

    // unsigned int i;
    DynamicalSystemsSet allDS; // the list of DS

    // --- DS: manipulator arm ---

    // The dof are angles between ground and arm and between differents parts of the arm. (See corresponding .pdf for more details)

    // Initial position (angles in radian)
    SimpleVector q0(nDof), v0(nDof);
    q0.zero();
    v0.zero();
    q0(0) = 0.9;
    q0(1) = -1.5;
    q0(2) = 0.9;
    q0(3) = -1.5;
    SiconosVector * z = new SimpleVector(nDof * 5);
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


    LagrangianDS * arm = new LagrangianDS(1, q0, v0);

    // external plug-in
    arm->setComputeMassFunction("Two-linkFlexiblePlugin.so", "mass");
    arm->setComputeNNLFunction("Two-linkFlexiblePlugin.so", "NNL");
    arm->setComputeJacobianNNLFunction(1, "Two-linkFlexiblePlugin.so", "jacobianVNNL");
    arm->setComputeJacobianNNLFunction(0, "Two-linkFlexiblePlugin.so", "jacobianQNNL");
    arm->setComputeFIntFunction("Two-linkFlexiblePlugin.so", "U");
    arm->setComputeJacobianFIntFunction(1, "Two-linkFlexiblePlugin.so", "jacobFintV");
    arm->setComputeJacobianFIntFunction(0, "Two-linkFlexiblePlugin.so", "jacobFintQ");
    arm->setZPtr(z);

    allDS.insert(arm);

    // -------------------
    // --- Interactions---
    // -------------------

    //  - one with Lagrangian non linear relation to define contact with ground
    //  Both with newton impact nslaw.

    InteractionsSet allInteractions;

    // -- relations --

    NonSmoothLaw * nslaw = new NewtonImpactNSL(e);
    Relation * relation = new LagrangianScleronomousR("Two-linkFlexiblePlugin:h0", "Two-linkFlexiblePlugin:G0");
    Interaction * inter = new Interaction("floor-arm", allDS, 0, 2, nslaw, relation);
    //  Relation * relation0 = new LagrangianScleronomousR("Two-linkFlexiblePlugin:h3","Two-linkFlexiblePlugin:G3");
    //  Interaction * inter0 = new Interaction("wall-arm", allDS,1,2, nslaw, relation0);


    allInteractions.insert(inter);
    // allInteractions.insert(inter0);

    // --------------------------------
    // --- NonSmoothDynamicalSystem ---
    // -------------------------------

    NonSmoothDynamicalSystem * nsds = new NonSmoothDynamicalSystem(allDS, allInteractions);

    // -------------
    // --- Model ---
    // -------------

    Model * Manipulator = new Model(t0, T);
    Manipulator->setNonSmoothDynamicalSystemPtr(nsds); // set NonSmoothDynamicalSystem of this model

    // ----------------
    // --- Simulation ---
    // ----------------

    // -- Time discretisation --
    TimeDiscretisation * t = new TimeDiscretisation(h, Manipulator);

    TimeStepping* s = new TimeStepping(t);

    // -- OneStepIntegrators --
    OneStepIntegrator * OSI =  new Moreau(arm, 0.500001, s);

    // -- OneStepNsProblem --
    IntParameters iparam(5);
    iparam[0] = 20000; // Max number of iteration
    DoubleParameters dparam(5);
    dparam[0] = 1e-8; // Tolerance
    string solverName = "Lemke" ;
    NonSmoothSolver * mySolver = new NonSmoothSolver(solverName, iparam, dparam);
    // -- OneStepNsProblem --
    OneStepNSProblem * osnspb = new LCP(s, mySolver);

    cout << "=== End of model loading === " << endl;

    // =========================== End of model definition ===========================


    // ================================= Computation
    // --- Simulation initialization ---



    s->initialize();
    cout << "End of simulation initialisation" << endl;

    int k = 0;
    int N = (int)((T - t0) / h); // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 14;
    SimpleMatrix dataPlot(N + 1, outputSize);
    // For the initial time step:
    // time

    SiconosVector * q = arm->getQPtr();
    SiconosVector * v = arm->getVelocityPtr();
    SiconosVector * p = arm->getPPtr(2);
    // EventsManager * eventsManager = s->getEventsManagerPtr();

    dataPlot(k, 0) =  Manipulator->getT0();
    dataPlot(k, 1) = (*q)(0);
    dataPlot(k, 2) = (*q)(1);
    dataPlot(k, 3) = (inter->getY(0))(1);
    dataPlot(k, 4) = (*v)(0);
    dataPlot(k, 5) = (*v)(1);
    dataPlot(k, 6) = (inter->getY(0))(0) - 2;
    dataPlot(k, 7) = nimpact; //(inter->getY(1))(1);
    dataPlot(k, 8) = (*z)(6);
    dataPlot(k, 9) = (*z)(4); //L
    dataPlot(k, 10) = test;
    dataPlot(k, 11) = (*p)(1);
    dataPlot(k, 12) = (*z)(14);
    dataPlot(k, 13) = (*z)(15);



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

      dataPlot(k, 0) =  s->getNextTime();
      dataPlot(k, 1) = (*q)(0);
      dataPlot(k, 2) = (*q)(1);
      dataPlot(k, 3) = (inter->getY(0))(1);
      dataPlot(k, 4) = (*v)(0);
      dataPlot(k, 5) = (*v)(1);
      dataPlot(k, 6) = (inter->getY(0))(0) - 2;
      dataPlot(k, 7) = nimpact; //(inter->getY(1))(1);
      dataPlot(k, 8) = (*z)(6);
      if (test == 3) dataPlot(k, 9) = (*z)(4) / h;
      else dataPlot(k, 9) = (*z)(4);
      dataPlot(k, 10) = test;
      dataPlot(k, 12) = (*z)(14);
      dataPlot(k, 13) = (*z)(15);


      s->newtonSolve(criterion, maxIter);
      dataPlot(k, 11) = (*p)(1);
      (*z)(4) = (inter->getLambda(1))(1);
      s->nextStep();

      //    controller during impacts accumulation phase before the first impact
      if ((dataPlot(k, 13) <= 0.05) && (test == 0) && (dataPlot(k, 6) < 0.3))
      {
        (*z)(8) = dataPlot(k, 0);
        (*z)(5) = (*z)(14);
        // (*z)(10)= dataPlot(k,3);
        (*z)(7) = (*z)(9);
        arm->setComputeFIntFunction("Two-linkFlexiblePlugin.so", "U1");
        test = 1;
      }

      //controller during impacts accumulation phase after the first impact
      if ((dataPlot(k - 1, 11) > 0) && (test == 1))
      {
        //(*z)(8) = dataPlot(k-1,0);
        arm->setComputeFIntFunction("Two-linkFlexiblePlugin.so", "U2");
        test = 2;
      }
      if ((dataPlot(k, 11) > 0) && (test == 2))
        nimpact = nimpact + 1;

      // controller during constraint-motion phase.
      if ((dataPlot(k, 11) > 0) && (test == 2) && (dataPlot(k, 7) - dataPlot(k - 1, 7) == 1)) //  && (fabs((inter0->getY(1))(0))<1e-6))
      {
        // L= dataPlot(k,0)-(*z)(8);
        (*z)(8) = dataPlot(k, 0);
        arm->setComputeFIntFunction("Two-linkFlexiblePlugin.so", "U3");
        test = 3;
        nimpact = 0;
      }

      // change of control law with a particular design of the desired trajectory that guarantee the take-off
      if ((trunc((dataPlot(k, 0) + h) / (*z)(11)) > trunc((dataPlot(k, 0)) / (*z)(11))) && (test == 3))
      {
        (*z)(8) = dataPlot(k, 0) + h;
        (*z)(10) = (*z)(12);
        arm->setComputeFIntFunction("Two-linkFlexiblePlugin.so", "U4");
        test = 4;
        // L = 0;
      }

      //  controller during free-motion phase
      if (((*z)(13) - 0.01 >= 0) && (test == 4))
      {
        arm->setComputeFIntFunction("Two-linkFlexiblePlugin.so", "U");
        test = 0;
        (*z)(13) = 0;
      }
    }
    cout << endl << "End of computation - Number of iterations done: " << k << endl;
    cout << "Computation Time " << time.elapsed()  << endl;
    // --- Output files ---
    ioMatrix out("result.dat", "ascii");
    out.write(dataPlot, "noDim");

    //     // --- Free memory ---
    delete osnspb;
    delete t;
    delete OSI;
    delete s;
    delete Manipulator;
    delete nsds;
    delete inter;
    delete relation;
    delete nslaw;
    delete z;
    delete arm;
  }

  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
    <<< <<< < .mine
  {cout << "Exception caught in TwolinkFlexManip" << endl;}
  == == == =
  {cout << "Exception caught in TwolinkFlexManipulator" << endl;}
}





















// SiconosVector * z = new SimpleVector(nDof*4);
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
//     arm->setComputeMassFunction("Two-linkFlexiblePlugin.so","mass");
//     arm->setComputeNNLFunction("Two-linkFlexiblePlugin.so","NNL");
//     arm->setComputeJacobianNNLFunction(1,"Two-linkFlexiblePlugin.so","jacobianVNNL");
//     arm->setComputeJacobianNNLFunction(0,"Two-linkFlexiblePlugin.so","jacobianQNNL");
//     arm->setComputeFExtFunction("Two-linkFlexiblePlugin.so","U");
//     arm->setZPtr(z);

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
// //     SimpleVector b1(2);
// //     H1.zero();
// //     H1(0,0) =-1;
// //     H1(1,0) =1;

// //     b1(0) = PI;
// //     b1(1) = 0;


// //     NonSmoothLaw * nslaw2 = new NewtonImpactNSL(e2);
// //     Relation * relation1 = new LagrangianLinearTIR(H1,b1);
// //     Interaction * inter1 =  new Interaction("floor-arm2", allDS,1,2, nslaw2, relation1);

// //     SimpleMatrix H2(2,4);
// //     SimpleVector b2(2);
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
>>> >>> > .r1364


}















