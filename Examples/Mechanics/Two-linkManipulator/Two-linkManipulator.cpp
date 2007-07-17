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

#include "SiconosKernel.h"
#include <math.h>

#define PI 3.14159265

using namespace std;


int main(int argc, char* argv[])
{
  try
  {

    // ================= Creation of the model =======================

    // User-defined main parameters
    unsigned int nDof = 2;           // degrees of freedom for robot arm
    double t0 = 0;                   // initial computation time
    double T = 30;                   // final computation time
    double h = 0.001;                // time step
    double criterion = 0.0001;
    unsigned int maxIter = 20000;
    double e = 0.7;                  // nslaw
    double e2 = 0.0;
    int test = 0;

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
    q0(1) = -1.6;
    SiconosVector * z = new SimpleVector(nDof * 7);
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
    (*z)(11) = 0;
    (*z)(12) = 0;
    (*z)(13) = 0;


    LagrangianDS * arm = new LagrangianDS(1, q0, v0);

    // external plug-in
    arm->setComputeMassFunction("Two-linkPlugin.so", "mass");
    arm->setComputeNNLFunction("Two-linkPlugin.so", "NNL");
    arm->setComputeJacobianNNLFunction(1, "Two-linkPlugin.so", "jacobianVNNL");
    arm->setComputeJacobianNNLFunction(0, "Two-linkPlugin.so", "jacobianQNNL");
    arm->setComputeFExtFunction("Two-linkPlugin.so", "U");
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
    Relation * relation = new LagrangianScleronomousR("Two-linkPlugin:h0", "Two-linkPlugin:G0");
    Interaction * inter = new Interaction("floor-arm", allDS, 0, 2, nslaw, relation);


    SimpleMatrix H1(2, 2);
    SimpleVector b1(2);
    H1.zero();
    H1(0, 0) = -1;
    H1(1, 0) = 1;

    b1(0) = PI;
    b1(1) = 0;

    NonSmoothLaw * nslaw2 = new NewtonImpactNSL(e2);
    Relation * relation1 = new LagrangianLinearR(H1, b1);
    Interaction * inter1 =  new Interaction("floor-arm2", allDS, 1, 2, nslaw2, relation1);

    SimpleMatrix H2(2, 2);
    SimpleVector b2(2);
    H2.zero();
    H2(0, 1) = -1;
    H2(1, 1) = 1;

    b2(0) = 0.0001;
    b2(1) = PI - 0.0001;


    Relation * relation2 = new LagrangianLinearR(H2, b2);
    Interaction * inter2 =  new Interaction("singular-points", allDS, 2, 2, nslaw2, relation2);

    //  NonSmoothLaw * nslaw = new NewtonImpactNSL(e);
    //      Relation * relation = new LagrangianScleronomousR("Two-linkPlugin:h0","Two-linkPlugin:G0");
    //      Interaction * inter = new Interaction("floor-arm", allDS,0,1, nslaw, relation);

    allInteractions.insert(inter);
    allInteractions.insert(inter1);
    allInteractions.insert(inter2);
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
    OneStepNSProblem * osnspb = new LCP(s, "name", "Lemke", 200001, 0.0001);

    cout << "=== End of model loading === " << endl;

    // =========================== End of model definition ===========================  dataPlot(k,7) = (inter->getY(0))(0);


    // ================================= Computation
    // --- Simulation initialization ---



    s->initialize();
    cout << "End of simulation initialisation" << endl;

    int k = 0;
    int N = t->getNSteps(); // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 11;
    SimpleMatrix dataPlot(N + 1, outputSize);
    // For the initial time step:
    // time

    SiconosVector * q = arm->getQPtr();
    SiconosVector * v = arm->getVelocityPtr();
    // EventsManager * eventsManager = s->getEventsManagerPtr();

    dataPlot(k, 0) =  Manipulator->getCurrentT();
    dataPlot(k, 1) = (*q)(0);
    dataPlot(k, 2) = (*q)(1);
    dataPlot(k, 3) = (inter->getY(0))(1);
    dataPlot(k, 4) = (*v)(0);
    dataPlot(k, 5) = (*v)(1);
    dataPlot(k, 6) = (inter->getY(0))(0) - 2;
    dataPlot(k, 7) = (inter->getY(1))(1);
    dataPlot(k, 8) = (*z)(6);
    dataPlot(k, 9) = (*z)(4);
    dataPlot(k, 10) = test;

    bool isNewtonConverge = false;
    unsigned int nbNewtonStep = 0; // number of Newton iterations


    while (s->hasNextEvent())
    {
      (*z)(0) = (*q)(0);
      (*z)(1) = (*q)(1);
      (*z)(2) = (*v)(0);
      (*z)(3) = (*v)(1);

      // get current time step
      k++;

      dataPlot(k, 0) =  Manipulator->getCurrentT();
      dataPlot(k, 1) = (*q)(0);
      dataPlot(k, 2) = (*q)(1);
      dataPlot(k, 3) = (inter->getY(0))(1);
      dataPlot(k, 4) = (*v)(0);
      dataPlot(k, 5) = (*v)(1);
      dataPlot(k, 6) = (inter->getY(0))(0) - 2;
      dataPlot(k, 7) = (inter->getY(1))(1);
      dataPlot(k, 8) = (*z)(6);
      dataPlot(k, 9) = (*z)(4);
      dataPlot(k, 10) = test;

      isNewtonConverge = false;
      nbNewtonStep = 0; // number of Newton iterations


      s->newtonSolve(criterion, maxIter);
      (*z)(4) = (inter->getLambdaOld(1))(1);
      //  controller during impacts accumulation phase before the first impact
      if ((- dataPlot(k, 0) + trunc(dataPlot(k, 0) / (*z)(11)) * (*z)(11) + (*z)(11) / 2 <= 0.1) &&
          (test == 0))
      {
        (*z)(8) = dataPlot(k, 0) + h;
        (*z)(5) =  0.65 + 0.1 * cos(2 * PI * (*z)(8) / (*z)(11));
        (*z)(10) =  0.1 * sin(2 * PI * (*z)(8) / (*z)(11));
        (*z)(13) = 2 * 0.1 * (PI / (*z)(11)) * cos(2 * PI * (*z)(8) / (*z)(11));
        (*z)(7) = (*z)(9);
        arm->setComputeFExtFunction("Two-linkPlugin.so", "U10");
        test = 1;
      }

      //  controller during impacts accumulation phase after the first impact
      if (((*z)(4) >= 1e-12) && (test == 1))
      {
        arm->setComputeFExtFunction("Two-linkPlugin.so", "U11");
        test = 2;
      }

      // controller during constraint-motion phase.
      if ((fabs((inter->getY(1))(1)) <= 1e-4) && (test == 2) && ((inter->getY(0))(1) <= 1e-6))
      {
        (*z)(8) = dataPlot(k, 0);
        arm->setComputeFExtFunction("Two-linkPlugin.so", "U2");
        test = 3;
      }

      // change of control law with a particular design of the desired trajectory that guarantee the take-off
      if ((trunc((dataPlot(k, 0) + h) / (*z)(11)) > trunc((dataPlot(k, 0)) / (*z)(11))) && (test == 3))
      {
        (*z)(10) = dataPlot(k, 0) + h;
        (*z)(8) = (*z)(12);
        arm->setComputeFExtFunction("Two-linkPlugin.so", "U3");
        test = 4;
      }

      // change of desired trajectory during free-motion phase
      if (((*z)(9) > 0.999) && (test == 4))
      {
        arm->setComputeFExtFunction("Two-linkPlugin.so", "U");
        test = 0;
        (*z)(9) = 0;
      }
    }

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
    delete inter1;
    delete inter2;
    delete relation;
    delete relation1;
    delete relation2;
    delete nslaw;
    delete nslaw2;
    delete z;
    delete arm;
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
