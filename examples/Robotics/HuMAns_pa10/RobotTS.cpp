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

using namespace std;

int main(int argc, char* argv[])
{
  try
  {

    // ================= Creation of the model =======================

    // User-defined main parameters
    unsigned int nDof = 3;           // degrees of freedom for robot arm
    double t0 = 0;                   // initial computation time
    double T = 1.9;                   // final computation time
    double h = 0.005;                // time step
    double criterion = 0.0005;
    unsigned int maxIter = 20000;
    double e = 0.9;                  // restit. coef. for impact on the ground.
    double e2 = 0.0;                 // restit. coef for angular stops impacts.

    // -> mind to set the initial conditions below.

    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    // unsigned int i;
    // --- DS: robot arm ---

    // The dof are angles between ground and arm and between differents parts of the arm. (See Robot.fig for more details)
    //


    // Initial position (angles in radian)
    SP::SiconosVector q0(new SiconosVector(nDof)), v0(new SiconosVector(nDof));
    (*q0)(0) = 0.05;
    (*q0)(1) = 0.05;

    SP::LagrangianDS arm(new LagrangianDS(q0, v0));

    // external plug-in
    arm->setComputeMassFunction("RobotPlugin", "mass");
    arm->setComputeFGyrFunction("RobotPlugin", "FGyr");
    arm->setComputeJacobianFGyrqDotFunction("RobotPlugin", "jacobianVFGyr");
    arm->setComputeJacobianFGyrqFunction("RobotPlugin", "jacobianFGyrq");

    // -------------------
    // --- Interactions---
    // -------------------

    // Two interactions:
    //  - one with Lagrangian non linear relation to define contact with ground
    //  - the other to define angles limitations (articular stops), with lagrangian linear relation
    //  Both with newton impact ns laws.

    // -- relations --

    // => arm-floor relation
    SP::NonSmoothLaw nslaw(new NewtonImpactNSL(e));
    string G = "RobotPlugin:G2";
    SP::Relation relation(new LagrangianScleronomousR("RobotPlugin:h2", G));
    SP::Interaction inter(new Interaction(nslaw, relation, 0));

    // => angular stops

    //     SimpleMatrix H(6,3);
    //     SiconosVector b(6);
    //     H.zero();
    //     H(0,0) =-1;
    //     H(1,0) =1;
    //     H(2,1) =-1;
    //     H(3,1) =1;
    //     H(4,2) =-1;
    //     H(5,2) =1;

    //     b(0) = 1.7;
    //     b(1) = 1.7;
    //     b(2) = 0.3;
    //     b(3) = 0.3;
    //     b(4) = 3.14;
    //     b(5) = 3.14;
    double lim0 = 1.6;
    double lim1 = 3.1;  // -lim <= q[1] <= lim
    SP::SimpleMatrix H(new SimpleMatrix(4, 3));
    SP::SiconosVector b(new SiconosVector(4));
    H->zero();

    (*H)(0, 0) = -1;
    (*H)(1, 0) = 1;
    (*H)(2, 1) = -1;
    (*H)(3, 1) = 1;

    (*b)(0) = lim0;
    (*b)(1) = lim0;
    (*b)(2) = lim1;
    (*b)(3) = lim1;

    SP::NonSmoothLaw nslaw2(new NewtonImpactNSL(e2));
    SP::Relation relation2(new LagrangianLinearTIR(H, b));
    SP::Interaction inter2(new Interaction(nslaw2, relation2, 1));

    // -------------
    // --- Model ---
    // -------------

    SP::Model Robot(new Model(t0, T));

    // add the dynamical system in the non smooth dynamical system
    Robot->nonSmoothDynamicalSystem()->insertDynamicalSystem(arm);

    // link the interactions and the dynamical system
    Robot->nonSmoothDynamicalSystem()->link(inter, arm);
    Robot->nonSmoothDynamicalSystem()->link(inter2, arm);

    // ----------------
    // --- Simulation ---
    // ----------------

    // -- Time discretisation --
    SP::TimeDiscretisation t(new TimeDiscretisation(t0, h));

    SP::TimeStepping s(new TimeStepping(t));

    // -- OneStepIntegrators --
    SP::OneStepIntegrator OSI(new MoreauJeanOSI(arm, 0.500001));
    s->insertIntegrator(OSI);

    SP::OneStepNSProblem osnspb(new LCP("PGS"));

    osnspb->numericsSolverOptions()->iparam[0] = 30001;
    osnspb->numericsSolverOptions()->dparam[0] = 0.005;


    s->insertNonSmoothProblem(osnspb);
    Robot->setSimulation(s);
    cout << "=== End of model loading === " << endl;

    // =========================== End of model definition ===========================  dataPlot(k,7) = (*inter->y(0))(0);


    // ================================= Computation =================================

    // --- Simulation initialization ---
    Robot->initialize();
    cout << "End of model initialisation" << endl;

    int k = 0;
    int N = ceil((T - t0) / h) + 1;

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 9;
    SimpleMatrix dataPlot(N + 1, outputSize);
    // For the initial time step:
    // time

    SP::SiconosVector q = arm->q();
    SP::SiconosVector vel = arm->velocity();
    SP::SiconosVector y = inter->y(0);

    dataPlot(k, 0) =  Robot->t0();
    dataPlot(k, 1) = (*q)(0);
    dataPlot(k, 2) = (*vel)(0);
    dataPlot(k, 3) = (*q)(1);
    dataPlot(k, 4) = (*vel)(1);
    dataPlot(k, 5) = (*q)(2);
    dataPlot(k, 6) = (*vel)(2);
    dataPlot(k, 7) = (*y)(0);
    dataPlot(k, 8) = (*y)(1);

    // --- Time loop ---
    cout << "Start computation ... " << endl;
    boost::timer boostTimer;
    boostTimer.restart();

    while (s->hasNextEvent())
    {
      // get current time step
      k++;
      s->newtonSolve(criterion, maxIter);
      dataPlot(k, 0) =  s->nextTime();
      dataPlot(k, 1) = (*q)(0);
      dataPlot(k, 2) = (*vel)(0);
      dataPlot(k, 3) = (*q)(1);
      dataPlot(k, 4) = (*vel)(1);
      dataPlot(k, 5) = (*q)(2);
      dataPlot(k, 6) = (*vel)(2);
      dataPlot(k, 7) = (*y)(0);
      dataPlot(k, 8) = (*y)(1);
      s->nextStep();
    }

    cout << "End of computation - Number of iterations done: " << k << endl;
    cout << "Computation Time: " << boostTimer.elapsed()  << endl;

    cout << endl << "Output writing ..." << endl;
    // --- Output files ---
    ioMatrix::write("result.dat", "ascii", dataPlot, "noDim");
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
