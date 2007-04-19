/* Siconos-sample version 2.0.1, Copyright INRIA 2005-2006.
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

using namespace std;

int main(int argc, char* argv[])
{
  boost::timer t;
  t.restart();
  try
  {

    // ================= Creation of the model =======================

    // User-defined main parameters
    unsigned int nDof = 3;           // degrees of freedom for robot arm
    double t0 = 0;                   // initial computation time
    double T = 3;                   // final computation time
    double h = 0.005;                // time step
    double criterion = 0.005;
    unsigned int maxIter = 2000;
    double e = 0.9;                  // nslaw
    double e2 = 0.0;

    // -> mind to set the initial conditions below.

    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    // unsigned int i;
    DynamicalSystemsSet allDS; // the list of DS

    // --- DS: robot arm ---

    // The dof are angles between ground and arm and between differents parts of the arm. (See corresponding .pdf for more details)

    // Initial position (angles in radian)
    SimpleVector q0(nDof), v0(nDof);
    q0.zero();
    v0.zero();
    q0(0) = 0.05;
    q0(1) = 0.05;

    LagrangianDS * arm = new LagrangianDS(1, q0, v0);

    // external plug-in
    arm->setComputeMassFunction("RobotPlugin.so", "mass");
    arm->setComputeNNLFunction("RobotPlugin.so", "NNL");
    arm->setComputeJacobianNNLFunction(1, "RobotPlugin.so", "jacobianVNNL");
    arm->setComputeJacobianNNLFunction(0, "RobotPlugin.so", "jacobianQNNL");

    allDS.insert(arm);

    // -------------------
    // --- Interactions---
    // -------------------

    // Two interactions:
    //  - one with Lagrangian non linear relation to define contact with ground
    //  - the other to define angles limitations (articular stops), with lagrangian linear relation
    //  Both with newton impact nslaw.

    InteractionsSet allInteractions;

    // -- relations --

    NonSmoothLaw * nslaw = new NewtonImpactNSL(e);
    string G = "RobotPlugin:G2";
    Relation * relation = new LagrangianScleronomousR("RobotPlugin:h2", G);
    Interaction * inter = new Interaction("floor-arm", allDS, 0, 2, nslaw, relation);

    //     SimpleMatrix H(6,3);
    //     SimpleVector b(6);
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
    SimpleMatrix H(2, 3);
    SimpleVector b(2);
    H.zero();
    H(0, 1) = -1;
    H(1, 1) = 1;

    b(0) = 0.3;
    b(1) = 0.3;

    NonSmoothLaw * nslaw2 = new NewtonImpactNSL(e2);
    Relation * relation2 = new LagrangianLinearR(H, b);
    //     Interaction * inter2 =  new Interaction("floor-arm2", allDS,1,6, nslaw2, relation2);
    Interaction * inter2 =  new Interaction("floor-arm2", allDS, 1, 2, nslaw2, relation2);


    allInteractions.insert(inter);
    allInteractions.insert(inter2);

    // --------------------------------
    // --- NonSmoothDynamicalSystem ---
    // --------------------------------

    bool isBVP = 0;
    NonSmoothDynamicalSystem * nsds = new NonSmoothDynamicalSystem(allDS, allInteractions, isBVP);

    // -------------
    // --- Model ---
    // -------------

    Model * Robot = new Model(t0, T);
    Robot->setNonSmoothDynamicalSystemPtr(nsds); // set NonSmoothDynamicalSystem of this model

    // ----------------
    // --- Simulation ---
    // ----------------

    // -- Time discretisation --
    TimeDiscretisation * t = new TimeDiscretisation(h, Robot);

    TimeStepping* s = new TimeStepping(t);

    // -- OneStepIntegrators --
    OneStepIntegrator * OSI =  new Moreau(arm, 0.500001, s);

    // -- OneStepNsProblem --
    OneStepNSProblem * osnspb = new LCP(s, "name", "PGS", 2001, 0.005);

    cout << "=== End of model loading === " << endl;

    // =========================== End of model definition ===========================  dataPlot(k,7) = (inter->getY(0))(0);


    // ================================= Computation =================================

    // --- Simulation initialization ---
    s->initialize();
    cout << "End of simulation initialisation" << endl;

    int k = 0;
    int N = t->getNSteps(); // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 9;
    SimpleMatrix dataPlot(N + 1, outputSize);
    // For the initial time step:
    // time
    dataPlot(k, 0) =  Robot->getCurrentT();
    dataPlot(k, 1) = arm->getQ()(0);
    dataPlot(k, 2) = arm->getVelocity()(0);
    dataPlot(k, 3) = arm->getQ()(1);
    dataPlot(k, 4) = arm->getVelocity()(1);
    dataPlot(k, 5) = arm->getQ()(2);
    dataPlot(k, 6) = arm->getVelocity()(2);
    dataPlot(k, 7) = (inter->getY(0))(0);
    dataPlot(k, 8) = (inter->getY(0))(1);

    //    EventsManager * eventsManager = s->getEventsManagerPtr();

    // --- Time loop ---
    cout << "Start computation ... " << endl;
    //     bool isNewtonConverge = false;
    //     unsigned int nbNewtonStep = 0; // number of Newton iterations
    //     while(k<N)
    //       {
    //  // get current time step
    //  k++;
    //  Robot->setCurrentT( Robot->getCurrentT()+t->getH() );
    //  s->saveInMemory();
    //  isNewtonConverge = false;
    //  nbNewtonStep = 0;
    //    while((!isNewtonConverge)&&(nbNewtonStep<=maxIter))
    //    {
    //      nbNewtonStep++;
    //      s->computeFreeState();
    //      s->updateIndexSets();
    //      s->computeOneStepNSProblem("timeStepping");
    //      s->update(1);
    //      isNewtonConverge = s->newtonCheckConvergence(criterion);
    //    }
    //  if (!isNewtonConverge)
    //    cout << "Newton process stopped: reach max step number" <<endl ;

    //  // time step increment
    //  Robot->setCurrentT(Robot->getCurrentT()+t->getH()  );



    //  //  cout << k << endl;

    //  //cout << "step " << k << endl;
    //  // solve ...
    //  //  s->newtonSolve(criterion,maxIter);



    //  dataPlot(k, 0) =  k*t->getH();//Robot->getCurrentT();

    //  dataPlot(k,1) = arm->getQ()(0);
    //  dataPlot(k,2) = arm->getVelocity()(0);
    //  dataPlot(k,3) = arm->getQ()(1);
    //  dataPlot(k,4) = arm->getVelocity()(1);
    //  dataPlot(k,5) = arm->getQ()(2);
    //  dataPlot(k,6) = arm->getVelocity()(2);
    //  dataPlot(k,7) = (inter->getY(0))(0);
    //  dataPlot(k,8) = (inter->getY(0))(1);
    //       }
    while (s->hasNextEvent())
    {
      // get current time step
      k++;
      s->newtonSolve(criterion, maxIter);
      dataPlot(k, 0) =  Robot->getCurrentT();

      dataPlot(k, 1) = arm->getQ()(0);
      dataPlot(k, 2) = arm->getVelocity()(0);
      dataPlot(k, 3) = arm->getQ()(1);
      dataPlot(k, 4) = arm->getVelocity()(1);
      dataPlot(k, 5) = arm->getQ()(2);
      dataPlot(k, 6) = arm->getVelocity()(2);
      dataPlot(k, 7) = (inter->getY(0))(0);
      dataPlot(k, 8) = (inter->getY(0))(1);
    }

    cout << "End of computation - Number of iterations done: " << k << endl;

    // --- Output files ---
    ioMatrix out("result.dat", "ascii");
    out.write(dataPlot, "noDim");

    // --- Free memory ---
    delete osnspb;
    delete t;
    delete OSI;
    delete s;
    delete Robot;
    delete nsds;
    delete inter;
    delete inter2;
    delete relation;
    delete nslaw;
    delete relation2;
    delete nslaw2;
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
