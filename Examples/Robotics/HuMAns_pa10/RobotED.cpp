/* Siconos-sample version 3.0.0, Copyright INRIA 2005-2008.
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
// Keywords: LagrangianDS, LagrangianLinear relation, EventDriven, LCP.
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
    double e = 0.9;                  // restit. coef. for impact on the ground.
    double e2 = 0.0;                 // restit. coef for angular stops impacts.

    // -> mind to set the initial conditions below.

    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    // unsigned int i;
    DynamicalSystemsSet allDS; // the list of DS

    // --- DS: robot arm ---

    // The dof are angles between ground and arm and between differents parts of the arm. (See Robot.fig for more details)
    //


    // Initial position (angles in radian)
    SimpleVector q0(nDof), v0(nDof);
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
    //  Both with newton impact ns laws.

    InteractionsSet allInteractions; // The set of all interactions.

    // -- relations --

    // => arm-floor relation
    NonSmoothLaw * nslaw = new NewtonImpactNSL(e);
    string G = "RobotPlugin:G2";
    Relation * relation = new LagrangianScleronomousR("RobotPlugin:h2", G);
    Interaction * inter = new Interaction("floor-arm", allDS, 0, 2, nslaw, relation);

    // => angular stops

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
    double lim0 = 1.6;
    double lim1 = 3.1;  // -lim_i <= q[i] <= lim_i
    SimpleMatrix H(4, 3);
    SimpleVector b(4);
    H.zero();

    H(0, 0) = -1;
    H(1, 0) = 1;
    H(2, 1) = -1;
    H(3, 1) = 1;

    b(0) = lim0;
    b(1) = lim0;
    b(2) = lim1;
    b(3) = lim1;

    NonSmoothLaw * nslaw2 = new NewtonImpactNSL(e2);
    Relation * relation2 = new LagrangianLinearR(H, b);
    Interaction * inter2 =  new Interaction("floor-arm2", allDS, 1, 4, nslaw2, relation2);

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

    EventDriven* s = new EventDriven(t);

    // -- OneStepIntegrators --
    Lsodar * OSI =  new Lsodar(arm, s);
    // -- OneStepNsProblem --
    IntParameters iparam(5);
    iparam[0] = 20001; // Max number of iteration
    DoubleParameters dparam(5);
    dparam[0] =  0.005; // Tolerance
    string solverName = "PGS" ;
    NonSmoothSolver * mySolver = new NonSmoothSolver(solverName, iparam, dparam);
    OneStepNSProblem * impact = new LCP(s, mySolver, "impact");
    OneStepNSProblem * acceleration = new LCP(s, mySolver, "acceleration");

    cout << "=== End of model loading === " << endl;

    // =========================== End of model definition ===========================  dataPlot(k,7) = (inter->getY(0))(0);


    // ================================= Computation =================================

    // --- Simulation initialization ---
    s->setPrintStat(true);
    s->initialize();
    cout << "End of simulation initialisation" << endl;

    int k = 0;
    int N = 10630;

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 11;
    SimpleMatrix dataPlot(N + 1, outputSize);
    // For the initial time step:
    // time

    SiconosVector * q = arm->getQPtr();
    SiconosVector * vel = arm->getVelocityPtr();
    SiconosVector * y = inter->getYPtr(0);
    SiconosVector * yDot = inter->getYPtr(1);
    // When a non-smooth event occurs, pre-impact values are saved in memory vectors at pos. 1:
    SiconosVector * qMem = arm->getQMemoryPtr()->getSiconosVector(1);
    SiconosVector * velMem = arm->getVelocityMemoryPtr()->getSiconosVector(1);
    SiconosVector * yMem = inter->getYOldPtr(0);
    SiconosVector * yDotMem = inter->getYOldPtr(1);

    dataPlot(k, 0) =  Robot->getT0();
    dataPlot(k, 1) = (*q)(0);
    dataPlot(k, 2) = (*vel)(0);
    dataPlot(k, 3) = (*q)(1);
    dataPlot(k, 4) = (*vel)(1);
    dataPlot(k, 5) = (*q)(2);
    dataPlot(k, 6) = (*vel)(2);
    dataPlot(k, 7) = (*y)(0);
    dataPlot(k, 8) = (*y)(1);
    dataPlot(k, 9) = (*yDot)(0);
    dataPlot(k, 10) = (*yDot)(1);

    // --- Time loop ---
    cout << "Start computation ... " << endl;
    boost::timer boostTimer;
    boostTimer.restart();

    unsigned int numberOfEvent = 0 ;
    EventsManager * eventsManager = s->getEventsManagerPtr();
    bool nonSmooth = false;
    while (s->getNextTime() < T)
    {
      // get current time step
      k++;
      s->advanceToEvent();
      if (eventsManager->getNextEventPtr()->getType() == 2)
        nonSmooth = true;

      s->processEvents();
      // If the treated event is non smooth, we get the pre-impact state.
      if (nonSmooth)
      {
        dataPlot(k, 0) =  s->getStartingTime();
        dataPlot(k, 1) = (*qMem)(0);
        dataPlot(k, 2) = (*velMem)(0);
        dataPlot(k, 3) = (*qMem)(1);
        dataPlot(k, 4) = (*velMem)(1);
        dataPlot(k, 5) = (*qMem)(2);
        dataPlot(k, 6) = (*velMem)(2);
        dataPlot(k, 7) = (*yMem)(0);
        dataPlot(k, 8) = (*yMem)(1);
        dataPlot(k, 9) = (*yDotMem)(0);
        dataPlot(k, 10) = (*yDotMem)(1);


        k++;
      }
      dataPlot(k, 0) =  s->getStartingTime();
      dataPlot(k, 1) = (*q)(0);
      dataPlot(k, 2) = (*vel)(0);
      dataPlot(k, 3) = (*q)(1);
      dataPlot(k, 4) = (*vel)(1);
      dataPlot(k, 5) = (*q)(2);
      dataPlot(k, 6) = (*vel)(2);
      dataPlot(k, 7) = (*y)(0);
      dataPlot(k, 8) = (*y)(1);
      dataPlot(k, 9) = (*yDot)(0);
      dataPlot(k, 10) = (*yDot)(1);
      numberOfEvent++;
      //  cout << k << endl;
      //  if (k==N) break;
    }

    cout << "===== End of Event Driven simulation. " << numberOfEvent << " events have been processed. ==== " << endl << endl;
    cout << "Computation Time: " << boostTimer.elapsed()  << endl;

    cout << endl << "Output writing ..." << endl;
    // --- Output files ---
    ioMatrix out("result.dat", "ascii");
    out.write(dataPlot, "noDim");

    // --- Free memory ---
    delete impact;
    delete acceleration;
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
