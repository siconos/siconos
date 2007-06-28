/* Siconos-sample version 2.1.0, Copyright INRIA 2005-2006.
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

#include "Model.h"

#include "LagrangianLinearTIDS.h"
#include "LCP.h"
#include "NewtonImpactNSL.h"
#include "LagrangianLinearR.h"
#include "EventDriven.h"
#include "Lsodar.h"
#include "Interaction.h"
#include <sys/time.h>
#include <math.h>
#include <iostream>
#include <sstream>

using namespace std;

int main(int argc, char* argv[])
{
  try
  {

    // ================= Creation of the model =======================

    // User-defined main parameters
    unsigned int nDof = 3;           // degrees of freedom for robot arm
    double t0 = 0;                   // initial computation time
    double T = 1.2;                   // final computation time
    double h = 0.005;                // time step
    double e = 0.9;                  // nslaw
    double e2 = 0.1;

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

    LagrangianDS * arm = new LagrangianDS(1, nDof, q0, v0);

    // external plug-in
    arm->setComputeMassFunction("RobotPlugin.so", "mass");
    arm->setComputeNNLFunction("RobotPlugin.so", "NNL");
    arm->setComputeJacobianVelocityNNLFunction("RobotPlugin.so", "jacobianVNNL");
    arm->setComputeJacobianQNNLFunction("RobotPlugin.so", "jacobianQNNL");

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
    vector<string> G;
    G.reserve(1);
    G.push_back("RobotPlugin:G2");
    Relation * relation = new LagrangianR("scleronomic", "RobotPlugin:h2", G);
    Interaction * inter = new Interaction("floor-arm", allDS, 0, 2, nslaw, relation);

    SimpleMatrix H(2, 3);
    SimpleVector b(2);
    H.zero();
    H(0, 1) = -1;
    H(1, 1) = 1;

    b(0) = 0.3;
    b(1) = 0.3;

    NonSmoothLaw * nslaw2 = new NewtonImpactNSL(e2);
    Relation * relation2 = new LagrangianLinearR(H, b);
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

    EventDriven* s = new EventDriven(Robot);

    // -- Time discretisation --
    TimeDiscretisation * t = new TimeDiscretisation(h, s);

    // -- OneStepIntegrators --
    OneStepIntegrator * OSI =  new Lsodar(arm, s);
    OSI->setSizeMem(2); // Necessary to save pre and post impact values.

    // -- OneStepNsProblem --
    OneStepNSProblem * impact = new LCP(s, "impact", "NSQP", 20001, 0.5);
    OneStepNSProblem * acceleration = new LCP(s, "acceleration", "NSQP", 2001, 0.5);

    cout << "=== End of model loading === " << endl;

    // =========================== End of model definition ===========================  dataPlot(k,7) = (inter->getY(0))(0);


    // ================================= Computation =================================

    // --- Simulation initialization ---
    s->initialize();
    cout << "End of simulation initialisation" << endl;

    int k = 0; // Current step
    int N = 322; // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 9;
    SimpleMatrix dataPlot(N + 1, outputSize);
    // For the initial time step:
    // time
    dataPlot(k, 0) = 0;

    dataPlot(k, 1) = arm->getQ()(0);
    dataPlot(k, 2) = arm->getVelocity()(0);
    dataPlot(k, 3) = arm->getQ()(1);
    dataPlot(k, 4) = arm->getVelocity()(1);
    dataPlot(k, 5) = arm->getQ()(2);
    dataPlot(k, 6) = arm->getVelocity()(2);
    dataPlot(k, 7) = (inter->getY(0))(0);
    dataPlot(k, 8) = (inter->getY(0))(1);

    // --- Compute elapsed time ---
    double t1, t2, elapsed;
    struct timeval tp;
    int rtn;
    clock_t start, end;
    double elapsed2;
    start = clock();
    rtn = gettimeofday(&tp, NULL);
    t1 = (double)tp.tv_sec + (1.e-6) * tp.tv_usec;

    OSI->display();

    // --- Time loop ---
    EventsManager * eventsManager = s->getEventsManagerPtr();
    cout << " ==== Start of Event Driven simulation - This may take a while ... ====" << endl;

    unsigned int numberOfEvent = 0 ;
    while (eventsManager->hasNextEvent())
    {
      k++;
      s->computeOneStep();
      cout << Robot->getCurrentT() << endl;
      if (eventsManager->getCurrentEventPtr()->getType() == "NonSmoothEvent")
      {
        dataPlot(k, 0) = Robot->getCurrentT();

        dataPlot(k, 1) = (*arm->getQMemoryPtr()->getSiconosVector(1))(0);
        dataPlot(k, 2) = (*arm->getVelocityMemoryPtr()->getSiconosVector(1))(0);
        dataPlot(k, 3) = (*arm->getQMemoryPtr()->getSiconosVector(1))(1);
        dataPlot(k, 4) = (*arm->getVelocityMemoryPtr()->getSiconosVector(1))(1);
        dataPlot(k, 5) = (*arm->getQMemoryPtr()->getSiconosVector(1))(2);
        dataPlot(k, 6) = (*arm->getVelocityMemoryPtr()->getSiconosVector(1))(2);
        dataPlot(k, 7) = (inter->getY(0))(0);
        dataPlot(k, 8) = (inter->getY(0))(1);
        k++;
      }
      dataPlot(k, 0) = Robot->getCurrentT();
      dataPlot(k, 1) = arm->getQ()(0);
      dataPlot(k, 2) = arm->getVelocity()(0);
      dataPlot(k, 3) = arm->getQ()(1);
      dataPlot(k, 4) = arm->getVelocity()(1);
      dataPlot(k, 5) = arm->getQ()(2);
      dataPlot(k, 6) = arm->getVelocity()(2);
      dataPlot(k, 7) = (inter->getY(0))(0);
      dataPlot(k, 8) = (inter->getY(0))(1);
      numberOfEvent++;
    }
    cout << "===== End of Event Driven simulation. " << numberOfEvent << " events have been processed. ==== " << endl;
    cout << "Number of impacts: " << k - numberOfEvent << endl;
    end = clock();
    rtn = gettimeofday(&tp, NULL);
    t2 = (double)tp.tv_sec + (1.e-6) * tp.tv_usec;
    elapsed = t2 - t1;
    elapsed2 = (end - start) / (double)CLOCKS_PER_SEC;
    cout << "time = " << elapsed << " --- cpu time " << elapsed2 << endl;

    // --- Output files ---
    ioMatrix io("result.dat", "ascii");
    io.write(dataPlot, "noDim");

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
