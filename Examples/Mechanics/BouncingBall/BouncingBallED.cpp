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

/*!\file BouncingBallED.cpp
\brief \ref EMBouncingBall - C++ input file, Event-Driven version - V. Acary, F. Perignon.

A Ball bouncing on the ground.
Direct description of the model without XML input.
Simulation with an Event-Driven scheme.
*/

#include "SiconosKernel.h"

using namespace std;

int main(int argc, char* argv[])
{
  try
  {

    // ================= Creation of the model =======================

    // User-defined main parameters
    unsigned int dsNumber = 1;       // the bouncing ball and the ground
    unsigned int nDof = 3;           // degrees of freedom for the ball
    double t0 = 0;                   // initial computation time
    double T = 10;                    // final computation time
    double h = 0.05;                // time step
    double position_init = 1.0;      // initial position for lowest bead.
    double velocity_init = 0.0;      // initial velocity for lowest bead.
    string solverName = "Lemke" ;
    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    DynamicalSystemsSet allDS; // the list of DS
    SiconosMatrix *Mass = new SimpleMatrix(nDof, nDof);
    Mass->eye();

    // -- Initial positions and velocities --
    vector<SimpleVector *> q0;
    vector<SimpleVector *> velocity0;
    q0.resize(dsNumber, NULL);
    velocity0.resize(dsNumber, NULL);
    q0[0] = new SimpleVector(nDof);
    velocity0[0] = new SimpleVector(nDof);
    (*(q0[0]))(0) = position_init;
    (*(velocity0[0]))(0) = velocity_init;
    LagrangianLinearTIDS* lds = new LagrangianLinearTIDS(0, *(q0[0]), *(velocity0[0]), *Mass);

    allDS.insert(lds);
    lds->setComputeFExtFunction("BallPlugin.so", "ballFExt");

    // --------------------
    // --- Interactions ---
    // --------------------

    InteractionsSet allInteractions;

    // -- nslaw --
    double e = 0.9;

    // Interaction ball-floor
    //
    SiconosMatrix *H = new SimpleMatrix(1, nDof);
    (*H)(0, 0) = 1.0;
    NonSmoothLaw * nslaw0 = new NewtonImpactNSL(e);
    Relation * relation0 = new LagrangianLinearR(*H);

    Interaction * inter = new Interaction("floor-ball", allDS, 0, 1, nslaw0, relation0);
    allInteractions.insert(inter);

    // --------------------------------
    // --- NonSmoothDynamicalSystem ---
    // --------------------------------
    bool isBVP = 0;
    NonSmoothDynamicalSystem * nsds = new NonSmoothDynamicalSystem(allDS, allInteractions, isBVP);

    // -------------
    // --- Model ---
    // -------------

    Model * ball = new Model(t0, T);
    ball->setNonSmoothDynamicalSystemPtr(nsds); // set NonSmoothDynamicalSystem of this model

    // ----------------
    // --- Simulation ---
    // ----------------

    // -- Time discretisation --
    TimeDiscretisation * t = new TimeDiscretisation(h, ball);

    EventDriven* s = new EventDriven(t);

    // -- OneStepIntegrators --
    OneStepIntegrator * OSI = new Lsodar(lds, s);

    // -- OneStepNsProblem --
    OneStepNSProblem * impact = new LCP(s, "impact", solverName, 101, 0.0001, "max", 0.6);
    OneStepNSProblem * acceleration = new LCP(s, "acceleration", solverName, 101, 0.0001, "max", 0.6);

    cout << "=== End of model loading === " << endl;
    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    // --- Simulation initialization ---
    s->initialize();

    cout << "End of simulation initialisation" << endl;

    int N = 10534; // Number of saved points: depends on the number of events ...

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 4;
    SimpleMatrix dataPlot(N + 1, outputSize);
    // For the initial time step:
    // time

    int k = 0;
    dataPlot(k, 0) = s->getCurrentTime();
    dataPlot(k, 1) = lds->getX()(0);
    dataPlot(k, 2) = lds->getX()(3);
    dataPlot(k, 3) = (ball->getNonSmoothDynamicalSystemPtr()->getInteractionPtr(0)->getLambda(1))(0);

    // --- Time loop ---
    cout << " ==== Start of Event Driven simulation - This may take a while ... ====" << endl;
    // --- Compute elapsed time ---
    double t1, t2, elapsed;
    struct timeval tp;
    int rtn;
    clock_t start, end;
    double elapsed2;
    start = clock();
    rtn = gettimeofday(&tp, NULL);
    t1 = (double)tp.tv_sec + (1.e-6) * tp.tv_usec;

    EventsManager * eventsManager = s->getEventsManagerPtr();
    unsigned int numberOfEvent = 0 ;
    while (s->hasNextEvent())
    {
      k++;
      s->advanceToEvent();

      s->processEvents();
      // If the treated event is non smooth, we save pre-impact state.
      if (eventsManager->getCurrentEventPtr()->getType() == "NonSmoothEvent")
      {
        dataPlot(k, 0) = s->getCurrentTime();
        dataPlot(k, 1) = (*lds->getQMemoryPtr()->getSiconosVector(1))(0);
        dataPlot(k, 2) = (*lds->getVelocityMemoryPtr()->getSiconosVector(1))(0);
        k++;
      }

      dataPlot(k, 0) = s->getCurrentTime();
      dataPlot(k, 1) = lds->getX()(0);
      dataPlot(k, 2) = lds->getX()(3);
      dataPlot(k, 3) = (ball->getNonSmoothDynamicalSystemPtr()->getInteractionPtr(0)->getLambda(1))(0);

      numberOfEvent++;
    }
    end = clock();
    rtn = gettimeofday(&tp, NULL);
    t2 = (double)tp.tv_sec + (1.e-6) * tp.tv_usec;
    elapsed = t2 - t1;
    elapsed2 = (end - start) / (double)CLOCKS_PER_SEC;
    cout << "time = " << elapsed << " --- cpu time " << elapsed2 << endl;
    // --- Output files ---
    ioMatrix io("result.dat", "ascii");
    io.write(dataPlot, "noDim");
    cout << "===== End of Event Driven simulation. " << numberOfEvent << " events have been processed. ==== " << endl;

    // --- Free memory ---

    delete impact;
    delete acceleration;
    delete t;
    delete s;
    delete ball;
    delete nsds;
    delete relation0;
    delete nslaw0;
    delete H;
    delete lds;
    delete q0[0];
    delete velocity0[0];
    delete Mass;
    delete OSI;
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
