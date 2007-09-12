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


// =============================== Double Pendulum Example ===============================
//
// Author: Vincent Acary
//
// Keywords: LagrangianDS, LagrangianLinear relation, Moreau TimeStepping, LCP.
//
// =============================================================================================

#include "SiconosKernel.h"
#include <stdlib.h>
using namespace std;


double gravity = 10.0;
double m1 = 1.0;
double l1 = 1.0 ;





int main(int argc, char* argv[])
{
  try
  {

    // ================= Creation of the model =======================

    // User-defined main parameters
    unsigned int nDof = 1;           // degrees of freedom for robot arm
    double t0 = 0;                   // initial computation time
    double T = 50.0;                   // final computation time
    double h = 0.0005;                // time step
    double criterion = 0.00005;
    unsigned int maxIter = 2000;
    double e = 0.9;                  // nslaw

    // -> mind to set the initial conditions below.

    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    // unsigned int i;
    DynamicalSystemsSet allDS; // the list of DS

    // --- DS: Double Pendulum ---

    // Initial position (angles in radian)
    SimpleVector q0(nDof), v0(nDof);
    q0.zero();
    v0.zero();
    q0(0) = 1.5;

    LagrangianDS * simplependulum = new LagrangianDS(1, q0, v0);

    SiconosMatrix *Mass = new SimpleMatrix(nDof, nDof);
    (*Mass)(0, 0) = m1 * l1;
    simplependulum->setMassPtr(Mass);


    // external plug-in
    //simplependulum->setComputeMassFunction("SimplePendulumPlugin.so","mass");


    simplependulum->setComputeFIntFunction("SimplePendulumPlugin.so", "FInt");
    simplependulum->setComputeJacobianFIntFunction(1, "SimplePendulumPlugin.so", "jacobianVFInt");
    simplependulum->setComputeJacobianFIntFunction(0, "SimplePendulumPlugin.so", "jacobianQFInt");

    allDS.insert(simplependulum);

    // -------------------
    // --- Interactions---
    // -------------------


    InteractionsSet allInteractions;

    // -- relations --


    //     SimpleMatrix H(1,2);
    //     SimpleVector b(1);
    //     H.zero();
    //     H(0,0) =1.0;
    //     H(0,1) =0.0;

    //     b(0) = 0.0;


    //     NonSmoothLaw * nslaw = new NewtonImpactNSL(e);
    //     Relation * relation = new LagrangianLinearR(H,b);
    //     Interaction * inter =  new Interaction("floor-mass1", allDS,1,1, nslaw, relation);


    string G = "DoublePendulumPlugin:G0";
    NonSmoothLaw * nslaw = new NewtonImpactNSL(e);
    Relation * relation = new LagrangianScleronomousR("SimplePendulumPlugin:h0", G);
    Interaction * inter = new Interaction("floor-mass1", allDS, 1, 1, nslaw, relation);


    allInteractions.insert(inter);


    // --------------------------------
    // --- NonSmoothDynamicalSystem ---
    // --------------------------------

    bool isBVP = false;
    NonSmoothDynamicalSystem * nsds = new NonSmoothDynamicalSystem(allDS, allInteractions, isBVP);

    // -------------
    // --- Model ---
    // -------------

    Model * Pendulum = new Model(t0, T);
    Pendulum->setNonSmoothDynamicalSystemPtr(nsds); // set NonSmoothDynamicalSystem of this model

    // ----------------
    // --- Simulation ---
    // ----------------

    // -- Time discretisation --
    TimeDiscretisation * t = new TimeDiscretisation(h, Pendulum);

    TimeStepping* s = new TimeStepping(t);

    // -- OneStepIntegrators --

    //double theta=0.500001;
    double theta = 0.500001;

    OneStepIntegrator * OSI =  new Moreau(simplependulum, theta, s);

    // -- OneStepNsProblem --
    OneStepNSProblem * osnspb = new LCP(s, "name", "Lemke", 2001, 0.005);

    cout << "=== End of model loading === " << endl;

    // =========================== End of model definition ===========================  dataPlot(k,7) = (inter->getY(0))(0);


    // ================================= Computation =================================


    // --- Simulation initialization ---
    s->initialize();
    cout << "End of simulation initialisation" << endl;

    int k = 0;
    int N = t->getNSteps(); // Number of time steps
    cout << "Number of time step" << N << endl;
    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 11;
    SimpleMatrix dataPlot(N + 1, outputSize);
    // For the initial time step:
    // time
    dataPlot(k, 0) =  Pendulum->getT0();
    dataPlot(k, 1) = simplependulum->getQ()(0);
    dataPlot(k, 2) = simplependulum->getVelocity()(0);
    dataPlot(k, 3) =  l1 * sin(simplependulum->getQ()(0));
    dataPlot(k, 4) = -l1 * cos(simplependulum->getQ()(0));
    dataPlot(k, 5) =  l1 * cos(simplependulum->getQ()(0)) * (simplependulum->getVelocity()(0));
    // --- Compute elapsed time ---
    double t1, t2, elapsed;
    struct timeval tp;
    int rtn;
    clock_t start, end;
    double elapsed2;
    start = clock();
    rtn = gettimeofday(&tp, NULL);
    t1 = (double)tp.tv_sec + (1.e-6) * tp.tv_usec;

    //    EventsManager * eventsManager = s->getEventsManagerPtr();

    // --- Time loop ---
    cout << "Start computation ... " << endl;
    cout << "Number of time step" << N << "\n";
    while (s->hasNextEvent())
    {
      k++;
      if (!(div(k, 1000).rem))  cout << "Step number " << k << "\n";

      // Solve problem
      s->newtonSolve(criterion, maxIter);

      // Data Output
      dataPlot(k, 0) =  s->getStartingTime();
      dataPlot(k, 1) = simplependulum->getQ()(0);
      dataPlot(k, 2) = simplependulum->getVelocity()(0);
      dataPlot(k, 3) =  l1 * sin(simplependulum->getQ()(0));
      dataPlot(k, 4) = -l1 * cos(simplependulum->getQ()(0));
      dataPlot(k, 5) =  l1 * cos(simplependulum->getQ()(0)) * (simplependulum->getVelocity()(0));
    }

    end = clock();
    rtn = gettimeofday(&tp, NULL);
    t2 = (double)tp.tv_sec + (1.e-6) * tp.tv_usec;
    elapsed = t2 - t1;
    elapsed2 = (end - start) / (double)CLOCKS_PER_SEC;
    cout << "time = " << elapsed << " --- cpu time " << elapsed2 << endl;
    cout << "End of computation - Number of iterations done: " << k << endl;

    // --- Output files ---
    ioMatrix out("SimplePendulumResult.dat", "ascii");
    out.write(dataPlot, "noDim");

    // --- Free memory ---
    delete osnspb;
    delete t;
    delete OSI;
    delete s;
    delete Pendulum;
    delete nsds;
    delete inter;
    //    delete relation;
    delete nslaw;
    delete simplependulum;
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
