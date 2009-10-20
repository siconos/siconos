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


// =============================== Double Pendulum Example ===============================
//
// Author: Vincent Acary
//
// Keywords: LagrangianDS, LagrangianLinear relation, Moreau TimeStepping, LCP.
//
// =============================================================================================

#include "SiconosKernel.hpp"
#include <stdlib.h>
using namespace std;


double gravity = 10.0;
double m1 = 1.0;
double m2 = 1.0 ;
double l1 = 1.0 ;
double l2 = 1.0 ;

int main(int argc, char* argv[])
{
  try
  {

    // ================= Creation of the model =======================

    // User-defined main parameters
    unsigned int nDof = 2;           // degrees of freedom for robot arm
    double t0 = 0;                   // initial computation time
    double T = 50.0;                   // final computation time
    double h = 0.005;                // time step
    double e = 1.0;                  // nslaw
    double e1 = 0.0;

    // -> mind to set the initial conditions below.

    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    // unsigned int i;
    DynamicalSystemsSet allDS; // the list of DS

    // --- DS: Double Pendulum ---

    // Initial position (angles in radian)
    SP::SimpleVector q0(new SimpleVector(nDof));
    SP::SimpleVector v0(new SimpleVector(nDof));
    q0->zero();
    v0->zero();
    (*q0)(0) = 1.5;
    (*q0)(1) = 1.5;

    SP::LagrangianDS doublependulum(new LagrangianDS(q0, v0, "DoublePendulumPlugin:mass"));

    // external plug-in
    doublependulum->setComputeNNLFunction("DoublePendulumPlugin.so", "NNL");
    doublependulum->setComputeJacobianNNLFunction(1, "DoublePendulumPlugin.so", "jacobianVNNL");
    doublependulum->setComputeJacobianNNLFunction(0, "DoublePendulumPlugin.so", "jacobianQNNL");
    doublependulum->setComputeFIntFunction("DoublePendulumPlugin.so", "FInt");
    doublependulum->setComputeJacobianFIntFunction(1, "DoublePendulumPlugin.so", "jacobianVFInt");
    doublependulum->setComputeJacobianFIntFunction(0, "DoublePendulumPlugin.so", "jacobianQFInt");

    allDS.insert(doublependulum);

    // -------------------
    // --- Interactions---
    // -------------------


    InteractionsSet allInteractions;

    // -- relations --
    string G = "DoublePendulumPlugin:G0";
    SP::NonSmoothLaw nslaw(new NewtonImpactNSL(e));
    SP::Relation relation(new LagrangianScleronomousR("DoublePendulumPlugin:h0", G));
    SP::Interaction inter(new Interaction("floor-mass1", allDS, 1, 1, nslaw, relation));

    string G1 = "DoublePendulumPlugin:G1";
    SP::NonSmoothLaw nslaw1(new NewtonImpactNSL(e1));
    SP::Relation relation1(new LagrangianScleronomousR("DoublePendulumPlugin:h1", G1));
    SP::Interaction inter1(new Interaction("floor-mass2", allDS, 2, 1, nslaw1, relation1));

    allInteractions.insert(inter);
    allInteractions.insert(inter1);

    // -------------
    // --- Model ---
    // -------------

    SP::Model Pendulum(new Model(t0, T, allDS, allInteractions));

    // ----------------
    // --- Simulation ---
    // ----------------

    // -- Time discretisation --
    SP::TimeDiscretisation t(new TimeDiscretisation(t0, h));

    SP::EventDriven s(new EventDriven(t));

    // -- OneStepIntegrators --
    SP::Lsodar OSI(new Lsodar(doublependulum));
    s->recordIntegrator(OSI);

    // -- OneStepNsProblem --
    IntParameters iparam(5);
    iparam[0] = 1010; // Max number of iteration
    DoubleParameters dparam(5);
    dparam[0] = 0.0001; // Tolerance
    string solverName = "Lemke" ;
    SP::NonSmoothSolver mySolver(new NonSmoothSolver(solverName, iparam, dparam));
    SP::OneStepNSProblem impact(new LCP(mySolver, "impact"));
    SP::OneStepNSProblem acceleration(new LCP(mySolver, "acceleration"));
    s->recordNonSmoothProblem(impact);
    s->recordNonSmoothProblem(acceleration);

    cout << "=== End of model loading === " << endl;

    // =========================== End of model definition ===========================  dataPlot(k,7) = (inter->getY(0))(0);


    // ================================= Computation =================================

    // --- Simulation initialization ---
    Pendulum->initialize(s);
    cout << "End of simulation initialisation" << endl;

    int k = 0;
    int N = 10363; // Number of time steps
    cout << "Number of time step" << N << endl;
    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 11;

    SimpleMatrix dataPlot(N, outputSize);
    // For the initial time step:
    // time
    SP::SiconosVector q = doublependulum->getQPtr();
    SP::SiconosVector v = doublependulum->getVelocityPtr();

    dataPlot(k, 0) =  t0;
    dataPlot(k, 1) = (*q)(0);
    dataPlot(k, 2) = (*v)(0);
    dataPlot(k, 3) = (*q)(1);
    dataPlot(k, 4) = (*v)(1);
    dataPlot(k, 5) =  l1 * sin((*q)(0));
    dataPlot(k, 6) = -l1 * cos((*q)(0));
    dataPlot(k, 7) =  l1 * sin((*q)(0)) + l2 * sin((*q)(1));
    dataPlot(k, 8) = -l1 * cos((*q)(0)) - l2 * cos((*q)(1));
    dataPlot(k, 9) =  l1 * cos((*q)(0)) * ((*v)(0));
    dataPlot(k, 10) = l1 * cos((*q)(0)) * ((*v)(0)) + l2 * cos((*q)(1)) * ((*v)(1));
    // --- Compute elapsed time ---
    boost::timer time;
    time.restart();
    // --- Time loop ---
    cout << "Start computation ... " << endl;
    bool nonSmooth = false;

    SP::EventsManager eventsManager = s->getEventsManagerPtr();
    while (s->getNextTime() < T)
    {
      k++;
      //  if (!(div(k,1000).rem))  cout <<"Step number "<< k << "\n";

      s->advanceToEvent();
      if (eventsManager->getNextEventPtr()->getType() == 2)
        nonSmooth = true;
      s->processEvents();

      if (nonSmooth)
      {
        dataPlot(k, 0) =  s->getStartingTime();
        dataPlot(k, 1) = (*q)(0);
        dataPlot(k, 2) = (*v)(0);
        dataPlot(k, 3) = (*q)(1);
        dataPlot(k, 4) = (*v)(1);
        dataPlot(k, 5) =  l1 * sin((*q)(0));
        dataPlot(k, 6) = -l1 * cos((*q)(0));
        dataPlot(k, 7) =  l1 * sin((*q)(0)) + l2 * sin((*q)(1));
        dataPlot(k, 8) = -l1 * cos((*q)(0)) - l2 * cos((*q)(1));
        dataPlot(k, 9) =  l1 * cos((*q)(0)) * ((*v)(0));
        dataPlot(k, 10) = l1 * cos((*q)(0)) * ((*v)(0)) + l2 * cos((*q)(1)) * ((*v)(1));
        k++;
      }

      dataPlot(k, 0) =  s->getStartingTime();
      dataPlot(k, 1) = (*q)(0);
      dataPlot(k, 2) = (*v)(0);
      dataPlot(k, 3) = (*q)(1);
      dataPlot(k, 4) = (*v)(1);
      dataPlot(k, 5) =  l1 * sin((*q)(0));
      dataPlot(k, 6) = -l1 * cos((*q)(0));
      dataPlot(k, 7) =  l1 * sin((*q)(0)) + l2 * sin((*q)(1));
      dataPlot(k, 8) = -l1 * cos((*q)(0)) - l2 * cos((*q)(1));
      dataPlot(k, 9) =  l1 * cos((*q)(0)) * ((*v)(0));
      dataPlot(k, 10) = l1 * cos((*q)(0)) * ((*v)(0)) + l2 * cos((*q)(1)) * ((*v)(1));

    }

    cout << "End of computation - Number of iterations done: " << k << endl;
    cout << "Computation Time " << time.elapsed()  << endl;

    // --- Output files ---
    ioMatrix out("DoublePendulumResult.dat", "ascii");
    out.write(dataPlot, "noDim");
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
