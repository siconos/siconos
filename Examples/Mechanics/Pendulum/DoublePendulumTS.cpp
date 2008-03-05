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

#include "SiconosKernel.h"
#include <stdlib.h>
using namespace std;

#include <boost/progress.hpp>

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
    double T = 40.0;                   // final computation time
    double h = 0.0005;                // time step
    double criterion = 0.05;
    unsigned int maxIter = 20000;
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
    SimpleVector q0(nDof), v0(nDof);
    q0.zero();
    v0.zero();
    q0(0) = 1.5;
    q0(1) = 1.5;

    LagrangianDS * doublependulum = new LagrangianDS(1, q0, v0);

    // external plug-in
    doublependulum->setComputeMassFunction("DoublePendulumPlugin.so", "mass");
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


    //     SimpleMatrix H(1,2);
    //     SimpleVector b(1);
    //     H.zero();
    //     H(0,0) =1.0;
    //     H(0,1) =0.0;

    //     b(0) = 0.0;


    //     NonSmoothLaw * nslaw = new NewtonImpactNSL(e);
    //     Relation * relation = new LagrangianLinearR(H,b);
    //     Interaction * inter =  new Interaction("floor-mass1", allDS,1,1, nslaw, relation);



    //     SimpleMatrix H2(1,2);
    //     SimpleVector b2(1);
    //     H2.zero();
    //     H2(0,0) =0.0;
    //     H2(0,1) =1;

    //     b2(0) = 0.0;

    //     NonSmoothLaw * nslaw1 = new NewtonImpactNSL(e1);
    //     Relation * relation1 = new LagrangianLinearR(H2,b2);
    //     Interaction * inter1 =  new Interaction("floor-mass2", allDS,2,1, nslaw1, relation1);





    string G = "DoublePendulumPlugin:G0";
    NonSmoothLaw * nslaw = new NewtonImpactNSL(e);
    Relation * relation = new LagrangianScleronomousR("DoublePendulumPlugin:h0", G);
    Interaction * inter = new Interaction("floor-mass1", allDS, 1, 1, nslaw, relation);

    string G1 = "DoublePendulumPlugin:G1";
    NonSmoothLaw * nslaw1 = new NewtonImpactNSL(e1);
    Relation * relation1 = new LagrangianScleronomousR("DoublePendulumPlugin:h1", G1);
    Interaction * inter1 = new Interaction("floor-mass2", allDS, 2, 1, nslaw1, relation1);

    allInteractions.insert(inter);
    allInteractions.insert(inter1);

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

    Moreau * OSI =  new Moreau(doublependulum, theta, s);

    // -- OneStepNsProblem --
    IntParameters iparam(5);
    iparam[0] = 2001; // Max number of iteration
    DoubleParameters dparam(5);
    dparam[0] = 0.0000005; // Tolerance
    string solverName = "Lemke" ;
    NonSmoothSolver * mySolver = new NonSmoothSolver(solverName, iparam, dparam);
    OneStepNSProblem * osnspb = new LCP(s, mySolver);

    cout << "=== End of model loading === " << endl;

    // =========================== End of model definition ===========================  dataPlot(k,7) = (inter->getY(0))(0);


    // ================================= Computation =================================

    // --- Simulation initialization ---
    s->initialize();
    cout << "End of simulation initialisation" << endl;

    int k = 0;
    int N = t->getNSteps(); // Number of time steps
    cout << "Number of time step   " << N << endl;
    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 12;
    SimpleMatrix dataPlot(N + 1, outputSize);
    // For the initial time step:
    // time
    SiconosVector * q = doublependulum->getQPtr();
    SiconosVector * v = doublependulum->getVelocityPtr();

    dataPlot(k, 0) =  Pendulum->getT0();
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

    boost::timer time;
    time.restart();

    // --- Time loop ---
    cout << "Start computation ... " << endl;

    boost::progress_display show_progress(N);

    while (s->hasNextEvent())
    {
      k++;
      ++show_progress;
      //  if (!(div(k,1000).rem))  cout <<"Step number "<< k << "\n";

      // Solve problem
      s->newtonSolve(criterion, maxIter);
      s->nextStep();
      // Data Output
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

    // --- Free memory ---
    delete osnspb;
    delete t;
    delete OSI;
    delete s;
    delete Pendulum;
    delete nsds;
    delete inter;
    delete inter1;
    delete relation;
    delete nslaw;
    delete relation1;
    delete nslaw1;
    delete doublependulum;
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
