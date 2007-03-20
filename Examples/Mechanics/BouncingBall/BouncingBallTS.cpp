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

/*!\file BouncingBallTS.cpp
\brief \ref EMBouncingBall - C++ input file, Time-Stepping version - V. Acary, F. Perignon.

A Ball bouncing on the ground.
Direct description of the model without XML input.
Simulation with a Time-Stepping scheme.
*/

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
    unsigned int nDof = 3;           // degrees of freedom for the ball
    double t0 = 0;                   // initial computation time
    double T = 10.0;                   // final computation time
    double h = 0.005;                // time step
    double position_init = 1.0;      // initial position for lowest bead.
    double velocity_init = 0.0;      // initial velocity for lowest bead.
    double theta = 0.5;              // theta for Moreau integrator
    string solverName = "Lemke" ;
    double R = 0.1; // Ball radius
    double m = 1; // Ball mass
    double g = 9.81; // Gravity
    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    DynamicalSystemsSet allDS; // the list of DS
    SiconosMatrix *Mass = new SimpleMatrix(nDof, nDof);
    (*Mass)(0, 0) = m;
    (*Mass)(1, 1) = m;
    (*Mass)(2, 2) = 3. / 5 * m * R * R;

    // -- Initial positions and velocities --
    SiconosVector * q0 = new SimpleVector(nDof);
    SiconosVector * v0 = new SimpleVector(nDof);
    (*q0)(0) = position_init;
    (*v0)(0) = velocity_init;
    // -- The dynamical system --
    LagrangianLinearTIDS* ball = new LagrangianLinearTIDS(0, *q0, *v0, *Mass);
    allDS.insert(ball);
    // -- Set external forces (weight) --
    SiconosVector * weight = new SimpleVector(nDof);
    (*weight)(0) = -m * g;
    ball->setFExtPtr(weight);

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

    Model * bouncingBall = new Model(t0, T);
    bouncingBall->setNonSmoothDynamicalSystemPtr(nsds); // set NonSmoothDynamicalSystem of this model

    // ------------------
    // --- Simulation ---
    // ------------------

    // -- Time discretisation --
    TimeDiscretisation * t = new TimeDiscretisation(h, bouncingBall);

    TimeStepping* s = new TimeStepping(t);

    // -- OneStepIntegrators --
    OneStepIntegrator * OSI = new Moreau(ball, theta, s);

    // -- OneStepNsProblem --
    OneStepNSProblem * osnspb = new LCP(s, "LCP", solverName, 101, 0.0001);

    cout << "=== End of model loading === " << endl;
    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    // --- Simulation initialization ---

    s->initialize();
    cout << "End of simulation initialisation" << endl;

    int N = t->getNSteps(); // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 4;
    SimpleMatrix dataPlot(N + 1, outputSize);

    SiconosVector * q = ball->getQPtr();
    SiconosVector * v = ball->getVelocityPtr();
    SiconosVector * lambda = bouncingBall->getNonSmoothDynamicalSystemPtr()->getInteractionPtr(0)->getLambdaPtr(1);

    dataPlot(0, 0) = bouncingBall->getT0();
    dataPlot(0, 1) = (*q)(0);
    dataPlot(0, 2) = (*v)(0);
    dataPlot(0, 3) = (*lambda)(0);
    cout << N << endl;
    // --- Time loop ---
    cout << "Start computation ... " << endl;
    // == Simulation, first possible writing ==
    int k = 0;
    for (k = 1 ; k < N + 1 ; ++k)
    {
      s->computeOneStep();
      // --- Get values to be plotted ---
      dataPlot(k, 0) =  bouncingBall->getCurrentT();
      dataPlot(k, 1) = (*q)(0);
      dataPlot(k, 2) = (*v)(0);
      dataPlot(k, 3) = (*lambda)(0);
      s->nextStep();
    }
    cout << "End of computation - Number of iterations done: " << k - 1 << endl;

    // == Simulation, second possible writing ==
    //     while(s->hasNextEvent())
    //       {
    //  s->computeOneStep();
    //  k++;
    //  dataPlot(k, 0) = bouncingBall->getCurrentT();
    //  dataPlot(k,1) = ball->getQ()(0);
    //  dataPlot(k,2) = ball->getVelocity()(0);
    //  dataPlot(k, 3) = (bouncingBall->getNonSmoothDynamicalSystemPtr()->getInteractionPtr(0)->getLambda(1))(0);
    //  s->nextStep();
    //  // --- Get values to be plotted ---
    //       }
    //     cout<<"End of computation - Number of iterations done: "<<k<<endl;

    // --- Output files ---
    ioMatrix io("result.dat", "ascii");
    io.write(dataPlot, "noDim");
    // --- Free memory ---
    delete osnspb;
    delete t;
    delete s;
    delete bouncingBall;
    delete nsds;
    delete relation0;
    delete nslaw0;
    delete H;
    delete ball;
    delete q0;
    delete v0;
    delete Mass;
    delete OSI;
  }

  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught in BouncingBallTS.cpp" << endl;
  }
  cout << "Computation Time " << t.elapsed()  << endl;


}
