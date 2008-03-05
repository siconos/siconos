/* Siconos-sample version 2.0.0, Copyright INRIA 2005-2008.
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
 *
 * 2 d.o.f system simulation - 3D frictionl contact problem in presence of a rigid foundations
 *
 * 21/06/2007- Authors: houari khenous
 *
 * Last modification 27/11/2007, H. Khenous
*/
// =============================== Two dof oscillator simulation ===============================
//
//                  v1 --->        <--- v2
//  |               m1                  m2             |
//  |                _     spring (k)   _              |
//  |_______________(_)xxxxxxxxxxxxxxxx(_)_____________|
//
//  0                                                  1
//  |--------------------------------------------------|------> x-axis
//
//
// =============================== 2dof oscillator simulation ===============================
//
// Keywords: LagrangianLinearTIDS relation, Moreau TimeStepping, NLGS, NLGSNEWTON.
//
// ======================================================================================================

#include "SiconosKernel.h"

using namespace std;


int main(int argc, char* argv[])
{
  boost::timer time;
  time.restart();
  try
  {

    // ================= Creation of the model =======================

    // User-defined main parameters
    unsigned int nDof = 2;           // degrees of freedom for the ball
    double t0 = 0;                   // initial computation time
    double T = 100.0;                 // final computation time
    double h = 0.005;                 // time step

    double pos1 = 0.2;               // initial position for m1.
    double v1   = 0.0;               // initial velocity for m1
    double pos2 = 0.8;               // initial position for m2
    double v2   = -0.1;               // initial velocity for m2

    double k = 0.2; // stiffness coefficient
    double L = 0.5; // initial lenth
    double m1 = 1.; //  m1
    double m2 = 1.; //  m2

    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    cout << "====> Model loading ..." << endl << endl;
    DynamicalSystemsSet allDS; // the list of DS
    SiconosMatrix *M = new SimpleMatrix(nDof, nDof);
    (*M)(0, 0) = m1;
    (*M)(1, 1) = m2;
    SiconosMatrix *K = new SimpleMatrix(nDof, nDof);
    (*K)(0, 0) = k;
    (*K)(0, 1) = -k;
    (*K)(1, 0) = -k;
    (*K)(1, 1) = k;
    SiconosMatrix *C = new SimpleMatrix(nDof, nDof);

    // -- Initial positions and velocities --
    SiconosVector * q0 = new SimpleVector(nDof);
    SiconosVector * v0 = new SimpleVector(nDof);
    (*q0)(0) = pos1;
    (*v0)(0) = v1;
    (*q0)(1) = pos2;
    (*v0)(1) = v2;


    // -- The dynamical system --
    LagrangianLinearTIDS* oscillator = new LagrangianLinearTIDS(0, *q0, *v0, *M, *K, *C);
    allDS.insert(oscillator);

    // -- Set external forces (weight) --
    SiconosVector * weight = new SimpleVector(nDof);
    double g = 9.81;
    (*weight)(0) = -m1 * g;
    oscillator->setFExtPtr(weight);

    // -- Set internal forces (stiffness) --
    SiconosVector * Stiff = new SimpleVector(nDof);
    (*Stiff)(0) = -k * L;
    (*Stiff)(1) =  k * L;
    oscillator->setFExtPtr(Stiff);


    // ==> at this point, all the required dynamical systems are saved in allDS.

    // --------------------
    // --- Interactions ---
    // --------------------

    InteractionsSet allInteractions;

    // -- nslaw --
    double mu = 0.1;
    double e  = 0.9;
    double R  = 0.05;

    // Interaction ball-floor
    //
    SiconosMatrix *H = new SimpleMatrix(3, nDof);
    SiconosVector *b = new SimpleVector(3);
    (*H)(0, 0) = 1.;
    (*b)(0)   = -R;

    NonSmoothLaw * nslaw = new NewtonImpactFrictionNSL(e, e, mu, 3);
    Relation * relation = new LagrangianLinearR(*H, *b);

    Interaction * inter = new Interaction("wall", allDS, 0, 3, nslaw, relation);
    allInteractions.insert(inter);

    // Interaction ball-Wall
    //
    SiconosMatrix *H1 = new SimpleMatrix(3, nDof);
    SiconosVector *b1 = new SimpleVector(3);
    (*H1)(0, 1) = -1.;
    (*b1)(0)   = 0.7 - R;

    Relation * relation1 = new LagrangianLinearR(*H1, *b1);
    Interaction * inter1 = new Interaction("wall2", allDS, 1, 3, nslaw, relation1);
    allInteractions.insert(inter1);

    // --------------------------------
    // --- NonSmoothDynamicalSystem ---
    // --------------------------------
    bool isBVP = 0;
    NonSmoothDynamicalSystem * nsds = new NonSmoothDynamicalSystem(allDS, allInteractions, isBVP);

    // -------------
    // --- Model ---
    // -------------

    Model * two_dof_oscillator = new Model(t0, T);
    two_dof_oscillator->setNonSmoothDynamicalSystemPtr(nsds); // set NonSmoothDynamicalSystem of this model

    // ------------------
    // --- Simulation ---
    // ------------------

    // -- Time discretisation --
    TimeDiscretisation * GLOB_T = new TimeDiscretisation(h, two_dof_oscillator);

    TimeStepping* GLOB_SIM = new TimeStepping(GLOB_T);

    // -- OneStepIntegrators --
    OneStepIntegrator * OSI;
    OSI = new Moreau(oscillator, 0.5000001 , GLOB_SIM);

    // -- OneStepNsProblem --
    string solverName = "NSGS";      // solver algorithm used for non-smooth problem
    IntParameters iparam(5);
    iparam[0] = 1010; // Max number of iteration
    iparam[4] = 0; // Solver/formulation  0: projection, 1: Newton/AlartCurnier, 2: Newton/Fischer-Burmeister

    DoubleParameters dparam(5);
    dparam[0] = 1e-3; // Tolerance
    NonSmoothSolver * Mysolver = new NonSmoothSolver(solverName, iparam, dparam);
    FrictionContact* osnspb = new FrictionContact(GLOB_SIM, 3, Mysolver);
    osnspb->setNumericsVerboseMode(0);

    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    // --- Simulation initialization ---

    cout << "====> Simulation initialisation ..." << endl << endl;

    GLOB_SIM->initialize();

    int N = GLOB_T->getNSteps(); // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 6;
    SimpleMatrix dataPlot(N + 1, outputSize);

    SiconosVector * q = oscillator->getQPtr();
    SiconosVector * v = oscillator->getVelocityPtr();
    SiconosVector * lambda = two_dof_oscillator->getNonSmoothDynamicalSystemPtr()->getInteractionPtr(0)->getLambdaPtr(1);

    dataPlot(0, 0) = two_dof_oscillator->getT0();
    dataPlot(0, 1) = (*q)(0);
    dataPlot(0, 2) = (*v)(0);
    dataPlot(0, 3) = (*lambda)(0);
    dataPlot(0, 4) = (*q)(1);
    dataPlot(0, 5) = (*v)(1);

    // --- Time loop ---
    cout << "====> Start computation ... " << endl << endl;
    // ==== Simulation loop - Writing without explicit event handling =====
    int i = 0;
    boost::progress_display show_progress(N);
    for (i = 1 ; i < N + 1 ; ++i)
    {
      GLOB_SIM->computeOneStep();
      // --- Get values to be plotted ---
      dataPlot(i, 0) = GLOB_SIM->getNextTime();
      dataPlot(i, 1) = (*q)(0);
      dataPlot(i, 2) = (*v)(0);
      dataPlot(i, 3) = (*lambda)(0);
      dataPlot(i, 4) = (*q)(1);
      dataPlot(i, 5) = (*v)(1);

      GLOB_SIM->nextStep();
      ++show_progress;
    }
    cout << "End of computation - Number of iterations done: " << i - 1 << endl;

    // --- Output files ---
    cout << "====> Output file writing ..." << endl;
    ioMatrix io("result.dat", "ascii");
    io.write(dataPlot, "noDim");
    // --- Free memory ---
    delete osnspb;
    delete GLOB_T;
    delete GLOB_SIM;
    delete OSI;
    delete two_dof_oscillator;
    delete nsds;
    delete inter;
    delete relation;
    delete nslaw;
    delete H;
    delete b;
    delete Stiff;
    delete oscillator;
    delete q0;
    delete v0;
    delete M;
    delete K;
    delete C;
  }

  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught in \'sample/Two_dof_oscillator\' " << endl;
  }
  cout << "Computation Time " << time.elapsed()  << endl;


}


