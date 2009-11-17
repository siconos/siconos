/* Siconos-Example version 3.0.0, Copyright INRIA 2005-2008.
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
 * Foundation, Inc., 51 Franklin St, Fifth FLOOR, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
 *
 */

#ifdef WithQGLViewer
#include "DrawUtilities.h"
#include <QGLViewer/qglviewer.h>
#endif
#include "BilliardModel.h"
#include "environment.h"
/**
   Billiard Table simulation - 3D frictionl contact problem in presence of a rigid foundation

   The Siconos Model is built and initialized thanks to initSiconos() function.
   All related functions are in src/BilliardModel.cpp

   The Simulation is saved in global variable GLOB_SIM

   The visualization is done thanks to drawstuff tools. For more details see the README file.

   16/05/2007- Authors: Houari Khenous, Roger Pissard, Franck Perignon

   Last modification 29/05/2008, F. Perignon


*/

using namespace std;

BilliardModel::BilliardModel(unsigned int n):
  numberOfFloors(n), numberOfSpheres(SUM(numberOfFloors)), nDof(6), iter_k(1)
{
  allSpheres.resize(numberOfSpheres);
}

BilliardModel::~BilliardModel()
{}

void BilliardModel::initialize()
{

  // initial computation time
  double t0 = 0.0;
  // final computation time
  double T = 2;
  // Default time step
  double h = 0.05;
  // ================= Creation of the model =======================

  // ---------------------------------------
  // --- Build the Dynamical systems set ---
  // ---------------------------------------
  DynamicalSystemsSet allDS; // the list of DS
  buildDynamicalSystems();

  for (DSLIST::iterator it = allSpheres.begin(); it != allSpheres.end(); ++it)
    allDS.insert(*it);

  // -------------------
  // --- Interactions---
  // -------------------
  InteractionsSet allInteractions;
  buildInteractions(allInteractions);

  // -------------
  // --- Model ---
  // -------------

  // initial computation time
  billiard.reset(new Model(t0, T, allDS, allInteractions));

  // ----------------
  // --- Simulation ---
  // ----------------

  // -- Time-discretisation and Simulation --
  SP::TimeDiscretisation t(new TimeDiscretisation(t0, h));
  SP::TimeStepping s(new TimeStepping(t));

  // -- OneStepIntegrators --
  SP::OneStepIntegrator OSI(new Moreau(allDS , 0.5000001));
  s->insertIntegrator(OSI);

  // -- OneStepNsProblem --
  string solverName = "NSGS";      // solver algorithm used for non-smooth problem
  IntParameters iparam(5);
  iparam[0] = 10100; // Max number of iteration
  // Solver/formulation
  // 0: projection, 1: Newton/AlartCurnier, 2: Newton/Fischer-Burmeister, 3: Path/Glocker
  iparam[4] = 1;
  DoubleParameters dparam(5);
  dparam[0] = 1e-6; // Tolerance
  dparam[2] = 1e-8; // Local Tolerance
  SP::NonSmoothSolver Mysolver(new NonSmoothSolver(solverName, iparam, dparam));
  SP::FrictionContact osnspb(new FrictionContact(3, Mysolver));
  s->insertNonSmoothProblem(osnspb);
  //osnspb->setNumericsVerboseMode(1);
  //  osnspb->setMStorageType(1);
  cout << "=== End of model loading === " << endl;

  // =========================== End of model definition

  // --- Simulation initialization ---

  cout << "..." << endl;
  billiard->initialize(s);
  cout << "=== End of simulation initialisation ===" << endl;

#ifndef WithQGLViewer
  unsigned int N = int((T - t0) / h);
  unsigned int outputSize = 1 + 3 * numberOfSpheres;
  // Output matrix
  dataPlot.reset(new SimpleMatrix(N + 1, outputSize));

  (*dataPlot)(0, 0) = t0; // time
  unsigned int i = 0;
  for (DSLIST::iterator it = allSpheres.begin(); it != allSpheres.end(); ++it)
  {
    (*dataPlot)(0, 1 + 3 * i) = (*it)->getQ(2); // q(2)
    (*dataPlot)(0, 2 + 3 * i) = (*it)->getVelocity(2); // v(2)
    (*dataPlot)(0, 3 + 3 * i) = (*it)->getQ(1); // v(2)
    // (*dataPlot)(0,2+2*i) = (*it)->getVelocity(2); // v(2)
    i++;
  }
#endif
}

void BilliardModel::compute()
{
  billiard->getSimulationPtr()->advanceToEvent();
#ifndef WithQGLViewer
  unsigned int i = 0;

  (*dataPlot)(iter_k, 0) =  billiard->getSimulationPtr()->getNextTime(); // time
  for (DSLIST::iterator it = allSpheres.begin(); it != allSpheres.end(); ++it)
  {
    (*dataPlot)(iter_k, 1 + 3 * i) = (*it)->getQ(2); // q(2)
    (*dataPlot)(iter_k, 2 + 3 * i) = (*it)->getVelocity(2); // v(2)
    (*dataPlot)(iter_k, 3 + 3 * i) = (*it)->getQ(1); // v(2)
    //     (*dataPlot)(iter_k,2+2*i) = (*it)->getVelocity(2); // v(2)
    i++;
  }
#endif
  billiard->getSimulationPtr()->processEvents();
  iter_k++;
}

bool BilliardModel::isSimulationFinished()
{
  return !(billiard->getSimulationPtr()->getNextTime() < billiard->getFinalT());
}

void BilliardModel::draw()
{
#ifdef WithQGLViewer

  DrawUtilities::drawHorizontalPlane(Ground);
  //   DrawUtilities::drawYVerticalPlane(Wall);
  //   DrawUtilities::drawYVerticalPlane(-Wall);
  //   DrawUtilities::drawXVerticalPlane(-Wall);
  //   DrawUtilities::drawXVerticalPlane(Wall);

  double x, y, z, r, color = 0.1;

  //  glLoadIdentity();
  for (DSLIST::iterator it = allSpheres.begin(); it != allSpheres.end(); ++it)
  {
    // (*it)->draw();

    x = (*it)->getQ(0);
    y = (*it)->getQ(1);
    z = (*it)->getQ(2);
    r = (*it)->getRadius();
    DrawUtilities::drawSphere(r, x, y, z, color);
    color += 0.2;
  }
#else
  cout << "Warning: can not call draw function for the model. You need to compile with option WithQGLViewer." << endl;
#endif
}

void BilliardModel::computeInitialPositions(Vectors q0, Vectors v0)
{
  // Balls position

  /* o
        o
     o     o
        o     o                     v0
     o     o     o                 <--   o
        o     o
     o     o
        o
     o
  */

  (*(q0[0]))(0) =  0.;
  (*(q0[0]))(1) =  0.;
  (*(q0[0]))(2) =  0.1;
  (*(q0[1]))(0) =  0.1;
  (*(q0[1]))(1) = -0.2;
  (*(q0[1]))(2) =  0.1;
  (*(q0[2]))(0) = -0.1;
  (*(q0[2]))(1) = -0.2;
  (*(q0[2]))(2) =  0.1;
  (*(q0[3]))(0) =  0.2;
  (*(q0[3]))(1) = -0.4;
  (*(q0[3]))(2) =  0.1;
  (*(q0[4]))(0) = -0.2;
  (*(q0[4]))(1) = -0.4;
  (*(q0[4]))(2) =  0.1;
  (*(q0[5]))(0) =  0.;
  (*(q0[5]))(1) = -0.4;
  (*(q0[5]))(2) =  0.1;
  (*(q0[6]))(0) =  0.1;
  (*(q0[6]))(1) = -0.6;
  (*(q0[6]))(2) =  0.1;
  (*(q0[7]))(0) = -0.1;
  (*(q0[7]))(1) = -0.6;
  (*(q0[7]))(2) =  0.1;

  (*(q0[8]))(0) =  0.;
  (*(q0[8]))(1) =  0.8;
  (*(q0[8]))(2) =  0.1;

  (*(v0[8]))(0) = -0.1;
  (*(v0[8]))(1) = -2.;

  (*(q0[9]))(0) =  0.3;
  (*(q0[9]))(1) = -0.6;
  (*(q0[9]))(2) =  0.1;
  (*(q0[10]))(0) = -0.3;
  (*(q0[10]))(1) = -0.6;
  (*(q0[10]))(2) =  0.1;
  (*(q0[11]))(0) =  0.2;
  (*(q0[11]))(1) = -0.8;
  (*(q0[11]))(2) =  0.1;
  (*(q0[12]))(0) = -0.2;
  (*(q0[12]))(1) = -0.8;
  (*(q0[12]))(2) =  0.1;
  (*(q0[13]))(0) =  0.;
  (*(q0[13]))(1) = -0.8;
  (*(q0[13]))(2) =  0.1;
  (*(q0[14]))(0) =  0.4;
  (*(q0[14]))(1) = -0.8;
  (*(q0[14]))(2) =  0.1;
  (*(q0[15]))(0) = -0.4;
  (*(q0[15]))(1) = -0.8;
  (*(q0[15]))(2) =  0.1;
}

void BilliardModel::buildDynamicalSystems()
{
  // Set the same radius and mass for all balls
  double Radius = DEFAULT_radius;
  double m = 1.0;

  // -- Initial positions and velocities --
  // q0[i] and v0[i] correspond to position and velocity of ball i.

  Vectors q0, v0;
  q0.resize(numberOfSpheres);
  v0.resize(numberOfSpheres);
  // Memory allocation for q0[i] and v0[i]
  for (unsigned int i = 0; i < numberOfSpheres; i++)
  {
    q0[i].reset(new SimpleVector(nDof));
    v0[i].reset(new SimpleVector(nDof));
  }

  // Computation of the position of all beads of the pyramid
  computeInitialPositions(q0, v0);

  // Build and insert the DS into allDS

  for (unsigned int i = 0; i < numberOfSpheres; i++)
    allSpheres[i].reset(new Sphere(Radius, m, *(q0[i]), *(v0[i])));

}


void BilliardModel::buildInteractions(InteractionsSet& allInteractions)
{
  unsigned int Fact = NBSpheres * (NBSpheres - 1) / 2;
  double e = 0.9;                  // nslaw
  double e2 = 0.9;                  // nslaw2
  double mu = 0.;
  DynamicalSystemsSet dsConcernedi, dsConcerned2 ;
  vector<SP::Relation> LLR(Fact);
  vector<SP::Relation> LLR1(NBSpheres);
  vector<SP::Relation> LLR1_(NBSpheres);
  vector<SP::Relation> LLR2(NBSpheres);
  vector<SP::Relation> LLR2_(NBSpheres);
  vector<SP::Relation> LLR3(NBSpheres);
  vector<SP::Relation> LLR3_(NBSpheres);

  SP::NonSmoothLaw nslaw1(new NewtonImpactFrictionNSL(e, e, mu, 3));

  // Interaction beads and plan1 (OXY)

  double R = DEFAULT_radius;

  SimpleVector b1(3);
  b1(0) = Ground - R;
  SimpleMatrix H1(3, nDof);
  H1(0, 2) = 1.0;
  H1(1, 0) = 1.0;
  H1(1, 4) = -R;
  H1(2, 1) = 1.0;
  H1(2, 3) =  R;

  SP::Interaction inter;

  for (int i = 0; i < NBSpheres; i++)
  {
    dsConcernedi.insert(allSpheres[i]);
    LLR1[i].reset(new LagrangianLinearTIR(H1, b1));
    inter.reset(new Interaction(dsConcernedi, i, 3, nslaw1, LLR1[i]));
    allInteractions.insert(inter);
    dsConcernedi.clear();
  }

  // Interaction beads and plan1 (-YOX)

  SimpleVector b1_(3);
  b1_(0) = Top - R;
  SimpleMatrix H1_(3, nDof);
  H1_(0, 2) = -1.0;
  H1_(1, 0) = 1.0;
  H1_(1, 4) = -R;
  H1_(2, 1) = 1.0;
  H1_(2, 3) =  R;

  for (int i = 0; i < NBSpheres; i++)
  {
    dsConcernedi.insert(allSpheres[i]);
    LLR1_[i].reset(new LagrangianLinearTIR(H1_, b1_));
    inter.reset(new Interaction(dsConcernedi, i, 3, nslaw1, LLR1_[i]));
    allInteractions.insert(inter);
    dsConcernedi.clear();
  }

  // Interaction beads and plan2 (OXZ)

  SimpleVector b2(3);
  b2(0) = Wall - R;
  SimpleMatrix H2(3, nDof);
  H2(0, 1) = 1.0;
  H2(1, 0) = 1.0;
  H2(1, 5) = -R;
  H2(2, 2) = 1.0;
  H2(2, 3) =  R;

  for (int i = 0; i < NBSpheres; i++)
  {
    dsConcernedi.insert(allSpheres[i]);
    LLR2[i].reset(new LagrangianLinearTIR(H2, b2));
    inter.reset(new Interaction(dsConcernedi, i, 3, nslaw1, LLR2[i]));
    allInteractions.insert(inter);
    dsConcernedi.clear();
  }

  // Interaction beads and plan2 (-ZOX)

  SimpleVector b2_(3);
  b2_(0) = Wall - R;
  SimpleMatrix H2_(3, nDof);
  H2_(0, 1) = -1.0;
  H2_(1, 0) = 1.0;
  H2_(1, 5) = -R;
  H2_(2, 2) = 1.0;
  H2_(2, 3) =  R;

  for (int i = 0; i < NBSpheres; i++)
  {
    dsConcernedi.insert(allSpheres[i]);
    LLR2_[i].reset(new LagrangianLinearTIR(H2_, b2_));
    inter.reset(new Interaction(dsConcernedi, i, 3, nslaw1, LLR2_[i]));
    allInteractions.insert(inter);
    dsConcernedi.clear();
  }

  // Interaction beads and plan3 (OYZ)

  SimpleVector b3(3);
  b3(0) = Wall - R;
  SimpleMatrix H3(3, nDof);
  H3(0, 0) = 1.0;
  H3(1, 1) = 1.0;
  H3(1, 5) = -R;
  H3(2, 2) = 1.0;
  H3(2, 4) =  R;

  for (int i = 0; i < NBSpheres; i++)
  {
    dsConcernedi.insert(allSpheres[i]);
    LLR3[i].reset(new LagrangianLinearTIR(H3, b3));
    inter.reset(new Interaction(dsConcernedi, i, 3, nslaw1, LLR3[i]));
    allInteractions.insert(inter);
    dsConcernedi.clear();
  }

  // Interaction beads and plan3 (-ZOY)

  SimpleVector b3_(3);
  b3_(0) = Wall - R;
  SimpleMatrix H3_(3, nDof);
  H3_(0, 0) = -1.0;
  H3_(1, 1) = 1.0;
  H3_(1, 5) = -R;
  H3_(2, 2) = 1.0;
  H3_(2, 4) =  R;

  for (int i = 0; i < NBSpheres; i++)
  {
    dsConcernedi.insert(allSpheres[i]);
    LLR3_[i].reset(new LagrangianLinearTIR(H3_, b3_));
    inter.reset(new Interaction(dsConcernedi, i, 3, nslaw1, LLR3_[i]));
    allInteractions.insert(inter);
    dsConcernedi.clear();
  }

  // Interaction between beads

  SP::NonSmoothLaw nslaw2(new NewtonImpactFrictionNSL(e2, e2, mu, 3));

  int l = 0;
  for (int i = 0; i < NBSpheres; i++)
  {
    dsConcerned2.insert(allSpheres[i]);
    for (int j = 0; j < NBSpheres; j++)
    {
      if (j > i)
      {
        dsConcerned2.insert(allSpheres[j]);
        LLR[l].reset(new LagrangianScleronomousR("3DDrawPlugin:h0", "3DDrawPlugin:G0"));
        inter.reset(new Interaction(dsConcerned2, l, 3, nslaw2, LLR[l]));
        allInteractions.insert(inter);
        dsConcerned2.erase(allSpheres[j]);
        l = l + 1;
      }
    }
    dsConcerned2.clear();
  }

}

void BilliardModel::end()
{
  cout << "End of computation - Number of iterations done: " << iter_k << endl;
#ifndef WithQGLViewer
  // --- Output files ---
  ioMatrix io("result.dat", "ascii");
  io.write(*dataPlot, "noDim");
#endif
}
