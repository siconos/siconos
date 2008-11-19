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
#include "CanonBallsModel.h"
#include "environment.h"
/**
   Pyramidal arrangement of Canon Balls simulation - 3D frictionl contact problem in presence of a rigid foundation

   The Siconos Model is built and initialized thanks to initSiconos() function.
   All related functions are in src/CanonModel.cpp

   The Simulation is saved in global variable GLOB_SIM

   The visualization is done thanks to drawstuff tools. For more details see the README file.

   16/05/2007- \authors: Houari Khenous, Roger Pissard, Franck Perignon

   Last modification 29/05/2008, F. Perignon


*/

using namespace std;

CanonBallsModel::CanonBallsModel(unsigned int n):
  numberOfFloors(n), numberOfSpheres(SUM(numberOfFloors)), nDof(6), iter_k(1)
{
  allSpheres.resize(numberOfSpheres);
}

CanonBallsModel::~CanonBallsModel()
{}

void CanonBallsModel::initialize()
{

  // initial computation time
  double t0 = 0.0;
  // final computation time
  double T = 1.0;
  // Default time step
  double h = 0.005;
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
  canonballs.reset(new Model(t0, T, allDS, allInteractions));

  // ----------------
  // --- Simulation ---
  // ----------------

  // -- Time-discretisation and Simulation --
  SP::TimeDiscretisation t(new TimeDiscretisation(t0, h));
  SP::TimeStepping s(new TimeStepping(t));

  // -- OneStepIntegrators --
  SP::OneStepIntegrator OSI(new Moreau(allDS , 0.5000001));
  s->recordIntegrator(OSI);

  // -- OneStepNsProblem --
  string solverName = "NSGS";      // solver algorithm used for non-smooth problem
  IntParameters iparam(5);
  iparam[0] = 10100; // Max number of iteration
  // Solver/formulation
  // 0: projection, 1: Newton/AlartCurnier, 2: Newton/Fischer-Burmeister, 3: Path/Glocker
  iparam[4] = 4;
  DoubleParameters dparam(5);
  dparam[0] = 1e-6; // Tolerance
  dparam[2] = 1e-8; // Local Tolerance
  SP::NonSmoothSolver Mysolver(new NonSmoothSolver(solverName, iparam, dparam));
  SP::FrictionContact osnspb(new FrictionContact(3, Mysolver));
  s->recordNonSmoothProblem(osnspb);
  //osnspb->setNumericsVerboseMode(1);
  //  osnspb->setMStorageType(1);
  cout << "=== End of model loading === " << endl;

  // =========================== End of model definition

  // --- Simulation initialization ---

  cout << "..." << endl;
  canonballs->initialize(s);
  cout << "=== End of model initialisation ===" << endl;

#ifndef WithQGLViewer
  unsigned int N = int((T - t0) / h);
  unsigned int outputSize = 1 + 3 * numberOfSpheres;
  // Output matrix
  dataPlot.reset(new SimpleMatrix(N + 2, outputSize));

  (*dataPlot)(0, 0) = t0; // time
  unsigned int i = 0;
  for (DSLIST::iterator it = allSpheres.begin(); it != allSpheres.end(); ++it)
  {
    (*dataPlot)(0, 1 + 3 * i) = (*it)->getQ(0); // q(2)
    (*dataPlot)(0, 2 + 3 * i) = (*it)->getQ(1); // v(2)
    (*dataPlot)(0, 3 + 3 * i) = (*it)->getQ(2); // v(2)
    // (*dataPlot)(0,2+2*i) = (*it)->getVelocity(2); // v(2)
    i++;
  }
#endif
}

void CanonBallsModel::compute()
{
  canonballs->getSimulationPtr()->advanceToEvent();
#ifndef WithQGLViewer
  unsigned int i = 0;

  (*dataPlot)(iter_k, 0) =  canonballs->getSimulationPtr()->getNextTime(); // time
  for (DSLIST::iterator it = allSpheres.begin(); it != allSpheres.end(); ++it)
  {
    (*dataPlot)(iter_k, 1 + 3 * i) = (*it)->getQ(0); // q(2)
    (*dataPlot)(iter_k, 2 + 3 * i) = (*it)->getQ(1); // v(2)
    (*dataPlot)(iter_k, 3 + 3 * i) = (*it)->getQ(2); // v(2)
    //     (*dataPlot)(iter_k,2+2*i) = (*it)->getVelocity(2); // v(2)
    i++;
  }
#endif
  canonballs->getSimulationPtr()->processEvents();
  iter_k++;
}

bool CanonBallsModel::isSimulationFinished()
{
  return !(canonballs->getSimulationPtr()->getNextTime() < canonballs->getFinalT());
}

void CanonBallsModel::draw()
{
#ifdef WithQGLViewer

  //  DrawUtilities::drawHorizontalPlane(GROUND);
  //   DrawUtilities::drawYVerticalPlane(WALL);
  //   DrawUtilities::drawYVerticalPlane(-WALL);
  //   DrawUtilities::drawXVerticalPlane(-WALL);
  //   DrawUtilities::drawXVerticalPlane(WALL);

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

void CanonBallsModel::computeInitialPositions(Vectors q0, Vectors v0, double Radius)
{

  double *X = (double*)malloc(3 * sizeof(double));
  unsigned int K;
  // Initial position of the first bead
  double RR = Radius * 1.00001; //  Radius + x% for "almost" contact
  (*(q0[0]))(0) =  0.;
  (*(q0[0]))(1) =  0.;
  (*(q0[0]))(2) = (numberOfFloors - 1) * 2.0 * sqrt(2.0 / 3.0) * RR + GROUND + RR;
  // Loop over all Floors
  for (unsigned int i = 2 ; i < numberOfFloors + 1 ; ++i)
  {
    for (unsigned int j = SUM(i - 2) ; j < SUM(i - 1) ; ++j)
    {
      for (unsigned int l = 0 ; l < 3 ; l++)
      {
        // Compute current ball position according to those of the floor above
        Q(i, j,  l , X, (*(q0[j]))(0), (*(q0[j]))(1), (*(q0[j]))(2), RR);

        if (l == 0) K = 0;
        else
          K = qlq(i, j);

        for (unsigned int k = 0 ; k < 3 ; k++)
          (*(q0[j + sum(i) - l - K]))(k) = X[k];
      }
    }
  }
  free(X);
}

void CanonBallsModel::buildDynamicalSystems()
{
  // Set the same radius and mass for all balls
  double Radius = 0.1;
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
  computeInitialPositions(q0, v0, Radius);

  // Build and insert the DS into allDS

  for (unsigned int i = 0; i < numberOfSpheres; i++)
    allSpheres[i].reset(new Sphere(Radius, m, *(q0[i]), *(v0[i])));
}


void CanonBallsModel::buildInteractions(InteractionsSet& allInteractions)
{
  // Definition of some obstacles

  // A ceiling (z=TOP)
  bool hasCeil = false;
  // A ground (z=Ground)
  bool hasGround = true;
  // Some walls
  bool obst_y_p = true;                    //  for y --> +
  bool obst_y_m = true;                    //  for y --> -
  bool obst_x_p = true;                    //  for x --> +
  bool obst_x_m = true;                    //  for x --> -

  int Fact = (numberOfSpheres) * (numberOfSpheres - 1) / 2;
  vector<SP::Relation> LLR(Fact);
  vector<SP::Relation> LLR1(numberOfSpheres);
  vector<SP::Relation> LLR1_(numberOfSpheres);
  vector<SP::Relation> LLR2(numberOfSpheres);
  vector<SP::Relation> LLR2_(numberOfSpheres);
  vector<SP::Relation> LLR3(numberOfSpheres);
  vector<SP::Relation> LLR3_(numberOfSpheres);

  double Radius = 0.1;
  double e  = 0.1;
  double mu = 0.4;
  SP::NonSmoothLaw nslaw1(new NewtonImpactFrictionNSL(e, e, mu, 3));

  SimpleVector bground(3), bceil(3), bwallYp(3), bwallYm(3), bwallXp(3), bwallXm(3);
  SimpleMatrix Hground(3, nDof), Hceil(3, nDof), HwallYp(3, nDof), HwallYm(3, nDof), HwallXp(3, nDof), HwallXm(3, nDof) ;
  if (hasGround)
  {
    bground(0) =  - GROUND - Radius;
    Hground(0, 2) = 1.0;
    Hground(1, 0) = 1.0;
    Hground(1, 4) = -Radius;
    Hground(2, 1) = 1.0;
    Hground(2, 3) =  Radius;
  }
  // Interaction beads and ceiling
  if (hasCeil)
  {
    bceil(0) = TOP - Radius;
    Hceil(0, 2) = -1.0;
    Hceil(1, 0) = 1.0;
    Hceil(1, 4) = -Radius;
    Hceil(2, 1) = 1.0;
    Hceil(2, 3) =  Radius;
  }
  if (obst_y_p)
  {
    bwallYp(0) = WALL - Radius;
    HwallYp(0, 1) = 1.0;
    HwallYp(1, 0) = 1.0;
    HwallYp(1, 5) = -Radius;
    HwallYp(2, 2) = 1.0;
    HwallYp(2, 3) =  Radius;
  }
  if (obst_y_m)
  {
    bwallYm(0) = WALL - Radius;
    HwallYm(0, 1) = -1.0;
    HwallYm(1, 0) = 1.0;
    HwallYm(1, 5) = -Radius;
    HwallYm(2, 2) = 1.0;
    HwallYm(2, 3) =  Radius;
  }
  if (obst_x_p)
  {
    bwallXp(0) = WALL - Radius;
    HwallXp(0, 0) = 1.0;
    HwallXp(1, 1) = 1.0;
    HwallXp(1, 5) = -Radius;
    HwallXp(2, 2) = 1.0;
    HwallXp(2, 4) =  Radius;
  }
  if (obst_x_m)
  {
    bwallXm(0) = WALL - Radius;
    HwallXm(0, 0) = -1.0;
    HwallXm(1, 1) = 1.0;
    HwallXm(1, 5) = -Radius;
    HwallXm(2, 2) = 1.0;
    HwallXm(2, 4) =  Radius;
  }


  // All Dynamical Systems can interact with ground/ceil/walls
  //  for (unsigned int i=0;i<numberOfSpheres;i++)
  //    {
  // Interaction beads and ground (z=0)
  SP::Interaction inter;
  for (unsigned int i = 0; i < numberOfSpheres; i++)
  {
    if (hasGround)
    {
      LLR1[i].reset(new LagrangianLinearTIR(Hground, bground));
      inter.reset(new Interaction(allSpheres[i], i, 3, nslaw1, LLR1[i]));
      allInteractions.insert(inter);
    }
    // Interaction beads and ceiling
    if (hasCeil)
    {
      LLR1_[i].reset(new LagrangianLinearTIR(Hceil, bceil));
      inter.reset(new Interaction(allSpheres[i], i, 3, nslaw1, LLR1_[i]));
      allInteractions.insert(inter);
    }
    // Interaction beads and plan2 (OXZ)
    if (obst_y_p)
    {
      LLR2[i].reset(new LagrangianLinearTIR(HwallYp, bwallYp));
      inter.reset(new Interaction(allSpheres[i], i, 3, nslaw1, LLR2[i]));
      allInteractions.insert(inter);

    }
    // Interaction beads and plan2 (-ZOX)
    if (obst_y_m)
    {
      LLR2_[i].reset(new LagrangianLinearTIR(HwallYm, bwallYm));
      inter.reset(new Interaction(allSpheres[i], i, 3, nslaw1, LLR2_[i]));
      allInteractions.insert(inter);
    }
    // Interaction beads and plan3 (OYZ)
    if (obst_x_p)
    {
      LLR3[i].reset(new LagrangianLinearTIR(HwallXp, bwallXp));
      inter.reset(new Interaction(allSpheres[i], i, 3, nslaw1, LLR3[i]));
      allInteractions.insert(inter);
    }
    // Interaction beads and plan3 (-ZOY)
    if (obst_x_m)
    {
      LLR3_[i].reset(new LagrangianLinearTIR(HwallXm, bwallXm));
      inter.reset(new Interaction(allSpheres[i], i, 3, nslaw1, LLR3_[i]));
      allInteractions.insert(inter);
    }
  }

  // Interaction between beads

  // frictional contact condition between beads
  double e2 = 0.9;
  SP::NonSmoothLaw nslaw2(new NewtonImpactFrictionNSL(e2, e2, mu, 3));
  unsigned int l = 0;
  DynamicalSystemsSet dsConcerned;
  for (unsigned int i = 0; i < numberOfSpheres; i++)
  {
    dsConcerned.insert(allSpheres[i]);
    for (unsigned int j = 0; j < numberOfSpheres; j++)
    {
      if (j > i)
      {
        dsConcerned.insert(allSpheres[j]);
        LLR[l].reset(new LagrangianScleronomousR("CanonPlugin:h0", "CanonPlugin:G0"));
        inter.reset(new Interaction(dsConcerned, l, 3, nslaw2, LLR[l]));
        allInteractions.insert(inter);
        dsConcerned.erase(allSpheres[j]);
        l++;
      }
    }
    dsConcerned.clear();
  }

}

void CanonBallsModel::end()
{
  cout << "End of computation - Number of iterations done: " << iter_k << endl;
#ifndef WithQGLViewer
  // --- Output files ---
  ioMatrix io("result.dat", "ascii");
  io.write(*dataPlot, "noDim");
#endif
}
