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

#include "BallsModel.h"
#include "environment.h"
#ifdef WithQGLViewer
#include "DrawUtilities.h"
#include <QGLViewer/qglviewer.h>
#endif

/**

 */
using namespace std;

BallsModel::BallsModel(unsigned int n):
  numberOfSpheres(n), nDof(6), iter_k(1)
{
  allSpheres.resize(numberOfSpheres);
}

BallsModel::~BallsModel()
{}

void BallsModel::initialize()
{
  // initial computation time
  double t0 = 0.0;
  // final computation time
  double T = 3.;
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

  // --------------------------------
  // --- NonSmoothDynamicalSystem ---
  // --------------------------------
  SP::NonSmoothDynamicalSystem nsds(new NonSmoothDynamicalSystem(allDS, allInteractions));

  // -------------
  // --- Model ---
  // -------------

  // initial computation time
  balls.reset(new Model(t0, T));
  balls->setNonSmoothDynamicalSystemPtr(nsds); // set NonSmoothDynamicalSystem of this model

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
  string solverName = "Lemke";      // solver algorithm used for non-smooth problem
  IntParameters iparam(5);
  iparam[0] = 1000; // Max number of iteration
  DoubleParameters dparam(5);
  dparam[0] = 1e-6; // Tolerance
  //dparam[2] = 1e-8; // Local Tolerance
  SP::NonSmoothSolver Mysolver(new NonSmoothSolver(solverName, iparam, dparam));
  SP::LCP osnspb(new LCP(Mysolver));
  s->recordNonSmoothProblem(osnspb);

  //osnspb->setNumericsVerboseMode(1);
  //  osnspb->setMStorageType(1);
  cout << "=== End of model loading === " << endl;

  // =========================== End of model definition

  // --- Simulation initialization ---

  cout << "..." << endl;
  balls->initialize(s);
  cout << "=== End of simulation initialisation ===" << endl;

  unsigned int N = int((T - t0) / h);
  unsigned int outputSize = 1 + 2 * numberOfSpheres;

#ifndef WithQGLViewer
  // Output matrix
  dataPlot.reset(new SimpleMatrix(N + 2, outputSize)); // Output data matrix
  (*dataPlot)(0, 0) = t0; // time
  unsigned int i = 0;
  for (DSLIST::iterator it = allSpheres.begin(); it != allSpheres.end(); ++it)
  {
    (*dataPlot)(0, 1 + 2 * i) = (*it)->getQ(2);
    (*dataPlot)(0, 2 + 2 * i) = (*it)->getVelocity(2);
    if ((*it)->number() == 9)
      break;
    i++;
  }
#endif
}

void BallsModel::compute()
{
  balls->simulation()->advanceToEvent();
  unsigned int i = 0;

#ifndef WithQGLViewer
  (*dataPlot)(iter_k, 0) =  balls->simulation()->getNextTime(); // time
  for (DSLIST::iterator it = allSpheres.begin(); it != allSpheres.end(); ++it)
  {
    (*dataPlot)(iter_k, 1 + 2 * i) = (*it)->getQ(2);
    (*dataPlot)(iter_k, 2 + 2 * i) = (*it)->getVelocity(2);
    if ((*it)->number() == 9)
      break;
    i++;
  }
#endif
  balls->simulation()->processEvents();
  iter_k++;
}

bool BallsModel::isSimulationFinished()
{
  return !(balls->simulation()->getNextTime() < balls->getFinalT());
}


void BallsModel::draw()
{
#ifdef WithQGLViewer


  DrawUtilities::drawHorizontalPlane(Ground);

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

void BallsModel::buildDynamicalSystems()
{
  unsigned int i;

  // Set the same radius and mass for all balls
  double Radius = DEFAULT_radius;
  double m = 1;   // beads mass

  // -- Initial positions and velocities --
  // q0[i] and v0[i] correspond to position and velocity of ball i.
  Vectors q0, v0;
  numberOfSpheres = NB_Spheres;
  q0.resize(numberOfSpheres);
  v0.resize(numberOfSpheres);

  double increment_position = 3;   // initial position increment from one DS to the following
  double increment_velocity = 0;   // initial velocity increment from one DS to the following
  double position_init = 4.0;     // initial position for lowest bead.
  double velocity_init = 0.0;      // initial velocity for lowest bead.

  for (i = 0; i < numberOfSpheres; i++)
  {
    // Memory allocation for q0[i] and v0[i]
    q0[i].reset(new SimpleVector(nDof));
    v0[i].reset(new SimpleVector(nDof));
    // set values
    (*(q0[i]))(2) = position_init;
    (*(v0[i]))(2) = velocity_init;
    // Create and insert in allDS a new Lagrangian Linear Dynamical System ...
    allSpheres[i].reset(new Sphere(Radius, m, *(q0[i]), *(v0[i])));
    position_init += increment_position;
    velocity_init += increment_velocity;
  }
}

void BallsModel::buildInteractions(InteractionsSet& allInteractions)
{
  // The total number of Interactions
  int interactionNumber = NB_Spheres;

  // Interaction first bead and floor
  // A set for the systems handles by the "current" Interaction
  DynamicalSystemsSet dsConcerned;
  // Only the "bottom" bead is concerned by this first Interaction,
  // therefore DynamicalSystem number 0.
  dsConcerned.insert(allSpheres[0]);

  // -- Newton impact law --
  double e = 0.9;

  // Lagrangian Relation
  unsigned int interactionSize = 1; // y vector size
  SP::SiconosMatrix H(new SimpleMatrix(interactionSize, nDof));
  (*H)(0, 2) = 1.0;
  SP::NonSmoothLaw nslaw0(new NewtonImpactNSL(e));
  SP::SiconosVector b(new SimpleVector(interactionSize));
  double R = DEFAULT_radius;
  (*b)(0) = -R - Ground;
  SP::Relation relation0(new LagrangianLinearTIR(*H, *b));
  unsigned int num = 0 ; // an id number for the Interaction
  SP::Interaction inter0(new Interaction(dsConcerned, num, interactionSize, nslaw0, relation0));
  allInteractions.insert(inter0);

  // Interactions ball-ball
  CheckInsertInteraction checkInter;
  // A vector that handles all the relations
  vector<SP::Relation> LLR(interactionNumber - 1);
  SP::SiconosMatrix H1(new SimpleMatrix(1, 2 * nDof));
  (*b)(0) = -2 * R;
  if (numberOfSpheres > 1)
  {
    (*H1)(0, 2) = -1.0;
    (*H1)(0, 8) = 1.0;
    for (int i = 1; i < interactionNumber; i++)
    {
      // The systems handled by the current Interaction ...
      dsConcerned.clear();
      dsConcerned.insert(allSpheres[i - 1]);
      dsConcerned.insert(allSpheres[i]);
      // The relations
      // Since Ri=Rj and h=0, we do not need to set b.
      LLR[i - 1].reset(new LagrangianLinearTIR(*H1, *b));
      SP::Interaction inter(new Interaction(dsConcerned, i, interactionSize, nslaw0, LLR[i - 1]));
      checkInter = allInteractions.insert(inter);
    }
  }
}

void BallsModel::end()
{
  cout << "End of computation - Number of iterations done: " << iter_k << endl;
#ifndef WithQGLViewer
  // --- Output files ---
  ioMatrix io("result.dat", "ascii");
  io.write(*dataPlot, "noDim");
#endif
}
