/* Siconos-sample version 2.0.0, Copyright INRIA 2005-2006.
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
 * Two beads 3D frictionl contact problem in presence of a rigid foundation
 * 27/06/2007- Authors: houari khenous

*/
// =============================== Two dof oscillator simulation ===============================
//
// ======================================================================================================

#include "SiconosKernel.h"
#include <sys/time.h>
#include <math.h>
#include <iostream>
#include <fstream>

using namespace std;

#include <drawstuff/drawstuff.h>

#define WALL 1       // Positions of walls

Simulation * GLOB_SIM;
TimeDiscretisation * GLOB_T;
LagrangianLinearTIDS * oscillator;
Model * two_dof_oscillator;
SiconosVector * GLOB_POS;
SiconosVector * GLOB_VELO;
SiconosVector * GLOB_LAMBDA;



int GLOB_COMPUTE;
int GLOB_STEP;
EventsManager * GLOB_EVT;

// Global variables for computation of CPU time, number of iterations and for curves ploting

#define PLOTMAX 2000
unsigned int outputSize = 4;
SimpleMatrix dataPlot(PLOTMAX + 1, outputSize);
int i = 0; // index for output.

extern void initSiconos();
extern void computeSiconos();


void Start()
{
  GLOB_COMPUTE = false;
  GLOB_STEP = false;
  initSiconos();
}

void Draw_two_dof_oscillator(LagrangianLinearTIDS *lds, float radius)
{
  float R[12];
  float pos1[3];
  float pos2[3];


  // Translation
  pos1[0] = 0.;
  pos1[1] = lds->getQ()(0);
  pos1[2] = radius;
  pos2[0] = 0.;
  pos2[1] = lds->getQ()(1);
  pos2[2] = radius;


  // Rotation
  R[0] = 1;
  R[1] = 0;
  R[2] = 0;
  R[4] = 0;
  R[5] = 1;
  R[6] = 0;
  R[8] = 0;
  R[9] = 0;
  R[10] = 1;

  R[3] = R[7] = R[11] = 0;

  dsSetTexture(DS_WOOD);
  dsSetColor(1, 0.8f, 0.6f);
  dsDrawSphere(pos1, R, radius);
  dsSetTexture(DS_NONE);
  dsSetColor(0.6f, 0.6f, 1);
  dsDrawSphere(pos2, R, radius);
}

void SimuLoop(int pause)
{
  float radius;

  if (GLOB_COMPUTE == true)
    computeSiconos();

  if (GLOB_STEP == true)
    GLOB_COMPUTE = false;

  // Ball Radius
  radius = 0.05;
  Draw_two_dof_oscillator(oscillator, radius);
}

void Command(int cmd)
{
  dsPrint("received command %d (`%c')\n", cmd, cmd);

  switch (cmd)
  {
  case 's':
    cout << "coucou" << "\n";
    break;
  case 'r':
    cout << "-- Run Simu" << endl;
    GLOB_COMPUTE = true;
    break;
  case 'p':
    cout << "-- Pause Simu" << endl;
    GLOB_COMPUTE = false;
    break;
  case 'f':
    cout << "-- Step mode" << endl;
    GLOB_STEP = !GLOB_STEP;
    break;
  default:
    cout << "-- Press " << endl;
    cout << " s - display DS status " << endl;
    cout << " r - run simu  " << endl;
    cout << " p - pause simu " << endl;
    cout << " f - step mode (toggle) " << endl;
  }
}

int main(int argc, char* argv[])
{

  // setup pointers to callback functions
  dsFunctions fn;
  fn.version = DS_VERSION;
  fn.start = &Start;
  fn.step = &SimuLoop;
  fn.command = Command;
  fn.stop = 0;
  fn.path_to_textures = "./textures/";
  // run simulation
  dsSimulationLoop(argc, argv, 400, 400, &fn);
  return 0;
}


void initSiconos()
{
  try
  {

    // ================= Creation of the model =======================

    // User-defined main parameters
    unsigned int nDof = 2;           // degrees of freedom for the ball
    double t0 = 0;                   // initial computation time
    double T = 10.0;                 // final computation time
    double h = 0.005;                 // time step

    double pos1 = 4.0;               // initial position for m1.
    double v1   = 0.0;               // initial velocity for m1
    double pos2 = 4.8;               // initial position for m2
    double v2   = 0.0;               // initial velocity for m2


    //string solverName = "NLGSNEWTON";      // solver algorithm used for non-smooth problem
    string solverName = "NLGS";      // solver algorithm used for non-smooth problem
    // string solverName = "Lemke" ;

    double k = 0.5; // stiffness coefficient
    double L = 0.5; // initial lenth
    double m1 = 1.; //  m1
    double m2 = 1.; //  m2
    //double g = 9.81; // Gravity

    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    cout << "====> Model loading ..." << endl << endl;
    DynamicalSystemsSet allDS; // the list of DS
    SiconosMatrix *M = new SimpleMatrix(nDof, nDof);
    (*M)(0, 0) = m1;
    (*M)(1, 1) = m2;
    SiconosMatrix *K = new SimpleMatrix(nDof, nDof);
    (*K)(0, 0) =  k;
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
    oscillator = new LagrangianLinearTIDS(0, *q0, *v0, *M, *K, *C);
    allDS.insert(oscillator);

    // // -- Set external forces (weight) --
    //     SiconosVector * weight = new SimpleVector(nDof);
    //     (*weight)(0) = -m*g;
    //     oscillator->setFExtPtr(weight);

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
    double e = 0.9;

    // Interaction ball-floor
    //
    SiconosMatrix *H = new SimpleMatrix(1, nDof);
    (*H)(0, 0) = 1.0;
    NonSmoothLaw * nslaw0 = new NewtonImpactNSL(e);
    Relation * relation0 = new LagrangianLinearR(*H);

    Interaction * inter = new Interaction("wall", allDS, 0, 1, nslaw0, relation0);
    allInteractions.insert(inter);

    // --------------------------------
    // --- NonSmoothDynamicalSystem ---
    // --------------------------------
    bool isBVP = 0;
    NonSmoothDynamicalSystem * nsds = new NonSmoothDynamicalSystem(allDS, allInteractions, isBVP);

    // -------------
    // --- Model ---
    // -------------

    two_dof_oscillator = new Model(t0, T);
    two_dof_oscillator->setNonSmoothDynamicalSystemPtr(nsds); // set NonSmoothDynamicalSystem of this model

    // ----------------
    // --- Simulation ---
    // ----------------

    GLOB_T = new TimeDiscretisation(h, two_dof_oscillator);
    GLOB_SIM = new TimeStepping(GLOB_T);

    // -- OneStepIntegrators --
    OneStepIntegrator * OSI = new Moreau(allDS , 0.5000001 , GLOB_SIM);

    // -- OneStepNsProblem --


    //OneStepNSProblem * osnspb = new LCP(GLOB_SIM ,"FrictionContact3D",solverName,101, 0.001);

    OneStepNSProblem * osnspb = new FrictionContact3D(GLOB_SIM , "FrictionContact3D", solverName, 1000001, 0.001);

    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    // --- Simulation initialization ---

    cout << "====> Simulation initialisation ..." << endl << endl;

    GLOB_SIM->initialize();

    GLOB_POS = oscillator->getQPtr();
    GLOB_VELO = oscillator->getVelocityPtr();
    GLOB_LAMBDA = two_dof_oscillator->getNonSmoothDynamicalSystemPtr()->getInteractionPtr(0)->getLambdaPtr(1);

    dataPlot(0, 0) = two_dof_oscillator->getT0();
    dataPlot(0, 1) = (*GLOB_POS)(0);
    dataPlot(0, 2) = (*GLOB_VELO)(0);
    dataPlot(0, 3) = (*GLOB_LAMBDA)(0);

    // ================================= Computation =================================

    // --- Time loop ---
    cout << "====> Start computation ... " << endl << endl;

    GLOB_EVT = GLOB_SIM->getEventsManagerPtr();
  }
  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught in \'sample/MultiBeadsColumn Init\'" << endl;
  }
}

void computeSiconos()
{
  try
  {
    // --- simulation solver ---
    if (GLOB_EVT->hasNextEvent())
    {
      GLOB_SIM->advanceToEvent();
      GLOB_SIM->processEvents();
      // --- Get values to be plotted ---
      i++;
      dataPlot(i, 0) =  two_dof_oscillator->getCurrentT();
      dataPlot(i, 1) = (*GLOB_POS)(0);
      dataPlot(i, 2) = (*GLOB_VELO)(0);
      dataPlot(i, 3) = (*GLOB_LAMBDA)(0);

      cout << "x1 =  " << (*GLOB_POS)(0) << endl;
      cout << "x2 =  " << (*GLOB_POS)(1) << endl;
    }

    cout << "End of computation - Number of iterations done: " << i << endl;
    // --- Output files ---
    cout << "====> Output file writing ..." << endl;
    ioMatrix io("result.dat", "ascii");
    io.write(dataPlot, "noDim");
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


