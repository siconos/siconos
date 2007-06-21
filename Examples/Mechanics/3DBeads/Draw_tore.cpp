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
 * Torus 3D frictionl contact problem in presence of a rigid foundations
 * 21/06/2007- Authors: houari khenous

*/
// =============================== Torus simulation ===============================
//
// Keywords: LagrangianLinearTIDS relation, Moreau TimeStepping, NLGS, NLGSNEWTON.
//
// ======================================================================================================

#include "SiconosKernel.h"
#include <sys/time.h>
#include <math.h>
#include <iostream>
#include <fstream>

using namespace std;

#include <drawstuff/drawstuff.h>

#define FEM  324      // the number of dof

#define cp  24      // the number of contact point


Simulation* GLOB_SIM;
TimeDiscretisation * GLOB_T;
LagrangianLinearTIDS * GLOB_tabLDS;

int GLOB_COMPUTE;
int GLOB_STEP;
EventsManager * GLOB_EVT;

// Global variables for computation of CPU time, number of iterations and for curves ploting

#define PLOTMAX 2000
unsigned int outputSize = 3;
SimpleMatrix dataPlot(PLOTMAX + 1, outputSize);
int k_iter = 0; // index for output.


// #define PLOTMAX 1000
// #define OUT_SIZE 4
// SimpleMatrix GLOB_DATAPLOT( PLOTMAX+1,OUT_SIZE);

extern void initSiconos();
extern void computeSiconos();

void Start()
{
  GLOB_COMPUTE = false;
  GLOB_STEP = false;
  initSiconos();
}

struct my_struct_c
{
  float V3[3];
};

void DrawFEM(LagrangianLinearTIDS *lds, float radius)
{
  float R[12];
  float pos[3];

  R[0] = 1;
  R[1] = 0;
  R[2] = 0;
  R[4] = 0;
  R[5] = 1;
  R[6] = 0;
  R[8] = 0;
  R[9] = 0;
  R[10] = 1;

  dsSetTexture(DS_NONE);
  dsSetColor(0.6f, 0.6f, 1);


  typedef struct my_struct_c my_struct;
  vector<my_struct *> my_vector;

  for (int i = 0; i < FEM / 3; i++)
  {
    my_struct *MS = new my_struct;
    MS->V3[0] = lds->getQ()(3 * i);
    MS->V3[1] = lds->getQ()(3 * i + 1);
    MS->V3[2] = lds->getQ()(3 * i + 2);

    my_vector.push_back(MS);
    cout << MS->V3[0] << endl;
    cout << MS->V3[1] << endl;
    cout << MS->V3[2] << endl;
  }

  vector<my_struct *>::iterator my_it;
  //deque<my_struct *>::iterator my_it;

  my_it = my_vector.begin();

  while (my_it != my_vector.end())
  {
    pos[0] = (*my_it)->V3[0];
    pos[1] = (*my_it)->V3[1];
    pos[2] = (*my_it)->V3[2];

    dsDrawSphere(pos, R, radius);
    my_it++;
  }

}

void SimuLoop(int pause)
{
  float radius;

  if (GLOB_COMPUTE == true)
    computeSiconos();
  if (GLOB_STEP == true)
    GLOB_COMPUTE = false;

  // Ball Radius
  radius = 0.3;

  DrawFEM(GLOB_tabLDS, radius);
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

    double t0 = 0;                    // initial computation time
    double T = 1;                    // final computation time
    double h = 0.005;                 // time step

    string solverName = "NLGSNEWTON";      // solver algorithm used for non-smooth problem
    //string solverName = "NLGS";      // solver algorithm used for non-smooth problem

    double e = 0.8;                  // nslaw
    double mu = 10.;
    double rho = 6e-6;
    double g = 9.81;

    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    unsigned int i;
    unsigned int j;

    DynamicalSystemsSet allDS; // the list of DS

    SiconosMatrix *M;
    SiconosMatrix *K;
    SiconosMatrix *C;
    SiconosMatrix *Position;

    // mass matrix
    M = new SimpleMatrix("mass", 1);

    // rigid matrix
    K = new SimpleMatrix("rigid", 1);

    // maillage matrix
    Position = new SimpleMatrix("position", 1);


    // amortissement matrix
    C = new SimpleMatrix(FEM, FEM);

    // -- Initial positions and velocities --
    // q0 and v0.

    SimpleVector* q0 = new SimpleVector(FEM);
    SimpleVector* v0 = new SimpleVector(FEM);

    double gap = 1.;

    // Memory allocation for q0[i] and v0[i]
    // add the gap with ground

    for (i = 0; i < FEM / 3; i++)
    {
      (*q0)(3 * i) = (*Position)(3 * i, 0);
      (*q0)(3 * i + 1) = (*Position)(3 * i + 1, 0);
      (*q0)(3 * i + 2) = (*Position)(3 * i + 2, 0) + gap;
    }
    GLOB_tabLDS = new LagrangianLinearTIDS(0, *q0, *v0, *M, *K, *C);

    allDS.insert(GLOB_tabLDS);


    // -- Set external forces (weight) --
    SiconosVector * weight = new SimpleVector(FEM);
    for (i = 0; i < FEM / 3; i++)
      (*weight)(3 * i + 2) = -rho * g;
    GLOB_tabLDS->setFExtPtr(weight);

    // ==> at this point, all the required dynamical systems are saved in allDS.


    // -------------------
    // --- Interactions---
    // -------------------
    InteractionsSet allInteractions;

    //cp = [0 1 2 6 7 8 12 13 14 18 19 20 24 25 26 30 31 32 54 55 56 57 58 59 66 67 68 72 73 74 75 76 77 78 79 80 81 82 83 87 88 89 99 100 101 105 106 107 198 199 200 201 202 203 204 205 206 213 214 215 288 289 290 291 292 293 294 295 296 303 304 305]

    SiconosMatrix *H = new SimpleMatrix(cp, FEM);
    SimpleVector *b = new SimpleVector(cp);


    std::vector<int> v;
    v.resize(cp / 3);

    v[0] = 0;
    v[1] = 6;
    v[2] = 12;
    v[3] = 18;
    v[4] = 24;
    v[5] = 30;
    v[6] = 54;
    v[7] = 57;
    v[8] = 66;
    v[9] = 72;
    v[10] = 75;
    v[11] = 78;
    v[12] = 81;
    v[13] = 87;
    v[14] = 99;
    v[15] = 105;
    v[16] = 198;
    v[17] = 201;
    v[18] = 204;
    v[19] = 213;
    v[20] = 288;
    v[21] = 291;
    v[22] = 294;
    v[23] = 303;


    j = 0;
    for (size_t m = 0, size = v.size(); m < size; ++m)
    {
      (*H)(3 * j + 2, v[m] + 2) = 1.;
      (*H)(3 * j, v[m]) = 1.;
      (*H)(3 * j + 1, v[m] + 1) = 1.;
      (*b)(3 * j + 2) = (*Position)(v[m] + 2, 0) + gap;
      ++j;
    }

    NonSmoothLaw* nslaw = new NewtonImpactFrictionNSL(e, e, mu, 3);
    Relation* relation = new LagrangianLinearR(*H, *b);
    Interaction * inter = new Interaction("bead1", allDS, 0, cp, nslaw, relation);


    allInteractions.insert(inter);

    // --------------------------------
    // --- NonSmoothDynamicalSystem ---
    // --------------------------------
    bool isBVP = 0;
    NonSmoothDynamicalSystem * nsds = new NonSmoothDynamicalSystem(allDS, allInteractions, isBVP);

    // -------------
    // --- Model ---
    // -------------

    Model * TORE = new Model(t0, T);
    TORE->setNonSmoothDynamicalSystemPtr(nsds); // set NonSmoothDynamicalSystem of this model

    // ----------------
    // --- Simulation ---
    // ----------------

    GLOB_T = new TimeDiscretisation(h, TORE);
    GLOB_SIM = new TimeStepping(GLOB_T);

    // -- OneStepIntegrators --
    OneStepIntegrator * OSI = new Moreau(allDS , 0.5000001 , GLOB_SIM);

    // -- OneStepNsProblem --
    //OneStepNSProblem * osnspb = new LCP(GLOB_SIM,"FrictionContact3D",solverName,101,0.001);

    OneStepNSProblem * osnspb = new FrictionContact3D(GLOB_SIM , "FrictionContact3D", solverName, 100, 0.001);

    cout << "=== End of model loading === " << endl;

    // =========================== End of model definition
    // ================================= Computation

    // ================================= Computation =================================

    // --- Simulation initialization ---

    GLOB_SIM->initialize();
    cout << "End of simulation initialisation" << endl;

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot

    //  unsigned int outputSize = 1+3*DSNUMBER;
    //     SimpleMatrix dataPlot(N+1,outputSize);
    //     dataPlot(k_iter,0) = k_iter*GLOB_T->getH();
    //     dataPlot(k_iter,1) = GLOB_tabLDS[0]->getQ()(2);
    //     dataPlot(k_iter,2) = GLOB_tabLDS[0]->getVelocity()(2);
    //     dataPlot(k_iter,3) = (TORE->getNonSmoothDynamicalSystemPtr()->getInteractionPtr(1)->getLambda(1))(0);

    //   dataPlot(k_iter,4) = GLOB_tabLDS[1]->getQ()(2);
    //    dataPlot(k_iter,5) = GLOB_tabLDS[1]->getVelocity()(2);
    //    dataPlot(k_iter,6) = (TORE->getNonSmoothDynamicalSystemPtr()->getInteractionPtr(1)->getLambda(1))(0);

    // --- Time loop ---
    cout << "Start computation ... " << endl;

    GLOB_EVT = GLOB_SIM->getEventsManagerPtr();
  }
  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught in \'sample/NLM Init\'" << endl;
  }
}

void computeSiconos()
{
  try
  {

    // --- Time loop ---
    cout << "Start computation ... " << endl;

    // --- simulation solver ---
    if (GLOB_EVT->hasNextEvent())
    {
      GLOB_SIM->advanceToEvent();
      GLOB_SIM->processEvents();

      // --- Get values to be plotted ---
      k_iter++;
      //  dataPlot(k_iter,0) = k_iter*GLOB_T->getH();
      //  dataPlot(k_iter,1) = GLOB_tabLDS[0]->getQ()(2);
      //  dataPlot(k_iter,2) = GLOB_tabLDS[0]->getVelocity()(2);
      //  dataPlot(k_iter,3) = (TORE->getNonSmoothDynamicalSystemPtr()->getInteractionPtr(0)->getLambda(1))(0);
      //  dataPlot(k_iter,4) = GLOB_tabLDS[1]->getQ()(2);
      //  dataPlot(k_iter,5) = GLOB_tabLDS[1]->getVelocity()(2);
      //  dataPlot(k_iter,6) = (TORE->getNonSmoothDynamicalSystemPtr()->getInteractionPtr(1)->getLambda(1))(0);
    }
    cout << "End of computation - Number of iterations done: " << k_iter << endl;

    // --- Output files ---
    //    ioMatrix io("result.dat", "ascii");
    //     io.write(dataPlot,"noDim");
    //    cout<<"End of computation - Number of iterations done: "<<k<<endl;

  }
  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught in \'sample/TORE Init\'" << endl;
  }

}





