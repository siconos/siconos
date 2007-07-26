/* Siconos-sample version 2.0.0, Copyright INRIA 2005-2007.
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
 * 26/05/2007- Authors: houari khenous

*/
// =============================== Multi bouncing beads couple simulation ===============================
//
// Keywords: LagrangianLinearDS, LagrangianDS relation, Moreau TimeStepping, newton method.
//
// ======================================================================================================

#include "SiconosKernel.h"
#include <sys/time.h>
#include <math.h>
#include <iostream>
#include <fstream>

using namespace std;

#include <drawstuff/drawstuff.h>


#define COUPLE   3      // the number of dynamical systems

#define WALL 1       // Positions of walls
#define TOP 1       // Positions of walls
#define GROUND 0       // Positions of walls


Simulation* GLOB_SIM;
TimeDiscretisation * GLOB_T;
LagrangianDS * GLOB_tabLDS[COUPLE];
//Model * BeadsCOUPLE;


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

void DrawCouple(LagrangianDS *lds, float radius)
{
  float R[12];
  float pos[3], pos11[3], pos22[3];
  float theta, phi, psi;
  // Translation
  pos[0] = lds->getQ()(0);
  pos[1] = lds->getQ()(1);
  pos[2] = lds->getQ()(2);
  // Rotation
  theta = lds->getQ()(3);
  phi = lds->getQ()(4);
  psi = lds->getQ()(5);

  // R[0] = cos(psi)*cos(phi) - cos(theta)*sin(psi)*sin(phi);
  //   R[1] = cos(psi)*sin(phi) + cos(theta)*cos(psi)*sin(phi);
  //   R[2] = sin(theta)*sin(psi);
  //   R[4] = -sin(psi)*cos(phi) - cos(theta)*sin(psi)*cos(phi);
  //   R[5] = -sin(psi)*sin(phi) + cos(theta)*cos(psi)*cos(phi);
  //   R[6] = sin(theta)*cos(psi);
  //   R[8] = sin(theta)*sin(phi);
  //   R[9] = -sin(theta)*cos(phi);
  //   R[10] = cos(theta);
  //   R[3]=R[7]=R[11]=0;

  R[0] = 1;
  R[1] = 0;
  R[2] = 0;
  R[4] = 0;
  R[5] = cos(theta);
  R[6] = -sin(theta);
  R[8] = 0;
  R[9] = sin(theta);
  R[10] = cos(theta);


  pos11[0] = pos[0] - radius * sin(theta) * sin(phi);
  pos11[1] = pos[1] +  radius * sin(theta) * cos(phi);
  pos11[2] = pos[2] - radius * cos(theta);

  pos22[0] = pos[0] + radius * sin(theta) * sin(phi);
  pos22[1] = pos[1] - radius * sin(theta) * cos(phi);
  pos22[2] = pos[2] + radius * cos(theta);

  dsSetTexture(DS_WOOD);
  dsSetColor(1, 0.8f, 0.6f);
  dsDrawSphere(pos11, R, radius);

  dsSetTexture(DS_NONE);
  dsSetColor(0.6f, 0.6f, 1);
  dsDrawSphere(pos22, R, radius);


}
void Drawwall()
{

  int Wall_z_p = 0;                    //  for z --> +
  int Wall_y_p = 1;                    //  for y --> +
  int Wall_y_m = 1;                    //  for y --> -
  int Wall_x_p = 1;                    //  for x --> +
  int Wall_x_m = 1;                    //  for x --> -

  double pos1[3];
  pos1[0] = pos1[1] = WALL;
  pos1[2] = GROUND;
  double pos2[3];
  pos2[0] = -WALL;
  pos2[1] = WALL;
  pos2[2] = GROUND;
  double pos3[3];
  pos3[0] = pos3[1] = -WALL;
  pos3[2] = GROUND;
  double pos4[3];
  pos4[0] = WALL;
  pos4[1] = -WALL;
  pos4[2] = GROUND;

  dsSetColor(1, 0.8f, 0.6f);
  int k;

  //  y-wall of the cube (positive direction)
  if (Wall_y_p)
  {
    pos1[2] = pos2[2] = pos3[2] = pos4[2] = GROUND;
    for (k = 0; k < 15; ++k)
    {
      dsDrawLineD(pos1, pos2);
      pos1[2] += k * 0.01;
      pos2[2] += k * 0.01;
    }
  }
  //  x-wall of the cube (negative direction)
  if (Wall_x_m)
  {
    pos1[2] = pos2[2] = GROUND;
    for (k = 0; k < 15; ++k)
    {
      dsDrawLineD(pos2, pos3);
      pos2[2] += k * 0.01;
      pos3[2] += k * 0.01;
    }
  }
  //  y-wall of the cube (negative direction)
  if (Wall_y_m)
  {
    pos3[2] = GROUND;
    for (k = 0; k < 15; ++k)
    {
      dsDrawLineD(pos3, pos4);
      pos3[2] += k * 0.01;
      pos4[2] += k * 0.01;
    }
  }
  //  x-wall of the cube (positive direction)
  if (Wall_x_p)
  {
    pos4[2] = GROUND;
    for (k = 0; k < 15; ++k)
    {
      dsDrawLineD(pos4, pos1);
      pos4[2] += k * 0.01;
      pos1[2] += k * 0.01;
    }
  }
  //  Top wall of the cube
  if (Wall_z_p)
  {
    pos1[2] = pos2[2] = pos3[2] = pos4[2] = TOP;
    dsDrawLineD(pos1, pos2);
    dsDrawLineD(pos2, pos3);
    dsDrawLineD(pos3, pos4);
    dsDrawLineD(pos4, pos1);
    for (k = 0; k < 20; ++k)
    {
      dsDrawLineD(pos3, pos4);
      pos3[1] += k * 0.01;
      pos4[1] += k * 0.01;
    }
  }
}

void SimuLoop(int pause)
{
  int i;
  float radius;

  if (GLOB_COMPUTE == true)
    computeSiconos();
  if (GLOB_STEP == true)
    GLOB_COMPUTE = false;

  // Ball Radius
  radius = 0.1;

  for (i = 0; i < COUPLE; i++)
  {
    DrawCouple(GLOB_tabLDS[i], radius);
  }
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

    unsigned int nDof = 6;            // degrees of freedom for beads

    double t0 = 0;                    // initial computation time
    double T = 10.;                    // final computation time
    double h = 0.005;                 // time step

    string solverName = "NLGSNEWTON";      // solver algorithm used for non-smooth problem
    //string solverName = "NLGS";      // solver algorithm used for non-smooth problem
    //string solverName = "Lemke";      // solver algorithm used for non-smooth problem

    double e  = 0.9;                  // nslaw
    double e2  = 0.5;
    double mu = 0.;
    double PI = 3.14;


    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    int Fact;
    Fact = (COUPLE) * (COUPLE - 1) / 2;

    unsigned int i;
    unsigned int j;
    unsigned int l;

    DynamicalSystemsSet allDS; // the list of DS
    CheckInsertDS checkDS;

    // -- Initial positions and velocities --

    vector<SimpleVector *> q0;
    vector<SimpleVector *> v0;
    q0.resize(COUPLE, NULL);
    v0.resize(COUPLE, NULL);

    // Memory allocation for q0[i] and v0[i]
    for (i = 0; i < COUPLE; i++)
    {
      q0[i] = new SimpleVector(nDof);
      v0[i] = new SimpleVector(nDof);
    }

    // set values

    (*(q0[0]))(0) = 0.;
    (*(q0[0]))(1) = 0.;
    (*(q0[0]))(2) =  0.15;
    (*(q0[0]))(3) =  PI / 2;
    (*(v0[0]))(3) =  4.;
    (*(q0[1]))(0) = 0.;
    (*(q0[1]))(1) = 0.;
    (*(q0[1]))(2) =  0.40;
    (*(q0[1]))(3) =  PI / 2;
    (*(v0[0]))(3) =  4.;
    (*(q0[2]))(0) = 0.;
    (*(q0[2]))(1) = 0.;
    (*(q0[2]))(2) =  0.65;
    (*(q0[2]))(3) =  PI / 3;
    (*(q0[2]))(4) =  PI / 3;
    (*(v0[2]))(3) =  4.;

    //  for (i=0;i<COUPLE;i++)
    //       {
    //  (*(q0[i]))(0) = 0.; (*(q0[i]))(1) = 0.5*(i+1.);  (*(q0[i]))(2) =  0.3; (*(q0[i]))(3) =  PI/2;
    //  //  (*(q0[i]))(1) = 0.5*i;
    //       }

    for (i = 0; i < COUPLE; i++)
    {
      GLOB_tabLDS[i] = new LagrangianDS(i, *(q0[i]), *(v0[i]));
      checkDS = allDS.insert(GLOB_tabLDS[i]);
      (static_cast<LagrangianDS*>(*(checkDS.first)))->setComputeFExtFunction("NLMPlugin.so", "gravity");
      (static_cast<LagrangianDS*>(*(checkDS.first)))->setComputeMassFunction("NLMPlugin.so", "mass");
      (static_cast<LagrangianDS*>(*(checkDS.first)))->setComputeNNLFunction("NLMPlugin.so", "NNL");
      (static_cast<LagrangianDS*>(*(checkDS.first)))->setComputeJacobianNNLFunction(0, "NLMPlugin.so", "jacobianQNNL");
      (static_cast<LagrangianDS*>(*(checkDS.first)))->setComputeJacobianNNLFunction(1, "NLMPlugin.so", "jacobianVNNL");
    }

    // ==> at this point, all the required dynamical systems are saved in allDS.

    // -------------------
    // --- Interactions---
    // -------------------
    InteractionsSet allInteractions;
    CheckInsertInteraction checkInter;

    // Interactions with walls

    // Interactions between beads couple
    vector<string> id22;
    id22.resize(COUPLE);
    vector<string> idd;
    idd.resize(Fact);

    DynamicalSystemsSet dsConcernedii;
    DynamicalSystemsSet dsConcerned22;

    vector<Relation*> LLRR(Fact);
    vector<Relation*> LLR11(COUPLE);
    vector<Relation*> LLR22(COUPLE);

    NonSmoothLaw * nslaw11 = new NewtonImpactFrictionNSL(e, e, mu, 3);

    for (i = 0; (int)i < COUPLE; i++)
    {
      dsConcernedii.insert(GLOB_tabLDS[i]);
      ostringstream ostr;
      ostr << i;
      id22[i] = ostr.str();
      LLR11[i] = new LagrangianScleronomousR("NLMPlugin:h1", "NLMPlugin:G1");
      checkInter = allInteractions.insert(new Interaction(id22[i], dsConcernedii, i, 3, nslaw11, LLR11[i]));
      dsConcernedii.clear();
    }

    for (i = 0; (int)i < COUPLE; i++)
    {
      dsConcernedii.insert(GLOB_tabLDS[i]);
      ostringstream ostr;
      ostr << i;
      id22[i] = ostr.str();
      LLR22[i] = new LagrangianScleronomousR("NLMPlugin:h2", "NLMPlugin:G2");
      checkInter = allInteractions.insert(new Interaction(id22[i], dsConcernedii, i, 3, nslaw11, LLR22[i]));
      dsConcernedii.clear();
    }


    // Interaction between beads

    NonSmoothLaw * nslaw12 = new NewtonImpactFrictionNSL(e2, e2, mu, 3);

    l = 0;
    for (i = 0; (int)i < COUPLE; i++)
    {
      dsConcerned22.insert(GLOB_tabLDS[i]);
      for (j = 0; (int)j < COUPLE; j++)
      {
        if (j > i)
        {
          dsConcerned22.insert(GLOB_tabLDS[j]);
          ostringstream ostr;
          ostr << l;
          idd[l] = ostr.str();
          LLRR[l] = new LagrangianScleronomousR("NLMPlugin:h22", "NLMPlugin:G22");
          checkInter = allInteractions.insert(new Interaction(idd[l], dsConcerned22, l, 3, nslaw12, LLRR[l]));
          dsConcerned22.erase(GLOB_tabLDS[j]);
          l = l + 1;
        }
      }
      dsConcerned22.clear();
    }

    l = 0;
    for (i = 0; (int)i < COUPLE; i++)
    {
      dsConcerned22.insert(GLOB_tabLDS[i]);
      for (j = 0; (int)j < COUPLE; j++)
      {
        if (j > i)
        {
          dsConcerned22.insert(GLOB_tabLDS[j]);
          ostringstream ostr;
          ostr << l;
          idd[l] = ostr.str();
          LLRR[l] = new LagrangianScleronomousR("NLMPlugin:h21", "NLMPlugin:G21");
          checkInter = allInteractions.insert(new Interaction(idd[l], dsConcerned22, l, 3, nslaw12, LLRR[l]));
          dsConcerned22.erase(GLOB_tabLDS[j]);
          l = l + 1;
        }
      }
      dsConcerned22.clear();
    }

    l = 0;
    for (i = 0; (int)i < COUPLE; i++)
    {
      dsConcerned22.insert(GLOB_tabLDS[i]);
      for (j = 0; (int)j < COUPLE; j++)
      {
        if (j > i)
        {
          dsConcerned22.insert(GLOB_tabLDS[j]);
          ostringstream ostr;
          ostr << l;
          idd[l] = ostr.str();
          LLRR[l] = new LagrangianScleronomousR("NLMPlugin:h12", "NLMPlugin:G12");
          checkInter = allInteractions.insert(new Interaction(idd[l], dsConcerned22, l, 3, nslaw12, LLRR[l]));
          dsConcerned22.erase(GLOB_tabLDS[j]);
          l = l + 1;
        }
      }
      dsConcerned22.clear();
    }

    l = 0;
    for (i = 0; (int)i < COUPLE; i++)
    {
      dsConcerned22.insert(GLOB_tabLDS[i]);
      for (j = 0; (int)j < COUPLE; j++)
      {
        if (j > i)
        {
          dsConcerned22.insert(GLOB_tabLDS[j]);
          ostringstream ostr;
          ostr << l;
          idd[l] = ostr.str();
          LLRR[l] = new LagrangianScleronomousR("NLMPlugin:h11", "NLMPlugin:G11");
          checkInter = allInteractions.insert(new Interaction(idd[l], dsConcerned22, l, 3, nslaw12, LLRR[l]));
          dsConcerned22.erase(GLOB_tabLDS[j]);
          l = l + 1;
        }
      }
      dsConcerned22.clear();
    }



    // --------------------------------
    // --- NonSmoothDynamicalSystem ---
    // --------------------------------
    bool isBVP = 0;
    NonSmoothDynamicalSystem * nsds = new NonSmoothDynamicalSystem(allDS, allInteractions, isBVP);

    // -------------
    // --- Model ---
    // -------------

    Model * BeadsCOUPLE = new Model(t0, T);
    BeadsCOUPLE->setNonSmoothDynamicalSystemPtr(nsds); // set NonSmoothDynamicalSystem of this model

    // ----------------
    // --- Simulation ---
    // ----------------

    GLOB_T = new TimeDiscretisation(h, BeadsCOUPLE);
    GLOB_SIM = new TimeStepping(GLOB_T);

    // -- OneStepIntegrators --
    OneStepIntegrator * OSI = new Moreau(allDS , 0.5000001 , GLOB_SIM);

    // -- OneStepNsProblem --
    OneStepNSProblem * osnspb = new FrictionContact3D(GLOB_SIM , "FrictionContact3D", solverName, 10000000, 0.001);
    //OneStepNSProblem * osnspb = new LCP(GLOB_SIM,"name",solverName,100, 0.001,1);

    cout << "=== End of model loading === " << endl;

    // =========================== End of model definition
    // ================================= Computation

    // ================================= Computation =================================

    // --- Simulation initialization ---

    GLOB_SIM->initialize();

    cout << "End of simulation initialisation" << endl;

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot

    //int N = GLOB_T->getNSteps(); // Number of time steps

    dataPlot(k_iter, 0) = k_iter * GLOB_T->getH();
    // dataPlot(k_iter,1) = (BeadsCOUPLE->getNonSmoothDynamicalSystemPtr()->getInteractionPtr(0)->getY(0))(0);
    dataPlot(k_iter, 1) = GLOB_tabLDS[0]->getQ()(2);
    dataPlot(k_iter, 2) = GLOB_tabLDS[0]->getVelocity()(2);

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
    //  cout << "Start computation ... " << endl;

    // --- simulation solver ---
    if (GLOB_EVT->hasNextEvent())
    {
      GLOB_SIM->advanceToEvent();
      GLOB_SIM->processEvents();
      // --- Get values to be plotted ---
      k_iter++;
      dataPlot(k_iter, 0) = k_iter * GLOB_T->getH();
      //dataPlot(k_iter,1) = (BeadsCOUPLE->getNonSmoothDynamicalSystemPtr()->getInteractionPtr(0)->getY(0))(0);
      dataPlot(k_iter, 1) = GLOB_tabLDS[0]->getQ()(2);
      dataPlot(k_iter, 2) = GLOB_tabLDS[0]->getVelocity()(2);

    }
    cout << "End of computation - Number of iterations done: " << k_iter << endl;


    // --- Output files ---
    ioMatrix io("result.dat", "ascii");
    io.write(dataPlot, "noDim");

  }

  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught in \'sample/CoupleBeads\'" << endl;
  }
}


