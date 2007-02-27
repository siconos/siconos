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
 * Two beads 3D frictionless contact problem in presence of a rigid foundation
 * 30/01/2007- Authors: houari khenous & Roger Pissard

*/
// =============================== Multi bouncing beads column simulation ===============================
//  N beads between a floor and a ceiling ...
// Keywords: LagrangianLinearDS, LagrangianLinear relation, Moreau TimeStepping, LCP.
//
// ======================================================================================================

#include "SiconosKernel.h"
#include <sys/time.h>
#include <math.h>
#include <iostream>
#include <fstream>

using namespace std;

#include <drawstuff/drawstuff.h>

#define DSNUMBER  20       // the number of dynamical systems

Simulation* GLOB_SIM;
TimeDiscretisation * GLOB_T;
LagrangianDS *GLOB_tabLDS[DSNUMBER];

int GLOB_COMPUTE;
int GLOB_STEP;
EventsManager * GLOB_EVT;

vector<Relation*> GLOB_LLR(3);


#define PLOTMAX 1000
#define OUT_SIZE 4
SimpleMatrix GLOB_DATAPLOT(PLOTMAX + 1, OUT_SIZE);

extern void initSiconos();
extern void computeSiconos();


void Start()
{
  GLOB_COMPUTE = false;
  GLOB_STEP = false;
  initSiconos();
}

void DrawBall(LagrangianDS *lds, float radius)
{
  float R[12];
  float pos[3];
  float theta, phi, psi;

  // Translation
  pos[0] = lds->getQ()(0);
  pos[1] = lds->getQ()(1);
  pos[2] = lds->getQ()(2);

  // Rotation
  theta = lds->getQ()(3);
  phi = lds->getQ()(4);
  psi = lds->getQ()(5);
  R[0] = 1;
  R[1] = 0;
  R[2] = 0;
  R[4] = 0;
  R[5] = cos(theta);
  R[6] = -sin(theta);
  R[8] = 0;
  R[9] = sin(theta);
  R[10] = cos(theta);

  R[3] = R[7] = R[11] = 0;

  dsSetTexture(DS_WOOD);
  dsSetColor(1, 0.8f, 0.6f);

  dsDrawSphere(pos, R, radius);

}

void Drawwall()
{
  double pos1[3];
  pos1[0] = pos1[1] = 1.;
  pos1[2] = 0.;
  double pos2[3];
  pos2[0] = -1.;
  pos2[1] = 1.;
  pos2[2] = 0.;
  double pos3[3];
  pos3[0] = pos3[1] = -1.;
  pos3[2] = 0.;
  double pos4[3];
  pos4[0] = 1.;
  pos4[1] = -1.;
  pos4[2] = 0.;

  dsSetColor(1, 0.8f, 0.6f);
  int k;
  pos1[2] = pos2[2] = pos3[2] = pos4[2] = 1.;
  dsDrawLineD(pos1, pos2);
  dsDrawLineD(pos2, pos3);
  dsDrawLineD(pos3, pos4);
  dsDrawLineD(pos4, pos1);

  pos1[2] = pos2[2] = pos3[2] = pos4[2] = 0;
  for (k = 0; k < 15; ++k)
  {
    dsDrawLineD(pos1, pos2);
    pos1[2] += k * 0.01;
    pos2[2] += k * 0.01;
  }
  pos1[2] = pos2[2] = 0;
  for (k = 0; k < 15; ++k)
  {
    dsDrawLineD(pos2, pos3);
    pos2[2] += k * 0.01;
    pos3[2] += k * 0.01;
  }
  pos3[2] = 0;
  for (k = 0; k < 15; ++k)
  {
    dsDrawLineD(pos3, pos4);
    pos3[2] += k * 0.01;
    pos4[2] += k * 0.01;
  }
  pos4[2] = 0;
  for (k = 0; k < 15; ++k)
  {
    dsDrawLineD(pos4, pos1);
    pos4[2] += k * 0.01;
    pos1[2] += k * 0.01;
  }
  pos1[2] = pos2[2] = pos3[2] = pos4[2] = 0;
  pos3[2] = pos4[2] = 1.;
  for (k = 0; k < 20; ++k)
  {
    dsDrawLineD(pos3, pos4);
    pos3[1] += k * 0.01;
    pos4[1] += k * 0.01;
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

  for (i = 0; i < DSNUMBER; i++)
  {
    DrawBall(GLOB_tabLDS[i], radius);
  }
  Drawwall();

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

    double m = 2.;                   // mass of balls
    double R = 0.1;                   // radius of balls

    double t0 = 0;                    // initial computation time
    double T = 10;                    // final computation time
    double h = 0.005;                 // time step

    string solverName = "Lemke";      // solver algorithm used for non-smooth problem

    double e = 0.5;                  // nslaw
    double e2 = 0.5;                  // nslaw2
    double mu = 10.;
    double mu2 = 0.5;



    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    int Fact;
    Fact = (DSNUMBER) * (DSNUMBER - 1) / 2;

    unsigned int i;
    unsigned int j;
    unsigned int l;
    DynamicalSystemsSet allDS; // the list of DS
    CheckInsertDS checkDS;

    SiconosMatrix *Mass = new SimpleMatrix(nDof, nDof);
    (*Mass)(0, 0) = (*Mass)(1, 1) = (*Mass)(2, 2) = m;    ;
    (*Mass)(3, 3) = (*Mass)(4, 4) = (*Mass)(5, 5) = 3. / 5 * R * R;

    // -- Initial positions and velocities --
    // q0[i] and v0[i] correspond to position and velocity of ball i.

    vector<SimpleVector *> q0;
    vector<SimpleVector *> v0;
    q0.resize(DSNUMBER, NULL);
    v0.resize(DSNUMBER, NULL);

    // Memory allocation for q0[i] and v0[i]
    for (i = 0; i < DSNUMBER; i++)
    {
      q0[i] = new SimpleVector(nDof);
      v0[i] = new SimpleVector(nDof);
    }

    // set values

    // billard

    //    (*(q0[0]))(0) =  0.;     (*(q0[0]))(1) =  0.;   (*(q0[0]))(2) =  0.1;
    //     (*(q0[1]))(0) =  0.1;    (*(q0[1]))(1) = -0.2;  (*(q0[1]))(2) =  0.1;
    //     (*(q0[2]))(0) = -0.1;    (*(q0[2]))(1) = -0.2;  (*(q0[2]))(2) =  0.1;
    //     (*(q0[3]))(0) =  0.2;    (*(q0[3]))(1) = -0.4;  (*(q0[3]))(2) =  0.1;
    //     (*(q0[4]))(0) = -0.2;    (*(q0[4]))(1) = -0.4;  (*(q0[4]))(2) =  0.1;
    //     (*(q0[5]))(0) =  0.;     (*(q0[5]))(1) = -0.4;  (*(q0[5]))(2) =  0.1;
    //     (*(q0[6]))(0) =  0.1;    (*(q0[6]))(1) = -0.6;  (*(q0[6]))(2) =  0.1;
    //     (*(q0[7]))(0) = -0.1;    (*(q0[7]))(1) = -0.6;  (*(q0[7]))(2) =  0.1;

    //     (*(q0[8]))(0) =  0.;     (*(q0[8]))(1) = 0.8;    (*(q0[8]))(2) =  0.1;

    //     (*(v0[8]))(0) =  -1; (*(v0[8]))(1) =  -20;

    //     (*(q0[9]))(0) =  0.3;    (*(q0[9]))(1) = -0.6;  (*(q0[9]))(2) =  0.1;
    //     (*(q0[10]))(0)= -0.3;    (*(q0[10]))(1)= -0.6;  (*(q0[10]))(2)=  0.1;
    //     (*(q0[11]))(0)=  0.2;    (*(q0[11]))(1)= -0.8;  (*(q0[11]))(2)=  0.1;
    //     (*(q0[12]))(0)= -0.2;    (*(q0[12]))(1)= -0.8;  (*(q0[12]))(2)=  0.1;
    //     (*(q0[13]))(0)=  0.;     (*(q0[13]))(1)= -0.8;  (*(q0[13]))(2)=  0.1;
    //     (*(q0[14]))(0)=  0.4;    (*(q0[14]))(1)= -0.8;  (*(q0[14]))(2)=  0.1;
    //     (*(q0[15]))(0)= -0.4;    (*(q0[15]))(1)= -0.8;  (*(q0[15]))(2)=  0.1;

    // Cube de billes


    (*(q0[0]))(0) =  0.2;
    (*(q0[0]))(1) = -0.2;
    (*(q0[0]))(2) =  0.2;
    (*(q0[1]))(0) =  0.2;
    (*(q0[1]))(1) =  0.2;
    (*(q0[1]))(2) =  0.2;
    (*(q0[2]))(0) = -0.2;
    (*(q0[2]))(1) =  0.2;
    (*(q0[2]))(2) =  0.2;
    (*(q0[3]))(0) = -0.2;
    (*(q0[3]))(1) = -0.2;
    (*(q0[3]))(2) =  0.2;
    (*(q0[4]))(0) =  0.25;
    (*(q0[4]))(1) = -0.2;
    (*(q0[4]))(2) =  0.4;
    (*(q0[5]))(0) =  0.25;
    (*(q0[5]))(1) =  0.2;
    (*(q0[5]))(2) =  0.4;
    (*(q0[6]))(0) = -0.25;
    (*(q0[6]))(1) =  0.2;
    (*(q0[6]))(2) =  0.4;
    (*(q0[7]))(0) = -0.25;
    (*(q0[7]))(1) = -0.2;
    (*(q0[7]))(2) =  0.4;

    (*(q0[8]))(0) =  0.;
    (*(q0[8]))(1) = 0.3;
    (*(q0[8]))(2) =  0.1;

    (*(v0[8]))(1) =  -10;

    (*(q0[9]))(0) =  0.2;
    (*(q0[9]))(1) =  0.2;
    (*(q0[9]))(2) =  0.6;
    (*(q0[10]))(0) = -0.2;
    (*(q0[10]))(1) =  0.2;
    (*(q0[10]))(2) =  0.6;
    (*(q0[11]))(0) = -0.2;
    (*(q0[11]))(1) = -0.2;
    (*(q0[11]))(2) =  0.6;
    (*(q0[12]))(0) =  0.25;
    (*(q0[12]))(1) = -0.2;
    (*(q0[12]))(2) =  0.8;
    (*(q0[13]))(0) =  0.25;
    (*(q0[13]))(1) =  0.2;
    (*(q0[13]))(2) =  0.8;
    (*(q0[14]))(0) = -0.25;
    (*(q0[14]))(1) =  0.2;
    (*(q0[14]))(2) =  0.9;
    (*(q0[15]))(0) = -0.25;
    (*(q0[15]))(1) = -0.2;
    (*(q0[15]))(2) =  0.9;
    (*(q0[16]))(0) =  0.35;
    (*(q0[16]))(1) = -0.35;
    (*(q0[16]))(2) =  0.1;
    (*(q0[17]))(0) =  0.35;
    (*(q0[17]))(1) = 0.35;
    (*(q0[17]))(2) =  0.1;
    (*(q0[18]))(0) = -0.35;
    (*(q0[18]))(1) = 0.35;
    (*(q0[18]))(2) =  0.1;
    (*(q0[19]))(0) = -0.35;
    (*(q0[19]))(1) = -0.35;
    (*(q0[19]))(2) =  0.1;


    for (i = 0; i < DSNUMBER; i++)
    {
      GLOB_tabLDS[i] = new LagrangianDS(i, *(q0[i]), *(v0[i]), *Mass);
      checkDS = allDS.insert(GLOB_tabLDS[i]);
      (static_cast<LagrangianDS*>(*(checkDS.first)))->setComputeFExtFunction("3DDrawPlugin.so", "gravity");
    }

    // ==> at this point, all the required dynamical systems are saved in allDS.

    // -------------------
    // --- Interactions---
    // -------------------
    InteractionsSet allInteractions;

    vector<string> id;
    vector<string> id2;
    id.resize(Fact);
    id2.resize(DSNUMBER);

    DynamicalSystemsSet dsConcernedi;
    DynamicalSystemsSet dsConcerned2 ;
    CheckInsertInteraction checkInter;
    vector<Relation*> LLR(Fact);
    vector<Relation*> LLR1(DSNUMBER);
    vector<Relation*> LLR1_(DSNUMBER);
    vector<Relation*> LLR2(DSNUMBER);
    vector<Relation*> LLR2_(DSNUMBER);
    vector<Relation*> LLR3(DSNUMBER);
    vector<Relation*> LLR3_(DSNUMBER);



    NonSmoothLaw * nslaw1 = new NewtonImpactNSL(e);

    // Interaction beads and plan1 (OXY)

    SiconosVector *b1 = new SimpleVector(1);
    (*b1)(0) = -R;
    SiconosMatrix *H1 = new SimpleMatrix(1, nDof);
    (*H1)(0, 2) = 1.0;

    Relation * relation1 = new LagrangianLinearR(*H1, *b1);

    for (i = 0; (int)i < DSNUMBER; i++)
    {
      dsConcernedi.insert(GLOB_tabLDS[i]);
      ostringstream ostr;
      ostr << i;
      id2[i] = ostr.str();
      LLR1[i] = new LagrangianLinearR(*relation1);
      checkInter = allInteractions.insert(new Interaction(id2[i], dsConcernedi, i, 1, nslaw1, LLR1[i]));
      dsConcernedi.clear();
    }

    // Interaction beads and plan1 (-YOX)

    SiconosVector *b1_ = new SimpleVector(1);
    (*b1_)(0) = 1. - R;
    SiconosMatrix *H1_ = new SimpleMatrix(1, nDof);
    (*H1_)(0, 2) = -1.0;

    Relation * relation1_ = new LagrangianLinearR(*H1_, *b1_);

    for (i = 0; (int)i < DSNUMBER; i++)
    {
      dsConcernedi.insert(GLOB_tabLDS[i]);
      ostringstream ostr;
      ostr << i;
      id2[i] = ostr.str();
      LLR1_[i] = new LagrangianLinearR(*relation1_);
      checkInter = allInteractions.insert(new Interaction(id2[i], dsConcernedi, i, 1, nslaw1, LLR1_[i]));
      dsConcernedi.clear();
    }


    // Interaction beads and plan2 (OXZ)

    SiconosVector *b2 = new SimpleVector(1);
    (*b2)(0) = 1. - R;
    SiconosMatrix *H2 = new SimpleMatrix(1, nDof);
    (*H2)(0, 1) = 1.0;

    Relation * relation2 = new LagrangianLinearR(*H2, *b2);

    for (i = 0; (int)i < DSNUMBER; i++)
    {
      dsConcernedi.insert(GLOB_tabLDS[i]);
      ostringstream ostr;
      ostr << i;
      id2[i] = ostr.str();
      LLR2[i] = new LagrangianLinearR(*relation2);
      checkInter = allInteractions.insert(new Interaction(id2[i], dsConcernedi, i, 1, nslaw1, LLR2[i]));
      dsConcernedi.clear();
    }


    // Interaction beads and plan2 (-ZOX)

    SiconosVector *b2_ = new SimpleVector(1);
    (*b2_)(0) = 1. - R;
    SiconosMatrix *H2_ = new SimpleMatrix(1, nDof);
    (*H2_)(0, 1) = -1.0;

    Relation * relation2_ = new LagrangianLinearR(*H2_, *b2_);

    for (i = 0; (int)i < DSNUMBER; i++)
    {
      dsConcernedi.insert(GLOB_tabLDS[i]);
      ostringstream ostr;
      ostr << i;
      id2[i] = ostr.str();
      LLR2_[i] = new LagrangianLinearR(*relation2_);
      checkInter = allInteractions.insert(new Interaction(id2[i], dsConcernedi, i, 1, nslaw1, LLR2_[i]));
      dsConcernedi.clear();
    }


    // Interaction beads and plan3 (OYZ)

    SiconosVector *b3 = new SimpleVector(1);
    (*b3)(0) = 1. - R;
    SiconosMatrix *H3 = new SimpleMatrix(1, nDof);
    (*H3)(0, 0) = 1.0;

    Relation * relation3 = new LagrangianLinearR(*H3, *b3);

    for (i = 0; (int)i < DSNUMBER; i++)
    {
      dsConcernedi.insert(GLOB_tabLDS[i]);
      ostringstream ostr;
      ostr << i;
      id2[i] = ostr.str();
      LLR3[i] = new LagrangianLinearR(*relation3);
      checkInter = allInteractions.insert(new Interaction(id2[i], dsConcernedi, i, 1, nslaw1, LLR3[i]));
      dsConcernedi.clear();
    }


    // Interaction beads and plan3 (-ZOY)

    SiconosVector *b3_ = new SimpleVector(1);
    (*b3_)(0) = 1. - R;
    SiconosMatrix *H3_ = new SimpleMatrix(1, nDof);
    (*H3_)(0, 0) = -1.0;

    Relation * relation3_ = new LagrangianLinearR(*H3_, *b3_);

    for (i = 0; (int)i < DSNUMBER; i++)
    {
      dsConcernedi.insert(GLOB_tabLDS[i]);
      ostringstream ostr;
      ostr << i;
      id2[i] = ostr.str();
      LLR3_[i] = new LagrangianLinearR(*relation3_);
      checkInter = allInteractions.insert(new Interaction(id2[i], dsConcernedi, i, 1, nslaw1, LLR3_[i]));
      dsConcernedi.clear();
    }


    // Interaction between beads

    NonSmoothLaw * nslaw2 = new NewtonImpactNSL(e2);

    vector<string> G;
    G.reserve(1);
    G.push_back("3DDrawPlugin:Gcontact");

    Relation * relation22 = new LagrangianR("scleronomic", "3DDrawPlugin:h0", G);

    l = 0;
    for (i = 0; (int)i < DSNUMBER; i++)
    {
      dsConcerned2.insert(GLOB_tabLDS[i]);
      for (j = 0; (int)j < DSNUMBER; j++)
      {
        if (j > i)
        {
          dsConcerned2.insert(GLOB_tabLDS[j]);

          ostringstream ostr;
          ostr << l;
          id[l] = ostr.str();
          LLR[l] = new LagrangianR(*relation22);
          checkInter = allInteractions.insert(new Interaction(id[l], dsConcerned2, l, 1, nslaw2, LLR[l]));
          dsConcerned2.erase(GLOB_tabLDS[j]);
          l = l + 1;
        }
      }
      dsConcerned2.clear();
    }

    // --------------------------------
    // --- NonSmoothDynamicalSystem ---
    // --------------------------------
    bool isBVP = 0;
    NonSmoothDynamicalSystem * nsds = new NonSmoothDynamicalSystem(allDS, allInteractions, isBVP);

    // -------------
    // --- Model ---
    // -------------

    Model * multiBeads = new Model(t0, T);
    multiBeads->setNonSmoothDynamicalSystemPtr(nsds); // set NonSmoothDynamicalSystem of this model

    // ----------------
    // --- Simulation ---
    // ----------------

    GLOB_T = new TimeDiscretisation(h, multiBeads);
    GLOB_SIM = new TimeStepping(GLOB_T);

    // -- OneStepIntegrators --
    OneStepIntegrator * OSI = new Moreau(allDS , 0.5000001 , GLOB_SIM);

    // -- OneStepNsProblem --
    OneStepNSProblem * osnspb = new LCP(GLOB_SIM , "FrictionContact3D", solverName, 101, 0.001);

    //OneStepNSProblem * osnspb = new FrictionContact3D(GLOB_SIM ,"FrictionContact3D",solverName,101, 0.001);

    cout << "=== End of model loading === " << endl;

    // =========================== End of model definition
    // ================================= Computation
    // --- Simulation initialization ---
    GLOB_SIM ->initialize();

    GLOB_EVT = GLOB_SIM->getEventsManagerPtr();

    cout << "End of simulation initialisation" << endl;

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
  unsigned int maxIter = 5000;
  double criterion = 0.001;

  try
  {

    // --- simulation solver ---
    if (GLOB_EVT->hasNextEvent())
    {
      GLOB_SIM->advanceToEvent();
      GLOB_SIM->processEvents();
    }
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


