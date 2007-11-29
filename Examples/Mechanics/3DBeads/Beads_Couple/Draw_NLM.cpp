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
 * Couple beads simulation - 3D frictionl contact problem in presence of a rigid foundation
 *
 * Each couple consists in two sticked beads where each one does not move with repsect to the other
 *
 * 26/05/2007- Authors: houari khenous
 *
 * For more details about visualization, see README
 *
 * For more details about computations, a document will be added to this directory soon
 *
 * Last modification 27/11/2007, H. Khenous

*/
// =============================== Multi contact beads couple simulation ===============================
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


#define COUPLE   2     // the number of dynamical systems

#define WALL 1       // Positions of walls
#define TOP 2.2       // Positions of walls
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


void DrawBox(float alpha)
{

  //  alpha = 0 signifie transparent, alpha = 1 signifie opaque
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glColor4f(0.5, 0.7, 0.9, alpha); // 3 premiers champs = RGB, quatri.ANhme champ = opacitNi

  glBegin(GL_QUADS);
  // Front Face
  glNormal3f(0.0f, 0.0f, 1.0f); // Normal Pointing Towards Viewer
  glVertex3f(-WALL, -WALL,  TOP);// Point 1 (Front)
  glVertex3f(WALL, -WALL,  TOP); // Point 2 (Front)
  glVertex3f(WALL,  WALL,  TOP); // Point 3 (Front)
  glVertex3f(-WALL,  WALL,  TOP);// Point 4 (Front)
  // Back Face
  glNormal3f(0.0f, 0.0f, -1.0f);// Normal Pointing Away From Viewer
  glVertex3f(-WALL, -WALL, GROUND);// Point 1 (Back)
  glVertex3f(-WALL,  WALL, GROUND);// Point 2 (Back)
  glVertex3f(WALL,  WALL, GROUND); // Point 3 (Back)
  glVertex3f(WALL, -WALL, GROUND); // Point 4 (Back)
  // Top Face
  glNormal3f(0.0f, 1.0f, GROUND); // Normal Pointing Up
  glVertex3f(-WALL,  WALL, GROUND);// Point 1 (Top)
  glVertex3f(-WALL,  WALL, TOP);// Point 2 (Top)
  glVertex3f(WALL,  WALL, TOP); // Point 3 (Top)
  glVertex3f(WALL,  WALL, GROUND); // Point 4 (Top)
  // Bottom Face
  glNormal3f(0.0f, -1.0f, GROUND);// Normal Pointing Down
  glVertex3f(-WALL, -WALL, GROUND);// Point 1 (Bottom)
  glVertex3f(WALL, -WALL, GROUND); // Point 2 (Bottom)
  glVertex3f(WALL, -WALL, TOP); // Point 3 (Bottom)
  glVertex3f(-WALL, -WALL, TOP);// Point 4 (Bottom)
  // Right face
  glNormal3f(1.0f, 0.0f, GROUND); // Normal Pointing Right
  glVertex3f(WALL, -WALL, GROUND); // Point 1 (Right)
  glVertex3f(WALL,  WALL, GROUND); // Point 2 (Right)
  glVertex3f(WALL,  WALL, TOP); // Point 3 (Right)
  glVertex3f(WALL, -WALL, TOP); // Point 4 (Right)
  // Left Face
  glNormal3f(-1.0f, 0.0f, GROUND);// Normal Pointing Left
  glVertex3f(-WALL, -WALL, GROUND);// Point 1 (Left)
  glVertex3f(-WALL, -WALL, TOP);// Point 2 (Left)
  glVertex3f(-WALL,  WALL, TOP);// Point 3 (Left)
  glVertex3f(-WALL,  WALL, GROUND);// Point 4 (Left)
  glEnd();// Done Drawing Quads

  glDisable(GL_BLEND);

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

  float alpha = 0.3; // alpha = 0 signifie transparent, alpha = 1 signifie opaque
  DrawBox(alpha);
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
  fn.path_to_textures = "../textures/";

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

    string solverName = "NSGS";      // solver algorithm used for non-smooth problem

    double e  = 0.9;                  // nslaw
    double e2 = 0.5;
    double mu = 0.;
    double PI = 3.14;
    double R  = 0.1;

    // 1 to take in account the obstacle and  0 no

    int obst_z_p = 0;                    //  for z --> +
    int obst_z_m = 1;                    //  for z --> -
    int obst_y_p = 1;                    //  for y --> +
    int obst_y_m = 1;                    //  for y --> -
    int obst_x_p = 1;                    //  for x --> +
    int obst_x_m = 1;                    //  for x --> -


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

    (*(q0[0]))(0) = 0.;
    (*(q0[0]))(1) = 0.;
    (*(q0[0]))(2) =  0.5;
    (*(q0[0]))(3) =  PI / 2;
    (*(q0[1]))(0) = 0.;
    (*(q0[1]))(1) = 0.;
    (*(q0[1]))(2) =  0.2;
    (*(q0[1]))(3) =  PI / 2;





    // set values

    //   /* Test of beads Couples in cube*/
    //     // 9*etage beads Couples in cube /* !!!!!! When you choose cote and etage you should verify that COUPLE = (cote^2)*etage + 1 (crazy beads couple)!!!!! */
    //     unsigned int cote  = 3;
    //     unsigned int etage  = 10;
    //     unsigned int k = 0;

    //     /* !!!!!!!!!!!!!!!!!!Very Important thing is that theta must be different from zero for stablity reason !!!!!!!!!!!!!!!!!!!!!!!*/
    //     for (j=0;j<COUPLE;j++)
    //       (*(q0[j]))(3) =  PI/2;


    //     for (j=0;j<etage;j++){
    //       for (k=0;k<cote;k++) {
    //  for (i=0;i<cote;i++) {
    //    if (j % 2 == 0){
    //      (*(q0[k*cote+i+j*cote*cote]))(0) = -0.5 + 5*k*R; (*(q0[k*cote+i+j*cote*cote]))(1) = -0.5 + 5*i*R; (*(q0[k*cote+i+j*cote*cote]))(2) = 0.12+(2*j)*R;
    //    }
    //    else{
    //      (*(q0[k*cote+i+j*cote*cote]))(0) = -0.4 + 5*k*R; (*(q0[k*cote+i+j*cote*cote]))(1) = -0.4 + 5*i*R+R/4;(*(q0[k*cote+i+j*cote*cote]))(2) = 0.12+(2*j)*R;
    //    }
    //  }
    //       }
    //     }

    //     /*  Crazy beads couple */

    //     (*(q0[COUPLE-1]))(0) = -0.8; (*(q0[COUPLE-1]))(1) = -0.8; (*(q0[COUPLE-1]))(2) =  0.7; (*(q0[COUPLE-1]))(3) =  PI/3; (*(q0[COUPLE-1]))(4) =  PI/3; (*(v0[COUPLE-1]))(2) =  1.; (*(v0[COUPLE-1]))(3) =  4.;


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
    vector<Relation*> LLR11z(COUPLE);
    vector<Relation*> LLR11z_(COUPLE);
    vector<Relation*> LLR22z(COUPLE);
    vector<Relation*> LLR22z_(COUPLE);
    vector<Relation*> LLR11y(COUPLE);
    vector<Relation*> LLR11y_(COUPLE);
    vector<Relation*> LLR22y(COUPLE);
    vector<Relation*> LLR22y_(COUPLE);
    vector<Relation*> LLR11x(COUPLE);
    vector<Relation*> LLR11x_(COUPLE);
    vector<Relation*> LLR22x(COUPLE);
    vector<Relation*> LLR22x_(COUPLE);


    NonSmoothLaw * nslaw11 = new NewtonImpactFrictionNSL(e, e, mu, 3);

    // Z axis
    if (obst_z_m)
    {
      for (i = 0; (int)i < COUPLE; i++)
      {
        dsConcernedii.insert(GLOB_tabLDS[i]);
        ostringstream ostr;
        ostr << i;
        id22[i] = ostr.str();
        LLR11z[i] = new LagrangianScleronomousR("NLMPlugin:h1z", "NLMPlugin:G1z");
        checkInter = allInteractions.insert(new Interaction(id22[i], dsConcernedii, i, 3, nslaw11, LLR11z[i]));
        dsConcernedii.clear();
      }

      for (i = 0; (int)i < COUPLE; i++)
      {
        dsConcernedii.insert(GLOB_tabLDS[i]);
        ostringstream ostr;
        ostr << i;
        id22[i] = ostr.str();
        LLR22z[i] = new LagrangianScleronomousR("NLMPlugin:h2z", "NLMPlugin:G2z");
        checkInter = allInteractions.insert(new Interaction(id22[i], dsConcernedii, i, 3, nslaw11, LLR22z[i]));
        dsConcernedii.clear();
      }
    }
    if (obst_z_p)
    {
      for (i = 0; (int)i < COUPLE; i++)
      {
        dsConcernedii.insert(GLOB_tabLDS[i]);
        ostringstream ostr;
        ostr << i;
        id22[i] = ostr.str();
        LLR11z_[i] = new LagrangianScleronomousR("NLMPlugin:h1z_", "NLMPlugin:G1z_");
        checkInter = allInteractions.insert(new Interaction(id22[i], dsConcernedii, i, 3, nslaw11, LLR11z_[i]));
        dsConcernedii.clear();
      }

      for (i = 0; (int)i < COUPLE; i++)
      {
        dsConcernedii.insert(GLOB_tabLDS[i]);
        ostringstream ostr;
        ostr << i;
        id22[i] = ostr.str();
        LLR22z_[i] = new LagrangianScleronomousR("NLMPlugin:h2z_", "NLMPlugin:G2z_");
        checkInter = allInteractions.insert(new Interaction(id22[i], dsConcernedii, i, 3, nslaw11, LLR22z_[i]));
        dsConcernedii.clear();
      }
    }
    // Y axis
    if (obst_y_m)
    {
      for (i = 0; (int)i < COUPLE; i++)
      {
        dsConcernedii.insert(GLOB_tabLDS[i]);
        ostringstream ostr;
        ostr << i;
        id22[i] = ostr.str();
        LLR11y[i] = new LagrangianScleronomousR("NLMPlugin:h1y", "NLMPlugin:G1y");
        checkInter = allInteractions.insert(new Interaction(id22[i], dsConcernedii, i, 3, nslaw11, LLR11y[i]));
        dsConcernedii.clear();
      }

      for (i = 0; (int)i < COUPLE; i++)
      {
        dsConcernedii.insert(GLOB_tabLDS[i]);
        ostringstream ostr;
        ostr << i;
        id22[i] = ostr.str();
        LLR22y[i] = new LagrangianScleronomousR("NLMPlugin:h2y", "NLMPlugin:G2y");
        checkInter = allInteractions.insert(new Interaction(id22[i], dsConcernedii, i, 3, nslaw11, LLR22y[i]));
        dsConcernedii.clear();
      }
    }
    if (obst_y_p)
    {
      for (i = 0; (int)i < COUPLE; i++)
      {
        dsConcernedii.insert(GLOB_tabLDS[i]);
        ostringstream ostr;
        ostr << i;
        id22[i] = ostr.str();
        LLR11y_[i] = new LagrangianScleronomousR("NLMPlugin:h1y_", "NLMPlugin:G1y_");
        checkInter = allInteractions.insert(new Interaction(id22[i], dsConcernedii, i, 3, nslaw11, LLR11y_[i]));
        dsConcernedii.clear();
      }

      for (i = 0; (int)i < COUPLE; i++)
      {
        dsConcernedii.insert(GLOB_tabLDS[i]);
        ostringstream ostr;
        ostr << i;
        id22[i] = ostr.str();
        LLR22y_[i] = new LagrangianScleronomousR("NLMPlugin:h2y_", "NLMPlugin:G2y_");
        checkInter = allInteractions.insert(new Interaction(id22[i], dsConcernedii, i, 3, nslaw11, LLR22y_[i]));
        dsConcernedii.clear();
      }
    }

    // X axis
    if (obst_x_m)
    {
      for (i = 0; (int)i < COUPLE; i++)
      {
        dsConcernedii.insert(GLOB_tabLDS[i]);
        ostringstream ostr;
        ostr << i;
        id22[i] = ostr.str();
        LLR11x[i] = new LagrangianScleronomousR("NLMPlugin:h1x", "NLMPlugin:G1x");
        checkInter = allInteractions.insert(new Interaction(id22[i], dsConcernedii, i, 3, nslaw11, LLR11x[i]));
        dsConcernedii.clear();
      }

      for (i = 0; (int)i < COUPLE; i++)
      {
        dsConcernedii.insert(GLOB_tabLDS[i]);
        ostringstream ostr;
        ostr << i;
        id22[i] = ostr.str();
        LLR22x[i] = new LagrangianScleronomousR("NLMPlugin:h2x", "NLMPlugin:G2x");
        checkInter = allInteractions.insert(new Interaction(id22[i], dsConcernedii, i, 3, nslaw11, LLR22x[i]));
        dsConcernedii.clear();
      }
    }
    if (obst_x_p)
    {
      for (i = 0; (int)i < COUPLE; i++)
      {
        dsConcernedii.insert(GLOB_tabLDS[i]);
        ostringstream ostr;
        ostr << i;
        id22[i] = ostr.str();
        LLR11x_[i] = new LagrangianScleronomousR("NLMPlugin:h1x_", "NLMPlugin:G1x_");
        checkInter = allInteractions.insert(new Interaction(id22[i], dsConcernedii, i, 3, nslaw11, LLR11x_[i]));
        dsConcernedii.clear();
      }

      for (i = 0; (int)i < COUPLE; i++)
      {
        dsConcernedii.insert(GLOB_tabLDS[i]);
        ostringstream ostr;
        ostr << i;
        id22[i] = ostr.str();
        LLR22x_[i] = new LagrangianScleronomousR("NLMPlugin:h2x_", "NLMPlugin:G2x_");
        checkInter = allInteractions.insert(new Interaction(id22[i], dsConcernedii, i, 3, nslaw11, LLR22x_[i]));
        dsConcernedii.clear();
      }
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

    cout << "Start computation ... " << endl;

    GLOB_EVT = GLOB_SIM->getEventsManagerPtr();
  }

  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught in \'sample/Beads Coulpe Init\'" << endl;
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
    cout << "Exception caught in \'sample/Beads Couple\'" << endl;
  }
}


