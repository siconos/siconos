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
 * 16/05/2007- Authors: houari khenous & Roger Pissard

*/
// =============================== Multi contact beads couple simulation ===============================
//
// ======================================================================================================

#include "SiconosKernel.h"
#include <sys/time.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <GL/gl.h>
#include <GL/glu.h>
#include <time.h>


//-----------------------------------------------------------------------------

#ifdef WIN32
#define NOMINMAX
#include <ctime>
#include <sys/timeb.h>
#include <time.h>
struct timeval
{
  long    tv_sec;         /* seconds */
  long    tv_usec;        /* and microseconds */
};

int gettimeofday(struct timeval *tp, void *tz)
{
  struct _timeb timebuffer;
  _ftime(&timebuffer);
  tp->tv_sec = static_cast<long>(timebuffer.time);
  tp->tv_usec = timebuffer.millitm * 1000;
  return 0;
}
#else
#include <sys/time.h>
#endif

//-----------------------------------------------------------------------------









using namespace std;

#include <drawstuff/drawstuff.h>

#define DSNUMBER   79      // the number of dynamical systems

#define WALL 7.5       // Positions of walls
#define TOP 2.2       // Positions of walls
#define GROUND 0       // Positions of walls


Simulation * GLOB_SIM;
TimeDiscretisation * GLOB_T;
LagrangianDS *GLOB_tabLDS[DSNUMBER];

int GLOB_COMPUTE;
int GLOB_STEP;
EventsManager * GLOB_EVT;

double T_max;
// Global variables for computation of CPU time, number of iterations and for curves ploting

// #define PLOTMAX 2000
// unsigned int outputSize = 5;
// SimpleMatrix dataPlot(PLOTMAX+1,outputSize);
int k_iter = 0; // index for output.


// #define PLOTMAX 1000
// #define OUT_SIZE 4
// SimpleMatrix GLOB_DATAPLOT( PLOTMAX+1,OUT_SIZE);

extern void initSiconos();
extern void computeSiconos();


void Start()
{

  float xyz[3], hpr[3];

  xyz[0] = 30.;
  xyz[1] = 7.;
  xyz[2] = 5.;
  hpr[0] = -170.;
  hpr[1] = hpr[2] = 0.;

  //  dsGetViewpoint (xyz,hpr);
  // hpr[0] += 7;
  dsSetViewpoint(xyz, hpr);

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

void DrawBox(float alpha)
{

  //  alpha = 0 signifie transparent, alpha = 1 signifie opaque
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glColor4f(0.5, 0.7, 0.9, alpha); // 3 premiers champs = RGB, quatri.ANhme champ = opacitNi

  glBegin(GL_QUADS);
  // Top Face
  glNormal3f(0.0f, 0.0f, 1.0f); // Normal Pointing Towards Viewer
  glVertex3f(-1.0f, -1.0f,  TOP);// Point 1 (Front)
  glVertex3f(1.0f, -1.0f,  TOP); // Point 2 (Front)
  glVertex3f(1.0f,  1.0f,  TOP); // Point 3 (Front)
  glVertex3f(-1.0f,  1.0f,  TOP);// Point 4 (Front)
  // Bottom Face
  glNormal3f(0.0f, 0.0f, -1.0f);// Normal Pointing Away From Viewer
  glVertex3f(-1.0f, -1.0f, GROUND);// Point 1 (Back)
  glVertex3f(-1.0f,  1.0f, GROUND);// Point 2 (Back)
  glVertex3f(1.0f,  1.0f, GROUND); // Point 3 (Back)
  glVertex3f(1.0f, -1.0f, GROUND); // Point 4 (Back)
  //Right Face
  glNormal3f(0.0f, 1.0f, 0.0f); // Normal Pointing Up
  glVertex3f(-1.0f,  1.0f, GROUND);// Point 1 (Top)
  glVertex3f(-1.0f,  1.0f, TOP);// Point 2 (Top)
  glVertex3f(1.0f,  1.0f, TOP); // Point 3 (Top)
  glVertex3f(1.0f,  1.0f, GROUND); // Point 4 (Top)
  // Left Face
  glNormal3f(0.0f, -1.0f, 0.0f);// Normal Pointing Down
  glVertex3f(-1.0f, -1.0f, GROUND);// Point 1 (Bottom)
  glVertex3f(1.0f, -1.0f, GROUND); // Point 2 (Bottom)
  glVertex3f(1.0f, -1.0f, TOP); // Point 3 (Bottom)
  glVertex3f(-1.0f, -1.0f, TOP);// Point 4 (Bottom)
  // Front Face
  glNormal3f(1.0f, 0.0f, 0.0f); // Normal Pointing Right
  glVertex3f(1.0f, -1.0f, GROUND); // Point 1 (Right)
  glVertex3f(1.0f,  1.0f, GROUND); // Point 2 (Right)
  glVertex3f(1.0f,  1.0f, TOP); // Point 3 (Right)
  glVertex3f(1.0f, -1.0f, TOP); // Point 4 (Right)
  // Back Face
  glNormal3f(-1.0f, 0.0f, 0.0f);// Normal Pointing Left
  glVertex3f(-1.0f, -1.0f, GROUND);// Point 1 (Left)
  glVertex3f(-1.0f, -1.0f, TOP);// Point 2 (Left)
  glVertex3f(-1.0f,  1.0f, TOP);// Point 3 (Left)
  glVertex3f(-1.0f,  1.0f, GROUND);// Point 4 (Left)
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
  radius = 1.;

  for (i = 0; i < DSNUMBER; i++)
  {
    DrawBall(GLOB_tabLDS[i], radius);
  }

  //  float alpha = 0.3; // alpha = 0 signifie transparent, alpha = 1 signifie opaque
  //   DrawBox(alpha);
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
  double T0, T1;
  timeval TIMER;

  gettimeofday(&TIMER, NULL);
  T0 = TIMER.tv_sec + (TIMER.tv_usec / 1000000.0);

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

  if (T_max == 2)
  {
    gettimeofday(&TIMER, NULL);
    T1 = TIMER.tv_sec + (TIMER.tv_usec / 1000000.0);
    cout << "Computation Time1: " << T1 - T0 << endl << endl;
  }
  return 0;

}

void initSiconos()
{
  try
  {

    // ================= Creation of the model =======================



    // User-defined main parameters

    unsigned int nDof = 6;            // degrees of freedom for beads

    double m = 1.;                   // mass of balls
    double R = 1.;                   // radius of balls

    double t0 = 0;                    // initial computation time
    double T = 2;                    // final computation time
    double h = 0.02;                 // time step


    string solverName = "NSGS";      // solver algorithm used for non-smooth problem
    //string solverName = "NLGS";      // solver algorithm used for non-smooth problem
    //string solverName = "PGS";      // solver algorithm used for non-smooth problem
    //string solverName = "Lemke";      // solver algorithm used for non-smooth problem

    double e  = 0.4;                  // nslaw
    double e2 = 0.4;                  // nslaw2
    double mu = 0.4;


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
    Fact = (DSNUMBER) * (DSNUMBER - 1) / 2;

    unsigned int i;
    unsigned int j;
    unsigned int l;

    DynamicalSystemsSet allDS; // the list of DS
    CheckInsertDS checkDS;

    SiconosMatrix *Mass = new SimpleMatrix(nDof, nDof);
    (*Mass)(0, 0) = (*Mass)(1, 1) = (*Mass)(2, 2) = m;    ;
    (*Mass)(3, 3) = (*Mass)(4, 4) = (*Mass)(5, 5) = 3. / 5 * m * R * R;

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

      //  (*(q0[i]))(2) = 1.2 + 3*i;
    }

    // set values


    unsigned int N1  = 32;

    for (j = 0; j < N1; j++)
    {
      if (j % 2 == 0)
      {
        (*(q0[j]))(0) = 1.75;
        (*(q0[j]))(1) = 1.25;
        (*(q0[j]))(2) = 6.25 + 3 * j;
      }
      else
      {
        (*(q0[j]))(0) = 1.75;
        (*(q0[j]))(1) = 1.75;
        (*(q0[j]))(2) = 6.25 + 3 * j;
      }
    }

    unsigned int N2  = 64;

    for (j = 32; j < N2; j++)
    {
      if (j % 2 == 0)
      {
        (*(q0[j]))(0) = 1.25;
        (*(q0[j]))(1) = 4.25;
        (*(q0[j]))(2) = 6.25 + 3 * (j - 32);
      }
      else
      {
        (*(q0[j]))(0) = 1.25;
        (*(q0[j]))(1) = 4.75;
        (*(q0[j]))(2) = 6.25 + 3 * (j - 32);
      }
    }


    unsigned int N3  = 79;

    for (j = 64; j < N3; j++)
    {
      if (j % 2 == 0)
      {
        (*(q0[j]))(0) = 1.75;
        (*(q0[j]))(1) = 7.25;
        (*(q0[j]))(2) = 6.25 + 3 * (j - 64);
      }
      else
      {
        (*(q0[j]))(0) = 1.75;
        (*(q0[j]))(1) = 7.75;
        (*(q0[j]))(2) = 6.25 + 3 * (j - 64);
      }
    }


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

    NonSmoothLaw * nslaw1 = new NewtonImpactFrictionNSL(e, e, mu, 3);

    // Interaction beads and plan1 (OXY)

    if (obst_z_m)
    {
      SiconosVector *b1 = new SimpleVector(3);
      (*b1)(0) = GROUND - R;
      SiconosMatrix *H1 = new SimpleMatrix(3, nDof);
      (*H1)(0, 2) = 1.0;
      (*H1)(1, 0) = 1.0;
      (*H1)(1, 4) = -R;
      (*H1)(2, 1) = 1.0;
      (*H1)(2, 3) =  R;

      for (i = 0; (int)i < DSNUMBER; i++)
      {
        dsConcernedi.insert(GLOB_tabLDS[i]);
        ostringstream ostr;
        ostr << i;
        id2[i] = ostr.str();
        LLR1[i] = new LagrangianLinearR(*H1, *b1);
        checkInter = allInteractions.insert(new Interaction(id2[i], dsConcernedi, i, 3, nslaw1, LLR1[i]));
        dsConcernedi.clear();
      }
    }

    // Interaction beads and plan1 (-YOX)

    if (obst_z_p)
    {
      SiconosVector *b1_ = new SimpleVector(3);
      (*b1_)(0) = TOP - R;
      SiconosMatrix *H1_ = new SimpleMatrix(3, nDof);
      (*H1_)(0, 2) = -1.0;
      (*H1_)(1, 0) = 1.0;
      (*H1_)(1, 4) = -R;
      (*H1_)(2, 1) = 1.0;
      (*H1_)(2, 3) =  R;

      for (i = 0; (int)i < DSNUMBER; i++)
      {
        dsConcernedi.insert(GLOB_tabLDS[i]);
        ostringstream ostr;
        ostr << i;
        id2[i] = ostr.str();
        LLR1_[i] = new LagrangianLinearR(*H1_, *b1_);
        checkInter = allInteractions.insert(new Interaction(id2[i], dsConcernedi, i, 3, nslaw1, LLR1_[i]));
        dsConcernedi.clear();
      }
    }


    // Interaction beads and plan2 (OXZ)

    if (obst_y_m)
    {
      SiconosVector *b2 = new SimpleVector(3);
      (*b2)(0) = GROUND - R;
      SiconosMatrix *H2 = new SimpleMatrix(3, nDof);
      (*H2)(0, 1) = 1.0;
      (*H2)(1, 0) = 1.0;
      (*H2)(1, 5) = -R;
      (*H2)(2, 2) = 1.0;
      (*H2)(2, 3) =  R;

      for (i = 0; (int)i < DSNUMBER; i++)
      {
        dsConcernedi.insert(GLOB_tabLDS[i]);
        ostringstream ostr;
        ostr << i;
        id2[i] = ostr.str();
        LLR2[i] = new LagrangianLinearR(*H2, *b2);
        checkInter = allInteractions.insert(new Interaction(id2[i], dsConcernedi, i, 3, nslaw1, LLR2[i]));
        dsConcernedi.clear();
      }
    }

    // Interaction beads and plan2 (-ZOX)

    if (obst_y_p)
    {
      SiconosVector *b2_ = new SimpleVector(3);
      (*b2_)(0) = 2 * WALL - R;
      SiconosMatrix *H2_ = new SimpleMatrix(3, nDof);
      (*H2_)(0, 1) = -1.0;
      (*H2_)(1, 0) = 1.0;
      (*H2_)(1, 5) = -R;
      (*H2_)(2, 2) = 1.0;
      (*H2_)(2, 3) =  R;

      for (i = 0; (int)i < DSNUMBER; i++)
      {
        dsConcernedi.insert(GLOB_tabLDS[i]);
        ostringstream ostr;
        ostr << i;
        id2[i] = ostr.str();
        LLR2_[i] = new LagrangianLinearR(*H2_, *b2_);
        checkInter = allInteractions.insert(new Interaction(id2[i], dsConcernedi, i, 3, nslaw1, LLR2_[i]));
        dsConcernedi.clear();
      }
    }

    // Interaction beads and plan3 (OYZ)

    if (obst_x_m)
    {
      SiconosVector *b3 = new SimpleVector(3);
      (*b3)(0) = GROUND - R;
      SiconosMatrix *H3 = new SimpleMatrix(3, nDof);
      (*H3)(0, 0) = 1.0;
      (*H3)(1, 1) = 1.0;
      (*H3)(1, 5) = -R;
      (*H3)(2, 2) = 1.0;
      (*H3)(2, 4) =  R;

      for (i = 0; (int)i < DSNUMBER; i++)
      {
        dsConcernedi.insert(GLOB_tabLDS[i]);
        ostringstream ostr;
        ostr << i;
        id2[i] = ostr.str();
        LLR3[i] = new LagrangianLinearR(*H3, *b3);
        checkInter = allInteractions.insert(new Interaction(id2[i], dsConcernedi, i, 3, nslaw1, LLR3[i]));
        dsConcernedi.clear();
      }
    }
    // Interaction beads and plan3 (-ZOY)

    if (obst_x_p)
    {
      SiconosVector *b3_ = new SimpleVector(3);
      (*b3_)(0) = 2 * WALL - R;
      SiconosMatrix *H3_ = new SimpleMatrix(3, nDof);
      (*H3_)(0, 0) = -1.0;
      (*H3_)(1, 1) = 1.0;
      (*H3_)(1, 5) = -R;
      (*H3_)(2, 2) = 1.0;
      (*H3_)(2, 4) =  R;

      for (i = 0; (int)i < DSNUMBER; i++)
      {
        dsConcernedi.insert(GLOB_tabLDS[i]);
        ostringstream ostr;
        ostr << i;
        id2[i] = ostr.str();
        LLR3_[i] = new LagrangianLinearR(*H3_, *b3_);
        checkInter = allInteractions.insert(new Interaction(id2[i], dsConcernedi, i, 3, nslaw1, LLR3_[i]));
        dsConcernedi.clear();
      }
    }

    // Interaction between beads

    // frictional contact condition between beads
    NonSmoothLaw * nslaw2 = new NewtonImpactFrictionNSL(e2, e2, mu, 3);

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
          LLR[l] = new LagrangianScleronomousR("3DDrawPlugin:h0", "3DDrawPlugin:G0");
          checkInter = allInteractions.insert(new Interaction(id[l], dsConcerned2, l, 3, nslaw2, LLR[l]));
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
    OneStepIntegrator * OSI = new Moreau(allDS , 0.5001 , GLOB_SIM);

    // -- OneStepNsProblem --


    //OneStepNSProblem * osnspb = new LCP(GLOB_SIM ,"FrictionContact3D",solverName,10000, 0.001);

    OneStepNSProblem * osnspb = new FrictionContact3D(GLOB_SIM , "FrictionContact3D", solverName, 10000, 0.0005);

    cout << "=== End of model loading === " << endl;

    // =========================== End of model definition

    // ================================= Computation =================================

    // --- Simulation initialization ---

    GLOB_SIM->initialize();

    cout << "End of simulation initialisation" << endl;
    //    N = GLOB_T->getNSteps(); // Number of time steps

    //     dataPlot(k_iter,0) = k_iter*GLOB_T->getH();
    //     dataPlot(k_iter,1) = GLOB_tabLDS[0]->getQ()(2);
    //     dataPlot(k_iter,2) = GLOB_tabLDS[0]->getVelocity()(2);
    //     // dataPlot(k_iter,3) = (multiBeads->getNonSmoothDynamicalSystemPtr()->getInteractionPtr(0)->getLambda(1))(0);
    //     dataPlot(k_iter,3) = GLOB_tabLDS[1]->getQ()(2);
    //     dataPlot(k_iter,4) = GLOB_tabLDS[1]->getVelocity()(2);
    //    dataPlot(k_iter,7) = (multiBeads->getNonSmoothDynamicalSystemPtr()->getInteractionPtr(1)->getLambda(1))(0);
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

      k_iter++;
      T_max = k_iter * GLOB_T->getH();
      //      dataPlot(k_iter,0) = k_iter*GLOB_T->getH();
      //       dataPlot(k_iter,1) = GLOB_tabLDS[0]->getQ()(2);
      //       dataPlot(k_iter,2) = GLOB_tabLDS[0]->getVelocity()(2);
      //       //    dataPlot(k_iter,3) = (multiBeads->getNonSmoothDynamicalSystemPtr()->getInteractionPtr(0)->getLambda(1))(0);
      //       dataPlot(k_iter,3) = GLOB_tabLDS[1]->getQ()(2);
      //       dataPlot(k_iter,4) = GLOB_tabLDS[1]->getVelocity()(2);
      //dataPlot(k_iter,7) = (multiBeads->getNonSmoothDynamicalSystemPtr()->getInteractionPtr(1)->getLambda(1))(0);

    }
    cout << "End of computation - Number of iterations done: " << k_iter << endl;



    // --- Output files ---
    //     ioMatrix io("result.dat", "ascii");
    //     io.write(dataPlot,"noDim");
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


