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
#include <drawstuff/drawstuff.h>

using namespace std;

Simulation* GLOB_SIM;
TimeDiscretisation * GLOB_T;
LagrangianDS * GLOB_LDS1;
LagrangianDS * GLOB_LDS2;
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

void SimuLoop(int pause)
{
  int k;
  float radius;

  if (GLOB_COMPUTE == true)
    computeSiconos();
  if (GLOB_STEP == true)
    GLOB_COMPUTE = false;


  // Ball Radius
  radius = 0.1;

  DrawBall(GLOB_LDS1, radius);
  DrawBall(GLOB_LDS2, radius);


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
    unsigned int dsNumber = 2;        // the number of dynamical systems
    unsigned int nDof = 6;            // degrees of freedom for beads

    double m1 = 1.;                   // mass of ball 1
    double R1 = 0.1;                   // radius of ball 1

    double m2 = 2.;                   // mass of ball 2
    double R2 = 0.1;                   // radius of ball 2

    double t0 = 0;                    // initial computation time
    double T = 10;                    // final computation time
    double h = 0.005;                 // time step

    string solverName = "Lemke";      // solver algorithm used for non-smooth problem

    double e = 0.5;                  // nslaw
    double e2 = 0.5;                  // nslaw2
    double mu = 100.;
    double mu2 = 0.5;

    // initial position ball 1
    double x1 = 0.;
    double yy1 = 0.;
    double z1 = 0.1;
    double theta1 = 0.;
    double phi1 = 0;
    double psi1 = 0;
    // initial velocity ball 1
    double dot_x1 = 0.;
    double dot_y1 = 0.;
    double dot_z1 = 0.;
    double dot_theta1 = 0.;
    double dot_phi1 = 0.;
    double dot_psi1 = 0.;
    // initial position ball 2
    double x2 = 0.5;
    double y2 = 0.5;
    double z2 = 0.2;
    double theta2 = 0.;
    double phi2 = 0;
    double psi2 = 0.;
    // initial velocity ball 2
    double dot_x2 = -2.;
    double dot_y2 = -2.;
    double dot_z2 = 2.;
    double dot_theta2 = 0.;
    double dot_phi2 = 0.;
    double dot_psi2 = 0.;

    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    unsigned int i;
    DynamicalSystemsSet allDS; // the list of DS
    CheckInsertDS checkDS;


    // -- Initial positions and velocities --
    vector<SimpleVector *> q0;
    vector<SimpleVector *> velocity0;
    q0.resize(dsNumber, NULL);
    velocity0.resize(dsNumber, NULL);

    q0[0] = new SimpleVector(nDof);
    velocity0[0] = new SimpleVector(nDof);

    (*(q0[0]))(0) = x1;
    (*(velocity0[0]))(0) = dot_x1;
    (*(q0[0]))(1) = yy1;
    (*(velocity0[0]))(1) = dot_y1;
    (*(q0[0]))(2) = z1;
    (*(velocity0[0]))(2) = dot_z1;
    (*(q0[0]))(3) = theta1;
    (*(velocity0[0]))(3) = dot_theta1;
    (*(q0[0]))(4) = phi1;
    (*(velocity0[0]))(4) = dot_phi1;
    (*(q0[0]))(5) = psi1;
    (*(velocity0[0]))(5) = dot_psi1;

    q0[1] = new SimpleVector(nDof);
    velocity0[1] = new SimpleVector(nDof);

    (*(q0[1]))(0) = x2;
    (*(velocity0[1]))(0) = dot_x2;
    (*(q0[1]))(1) = y2;
    (*(velocity0[1]))(1) = dot_y2;
    (*(q0[1]))(2) = z2;
    (*(velocity0[1]))(2) = dot_z2;
    (*(q0[1]))(3) = theta2;
    (*(velocity0[1]))(3) = dot_theta2;
    (*(q0[1]))(4) = phi2;
    (*(velocity0[1]))(4) = dot_phi2;
    (*(q0[1]))(5) = psi2;
    (*(velocity0[1]))(5) = dot_psi2;

    GLOB_LDS1 = new LagrangianDS(0, *(q0[0]), *(velocity0[0]));
    GLOB_LDS2 = new LagrangianDS(1, *(q0[1]), *(velocity0[1]));

    // weight of beads as internal forces

    SimpleVector * poids1 = new SimpleVector(nDof);
    SimpleVector * poids2 = new SimpleVector(nDof);
    double g1 = 9.81;
    double g2 = 9.81;

    (*poids1)(2) =  -m1 * g1;
    (*poids2)(2) =  -m2 * g2;

    GLOB_LDS1->setFExtPtr(poids1);
    GLOB_LDS2->setFExtPtr(poids2);


    // external forces plug-in

    GLOB_LDS1->setComputeMassFunction("3DBeadsPlugin.so", "Mass1");
    GLOB_LDS2->setComputeMassFunction("3DBeadsPlugin.so", "Mass2");

    allDS.insert(GLOB_LDS1);
    allDS.insert(GLOB_LDS2);

    // ==> at this point, all the required dynamical systems are saved in allDS.

    // -------------------
    // --- Interactions---
    // -------------------

    InteractionsSet allInteractions;
    int interactionNumber = dsNumber + 1;
    vector<string> id;
    id.resize(interactionNumber - 2);
    unsigned int nInter = 3; // number of relations in each interaction
    DynamicalSystemsSet dsConcerned0;
    DynamicalSystemsSet dsConcerned1;

    // Interaction beads and floor
    SiconosVector *b0 = new SimpleVector(3);
    SiconosVector *b1 = new SimpleVector(3);
    (*b0)(0) = -R1;
    (*b1)(0) = -R2;
    SiconosMatrix *H1 = new SimpleMatrix(3, nDof);
    (*H1)(0, 2) = 1.0;
    (*H1)(1, 0) = 1.0;
    (*H1)(1, 4) = -R1;
    (*H1)(2, 1) = 1.0;
    (*H1)(2, 3) =  R1;
    SiconosMatrix *H2 = new SimpleMatrix(3, nDof);
    (*H2)(0, 2) = 1.0;
    (*H2)(1, 0) = 1.0;
    (*H2)(1, 4) = -R2;
    (*H2)(2, 1) = 1.0;
    (*H2)(2, 3) =  R2;


    NonSmoothLaw * nslaw = new NewtonImpactFrictionNSL(e, e, mu, 3); // new NewtonImpactNSL(e); //

    Relation * relation0 = new LagrangianLinearR(*H1, *b0);
    Relation * relation1 = new LagrangianLinearR(*H2, *b1);
    dsConcerned0.insert(GLOB_LDS1);
    Interaction * inter0 = new Interaction("floor_bead1", dsConcerned0, 0, 3, nslaw, relation0);
    dsConcerned1.insert(GLOB_LDS2);
    Interaction * inter1 = new Interaction("floor_bead2", dsConcerned1, 1, 3, nslaw, relation1);

    allInteractions.insert(inter0);
    allInteractions.insert(inter1);

    dsConcerned0.clear();
    dsConcerned1.clear();

    // Interaction between beads

    DynamicalSystemsSet dsConcerned2 ;
    CheckInsertInteraction checkInter;


    NonSmoothLaw * nslaw2 = new NewtonImpactFrictionNSL(e2, e2, mu2, 3);// new NewtonImpactNSL(e2); //

    string G = "3DBeadsPlugin:G0";
    Relation * relation2 = new LagrangianScleronomousR("3DBeadsPlugin:h0", G);

    for (i = 1; (int)i < interactionNumber - 1; i++)
    {
      dsConcerned2.insert(allDS.getDynamicalSystemPtr(i - 1));
      dsConcerned2.insert(allDS.getDynamicalSystemPtr(i));
      ostringstream ostr;
      ostr << i;
      id[i - 1] = ostr.str();
      GLOB_LLR[i - 1] = new LagrangianR("3DBeadsPlugin:h0", G);
      checkInter = allInteractions.insert(new Interaction(id[i - 1], dsConcerned2, i, nInter, nslaw2, GLOB_LLR[i - 1]));
      dsConcerned2.clear();
    }

    delete relation2;

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


