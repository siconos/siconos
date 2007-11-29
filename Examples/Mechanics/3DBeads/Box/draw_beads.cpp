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
 *
 * 16/05/2007- Authors: houari khenous & Roger Pissard
 *
 * Last modification 27/11/2007, H. Khenous

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

using namespace std;

#include <drawstuff/drawstuff.h>

#define DSNUMBER  16     // the number of dynamical systems

#define WALL 1       // Positions of walls
#define TOP 2.2       // Positions of walls
#define GROUND 0       // Positions of walls

Simulation * GLOB_SIM;
TimeDiscretisation * GLOB_T;
LagrangianDS *GLOB_tabLDS[DSNUMBER];

int GLOB_COMPUTE;
int GLOB_STEP;
EventsManager * GLOB_EVT;

// Global variables for computation of CPU time, number of iterations and for curves ploting

#define PLOTMAX 2000
unsigned int outputSize = 5;
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
  glNormal3f(0.0f, 1.0f, 0.0f); // Normal Pointing Up
  glVertex3f(-WALL,  WALL, GROUND);// Point 1 (Top)
  glVertex3f(-WALL,  WALL, TOP);// Point 2 (Top)
  glVertex3f(WALL,  WALL, TOP); // Point 3 (Top)
  glVertex3f(WALL,  WALL, GROUND); // Point 4 (Top)
  // Bottom Face
  glNormal3f(0.0f, -1.0f, 0.0f);// Normal Pointing Down
  glVertex3f(-WALL, -WALL, GROUND);// Point 1 (Bottom)
  glVertex3f(WALL, -WALL, GROUND); // Point 2 (Bottom)
  glVertex3f(WALL, -WALL, TOP); // Point 3 (Bottom)
  glVertex3f(-WALL, -WALL, TOP);// Point 4 (Bottom)
  // Right face
  glNormal3f(1.0f, 0.0f, 0.0f); // Normal Pointing Right
  glVertex3f(WALL, -WALL, GROUND); // Point 1 (Right)
  glVertex3f(WALL,  WALL, GROUND); // Point 2 (Right)
  glVertex3f(WALL,  WALL, TOP); // Point 3 (Right)
  glVertex3f(WALL, -WALL, TOP); // Point 4 (Right)
  // Left Face
  glNormal3f(-1.0f, 0.0f, 0.0f);// Normal Pointing Left
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

  for (i = 0; i < DSNUMBER; i++)
  {
    DrawBall(GLOB_tabLDS[i], radius);
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

    double m = 1.;                   // mass of balls
    double R = 0.1;                   // radius of balls

    double t0 = 0;                    // initial computation time
    double T = 10.;                    // final computation time
    double h = 0.005;                 // time step


    string solverName = "NSGS";      // solver algorithm used for non-smooth problem
    //string solverName = "NLGS";      // solver algorithm used for non-smooth problem
    //string solverName = "PGS";      // solver algorithm used for non-smooth problem
    //string solverName = "Lemke";      // solver algorithm used for non-smooth problem

    double e  = 0.9;                  // nslaw
    double e2 = 0.9;                  // nslaw2
    double mu = 0.1;


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

      //   (*(q0[i]))(2) =  0.12 + 0.22*i;
      //     (*(q0[i]))(1) =  0.12 + 0.1*i;

    }


    //    // set values
    //   (*(q0[0]))(0) =  0.0;    (*(q0[0]))(1) =  0.3;  (*(q0[0]))(2) =  0.35;
    //   (*(q0[1]))(0) =  0.0;    (*(q0[1]))(1) =  0.3;  (*(q0[1]))(2) =  0.12;

    //   (*(q0[0]))(0) =  0.0;    (*(q0[0]))(1) =  0.;  (*(q0[0]))(2) =  0.12;
    //     (*(q0[1]))(0) =  0.0;    (*(q0[1]))(1) =  0.3;  (*(q0[1]))(2) =  0.12;

    //   (*(q0[0]))(0) =  0.0;    (*(q0[0]))(1) =  0.4;  (*(q0[0]))(2) =  0.1;
    //     (*(v0[0]))(0) =  0.;    (*(v0[0]))(1) =  -1.;  (*(v0[0]))(2) =  0.;

    //     (*(q0[1]))(0) =  0.0;    (*(q0[1]))(1) =  0.;  (*(q0[1]))(2) =  0.1;
    //     (*(v0[1]))(0) =  0.;    (*(v0[1]))(1) =  1.;  (*(v0[1]))(2) =  0.;

    //     (*(q0[2]))(0) =  0.;    (*(q0[2]))(1) =  0.8;  (*(q0[2]))(2) =  0.1;
    //     (*(q0[3]))(0) =  0.;    (*(q0[3]))(1) =  -0.5;  (*(q0[3]))(2) =  0.1;


    //   (*(q0[0]))(0) =  0.0;    (*(q0[0]))(1) =  0.0;  (*(q0[0]))(2) =  0.2;
    //      (*(q0[1]))(0) =  0.0;    (*(q0[1]))(1) =  0.;  (*(q0[1]))(2) =  0.5;
    //      (*(v0[1]))(3) =  10.;
    //      (*(q0[2]))(0) =  0.0;    (*(q0[2]))(1) =  0.;  (*(q0[2]))(2) =  0.6;


    // carreau
    //  (*(q0[0]))(0) =  0.0;    (*(q0[0]))(1) =  0.0;  (*(q0[0]))(2) =  0.1;
    //       (*(q0[1]))(0) =  0.5;    (*(q0[1]))(1) =  0.5;  (*(q0[1]))(2) =  0.3;
    //       (*(v0[1]))(0) = -1.;    (*(v0[1]))(1) = -1.;  (*(v0[1]))(2) =  1.3;


    //    // billard

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

    (*(v0[8]))(0) = -1;
    (*(v0[8]))(1) = -20.;

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


    //Cube de billes

    //  (*(q0[0]))(0) =  0.2;    (*(q0[0]))(1) = -0.2;  (*(q0[0]))(2) =  0.2;
    //     (*(q0[1]))(0) =  0.2;    (*(q0[1]))(1) =  0.2;  (*(q0[1]))(2) =  0.2;
    //     (*(q0[2]))(0) = -0.2;    (*(q0[2]))(1) =  0.2;  (*(q0[2]))(2) =  0.2;
    //     (*(q0[3]))(0) = -0.2;    (*(q0[3]))(1) = -0.2;  (*(q0[3]))(2) =  0.2;
    //     (*(q0[4]))(0) =  0.25;   (*(q0[4]))(1) = -0.2;  (*(q0[4]))(2) =  0.4;
    //     (*(q0[5]))(0) =  0.25;   (*(q0[5]))(1) =  0.2;  (*(q0[5]))(2) =  0.4;
    //     (*(q0[6]))(0) = -0.25;   (*(q0[6]))(1) =  0.2;  (*(q0[6]))(2) =  0.4;
    //     (*(q0[7]))(0) = -0.25;   (*(q0[7]))(1) = -0.2;  (*(q0[7]))(2) =  0.4;

    //     (*(q0[8]))(0) =  0.;     (*(q0[8]))(1) = 0.3;    (*(q0[8]))(2) =  0.1;

    //     (*(v0[8]))(1) =  -10;

    //     (*(q0[9]))(0) =  0.2;    (*(q0[9]))(1) =  0.2;  (*(q0[9]))(2) =  0.6;
    //     (*(q0[10]))(0)= -0.2;    (*(q0[10]))(1)=  0.2;  (*(q0[10]))(2)=  0.6;
    //     (*(q0[11]))(0)= -0.2;    (*(q0[11]))(1)= -0.2;  (*(q0[11]))(2)=  0.6;
    //     (*(q0[12]))(0)=  0.25;   (*(q0[12]))(1)= -0.2;  (*(q0[12]))(2)=  0.8;
    //     (*(q0[13]))(0)=  0.25;   (*(q0[13]))(1)=  0.2;  (*(q0[13]))(2)=  0.8;
    //     (*(q0[14]))(0)= -0.25;   (*(q0[14]))(1)=  0.2;  (*(q0[14]))(2)=  0.9;
    //     (*(q0[15]))(0)= -0.25;   (*(q0[15]))(1)= -0.2;  (*(q0[15]))(2)=  0.9;
    //     (*(q0[16]))(0)=  0.35;   (*(q0[16]))(1)=-0.35;  (*(q0[16]))(2)=  0.1;
    //     (*(q0[17]))(0)=  0.35;   (*(q0[17]))(1)= 0.35;  (*(q0[17]))(2)=  0.1;
    //     (*(q0[18]))(0)= -0.35;   (*(q0[18]))(1)= 0.35;  (*(q0[18]))(2)=  0.1;
    //     (*(q0[19]))(0)= -0.35;   (*(q0[19]))(1)=-0.35;  (*(q0[19]))(2)=  0.1;

    //    (*(q0[0]))(0) =  0.0;    (*(q0[0]))(1) =  0.3;  (*(q0[0]))(2) =  0.1;
    //     (*(v0[0]))(0) =  0.;    (*(v0[0]))(1) =  -1.;  (*(v0[0]))(2) =  0.12;

    //     (*(q0[1]))(0) =  0.0;    (*(q0[1]))(1) =  0.;  (*(q0[1]))(2) =  0.1;
    //     (*(v0[1]))(0) =  0.;    (*(v0[1]))(1) =  1.;  (*(v0[1]))(2) =  0.;

    //    // 25*etage beads in cube
    //     unsigned int cote  = 5;
    //     unsigned int etage = 4;
    //     unsigned int k = 0;
    //     for (j=0;j<etage;j++){
    //       for (k=0;k<cote;k++) {
    //  for (i=0;i<cote;i++) {
    //    if (j % 2 == 0){
    //      (*(q0[k*cote+i+j*cote*cote]))(0) = -0.6 + 3*k*R; (*(q0[k*cote+i+j*cote*cote]))(1) = -0.6 + 3*i*R; (*(q0[k*cote+i+j*cote*cote]))(2) = 0.11+(2*j)*R;
    //    }
    //    else{
    //      (*(q0[k*cote+i+j*cote*cote]))(0) = -0.6 + 3*k*R; (*(q0[k*cote+i+j*cote*cote]))(1) = -0.6 + 3*i*R+R/4;(*(q0[k*cote+i+j*cote*cote]))(2) = 0.11+(2*j)*R;
    //    }
    //  }
    //       }
    //     }

    //   (*(q0[DSNUMBER-2]))(0)= -0.8;   (*(q0[DSNUMBER-2]))(1)= -0.8;   (*(q0[DSNUMBER-2]))(2)=  (*(q0[DSNUMBER-3]))(2);    (*(v0[DSNUMBER-2]))(1) =  -1;(*(v0[DSNUMBER-2]))(2) =  -2;
    //     (*(q0[DSNUMBER-1]))(0)= -0.8;  (*(q0[DSNUMBER-1]))(1)= -0.8;  (*(q0[DSNUMBER-1]))(2)=  0.1;  (*(v0[DSNUMBER-1]))(0) =  10;(*(v0[DSNUMBER-1]))(1) =  10;


    // un grand merci

    //      for (i=0;i<DSNUMBER;i++)
    //  (*(q0[i]))(0) = 0.;

    //       // le M
    //       (*(q0[0]))(1) = -0.9;    (*(q0[0]))(2) =  0.15;
    //       (*(q0[1]))(1) = -0.9;    (*(q0[1]))(2) =  0.25;
    //       (*(q0[2]))(1) = -0.9;    (*(q0[2]))(2) =  0.35;
    //       (*(q0[3]))(1) = -0.9;    (*(q0[3]))(2) =  0.45;
    //       (*(q0[4]))(1) = -0.9;    (*(q0[4]))(2) =  0.55;
    //       //pas bouger
    //       (*(q0[31]))(1) = -0.8;    (*(q0[31]))(2) =  0.40;
    //       (*(q0[32]))(1) = -0.7;    (*(q0[32]))(2) =  0.32;
    //       (*(q0[33]))(1) = -0.6;    (*(q0[33]))(2) =  0.40;

    //       (*(q0[5]))(1) = -0.49;   (*(q0[5]))(2) =  0.15;
    //       (*(q0[6]))(1) = -0.49;   (*(q0[6]))(2) =  0.25;
    //       (*(q0[7]))(1) = -0.49;   (*(q0[7]))(2) =  0.35;
    //       (*(q0[8]))(1) = -0.49;   (*(q0[8]))(2) =  0.45;
    //       (*(q0[9]))(1) = -0.49;   (*(q0[9]))(2) =  0.55;

    //       // le E
    //       (*(q0[10]))(1)= -0.31;   (*(q0[10]))(2) =  0.25;
    //       (*(q0[11]))(1)= -0.31;   (*(q0[11]))(2) =  0.35;
    //       (*(q0[12]))(1)= -0.31;   (*(q0[12]))(2) =  0.45;
    //       (*(q0[13]))(1)= -0.31;   (*(q0[13]))(2) =  0.55;
    //       (*(q0[14]))(1)= -0.31;   (*(q0[14]))(2) =  0.65;
    //       //pas bouger
    //       (*(q0[34]))(1)= -0.2;    (*(q0[34]))(2) =  0.45;
    //       (*(q0[35]))(1)= -0.1;    (*(q0[35]))(2) =  0.45;
    //       (*(q0[36]))(1)= -0.2;    (*(q0[36]))(2) =  0.25;
    //       (*(q0[37]))(1)= -0.2;    (*(q0[37]))(2) =  0.05;
    //       (*(q0[38]))(1)= -0.1;    (*(q0[38]))(2) =  0.05;

    //       // le R
    //       (*(q0[15]))(1)=  0.09;   (*(q0[15]))(2) =  0.25;
    //       (*(q0[16]))(1)=  0.09;   (*(q0[16]))(2) =  0.35;
    //       (*(q0[17]))(1)=  0.09;   (*(q0[17]))(2) =  0.45;
    //       (*(q0[18]))(1)=  0.09;   (*(q0[18]))(2) =  0.55;
    //       (*(q0[19]))(1)=  0.09;   (*(q0[19]))(2) =  0.65;
    //       //pas bouger
    //       (*(q0[39]))(1)=  0.2;    (*(q0[39]))(2) =  0.45;
    //       (*(q0[40]))(1)=  0.3;    (*(q0[40]))(2) =  0.45;
    //       (*(q0[41]))(1)=  0.3;    (*(q0[41]))(2) =  0.35;
    //       (*(q0[42]))(1)=  0.2;    (*(q0[42]))(2) =  0.25;
    //       (*(q0[43]))(1)=  0.3;    (*(q0[43]))(2) =  0.15;
    //       (*(q0[44]))(1)=  0.3;    (*(q0[44]))(2) =  0.05;

    //       // le C
    //       (*(q0[20]))(1)=  0.49;    (*(q0[20]))(2) =  0.25;
    //       (*(q0[21]))(1)=  0.49;    (*(q0[21]))(2) =  0.35;
    //       (*(q0[22]))(1)=  0.49;    (*(q0[22]))(2) =  0.45;
    //       (*(q0[23]))(1)=  0.49;    (*(q0[23]))(2) =  0.55;
    //       (*(q0[24]))(1)=  0.49;    (*(q0[24]))(2) =  0.65;
    //       // pas bouger
    //       (*(q0[45]))(1)=  0.6;    (*(q0[45]))(2) =  0.45;
    //       (*(q0[46]))(1)=  0.7;    (*(q0[46]))(2) =  0.45;
    //       (*(q0[47]))(1)=  0.6;    (*(q0[47]))(2) =  0.05;
    //       (*(q0[48]))(1)=  0.7;    (*(q0[48]))(2) =  0.05;

    //       // le I
    //       (*(q0[25]))(1)=  0.9;    (*(q0[25]))(2) =  0.25;
    //       (*(q0[26]))(1)=  0.9;    (*(q0[26]))(2) =  0.35;
    //       (*(q0[27]))(1)=  0.9;    (*(q0[27]))(2) =  0.45;
    //       (*(q0[28]))(1)=  0.9;    (*(q0[28]))(2) =  0.55;
    //       (*(q0[29]))(1)=  0.9;    (*(q0[29]))(2) =  0.65;

    //       // le point du i

    //       (*(q0[30]))(1)=  0.9;    (*(q0[30]))(2) =  0.9;

    // // un autre merci

    //       for (i=0;i<DSNUMBER;i++)
    //  (*(q0[i]))(0) = 0.;

    //       // le M
    //       (*(q0[0]))(1) = -0.9;    (*(q0[0]))(2) =  0.05;
    //       (*(q0[1]))(1) = -0.9;    (*(q0[1]))(2) =  0.15;
    //       (*(q0[2]))(1) = -0.9;    (*(q0[2]))(2) =  0.25;
    //       (*(q0[3]))(1) = -0.9;    (*(q0[3]))(2) =  0.35;
    //       (*(q0[4]))(1) = -0.9;    (*(q0[4]))(2) =  0.45;

    //       (*(q0[5]))(1) = -0.8;    (*(q0[5]))(2) =  0.4;
    //       (*(q0[6]))(1) = -0.7;    (*(q0[6]))(2) =  0.3;
    //       (*(q0[7]))(1) = -0.6;    (*(q0[7]))(2) =  0.4;

    //       (*(q0[8]))(1) = -0.5;    (*(q0[8]))(2) =  0.05;
    //       (*(q0[9]))(1) = -0.5;    (*(q0[9]))(2) =  0.15;
    //       (*(q0[10]))(1)= -0.5;    (*(q0[10]))(2)=  0.25;
    //       (*(q0[11]))(1)= -0.5;    (*(q0[11]))(2)=  0.35;
    //       (*(q0[12]))(1)= -0.5;    (*(q0[12]))(2)=  0.45;

    //       // le E
    //       (*(q0[13]))(1)= -0.3;    (*(q0[13]))(2) =  0.05;
    //       (*(q0[14]))(1)= -0.3;    (*(q0[14]))(2) =  0.15;
    //       (*(q0[15]))(1)= -0.3;    (*(q0[15]))(2) =  0.25;
    //       (*(q0[16]))(1)= -0.3;    (*(q0[16]))(2) =  0.35;
    //       (*(q0[17]))(1)= -0.3;    (*(q0[17]))(2) =  0.45;

    //       (*(q0[18]))(1)= -0.2;    (*(q0[18]))(2) =  0.45;
    //       (*(q0[19]))(1)= -0.1;    (*(q0[19]))(2) =  0.45;
    //       (*(q0[20]))(1)= -0.2;    (*(q0[20]))(2) =  0.25;
    //       (*(q0[21]))(1)= -0.2;    (*(q0[21]))(2) =  0.05;
    //       (*(q0[22]))(1)= -0.1;    (*(q0[22]))(2) =  0.05;

    //       // le R
    //       (*(q0[23]))(1)=  0.1;    (*(q0[23]))(2) =  0.05;
    //       (*(q0[24]))(1)=  0.1;    (*(q0[24]))(2) =  0.15;
    //       (*(q0[25]))(1)=  0.1;    (*(q0[25]))(2) =  0.25;
    //       (*(q0[26]))(1)=  0.1;    (*(q0[26]))(2) =  0.35;
    //       (*(q0[27]))(1)=  0.1;    (*(q0[27]))(2) =  0.45;

    //       (*(q0[28]))(1)=  0.2;    (*(q0[28]))(2) =  0.45;
    //       (*(q0[29]))(1)=  0.3;    (*(q0[29]))(2) =  0.45;
    //       (*(q0[30]))(1)=  0.3;    (*(q0[30]))(2) =  0.35;
    //       (*(q0[31]))(1)=  0.2;    (*(q0[31]))(2) =  0.25;
    //       (*(q0[32]))(1)=  0.3;    (*(q0[32]))(2) =  0.15;
    //       (*(q0[33]))(1)=  0.3;    (*(q0[33]))(2) =  0.05;

    //       // le C
    //       (*(q0[34]))(1)=  0.5;    (*(q0[34]))(2) =  0.05;
    //       (*(q0[35]))(1)=  0.5;    (*(q0[35]))(2) =  0.15;
    //       (*(q0[36]))(1)=  0.5;    (*(q0[36]))(2) =  0.25;
    //       (*(q0[37]))(1)=  0.5;    (*(q0[37]))(2) =  0.35;
    //       (*(q0[38]))(1)=  0.5;    (*(q0[38]))(2) =  0.45;

    //       (*(q0[39]))(1)=  0.6;    (*(q0[39]))(2) =  0.45;
    //       (*(q0[40]))(1)=  0.7;    (*(q0[40]))(2) =  0.45;
    //       (*(q0[41]))(1)=  0.6;    (*(q0[41]))(2) =  0.05;
    //       (*(q0[42]))(1)=  0.7;    (*(q0[42]))(2) =  0.05;

    //       // le I
    //       (*(q0[43]))(1)=  0.9;    (*(q0[43]))(2) =  0.05;
    //       (*(q0[44]))(1)=  0.9;    (*(q0[44]))(2) =  0.15;
    //       (*(q0[45]))(1)=  0.9;    (*(q0[45]))(2) =  0.25;
    //       (*(q0[46]))(1)=  0.9;    (*(q0[46]))(2) =  0.35;
    //       (*(q0[47]))(1)=  0.9;    (*(q0[47]))(2) =  0.45;

    //       // le point du i

    //       (*(q0[48]))(1)=  0.9;    (*(q0[48]))(2) =  0.9;

    // pour le merci

    //   for (i=0;i<DSNUMBER;i++)
    //       {
    //  GLOB_tabLDS[i]=new LagrangianDS(i,*(q0[i]),*(v0[i]),*Mass);
    //  checkDS = allDS.insert(GLOB_tabLDS[i]);
    //  if (i < 31){
    //  (static_cast<LagrangianDS*>(*(checkDS.first)))->setComputeFExtFunction("3DDrawPlugin.so", "gravity");
    //  }
    //       }

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

    if (obst_y_p)
    {
      SiconosVector *b2 = new SimpleVector(3);
      (*b2)(0) = WALL - R;
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

    if (obst_y_m)
    {
      SiconosVector *b2_ = new SimpleVector(3);
      (*b2_)(0) = WALL - R;
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

    if (obst_x_p)
    {
      SiconosVector *b3 = new SimpleVector(3);
      (*b3)(0) = WALL - R;
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

    if (obst_x_m)
    {
      SiconosVector *b3_ = new SimpleVector(3);
      (*b3_)(0) = WALL - R;
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
    //      NonSmoothLaw * nslaw2 = new NewtonImpactFrictionNSL(e2, e2, mu, 3);

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
    OneStepIntegrator * OSI = new Moreau(allDS , 0.5000001 , GLOB_SIM);

    // -- OneStepNsProblem --


    //OneStepNSProblem * osnspb = new LCP(GLOB_SIM ,"FrictionContact3D",solverName,101, 0.001);

    OneStepNSProblem * osnspb = new FrictionContact3D(GLOB_SIM , "FrictionContact3D", solverName, 1000001, 0.001);

    cout << "=== End of model loading === " << endl;

    // =========================== End of model definition

    // ================================= Computation =================================

    // --- Simulation initialization ---

    GLOB_SIM->initialize();

    cout << "End of simulation initialisation" << endl;
    //    N = GLOB_T->getNSteps(); // Number of time steps

    dataPlot(k_iter, 0) = k_iter * GLOB_T->getH();
    dataPlot(k_iter, 1) = GLOB_tabLDS[0]->getQ()(2);
    dataPlot(k_iter, 2) = GLOB_tabLDS[0]->getVelocity()(2);
    // dataPlot(k_iter,3) = (multiBeads->getNonSmoothDynamicalSystemPtr()->getInteractionPtr(0)->getLambda(1))(0);
    dataPlot(k_iter, 3) = GLOB_tabLDS[1]->getQ()(2);
    dataPlot(k_iter, 4) = GLOB_tabLDS[1]->getVelocity()(2);
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
      dataPlot(k_iter, 0) = k_iter * GLOB_T->getH();
      dataPlot(k_iter, 1) = GLOB_tabLDS[0]->getQ()(2);
      dataPlot(k_iter, 2) = GLOB_tabLDS[0]->getVelocity()(2);
      //    dataPlot(k_iter,3) = (multiBeads->getNonSmoothDynamicalSystemPtr()->getInteractionPtr(0)->getLambda(1))(0);
      dataPlot(k_iter, 3) = GLOB_tabLDS[1]->getQ()(2);
      dataPlot(k_iter, 4) = GLOB_tabLDS[1]->getVelocity()(2);
      //dataPlot(k_iter,7) = (multiBeads->getNonSmoothDynamicalSystemPtr()->getInteractionPtr(1)->getLambda(1))(0);
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
    cout << "Exception caught in \'sample/MultiBeadsColumn\'" << endl;
  }



}


