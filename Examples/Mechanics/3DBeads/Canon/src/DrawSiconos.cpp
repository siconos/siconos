/* Siconos-sample version 2.0.0, Copyright INRIA 2005-2008.
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
 */

// =============================== Multi contact beads simulation ===============================
//
// Keywords: LagrangianLinearDS, LagrangianLinear relation, Moreau TimeStepping, Non-smooth Newton method.
//
// ======================================================================================================

#include "CanonBallsModel.h"
#include "environment.h"
#include "drawstuff.h"

CanonBallsModel* canon;
bool GLOB_COMPUTE;
bool GLOB_STEP;

using namespace std;

void Start()
{
  GLOB_STEP = false;
  GLOB_COMPUTE = false;
  canon = new CanonBallsModel(2);
  canon->initialize();
}

void DrawBox(float alpha)
{
  //  alpha = 0 signifie transparent, alpha = 1 signifie opaque
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glColor4f(0.5, 0.7, 0.9, alpha); // 3 premiers champs = RGB, quatriNhme champ = opacitNi

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
void computeSiconos()
{
  try
  {

    // --- simulation solver ---
    if (!canon->isSimulationFinished())
      canon->compute();
    else
      dsStop();
  }
  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught in \'sample/Canon_beads\'" << endl;
  }
}

void SimuLoop(int pause)
{
  if (GLOB_COMPUTE == true)
    computeSiconos();
  if (GLOB_STEP == true)
    GLOB_COMPUTE = false;

  canon->draw();

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

void End()
{
  canon->end();
}


void run(int argc, char* argv[])
{
  // setup pointers to callback functions
  dsFunctions fn;
  fn.version = DS_VERSION;
  fn.start = &Start;
  fn.step = &SimuLoop;
  fn.command = &Command;
  fn.stop = &End;
  fn.path_to_textures = "./textures/";

  // run simulation
  dsSimulationLoop(argc, argv, 400, 400, &fn);
}

