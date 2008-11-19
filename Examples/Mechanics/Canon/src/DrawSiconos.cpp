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
#ifdef WithQGLViewer
#include "DrawSiconos.h"
#include "environment.h"

using namespace std;
using namespace qglviewer;

// The siconos model
SP::CanonBallsModel model;

void Viewer::init()
{
  restoreStateFromFile();
  setSceneRadius(NBFloors * DEFAULT_radius * 2);
  showEntireScene();
  // color
  setBackgroundColor(QColor(200, 200, 200));
  setForegroundColor(QColor(0, 0, 0));

  camera()->setPosition(Vec(0.0, -1.0, 0.1));
  camera()->lookAt(sceneCenter());
  camera()->setType(Camera::ORTHOGRAPHIC);
  camera()->showEntireScene();

  // Model construction and initialisation
  model.reset(new CanonBallsModel(NBFloors));
  model->initialize();

  glPointSize(1.0);
  //setGridIsDrawn();
  help();
  // uncomment next line to start automatically the animation
  //startAnimation();
}

void Viewer::draw()
{
  model->draw();
}

void Viewer::animate()
{
  try
  {
    // --- simulation solver ---
    if (!model->isSimulationFinished())
      model->compute();
    else
    {
      stopAnimation();
      model->end();
    }
  }
  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught in DrawSiconos for beads model." << endl;
  }
}

QString Viewer::helpString() const
{
  QString text("<h2>Model of beads</h2>");
  text += "Press:<br>";
  text += " <b>Return</b> to start/stop the animation.<br>";
  text += " <b>A</b> to display axes.<br>";
  text += " <b>Esc</b> to quit.<br>";
  text += "Use the mouse to move the camera around the object. ";
  text += "You can respectively revolve around, zoom and translate with the three mouse buttons. ";
  text += "Left and middle buttons pressed together rotate around the camera view direction axis<br><br>";
  text += "Pressing <b>Alt</b> and one of the function keys (<b>F1</b>..<b>F12</b>) defines a camera keyFrame. ";
  text += "Simply press the function key again to restore it. Several keyFrames define a ";
  text += "camera path. Paths are saved when you quit the application and restored at next start.<br><br>";
  text += "Press <b>F</b> to display the frame rate, <b>A</b> for the world axis, ";
  text += "<b>Alt+Return</b> for full screen mode and <b>Control+S</b> to save a snapshot. ";
  text += "See the <b>Keyboard</b> tab in this window for a complete shortcut list.<br><br>";
  text += "Double clicks automates single click actions: A left button double click aligns the closer axis with the camera (if close enough). ";
  text += "A middle button double click fits the zoom of the camera and the right button re-centers the scene.<br><br>";
  text += "A left button double click while holding right button pressed defines the camera <i>Revolve Around Point</i>. ";
  text += "See the <b>Mouse</b> tab and the documentation of QGLViewer on web pages for details.<br><br>";
  return text;
}


#endif
