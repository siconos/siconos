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

/**
   Billiard table example using QGLViewer for drawing.

   Use:
   siconos -DWithQGLViewer main.cpp



 */
#ifdef WithQGLViewer
#include "DrawSiconos.h"
#include <qapplication.h>
#else
#include "BilliardModel.h"
#include "environment.h"
#endif

using namespace std;

int main(int argc, char* argv[])
{

  // 3D display using QGLViewer
#ifdef WithQGLViewer
  QApplication application(argc, argv);

  Viewer viewer;

#if QT_VERSION < 0x040000
  application.setMainWidget(&viewer);
#else
  viewer.setWindowTitle("animation");
#endif

  viewer.show();

  return application.exec();

#else // No viewer

  try
  {

    SP::BilliardModel column(new BilliardModel(4));
    column->initialize();
    // --- simulation solver ---
    while (!column->isSimulationFinished())
      column->compute();
    column->end();

  }
  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught in DrawSiconos for beads column." << endl;
  }


  return 0;
#endif

}
