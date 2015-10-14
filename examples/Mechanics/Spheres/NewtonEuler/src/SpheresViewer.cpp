/* Siconos-sample version 3.1.0, Copyright INRIA 2005-2009.
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

#include "SpheresViewer.hpp"

#include <Model.hpp>
#include <NonSmoothDynamicalSystem.hpp>

using namespace qglviewer;

void SpheresViewer::init()
{
  Siconos_.reset(new Spheres());

  BodiesViewer::init();

  setSceneRadius(500.);

  camera()->setPosition(qglviewer::Vec(100.0, 100.0, 100.0));
  camera()->setOrientation(1., 1.);

}


void SpheresViewer::draw()
{

  char qs[6];

  DSIterator itDS;
  SP::DynamicalSystemsSet involvedDS;
  SP::Interaction interaction;
  SP::Relation relation;

  for (unsigned int i = 0; i < GETNDS(Siconos_); i++)
  {
    if (shapes_[i]->selected())
    {
      drawSelectedQGLShape(*shapes_[i]);
    }
    else
    {
      drawQGLShape(*shapes_[i]);
    }
  }

  glColor3f(.45, .45, .45);
  glLineWidth(1.);
  drawGrid(100, 200);
  setGridIsDrawn();

  glColor3f(.1, .1, .3);
  drawVec(-100, 0, 100, 0);
  drawVec(0, -100, 0, 100);

  glColor3f(0, 0, 1);
  glLineWidth(4.);

  for (unsigned int i = 0 ; i < Siconos_->plans()->size(0) ; ++i)
  {
    double A = (*Siconos_->plans())(i, 0);
    double B = (*Siconos_->plans())(i, 1);
    double C = (*Siconos_->plans())(i, 2);
    double D = (*Siconos_->plans())(i, 3);

    double N2 = A * A + B * B + C * C;

    glPushMatrix();
    glLineWidth(1.);
    glTranslatef(-D * A / N2, -D * B / N2, -D * C / N2);
    glMultMatrixd(Quaternion(Vec(0., 0., 1.), Vec(A, B, C)).matrix());
    drawGrid(100, 200);
    glPopMatrix();

  }

  glColor3f(.1, .1, .1);
  glLineWidth(4.);
  QGLViewer::drawArrow(qglviewer::Vec(0, 0, .1), qglviewer::Vec(1, 0, .1), .01, 3);
  QGLViewer::drawArrow(qglviewer::Vec(0, 0, .1), qglviewer::Vec(0, 1, .1), .01, 3);

  glLineWidth(1.);
  for (unsigned int i = -100; i <= 100; i += 5)
  {
    sprintf(qs, "%d", i);
    //    print((float)i,-.8,qs,small_text);
    //print(-.8,(float)i,qs,small_text);
    drawVec((float)i, -.2, (float)i, .2);
    drawVec(-.2, (float)i, .2, (float)i);
  }
  for (unsigned int i = -100; i <= 100; i++)
  {
    drawVec((float)i, -.1, (float)i, .1);
    drawVec(-.1, (float)i, .1, (float)i);
  }
}



