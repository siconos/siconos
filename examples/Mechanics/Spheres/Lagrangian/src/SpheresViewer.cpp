/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
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



