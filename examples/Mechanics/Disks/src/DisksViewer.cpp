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

#include "DisksViewer.hpp"

//#include <Siconos/io/Register.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>

// FMatrix is forwarded
#include <SpaceFilter_impl.hpp>

using namespace qglviewer;

void DisksViewer::init()
{

  // siconos setup
  Siconos_.reset(new Disks());

  BodiesViewer::init();

  // 2D setup
  glDisable(GL_LIGHTING);
  camera()->setPosition(Vec(0.0, 0.0, 1.0));
  camera()->lookAt(sceneCenter());
  camera()->setType(Camera::ORTHOGRAPHIC);
  camera()->showEntireScene();
  constraint_.reset(new WorldConstraint());
  constraint_->setRotationConstraintType(AxisPlaneConstraint::FORBIDDEN);
  camera()->frame()->setConstraint(&*constraint_);

}




void DisksViewer::draw()
{

  int i;

  char qs[6];

  float lbd, w;

  float lbdmax = 0.;

  DSIterator itDS;
  SP::InteractionsGraph I1;
  SP::Interaction interaction;
  SP::Relation relation;

  if (Siconos_->model()->nonSmoothDynamicalSystem()->topology()->numberOfIndexSet() > 1)
  {
    I1 = Siconos_->model()->simulation()->indexSet(1);

    // calibration
    InteractionsGraph::VIterator ui, uiend;
    for (boost::tie(ui, uiend) = I1->vertices(); ui != uiend; ++ui)
    {
      lbdmax = fmax(I1->bundle(*ui)->lambdaOld(1)->getValue(0), lbdmax);
    }

    for (boost::tie(ui, uiend) = I1->vertices(); ui != uiend; ++ui)
    {
      interaction = I1->bundle(*ui);
      relation = interaction->relation();
      
      lbd = interaction->lambdaOld(1)->getValue(0);

      // screen width of interaction
      w = lbd / (2 * fmax(lbdmax, 1.)) + .03;
   
      // disk/disk
      
      SP::DynamicalSystem d1 = I1->properties(*ui).source;
      SP::DynamicalSystem d2 = I1->properties(*ui).target;

      SP::SiconosVector q1 = ask<ForPosition>(*d1);

      float x1 = (*q1)(0);
      float y1 = (*q1)(1);
      float r1 = ask<ForRadius>(*d1);


      if (d1 != d2)
      {
        SP::SiconosVector q2 = ask<ForPosition>(*d2);
        float x2 = (*q2)(0);
        float y2 = (*q2)(1);
        float r2 = ask<ForRadius>(*d2);

        float d = hypotf(x1 - x2, y1 - y2);

        glPushMatrix();

        glColor3f(.0f, .0f, .0f);
        drawRec(x1, y1, x1 + (x2 - x1)*r1 / d, y1 + (y2 - y1)*r1 / d, w);
        drawRec(x2, y2, x2 + (x1 - x2)*r2 / d, y2 + (y1 - y2)*r2 / d, w);

        glPopMatrix();
      }

      else
      {
        SP::SiconosMatrix jachq = ask<ForJachq>(*relation);
        double jx = jachq->getValue(0, 0);
        double jy = jachq->getValue(0, 1);
        double dj = hypot(jx, jy);

        glPushMatrix();

        glColor3f(.0f, .0f, .0f);
        drawRec(x1, y1, x1 - r1 * jx / dj, y1 - r1 * jy / dj, w);
        glPopMatrix();
      }
    }
  }


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

  if (Siconos_->plans())
  {
    for (unsigned int i = 0 ; i < Siconos_->plans()->size(0) ; ++i)
    {
      double A = (*Siconos_->plans())(i, 0);
      double B = (*Siconos_->plans())(i, 1);
      //double C = (*Siconos_->plans())(i,2);
      double xc = (*Siconos_->plans())(i, 3);
      double yc = (*Siconos_->plans())(i, 4);
      double w = fmin(1e10, (*Siconos_->plans())(i, 5));
      double H = hypot(A, B);

      if (w == 0) w = 1e10;

      //      assert ( fabs(A*xc + B*yc + C) <= std::numeric_limits<double>::epsilon() );

      drawVec(xc, yc, xc - 0.5 * w * B / H, yc + 0.5 * w * A / H);
      drawVec(xc, yc, xc + 0.5 * w * B / H, yc - 0.5 * w * A / H);
    }
  }


  if (Siconos_->movingPlans())
  {
    double time = Siconos_->model()->currentTime();
    for (unsigned int i = 0 ; i < Siconos_->movingPlans()->size1() ; ++i)
    {
      double A = (*Siconos_->movingPlans())(i, 0)(time);
      double B = (*Siconos_->movingPlans())(i, 1)(time);
      double C = (*Siconos_->movingPlans())(i, 2)(time);
      double w = 1e10;
      double H = hypot(A, B);
      double xc = 0.;
      double yc = 0.;

      if (fabs(C) > std::numeric_limits<double>::epsilon())
      {
        if (A == 0)
          // By+C=0
        {
          yc = -C / B;
        }
        else if (B == 0)
          // Ax+C=0
        {
          xc = -C / A;
        }
        else
          // Ax+By+C=0
        {
          if (xc != 0)
            yc = - (A * xc + C) / B;
          else
            xc = - (B * yc + C) / A;
        }
      }

      drawVec(xc, yc, xc - 0.5 * w * B / H, yc + 0.5 * w * A / H);
      drawVec(xc, yc, xc + 0.5 * w * B / H, yc - 0.5 * w * A / H);
    }
  }
  glColor3f(.1, .1, .1);
  glLineWidth(4.);
  QGLViewer::drawArrow(qglviewer::Vec(0, 0, .1), qglviewer::Vec(1, 0, .1), .01, 3);
  QGLViewer::drawArrow(qglviewer::Vec(0, 0, .1), qglviewer::Vec(0, 1, .1), .01, 3);

  glLineWidth(1.);
  for (i = -100; i <= 100; i += 5)
  {
    sprintf(qs, "%d", i);
    //    print((float)i,-.8,qs,small_text);
    //print(-.8,(float)i,qs,small_text);
    drawVec((float)i, -.2, (float)i, .2);
    drawVec(-.2, (float)i, .2, (float)i);
  }
  for (i = -100; i <= 100; i++)
  {
    drawVec((float)i, -.1, (float)i, .1);
    drawVec(-.1, (float)i, .1, (float)i);
  }
}


