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

#include "BulletViewer.hpp"
#include "BulletBodies.hpp"

#include "BulletSpaceFilter.hpp"

using namespace qglviewer;

btVector3 normalize(const btVector3& v)
{
  return (v / sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]));
};


void BulletViewer::init()
{
  // viewer and scene state
  restoreStateFromFile();
  setSceneRadius(20.);
  showEntireScene();

  // color
  setBackgroundColor(QColor(200, 200, 200));
  setForegroundColor(QColor(0, 0, 0));

  // 2D setup
  glDisable(GL_LIGHTING);
  camera()->setPosition(Vec(0.0, 0.0, 1.0));
  //  camera()->lookAt(sceneCenter());
  //  camera()->setType(Camera::ORTHOGRAPHIC);
  camera()->showEntireScene();
  setManipulatedFrame(camera()->frame());

  //  constraint_.reset(new WorldConstraint());
  //  constraint_->setRotationConstraintType(AxisPlaneConstraint::FORBIDDEN);
  //  camera()->frame()->setConstraint(&*constraint_);

  // help screen
  help();

  // siconos setup
  Siconos_.reset(new BulletBodies());

  Siconos_->init();

  ask<ForCollisionWorld>(*Siconos_->spaceFilter())->setDebugDrawer(&_debugDrawer);


  int i;
  DSIterator itDS;

  DynamicalSystemsGraph::VIterator dsi, dsend;
  boost::tie(dsi, dsend) = GETALLDS(Siconos_)->vertices();
  for (i = 0; dsi != dsend; ++i, ++dsi)
  {
    SP::DynamicalSystem ds = GETALLDS(Siconos_)->bundle(*dsi);
    insertQGLShape(ask<ForShape>(*ds), ds);
  }

  lastSelected_ = -1;

  initUCircle();

  myMouseBehavior_ = false;

  startAnimation();

}




void BulletViewer::draw()
{

  int i, mrow;
  clock_t lcptime;

  char qs[6];

  float lbd, lbdmax, w, x1, y1, x2, y2, d, r;

  DSIterator itDS;
  SP::DynamicalSystemsSet involvedDS;
  SP::UnitaryRelationsGraph I1 = Siconos_->model()->simulation()->indexSet(1);
  SP::Interaction interaction;
  SP::Relation relation;


  // calibration
  UnitaryRelationsGraph::VIterator ui, uiend;
  for (boost::tie(ui, uiend) = I1->vertices(); ui != uiend; ++ui)
  {
    lbdmax = fmax(I1->bundle(*ui)->interaction()->lambdaOld(1)->getValue(0), lbdmax);
  }

  for (boost::tie(ui, uiend) = I1->vertices(); ui != uiend; ++ui)
  {
    interaction = I1->bundle(*ui)->interaction();
    relation = interaction->relation();

    lbd = interaction->lambdaOld(1)->getValue(0);

    // screen width of interaction
    w = lbd / (2 * fmax(lbdmax, 1.)) + .03;

    involvedDS = interaction->dynamicalSystems();

    // disk/disk
    itDS = involvedDS->begin();

    SP::DynamicalSystem d1 = *itDS;

    SP::SiconosVector q1 = ask<ForPosition>(*d1);

    float x1 = (*q1)(0);
    float y1 = (*q1)(1);
    float r1 = ask<ForRadius>(*d1);


    if (involvedDS->size() == 2)
    {
      SimpleVector& cf = *ask<ForContactForce>(*relation);

      btManifoldPoint& cpoint = *ask<ForContactPoints>(*relation);

      btVector3 posa = cpoint.getPositionWorldOnA();
      btVector3 posb = cpoint.getPositionWorldOnB();
      btVector3 dirf(cf(0), cf(1), cf(2));
      btVector3 endf = posa + normalize(dirf) / 2.;
      btVector3 cnB = posb + normalize(cpoint.m_normalWorldOnB) / 2.;


      glPushMatrix();
      glColor3f(.80, 0, 0);
      glLineWidth(2.);
      QGLViewer::drawArrow(qglviewer::Vec(posa[0], posa[1], posa[2]), qglviewer::Vec(endf[0], endf[1], endf[2]), .05, 10.);
      glPopMatrix();

      glPushMatrix();
      glColor3f(0, .80, 0);
      glLineWidth(2.);
      QGLViewer::drawArrow(qglviewer::Vec(posb[0], posb[1], posb[2]), qglviewer::Vec(cnB[0], cnB[1], cnB[2]), .05, 10.);
      glPopMatrix();

    }

    else
    {
      SimpleVector& cf = *ask<ForContactForce>(*relation);

      btManifoldPoint& cpoint = *ask<ForContactPoints>(*relation);

      btVector3 posa = cpoint.getPositionWorldOnA();
      btVector3 posb = cpoint.getPositionWorldOnB();
      btVector3 dirf(cf(0), cf(1), cf(2));
      btVector3 endf = posa + normalize(dirf) / 2.;
      btVector3 cnB = posb + normalize(cpoint.m_normalWorldOnB) / 2.;

      glPushMatrix();
      glColor3f(.80, 0, 0);
      glLineWidth(2.);
      QGLViewer::drawArrow(qglviewer::Vec(posa[0], posa[1], posa[2]), qglviewer::Vec(endf[0], endf[1], endf[2]), .05, 10.);
      glPopMatrix();


      glPushMatrix();
      glColor3f(0, .80, 0);
      glLineWidth(2.);
      QGLViewer::drawArrow(qglviewer::Vec(posb[0], posb[1], posb[2]), qglviewer::Vec(cnB[0], cnB[1], cnB[2]), .05, 10.);
      glPopMatrix();
    }
  }

  for (i = 0; i < GETNDS(Siconos_); i++)
  {
    //    if (shapes_[i]->selected())
    //    {
    //      drawSelectedQGLShape(*shapes_[i]);
    //    }
    //    else
    //    {
    drawQGLShape(*shapes_[i]);
    //   }
  }

  // draw static objects
  std::vector<SP::btCollisionObject>& staticObjects = *(ask<ForStaticObjects>(*Siconos_->spaceFilter()));

  for (std::vector<SP::btCollisionObject>::iterator iso = staticObjects.begin();
       iso != staticObjects.end(); ++iso)
  {
    ask<ForCollisionWorld>(*Siconos_->spaceFilter())
    ->debugDrawObject((*iso)->getWorldTransform(),
                      (*iso)->getCollisionShape(),
                      btVector3(1, 1, 0));
  };
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
      double C = (*Siconos_->plans())(i, 2);
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
      double xc, yc;

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


void BulletViewer::animate()
{

  Siconos_->compute();
  //Siconos_->compute();
  //saveSnapshot();
}


void BulletViewer::mousePressEvent(QMouseEvent* e)
{
  if ((e->button() == Qt::LeftButton) && (e->modifiers() == Qt::NoButton))
    myMouseBehavior_ = true;
  else
    QGLViewer::mousePressEvent(e);
}

void BulletViewer::mouseMoveEvent(QMouseEvent *e)
{

  float coor[3];

  if (selectedName() >= 0)
  {
    lastSelected_ = selectedName();
    coor[0] = (float) e->x();
    coor[1] = (float) e->y();
    coor[2] = 0;

    camera()->getUnprojectedCoordinatesOf(coor, coor, NULL);

    SP::SiconosVector fext =
      ask<ForFExt>(*shapes_[selectedName()]->DS());

    SP::SiconosVector q =
      ask<ForPosition>(*shapes_[selectedName()]->DS());

    double massValue =
      ask<ForMassValue>(*shapes_[selectedName()]->DS());

    fext->setValue(0, ((float)coor[0] - q->getValue(0))*massValue * 30);

    fext->setValue(1, ((float)coor[1] - q->getValue(1))*massValue * 30);
  }

  QGLViewer::mouseMoveEvent(e);
}

void BulletViewer::mouseReleaseEvent(QMouseEvent* e)
{
  if (myMouseBehavior_)
    myMouseBehavior_ = false;

  if (lastSelected_ >= 0)
  {
    shapes_[lastSelected_]->restoreFExt();
    shapes_[lastSelected_]->nextSelection();
  };

  QGLViewer::mouseReleaseEvent(e);
}


QString BulletViewer::helpString() const
{
  QString text("<h2>S i m p l e V i e w e r</h2>");
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
  text += "See the <b>Mouse</b> tab and the documentation web pages for details.<br><br>";
  text += "Press <b>Escape</b> to exit the viewer.";
  return text;
}


