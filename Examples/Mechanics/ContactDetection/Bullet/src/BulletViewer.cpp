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

void BulletViewer::init()
{
  // viewer and scene state
  restoreStateFromFile();
  setSceneRadius(20.);
  showEntireScene();

  // color
  setBackgroundColor(QColor(200, 200, 200));
  setForegroundColor(QColor(0, 0, 0));

  camera()->setPosition(Vec(0.0, 0.0, 1.0));
  camera()->showEntireScene();
  setManipulatedFrame(camera()->frame());

  // shapeDrawer
  _shapeDrawer.enableTexture(true);
  _transparency = 1.f;

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

  btVector3 c1(.91, .1, .21);
  btVector3 c2(.21, .91, .1);
  btVector3 c3(.1, .21, .91);
  btVector3 c4(.21, .1, .91);
  btVector3 c5(.91, .21, .1);
  btVector3 c6(.7, .7, .7);
  btVector3 c7(.99, .99, .99);

  _colors.push_back(c1);
  _colors.push_back(c2);
  _colors.push_back(c3);
  _colors.push_back(c4);
  _colors.push_back(c5);
  _colors.push_back(c6);
  _colors.push_back(c7);

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

void BulletViewer::drawQGLShape(const QGLShape& fig)
{

  DynamicalSystem * ds = fig.DS().get();

  assert(ask<ForShape>(*ds) == BULLET);


  glClear(GL_STENCIL_BUFFER_BIT);
  glEnable(GL_CULL_FACE);


  GLfloat light_ambient[] = { btScalar(0.2), btScalar(0.2), btScalar(0.2), btScalar(1.0) };
  GLfloat light_diffuse[] = { btScalar(1.0), btScalar(1.0), btScalar(1.0), btScalar(1.0) };
  GLfloat light_specular[] = { btScalar(1.0), btScalar(1.0), btScalar(1.0), btScalar(1.0)};
  /*  light_position is NOT default value */
  GLfloat light_position0[] = { btScalar(1.0), btScalar(10.0), btScalar(1.0), btScalar(0.0)};
  GLfloat light_position1[] = { btScalar(-1.0), btScalar(-10.0), btScalar(-1.0), btScalar(0.0) };

  glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
  glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
  glLightfv(GL_LIGHT0, GL_POSITION, light_position0);

  glLightfv(GL_LIGHT1, GL_AMBIENT, light_ambient);
  glLightfv(GL_LIGHT1, GL_DIFFUSE, light_diffuse);
  glLightfv(GL_LIGHT1, GL_SPECULAR, light_specular);
  glLightfv(GL_LIGHT1, GL_POSITION, light_position1);

  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_LIGHT1);


  glShadeModel(GL_SMOOTH);
  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LESS);

  glClearColor(btScalar(0.7), btScalar(0.7), btScalar(0.7), btScalar(0));

  btCollisionObject* co = &*ask<ForCollisionObject>(*ds);
  btCollisionWorld* cow = &*ask<ForCollisionWorld>(*Siconos_->spaceFilter());


  //    ask<ForCollisionWorld>(*Siconos_->spaceFilter())->debugDrawObject(co->getWorldTransform(), co->getCollisionShape(), btVector3(1,1,0));


  btScalar m[16];
  btMatrix3x3 rot;
  co->getWorldTransform().getOpenGLMatrix(m);
  rot = co->getWorldTransform().getBasis();

  btVector3 aabbMin, aabbMax;
  cow->getBroadphase()->getBroadphaseAabb(aabbMin, aabbMax);

  aabbMin -= btVector3(BT_LARGE_FLOAT, BT_LARGE_FLOAT, BT_LARGE_FLOAT);
  aabbMax += btVector3(BT_LARGE_FLOAT, BT_LARGE_FLOAT, BT_LARGE_FLOAT);

  _shapeDrawer.drawOpenGL(_transparency, m, co->getCollisionShape(), _colors[ds->number() % _colors.size()], 0, aabbMin, aabbMax);

  unsigned int numManifolds =
    cow->getDispatcher()->getNumManifolds();

  for (unsigned int i = 0; i < numManifolds; ++i)
  {
    btPersistentManifold* contactManifold =
      cow->getDispatcher()->getManifoldByIndexInternal(i);

    unsigned int numContacts = contactManifold->getNumContacts();

    for (unsigned int j = 0; j < numContacts; ++j)
    {
      btManifoldPoint& pt = contactManifold->getContactPoint(j);
      glBegin(GL_LINES);
      glLineWidth(10.);
      glColor3f(0, 0, 0);

      btVector3 ptA = pt.getPositionWorldOnA();
      btVector3 ptB = pt.getPositionWorldOnB();

      glVertex3d(ptA.x(), ptA.y(), ptA.z());
      glVertex3d(ptB.x(), ptB.y(), ptB.z());
      glEnd();
    }
  }
}


void BulletViewer::draw()
{

  char qs[6];

  float lbdmax, w;

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

    involvedDS = interaction->dynamicalSystems();

    itDS = involvedDS->begin();

    SP::DynamicalSystem d1 = *itDS;

    SP::SiconosVector q1 = ask<ForPosition>(*d1);

    {
      SimpleVector& cf = *ask<ForContactForce>(*relation);
      double cfn = cf.norm2();

      w = fmax(.3, cfn / fmax(lbdmax, cfn));


      btManifoldPoint& cpoint = *ask<ForContactPoints>(*relation);

      btVector3 posa = cpoint.getPositionWorldOnA();
      btVector3 posb = cpoint.getPositionWorldOnB();
      btVector3 dirf(cf(0), cf(1), cf(2));
      btVector3 endf = posa + dirf.normalize() * w;
      btVector3 cnB = posb + cpoint.m_normalWorldOnB.normalize() / 4.;


      glPushMatrix();
      glColor3f(.80, 0, 0);
      glLineWidth(2.);
      QGLViewer::drawArrow(qglviewer::Vec(posa[0], posa[1], posa[2]), qglviewer::Vec(endf[0], endf[1], endf[2]), .05, 10.);
      glPopMatrix();

      glPushMatrix();
      glColor3f(0, .80, 0);
      glLineWidth(1.);
      QGLViewer::drawArrow(qglviewer::Vec(posb[0], posb[1], posb[2]), qglviewer::Vec(cnB[0], cnB[1], cnB[2]), .03, 10.);
      glPopMatrix();

    }
  }

  for (unsigned int i = 0; i < GETNDS(Siconos_); i++)
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

    btCollisionObject* co = &**iso;
    btCollisionWorld* cow = &*ask<ForCollisionWorld>(*Siconos_->spaceFilter());
    btScalar m[16];
    btMatrix3x3 rot;
    co->getWorldTransform().getOpenGLMatrix(m);
    rot = co->getWorldTransform().getBasis();

    btVector3 aabbMin, aabbMax;
    cow->getBroadphase()->getBroadphaseAabb(aabbMin, aabbMax);

    aabbMin -= btVector3(BT_LARGE_FLOAT, BT_LARGE_FLOAT, BT_LARGE_FLOAT);
    aabbMax += btVector3(BT_LARGE_FLOAT, BT_LARGE_FLOAT, BT_LARGE_FLOAT);

    btVector3 color = btVector3(0.81f, 0.81f, 0.81f);

    _shapeDrawer.drawOpenGL(1.f, m, co->getCollisionShape(), color, 0, aabbMin, aabbMax);

    //    ask<ForCollisionWorld>(*Siconos_->spaceFilter())
    //     ->debugDrawObject((*iso)->getWorldTransform(),
    //                   (*iso)->getCollisionShape(),
    //                   btVector3(1,1,0));
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
  {
    myMouseBehavior_ = false;
    QGLViewer::mousePressEvent(e);
  }
}

void BulletViewer::mouseMoveEvent(QMouseEvent *e)
{

  if (myMouseBehavior_)
  {
    qglviewer::Vec cpos = camera()->position();

    if (selectedName() >= 0)
    {
      if (!shapes_[selectedName()]->selected())
      {
        shapes_[selectedName()]->setSelection(true);
        shapes_[selectedName()]->saveFExt();

        SP::SimpleVector fext(new SimpleVector());
        fext->resize(ask<ForFExt>(*shapes_[selectedName()]->DS())->size());

        switch (Type::value(*shapes_[selectedName()]->DS()))
        {
        case Type::NewtonEulerDS :
        {
          boost::static_pointer_cast<NewtonEulerDS>(shapes_[selectedName()]->DS())->setFExtPtr(fext);
          break;
        };
        case Type::LagrangianDS :
        {
          boost::static_pointer_cast<LagrangianDS>(shapes_[selectedName()]->DS())->setFExtPtr(fext);
          break;
        };
        default :
        {
          assert(0);
        };
        }
      }

      SP::SiconosVector fext = ask<ForFExt>(*shapes_[selectedName()]->DS());

      lastSelected_ = selectedName();

      SP::SiconosVector q =
        ask<ForPosition>(*shapes_[selectedName()]->DS());

      double massValue =
        ask<ForMassValue>(*shapes_[selectedName()]->DS());

      qglviewer::ManipulatedFrame* frame =  shapes_[selectedName()]->frame();

      //      need ref frame -> ground
      //      et update frames pos with setTranslation & setRotation
      //      Quaternion r = frame->rotation();
      //      Vec t = frame->translation();
      //      printf("translation : %g,%g,%g ; rotation : %g,%g,%g\n", t[0],t[1],t[2],r[0],r[2],r[3]);

      fext->setValue(0, -((float)cpos[0] - q->getValue(0))*massValue * 30);
      fext->setValue(1, -((float)cpos[1] - q->getValue(1))*massValue * 30);
      fext->setValue(2, -((float)cpos[2] - q->getValue(2))*massValue * 30);

    }


  }
  else
  {
    QGLViewer::mouseMoveEvent(e);
  }
}

void BulletViewer::mouseReleaseEvent(QMouseEvent* e)
{
  if (myMouseBehavior_)
  {
    myMouseBehavior_ = false;

    if (lastSelected_ >= 0)
    {
      shapes_[lastSelected_]->restoreFExt();
      shapes_[lastSelected_]->setSelection(false);
      setManipulatedFrame(camera()->frame());
      lastSelected_ = -1;
      setSelectedName(-1);
    };
  }
  else
  {
    QGLViewer::mouseReleaseEvent(e);
  }
}

void BulletViewer::keyPressEvent(QKeyEvent* e)
{

  if ((e->key() == Qt::Key_T) && (e->modifiers() == Qt::NoButton))
  {
    _transparency -= .01;
    _transparency = fmax(0., _transparency);
    printf("transparency: %f\n", _transparency);

  }

  if ((e->key() == Qt::Key_R) && (e->modifiers() == Qt::NoButton))
  {
    _transparency += .01;
    _transparency = fmin(1., _transparency);
    printf("transparency: %f\n", _transparency);
  }

  QGLViewer::keyPressEvent(e);
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


