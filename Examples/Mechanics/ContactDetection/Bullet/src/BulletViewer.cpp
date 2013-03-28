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

#include <stdio.h>
using namespace qglviewer;

void BulletViewer::init()
{
  // siconos setup
  Siconos_.reset(new BulletBodies());

  BodiesViewer::init();

  // shapeDrawer
  _shapeDrawer.enableTexture(true);
  _transparency = 1.;

  // help screen
  help();


  ask<ForCollisionWorld>(*Siconos_->spaceFilter())->setDebugDrawer(&_debugDrawer);



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

  //  ask<ForCollisionWorld>(*Siconos_->spaceFilter())->getDebugDrawer()->setDebugMode(btIDebugDraw::DBG_DrawAabb | btIDebugDraw::DBG_DrawContactPoints);
  //  ask<ForCollisionWorld>(*Siconos_->spaceFilter())->debugDrawWorld();

  //ask<ForCollisionWorld>(*Siconos_->spaceFilter())->debugDrawObject(co->getWorldTransform(), co->getCollisionShape(), btVector3(1,1,0));


  btScalar m[16];
  btMatrix3x3 rot;
  co->getWorldTransform().getOpenGLMatrix(m);
  rot = co->getWorldTransform().getBasis();

  btVector3 aabbMin, aabbMax;
  cow->getBroadphase()->getBroadphaseAabb(aabbMin, aabbMax);

  aabbMin -= btVector3(BT_LARGE_FLOAT, BT_LARGE_FLOAT, BT_LARGE_FLOAT);
  aabbMax += btVector3(BT_LARGE_FLOAT, BT_LARGE_FLOAT, BT_LARGE_FLOAT);

  _shapeDrawer.drawOpenGL((btScalar*) &_transparency, (const btCollisionShape*) co->getCollisionShape(), (const btVector3&) _colors[ds->number() % _colors.size()], 2, (const btVector3&) aabbMin, (const btVector3&) aabbMax);

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

  //  ask<ForCollisionWorld>(*Siconos_->spaceFilter())->debugDrawWorld();


}


void BulletViewer::draw()
{

  char qs[6];

  float lbdmax = 0., w;

  DSIterator itDS;
  SP::DynamicalSystemsSet involvedDS;
  SP::Interaction interaction;
  SP::Relation relation;


  /*  static int file_counter=0;
    char filename[32];
    sprintf(filename,"data%d.txt",file_counter++);
    FILE* posf = fopen(filename,"w");


    fprintf(posf, "%g\n", Siconos_->model()->currentTime());

    for (unsigned int i=0;i<GETNDS(Siconos_);i++)
    {

      SP::DynamicalSystem d1 = shapes_[i]->DS();
      SP::SiconosVector q1 = ask<ForPosition>(*d1);

      qglviewer::Quaternion q = qglviewer::Quaternion(GETA2(d1),GETA3(d1),GETA4(d1),GETA1(d1));
      qglviewer::Vec axis;
      float angle;
      q.getAxisAngle(axis,angle);

      fprintf(posf,
              "%g %g %g %f %f %f %f file:///tmp/pyramid.vtu 0. 0.\n",
              (*q1)(0), (*q1)(1), (*q1)(2),
              angle, axis[0], axis[1], axis[2]);
              }*/


  // calibration
  if (Siconos_->model()->nonSmoothDynamicalSystem()->topology()->numberOfIndexSet() >= 2)
  {
    SP::InteractionsGraph I1 = Siconos_->model()->simulation()->indexSet(1);

    InteractionsGraph::VIterator ui, uiend;
    for (boost::tie(ui, uiend) = I1->vertices(); ui != uiend; ++ui)
    {
      lbdmax = fmax(I1->bundle(*ui)->lambdaOld(1)->getValue(0), lbdmax);
    }

    for (boost::tie(ui, uiend) = I1->vertices(); ui != uiend; ++ui)
    {
      interaction = I1->bundle(*ui);
      relation = interaction->relation();

      involvedDS = interaction->dynamicalSystems();

      itDS = involvedDS->begin();

      SP::DynamicalSystem d1 = *itDS;

      SP::SiconosVector q1 = ask<ForPosition>(*d1);

      {
        SiconosVector& cf = *ask<ForContactForce>(*relation);
        double cfn = cf.norm2();

        w = fmax(.3, cfn / fmax(lbdmax, cfn));


        btManifoldPoint& cpoint = *ask<ForContactPoint>(*relation);

        //      printf("cpoint.getDistance():%g\n",cpoint.getDistance());

        btVector3 posa = cpoint.getPositionWorldOnA();
        btVector3 posb = cpoint.getPositionWorldOnB();
        btVector3 dirf(cf(0), cf(1), cf(2));
        btVector3 endf = posa + dirf.normalize() * w;
        btVector3 cnB = posb + cpoint.m_normalWorldOnB.normalize() / 4.;

        //      printf("dist(posb,posa):%g\n",sqrt((posb-posa).dot(posb-posa)));

        glPushMatrix();
        glColor3f(.80, 0, 0);
        glLineWidth(2.);
        QGLViewer::drawArrow(qglviewer::Vec(posa[0], posa[1], posa[2]), qglviewer::Vec(endf[0], endf[1], endf[2]), .05, 10.);
        glPopMatrix();

        /*        fprintf(posf,
                        "%g %g %g %g %g %g %g Arrow %g %g\n",
                        posa[0], posa[1], posa[2],
                        0., cf(0), cf(1), cf(2), .05, w);*/

        glPushMatrix();
        glColor3f(0, .80, 0);
        glLineWidth(1.);
        QGLViewer::drawArrow(qglviewer::Vec(posb[0], posb[1], posb[2]), qglviewer::Vec(cnB[0], cnB[1], cnB[2]), .03, 10.);
        glPopMatrix();

        qglviewer::Vec cpo = qglviewer::Vec(posb[0], posb[1], posb[2]) - qglviewer::Vec(cnB[0], cnB[1], cnB[2]);

        /*        fprintf(posf,
                        "%g %g %g %g %g %g %g Arrow %g %g\n",
                        posa[0], posa[1], posa[2],
                        0., cpo[0], cpo[1], cpo[2], .05, .3);*/

      }
    }
    //    fclose(posf);
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

    btScalar transp = 1.;

    _shapeDrawer.drawOpenGL(&transp, (const btCollisionShape*) co->getCollisionShape(), (const btVector3&) color, 0, (const btVector3&) aabbMin, (const btVector3&) aabbMax);

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
      double xc = 0., yc = 0.;

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






