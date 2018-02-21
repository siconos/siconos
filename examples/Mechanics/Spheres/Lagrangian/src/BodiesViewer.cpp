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

#include "BodiesViewer.hpp"
#include <Model.hpp>
#include <Simulation.hpp>

#include <GL/gl.h>
#include <GL/glu.h>
#include <op3x3.h>


using namespace std;


void BodiesViewer::init()
{
  // viewer and scene state
  restoreStateFromFile();
  setSceneRadius(20.);
  showEntireScene();

  // color
  setBackgroundColor(QColor(200, 200, 200));
  setForegroundColor(QColor(0, 0, 0));

  camera()->setPosition(qglviewer::Vec(0.0, 0.0, 1.0));
  camera()->showEntireScene();
  setManipulatedFrame(camera()->frame());

  assert(Siconos_);
  Siconos_->init();

  //  Siconos_->model()->simulation()->setStaticLevels(true);

  stepSimulation_ = false;
  stepNow_ = false;
  //  setAnimationPeriod(0.);
  setAnimationPeriod(Siconos_->simulation()->timeStep() * 1000);

  DynamicalSystemsGraph::VIterator dsi, dsend;
  boost::tie(dsi, dsend) = GETALLDS(Siconos_)->vertices();
  for (unsigned int i = 0; dsi != dsend; ++i, ++dsi)
  {
    SP::DynamicalSystem ds = GETALLDS(Siconos_)->bundle(*dsi);
    insertQGLShape(ask<ForShape>(*ds), ds);
  }

  lastSelected_ = -1;

  initUCircle();

  myMouseBehavior_ = false;

  startAnimation();
}


static const int Ne = 64;

static float circcoords[Ne][2];

void BodiesViewer::initUCircle()
{
  int i;

  GLUquadricObj * q = gluNewQuadric();
  gluQuadricNormals(q, GLU_TRUE);
  gluQuadricTexture(q, GLU_TRUE);

  for (i = 0; i < Ne; i++)
  {
    circcoords[i][0] = cos(i * 2 * M_PI / Ne);
    circcoords[i][1] = sin(i * 2 * M_PI / Ne);
  }
}

void  BodiesViewer::drawUCircle()
{
  int i;

  glBegin(GL_LINE_LOOP);
  for (i = 0; i < Ne; i++)
    glVertex2fv(&circcoords[i][0]);
  glEnd();

  //   gluPartialDisk(q,.9,1.,20,64,double(0.),double(360.));

}

void BodiesViewer::drawUDisk()
{
  int i;

  glBegin(GL_POLYGON);
  for (i = 0; i < Ne; i++)
    glVertex2fv(&circcoords[i][0]);
  glEnd();

  //   gluPartialDisk(q,.9,1.,20,64,double(0.),double(360.));

}

/* Ticks inside a circle */
void BodiesViewer::drawUCircleTicks(float a)
{
  int i;

  glRotatef(a * 180. / M_PI, 0, 0, 1);
  glPointSize(4.);
  glBegin(GL_POINTS);

  for (i = 0; i < Ne; i += 8)
  {
    glVertex2f(circcoords[i][0], circcoords[i][1]);
  }
  glEnd();

  //   gluPartialDisk(q,.9,1.,20,64,double(0.),double(360.));

}

/* Ticks inside a disk */
void BodiesViewer::drawUDiskTicks(float a)
{
  int i;

  glRotatef(a * 180. / M_PI, 0, 0, 1);
  glBegin(GL_LINE_LOOP);
  for (i = 0; i < Ne; i += 8)
  {
    glVertex2f(0, 0);
    glVertex2f(circcoords[i][0], circcoords[i][1]);
  }
  glEnd();

  //   gluPartialDisk(q,.9,1.,20,64,double(0.),double(360.));

}

void BodiesViewer::drawUTriangle(float depth)
{
  glBegin(GL_TRIANGLES);
  glVertex3f(0.0f, 1.0f, depth);
  glVertex3f(-1.0f, -1.0f, depth);
  glVertex3f(1.0f, -1.0f, depth);
  glEnd();
}

void BodiesViewer::drawVec(float x1, float y1, float x2, float y2)
{
  glPushMatrix();
  glLineWidth(2.);
  glBegin(GL_LINES);
  glVertex2f(x1, y1);
  glVertex2f(x2, y2);
  glEnd();
  glPopMatrix();
}

void BodiesViewer::drawRec(float x1, float y1, float x2, float y2, float w, float z)
{

  float dx = x2 - x1;
  float dy = y2 - y1;

  float d = sqrt(dx * dx + dy * dy);

  float ddx = w * dx / d;
  float ddy = w * dy / d;

  glPushMatrix();
  glLineWidth(2.);
  glBegin(GL_POLYGON);
  glVertex3f(x1 + ddy, y1 - ddx, z);
  glVertex3f(x1 - ddy, y1 + ddx, z);

  glVertex3f(x2 - ddy, y2 + ddx, z);
  glVertex3f(x2 + ddy, y2 - ddx, z);

  glEnd();
  glPopMatrix();
}


/* always visible */
void BodiesViewer::drawArrow(float x1, float y1, float x2, float y2, float w)
{
  float th, dx, dy;
  drawRec(x1, y1, x2, y2, w, .01);

  glPushMatrix();
  glTranslatef(x2, y2, 0.1);
  glScalef(2.*w, 2.*w, 0);
  dx = x2 - x1;
  dy = y2 - y1;
  th = (y2 > y1 ? 1 : -1) * acos(dx / sqrt(dx * dx + dy * dy));
  glRotatef(-90. + 180.*th / M_PI, 0, 0, 1);
  glTranslatef(0., .5, 0.);
  drawUTriangle(0.1);
  glPopMatrix();
}

void BodiesViewer::drawPar(float x1, float y1, float z1, float x2, float y2, float z2, float w)
{

  float dx = x2 - x1;
  float dy = y2 - y1;
  float dz = z2 - z1;

  float d = sqrt(dx * dx + dy * dy + dz * dz);

  float ddx = w * dx / d;
  float ddy = w * dy / d;
  float ddz = w * dz / d;

  glPushMatrix();
  glLineWidth(2.);
  glBegin(GL_POLYGON);
  glVertex3f(x1 + ddy, y1 - ddx, z1 - ddz);
  glVertex3f(x1 - ddy, y1 + ddx, z1 + ddz);

  glVertex3f(x2 - ddy, y2 + ddx, z2 - ddz);
  glVertex3f(x2 + ddy, y2 - ddx, z2 + ddz);

  glEnd();
  glPopMatrix();
}

/* always visible */
void BodiesViewer::drawArrow(float x1, float y1, float z1, float x2, float y2, float z2, float w)
{
  float th, dx, dy, dz;
  drawPar(x1, y1, z1, x2, y2, z2, w);

  glPushMatrix();
  glTranslatef(x2, y2, z2);
  glScalef(2.*w, 2.*w, 2.*w);
  dx = x2 - x1;
  dy = y2 - y1;
  dz = z2 - z1;
  th = (y2 > y1 ? 1 : -1) * acos(dx / sqrt(dx * dx + dy * dy + dz * dz));
  glRotatef(-90. + 180.*th / M_PI, 0, 0, 1);
  glTranslatef(0., .5, 0.);
  drawUTriangle(0.1);
  glPopMatrix();
}

void BodiesViewer::drawCircleTicks(float x, float y, float a, float r)
{
  glPushMatrix();
  glTranslatef(x, y, 0);
  glScalef(r, r, 0);
  drawUCircleTicks(a);
  glPopMatrix();
}

void BodiesViewer::drawDiskTicks(float x, float y, float a, float r, float *c)
{

  float cl[3];
  cl[0] = c[0] * .4;
  cl[1] = c[0] * .4;
  cl[2] = c[0] * .4;

  glPushMatrix();
  glColor3fv(cl);
  glLineWidth(2.);
  glTranslatef(x, y, 0);
  glScalef(r, r, 0);
  drawUDiskTicks(a);
  glPopMatrix();
}

void BodiesViewer::drawSimpleCircle(float x, float y, float r)
{
  glPushMatrix();
  glTranslatef(x, y, 0);
  glScalef(r, r, 0);
  drawUCircle();
  glPopMatrix();
}


void BodiesViewer::drawCircle(float x, float y, float a, float r, float *c)
{
  glPushMatrix();
  glColor3f(0., 0., 0.);
  drawCircleTicks(x, y, a, r);
  glPopMatrix();
  glPushMatrix();
  glColor3fv(c);
  glTranslatef(x, y, 0);
  glScalef(r, r, 0);
  drawUCircle();
  glPopMatrix();
}

void BodiesViewer::drawDisk(float x, float y, float a, float r, float *c)
{
  float cl[3];

  glColor3fv(c);
  drawSimpleCircle(x, y, r);
  drawDiskTicks(x, y, a, r, c);

  cl[0] = c[0] * 1.5;
  cl[1] = c[0];
  cl[2] = c[0];
  glColor3fv(cl);
  glPushMatrix();
  glTranslatef(x, y, 0);
  glScalef(r, r, 0);
  drawUDisk();
  glPopMatrix();
}

void BodiesViewer::rotate(const float R[12])
{
  GLfloat matrix[16];
  matrix[0] = R[0];
  matrix[1] = R[4];
  matrix[2] = R[8];
  matrix[3] = 0;
  matrix[4] = R[1];
  matrix[5] = R[5];
  matrix[6] = R[9];
  matrix[7] = 0;
  matrix[8] = R[2];
  matrix[9] = R[6];
  matrix[10] = R[10];
  matrix[11] = 0;
  matrix[12] = 0;
  matrix[13] = 0;
  matrix[14] = 0;
  matrix[15] = 1;
  glMultMatrixf(matrix);
}

void BodiesViewer::drawSphere(float x, float y, float z, float theta,
                              float phi, float psi, float r, float *c)
{

  float cl[3] = { 0, 0, 0 };

  float clg1[3] = { 1, 0, 0 };
  float clg2[3] = { 0, 1, 0 };
  float clg3[3] = { 1, 1, 1 };

  GLUquadricObj *Disk1 = gluNewQuadric();
  GLUquadricObj *Disk2 = gluNewQuadric();

  GLUquadricObj *Diskphi = gluNewQuadric();
  GLUquadricObj *Disktheta = gluNewQuadric();
  GLUquadricObj *Diskpsi = gluNewQuadric();

  GLUquadricObj *Sphere = gluNewQuadric();


  GLint slices = 10; // number of slices around z-axis
  GLint stacks = 10; // number of stacks along z-axis
  GLint loops = 10; // number of stacks along z-axis

  glPushMatrix();

  glTranslatef(x, y, z);

  glRotatef(phi * 180. / M_PI, 0, 0, 1);

  glPushMatrix();
  glColor3fv(clg1);
  glRotatef(90, 1, 0, 0);
  gluDisk(Diskphi, r, 1.1 * r, slices, loops);
  glPopMatrix();


  glRotatef(theta * 180. / M_PI, 1, 0, 0);

  glPushMatrix();
  glColor3fv(clg2);
  glRotatef(90, 0, 1, 0);
  gluDisk(Disktheta, r, 1.1 * r, slices, loops);
  glPopMatrix();


  glRotatef(psi * 180. / M_PI, 0, 0, 1);

  glPushMatrix();
  glColor3fv(clg3);
  glRotatef(90, 0, 1, 0);
  gluDisk(Diskpsi, r, 1.1 * r, slices, loops);
  glPopMatrix();


  glColor3fv(c);
  gluSphere(Sphere, r, slices, stacks);

  glColor3fv(cl);

  glPushMatrix();
  glRotatef(90, 0, 1, 0);
  gluDisk(Disk1, r, 1.05 * r, slices, loops);
  glPopMatrix();

  glPushMatrix();
  glRotatef(90, 1, 0, 0);
  gluDisk(Disk2, r, 1.05 * r, slices, loops);
  glPopMatrix();

  glPopMatrix();
}



void BodiesViewer::drawSphere(float x, float y, float z, float a,
                              float b, float c, float d, float r, float *color)
{

  GLUquadricObj *Sphere = gluNewQuadric();
  GLUquadricObj *Disk1 = gluNewQuadric();

  GLint slices = 10; // number of slices around z-axis
  GLint stacks = 10; // number of stacks along z-axis
  GLint loops = 10; // number of stacks along z-axis


  qglviewer::Quaternion q = qglviewer::Quaternion(b, c, d, a);

  glPushMatrix();

  glTranslatef(x, y, z);

  glMultMatrixd(q.matrix());

  glColor3fv(color);
  gluSphere(Sphere, r, slices, stacks);

  glPushMatrix();
  glColor3f(.7, .7, .7);
  glRotatef(90, 0, 1, 0);
  gluDisk(Disk1, r, 1.1 * r, slices, loops);
  glPopMatrix();

  glPushMatrix();
  glColor3f(.7, .7, .7);
  glRotatef(90, 1, 0, 0);
  gluDisk(Disk1, r, 1.1 * r, slices, loops);
  glPopMatrix();

  glPopMatrix();
}

void BodiesViewer::drawPolyg(unsigned int nvertices, double *coor, float *c)
{

  float cl[3];
  cl[0] = c[0] * 1.5;
  cl[1] = c[0];
  cl[2] = c[0];
  double *pcoor = coor;

  glPushMatrix();

  glLineWidth(2.);

  glColor3fv(cl);

  glBegin(GL_LINE_LOOP);
  for (unsigned int i = 0; i < nvertices; ++i, ++pcoor, ++pcoor)
  {
    glVertex2dv(pcoor);
  }
  glEnd();

  glPopMatrix();

  glPushMatrix();

  glColor3fv(c);

  pcoor = coor;

  glBegin(GL_POLYGON);
  for (unsigned int i = 0; i < nvertices; ++i, ++pcoor, ++pcoor)
  {
    glVertex2dv(pcoor);
  }
  glEnd();

  glPopMatrix();
}

void BodiesViewer::print(float x, float y, const char *s, int size)
{
  qreal coor[3];
  int i, j, ps1, ps2;
  coor[0] = x;
  coor[1] = y;
  coor[2] = 0;

  camera()->getProjectedCoordinatesOf(coor, coor, NULL);
  i = (int) coor[0];
  j = (int) coor[1];

  QFont F = QApplication::font();

  //text size calibration
  ps1 = min(max((int)(camera()->position())[2], 600), size);
  ps2 = min(max((int) 400 / ps1, 600), size);

  F.setPointSize(ps2);
  drawText(i, j, QString(s), F);
}

void BodiesViewer::drawWithNames()
{
  for (unsigned int i = 0; i < GETNDS(Siconos_); i++)
  {

    glPushName(i);
    drawQGLShape(*shapes_[i]);
    glPopName();
  }
}

void BodiesViewer::postSelection(const QPoint& point)
{

  if (selectedName() == -1)
  {
    setManipulatedFrame(camera()->frame());
  }
  else
  {
    lastSelected_ = selectedName();
    setManipulatedFrame(shapes_[selectedName()]->frame());
  }

}

void BodiesViewer::drawQGLShape(const QGLShape& fig)
{

  DynamicalSystem *ds = fig.DS().get();
  float c[3];


  switch (fig.shape())
  {
  case CIRCLE :
    c[0] = .9 ;
    c[1] = .1;
    c[2] = .1;
    drawCircle(GETX(ds),
               GETY(ds),
               GETZ(ds),
               GETRADIUS(ds), c);
    break;

  case DISK :
    c[0] = 1 ;
    c[1] = .0;
    c[2] = 0;
    glLineWidth(2.);
    drawDisk(GETX(ds),
             GETY(ds),
             GETZ(ds),
             GETRADIUS(ds), c);
    break;

  case SPHERE:
    c[0] = .9 ;
    c[1] = .9;
    c[2] = .9;
    glLineWidth(2.);
    switch (fig.type())
    {
    case Type::LagrangianDS :
      drawSphere(GETX(ds),
                 GETY(ds),
                 GETZ(ds),
                 GETA1(ds),
                 GETA2(ds),
                 GETA3(ds),
                 GETRADIUS(ds), c);
      break;
    case Type::NewtonEulerDS :
      drawSphere(GETX(ds),
                 GETY(ds),
                 GETZ(ds),
                 GETA1(ds),
                 GETA2(ds),
                 GETA3(ds),
                 GETA4(ds),
                 GETRADIUS(ds), c);
      break;

    default :
    {};
    }
    break;

  default:
  {};

  };




};

void BodiesViewer::drawSelectedQGLShape(const QGLShape& fig)
{


  DynamicalSystem * ds = fig.DS().get();

  float c[3], dFe, Cal;



  switch (fig.shape())
  {
  case CIRCLE :
    c[0] = .7 ;
    c[1] = .5;
    c[2] = .5;
    glLineWidth(3.);
    drawCircle(GETX(ds), GETY(ds), GETZ(ds), GETRADIUS(ds), c);
    glColor3f(1., 0., 0.);
    dFe = hypot(GETXFE(ds), GETYFE(ds));
    Cal = log(dFe);

    drawArrow(GETX(ds), GETY(ds),
              GETX(ds) + Cal * GETXFE(ds) / dFe,
              GETY(ds) + Cal * GETYFE(ds) / dFe, Cal / 10.);
    break;

  case DISK :
    c[0] = .6 ;
    c[1] = .1;
    c[2] = .5;
    glLineWidth(3.);
    drawDisk(GETX(ds), GETY(ds), GETZ(ds), GETRADIUS(ds), c);
    glColor3f(1., 0., 0.);
    dFe = hypot(GETXFE(ds), GETYFE(ds));
    Cal = log(dFe);

    drawArrow(GETX(ds), GETY(ds),
              GETX(ds) + Cal * GETXFE(ds) / dFe,
              GETY(ds) + Cal * GETYFE(ds) / dFe, Cal / 10.);
    break;

  case SPHERE:
    c[0] = .6 ;
    c[1] = .1;
    c[2] = .5;
    glLineWidth(3.);
    switch (fig.type())
    {
    case Type::LagrangianDS :

      drawSphere(GETX(ds), GETY(ds), GETZ(ds),
                 GETA1(ds), GETA2(ds), GETA3(ds),
                 GETRADIUS(ds), c);
      break;
    case Type::NewtonEulerDS :
      drawSphere(GETX(ds), GETY(ds), GETZ(ds),
                 GETA1(ds), GETA2(ds), GETA3(ds), GETA4(ds),
                 GETRADIUS(ds), c);
      break;

    default :
    {};
    }
    glColor3f(1., 0., 0.);
    dFe = hypot3(ask<ForFExt>(*ds)->getArray());
    Cal = log(dFe);

    drawArrow(GETX(ds), GETY(ds), GETZ(ds),
              GETX(ds) + Cal * GETXFE(ds) / dFe,
              GETY(ds) + Cal * GETYFE(ds) / dFe,
              GETZ(ds) + Cal * GETZFE(ds), Cal / 10.);
    break;

  default:
  {};

  }
};

void BodiesViewer::insertQGLShape(SHAPE f, SP::DynamicalSystem ds)
{

  SP::QGLShape fig(new QGLShape(f, ds, camera()->frame()));

  shapes_.push_back(fig);
}

void BodiesViewer::mouseMoveEvent(QMouseEvent *e)
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

        SP::SiconosVector fext(new SiconosVector());
        fext->resize(ask<ForFExt>(*shapes_[selectedName()]->DS())->size());

        switch (Type::value(*shapes_[selectedName()]->DS()))
        {
        case Type::NewtonEulerDS :
        {
          std11::static_pointer_cast<NewtonEulerDS>(shapes_[selectedName()]->DS())->setFExtPtr(fext);
          break;
        };
        case Type::LagrangianDS :
        {
          std11::static_pointer_cast<LagrangianDS>(shapes_[selectedName()]->DS())->setFExtPtr(fext);
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

      //      need ref frame -> ground
      //      et update frames pos with setTranslation & setRotation
      //      Quaternion r = frame->rotation();
      //      Vec t = frame->translation();
      //      printf("translation : %g,%g,%g ; rotation : %g,%g,%g\n", t[0],t[1],t[2],r[0],r[2],r[3]);

      fext->setValue(0, -((float)cpos[0] - q->getValue(0))*massValue * 30);
      fext->setValue(1, -((float)cpos[1] - q->getValue(1))*massValue * 30);
      if (ask<ForNdof>(*shapes_[selectedName()]->DS()) > 3)
      {
        fext->setValue(2, -((float)cpos[2] - q->getValue(2))*massValue * 30);
      }

    }


  }
  else
  {
    QGLViewer::mouseMoveEvent(e);
  }
}

void BodiesViewer::mousePressEvent(QMouseEvent* e)
{
  if ((e->button() == Qt::LeftButton) && (e->modifiers() == Qt::NoButton))
    myMouseBehavior_ = true;
  else
  {
    myMouseBehavior_ = false;
    QGLViewer::mousePressEvent(e);
  }
}


void BodiesViewer::mouseReleaseEvent(QMouseEvent* e)
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


void BodiesViewer::keyPressEvent(QKeyEvent* e)
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

  if ((e->key() == (Qt::Key) 'S'))
  {
    stepSimulation_ = not stepSimulation_;
    stepNow_ = false;
    printf("step simulation mode\n");
  }

  if ((e->key() == (Qt::Key) '.'))
  {
    stepNow_ = true;
    printf("step now\n");
  }



  QGLViewer::keyPressEvent(e);
}

void BodiesViewer::animate()
{

  if (!stepSimulation_)
  {
    Siconos_->compute();
  }
  else if (stepNow_)
  {
    stepNow_ = false;
    Siconos_->compute();
  }
  //  saveSnapshot();
}

QString BodiesViewer::helpString() const
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

#ifdef QT_INTERFACE
BodiesViewer::BodiesViewer(QWidget *parent)
  : QGLViewer(parent)
{
};
#endif

