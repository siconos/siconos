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

#include "BodiesViewer.hpp"

#include "BulletSpaceFilter.hpp"
#include "BulletBodies.hpp"

#include "GL_ShapeDrawer.h"
#include <bullet/BulletDynamics/Dynamics/btRigidBody.h>

#include <GL/gl.h>
#include <GL/glu.h>
#include <op3x3.h>


using namespace std;

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
  gluDisk(Diskphi, r, 1.5 * r, slices, loops);
  glPopMatrix();


  glRotatef(theta * 180. / M_PI, 1, 0, 0);

  glPushMatrix();
  glColor3fv(clg2);
  glRotatef(90, 0, 1, 0);
  gluDisk(Disktheta, r, 1.5 * r, slices, loops);
  glPopMatrix();


  glRotatef(psi * 180. / M_PI, 0, 0, 1);

  glPushMatrix();
  glColor3fv(clg3);
  glRotatef(90, 0, 1, 0);
  gluDisk(Diskpsi, r, 1.5 * r, slices, loops);
  glPopMatrix();


  glColor3fv(c);
  gluSphere(Sphere, r, slices, stacks);

  glColor3fv(cl);

  glPushMatrix();
  glRotatef(90, 0, 1, 0);
  gluDisk(Disk1, r, 1.1 * r, slices, loops);
  glPopMatrix();

  glPushMatrix();
  glRotatef(90, 1, 0, 0);
  gluDisk(Disk2, r, 1.1 * r, slices, loops);
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
  float coor[3];
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
    setManipulatedFrame(shapes_[selectedName()]->getFrame());
    shapes_[selectedName()]->nextSelection();
  }

}

void BodiesViewer::drawQGLShape(const QGLShape& fig)
{

  DynamicalSystem *ds = fig.DS().get();
  float c[3];


  switch (fig.kind())
  {
  case POLYG :
    c[0] = .3 ;
    c[1] = .7;
    c[2] = .3;
    drawPolyg(ask<ForVerticesNumber>(*ds),
              ask<ForVertices>(*ds)->getArray(), c);
    break;

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

  case SPHERELDS :
    c[0] = 1 ;
    c[1] = .0;
    c[2] = 0;
    glLineWidth(2.);
    drawSphere(GETX(ds),
               GETY(ds),
               GETZ(ds),
               GETA1(ds),
               GETA2(ds),
               GETA3(ds),
               GETRADIUS(ds), c);
    break;
  case SPHERENEDS :
    c[0] = 1 ;
    c[1] = .0;
    c[2] = 0;
    glLineWidth(2.);
    drawSphere(GETX(ds),
               GETY(ds),
               GETZ(ds),
               GETA1(ds),
               GETA2(ds),
               GETA3(ds),
               GETA4(ds),
               GETRADIUS(ds), c);
    break;
  };

};

void BodiesViewer::drawSelectedQGLShape(const QGLShape& fig)
{


  DynamicalSystem * ds = fig.DS().get();

  float c[3], dFe, Cal;



  switch (fig.kind())
  {
  case POLYG :
    c[0] = .4 ;
    c[1] = .4;
    c[2] = .7;
    drawPolyg(ask<ForVerticesNumber>(*ds),
              ask<ForVertices>(*ds)->getArray(), c);

    glColor3f(1., 0., 0.);
    dFe = hypot(GETXFE(ds), GETYFE(ds));
    Cal = log(dFe);

    drawArrow(GETX(ds), GETY(ds),
              GETX(ds) + Cal * GETXFE(ds) / dFe,
              GETY(ds) + Cal * GETYFE(ds) / dFe, Cal / 10.);
    break;

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

  case SPHERELDS :
    c[0] = .6 ;
    c[1] = .1;
    c[2] = .5;
    glLineWidth(3.);
    drawSphere(GETX(ds), GETY(ds), GETZ(ds),
               GETA1(ds), GETA2(ds), GETA3(ds),
               GETRADIUS(ds), c);
    glColor3f(1., 0., 0.);
    dFe = hypot3(ask<ForFExt>(*ds)->getArray());
    Cal = log(dFe);

    drawArrow(GETX(ds), GETY(ds), GETZ(ds),
              GETX(ds) + Cal * GETXFE(ds) / dFe,
              GETY(ds) + Cal * GETYFE(ds) / dFe,
              GETZ(ds) + Cal * GETZFE(ds), Cal / 10.);
    break;


  case SPHERENEDS :
    c[0] = .6 ;
    c[1] = .1;
    c[2] = .5;
    glLineWidth(3.);
    drawSphere(GETX(ds), GETY(ds), GETZ(ds),
               GETA1(ds), GETA2(ds), GETA3(ds), GETA4(ds),
               GETRADIUS(ds), c);
    glColor3f(1., 0., 0.);
    dFe = hypot3(ask<ForFExt>(*ds)->getArray());
    Cal = log(dFe);
    break;
  }
};

void BodiesViewer::insertQGLShape(SHAPE f, SP::DynamicalSystem ds)
{

  SP::QGLShape fig(new QGLShape(f, ds));

  shapes_.push_back(fig);
}

#ifdef QT_INTERFACE
BodiesViewer::BodiesViewer(QWidget *parent)
  : QGLViewer(parent)
{
};
#endif

