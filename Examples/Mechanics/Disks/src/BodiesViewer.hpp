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

#ifndef BodiesViewer_hpp
#define BodiesViewer_hpp

/* Siconos */
#include <SiconosKernel.hpp>

#include <SpaceFilter.hpp>
#include "SiconosBodies.hpp"



/* QGLViewer */
#include <QGLViewer/qglviewer.h>
#include <qapplication.h>

#include <qevent.h>
#include <qmessagebox.h>

#ifdef QT_INTERFACE
#include "ui_viewerInterface.Qt4.h"
class ViewerInterface : public QDialog, public Ui::Dialog
{
public:
  ViewerInterface()
  {
    setupUi(this);
  }
};
#endif


/* needed shared pointers */
DEFINE_SPTR(Drawing);
DEFINE_SPTR(BodyDraw);
DEFINE_SPTR(SelectedBodyDraw);

/* convenient macros */
#define GETX(C) C->q()->getValue(0)
#define GETY(C) C->q()->getValue(1)
#define GETZ(C) C->q()->getValue(2)
#define GETA1(C) C->q()->getValue(3)
#define GETA2(C) C->q()->getValue(4)
#define GETA3(C) C->q()->getValue(5)
#define GETXFE(C) C->fExt()->getValue(0)
#define GETYFE(C) C->fExt()->getValue(1)
#define GETZFE(C) C->fExt()->getValue(2)

#define GETVX(C) C->velocity()->getValue(0)
#define GETVY(C) C->velocity()->getValue(1)
#define GETVZ(C) C->velocity()->getValue(2)

#define GETVA1(C) C->velocity()->getValue(3)
#define GETVA2(C) C->velocity()->getValue(4)
#define GETVA3(C) C->velocity()->getValue(5)

#define GETALLDS(M) M->model()->nonSmoothDynamicalSystem()->topology()->dSG(0)
#define GETNDS(M) GETALLDS(M)->size()

/* a drawing of a Siconos Lagrangian DS */
class Drawing
{

public:

  /* construction from a LagrangianDS */
  Drawing(SP::LagrangianDS D)
  {
    DS_ = D;
    frame_.reset(new qglviewer::ManipulatedFrame());
    savedFExt_.reset(new SimpleVector(getDS()->getNdof()));
    selected_ = false;
  };

  ~Drawing() {};

  /* pointer to DS */
  SP::LagrangianDS getDS()
  {
    return DS_;
  };

  /* selection with mouse */
  bool selected()
  {
    return selected_ ;
  };
  void nextSelection()
  {
    selected_ = !selected_ ;
  };

  /* identifiant */
  int getID()
  {
    return id_ ;
  };
  void setID(int i)
  {
    id_ = i;
  };


  /* External force set from mouse and restore */
  void saveFExt()
  {
    *savedFExt_ = *(getDS()->fExt()) ;
  };
  void restoreFExt()
  {
    *(getDS()->fExt()) = *savedFExt_;
  };

  /* DS frame */
  qglviewer::ManipulatedFrame * getFrame()
  {
    return frame_.get();
  };

protected:
  int id_;
  bool selected_;
  SP::SiconosVector savedFExt_;
  boost::shared_ptr<qglviewer::ManipulatedFrame> frame_;
  SP::LagrangianDS DS_;
};



/* QGLViewer main object */
class BodiesViewer : public QGLViewer
{

public:

#ifdef QT_INTERFACE
  Viewer(QWidget *parent);
#endif
  std::vector<SP::Drawing> drawings_;

public:
  virtual void draw() = 0;
  void drawWithNames();
  void initUCircle();
  void drawUCircle();
  void drawUDisk();
  void drawUCircleTicks(float a);
  void drawUDiskTicks(float a) ;
  void drawUTriangle(float depth = 0.) ;
  void drawVec(float x1, float y1, float x2, float y2) ;
  void drawRec(float x1, float y1, float x2, float y2, float w, float z = 0.);
  void drawArrow(float x1, float y1, float x2, float y2, float w);
  void drawPar(float x1, float y1, float z1, float x2, float y2, float z2, float w);
  void drawArrow(float x1, float y1, float z1, float x2, float y2, float z2, float w);
  void drawCircleTicks(float x, float y, float a, float r);
  void drawDiskTicks(float x, float y, float a, float r, float *c);
  void drawSimpleCircle(float x, float y, float r);
  void drawCircle(float x, float y, float a, float r, float *c);
  void drawDisk(float x, float y, float a, float r, float *c);
  void rotate(const float R[12]);
  void drawSphere(float x, float y, float z, float theta, float phi, float psi, float r, float *c);


protected :

  void postSelection(const QPoint& point);
  virtual void init() = 0;
  virtual void animate() = 0;
  virtual void mousePressEvent(QMouseEvent *) = 0;
  virtual void mouseMoveEvent(QMouseEvent *) = 0;
  virtual void mouseReleaseEvent(QMouseEvent *) = 0;

  virtual QString helpString() const = 0;
  boost::shared_ptr<qglviewer::WorldConstraint> constraint_;

  SP::SiconosBodies Siconos_;

  void print(float x, float y, const char *s, int size);

  int NDS_;


  SP::BodyDraw bodydraw_;
  SP::SelectedBodyDraw sbodydraw_;

  SP::OneStepNSProblem oneStepNSP;

  qglviewer::Vec selectedPoint_;

  int lastSelected_;

  bool myMouseBehavior_;

  long timeSiconos_;

  long timeGlob_;

};

class SetDrawing : public SiconosVisitor
{
private:
  BodiesViewer*  viewer_;
  int i_;

public:
  SetDrawing(BodiesViewer& viewer, int i) : i_(i)
  {
    viewer_ = &viewer ;
  };

  void visit(SP::Disk d)
  {
    SP::Drawing drawing(new Drawing(d));
    drawing->setID(i_);
    drawing->saveFExt();
    viewer_->drawings_.push_back(drawing);
  }

  void visit(SP::Circle d)
  {
    SP::Drawing drawing(new Drawing(d));
    drawing->setID(i_);
    drawing->saveFExt();
    viewer_->drawings_.push_back(drawing);
  }

  void visit(SP::SphereLDS d)
  {
    SP::Drawing drawing(new Drawing(d));
    drawing->setID(i_);
    drawing->saveFExt();
    viewer_->drawings_.push_back(drawing);
  }

};

/*
 * BodyQGL
 *
 */

class BodyDraw : public SiconosVisitor
{
private:
  BodiesViewer* viewer_;

public:
  BodyDraw(BodiesViewer& viewer)
  {
    viewer_ = &viewer ;
  };
  void visit(SP::Circle circle)
  {
    float c[3];
    c[0] = .9 ;
    c[1] = .1;
    c[2] = .1;
    viewer_->drawCircle(GETX(circle),
                        GETY(circle),
                        GETZ(circle),
                        circle->getRadius(), c);
  }
  void visit(SP::Disk disk)
  {
    float c[3];
    c[0] = 1 ;
    c[1] = .0;
    c[2] = 0;
    glLineWidth(2.);
    viewer_->drawDisk(GETX(disk),
                      GETY(disk),
                      GETZ(disk),
                      disk->getRadius(), c);
  }

  void visit(SP::SphereLDS sphere)
  {
    float c[3];
    c[0] = 1 ;
    c[1] = .0;
    c[2] = 0;
    glLineWidth(2.);
    viewer_->drawSphere(GETX(sphere),
                        GETY(sphere),
                        GETZ(sphere),
                        GETA1(sphere),
                        GETA2(sphere),
                        GETA3(sphere),
                        sphere->getRadius(), c);
  }

};


class SelectedBodyDraw : public SiconosVisitor
{
private:

  BodiesViewer* viewer_;

  inline double hypot3d(double x, double y, double z)
  {
    return sqrt(x * x + y * y + z * z);
  };

  inline float hypot3d(float x, float y, float z)
  {
    return sqrt(x * x + y * y + z * z);
  };

public:
  SelectedBodyDraw(BodiesViewer& viewer)
  {
    viewer_ = &viewer ;
  };
  void visit(SP::Circle circle)
  {
    float c[3], dFe, Cal;
    c[0] = .7 ;
    c[1] = .5;
    c[2] = .5;
    glLineWidth(3.);
    viewer_->drawCircle(GETX(circle), GETY(circle), GETZ(circle), circle->getRadius(), c);
    glColor3f(1., 0., 0.);
    dFe = hypot(GETXFE(circle), GETYFE(circle));
    Cal = log(dFe);

    viewer_->drawArrow(GETX(circle), GETY(circle),
                       GETX(circle) + Cal * GETXFE(circle) / dFe,
                       GETY(circle) + Cal * GETYFE(circle) / dFe, Cal / 10.);
  }
  void visit(SP::Disk disk)
  {
    float c[3], dFe, Cal;
    c[0] = .6 ;
    c[1] = .1;
    c[2] = .5;
    glLineWidth(3.);
    viewer_->drawDisk(GETX(disk), GETY(disk), GETZ(disk), disk->getRadius(), c);
    glColor3f(1., 0., 0.);
    dFe = hypot(GETXFE(disk), GETYFE(disk));
    Cal = log(dFe);

    viewer_->drawArrow(GETX(disk), GETY(disk),
                       GETX(disk) + Cal * GETXFE(disk) / dFe,
                       GETY(disk) + Cal * GETYFE(disk) / dFe, Cal / 10.);
  }

  void visit(SP::SphereLDS sphere)
  {
    float c[3], dFe, Cal;
    c[0] = .6 ;
    c[1] = .1;
    c[2] = .5;
    glLineWidth(3.);
    viewer_->drawSphere(GETX(sphere), GETY(sphere), GETZ(sphere),
                        GETA1(sphere), GETA2(sphere), GETA3(sphere),
                        sphere->getRadius(), c);
    glColor3f(1., 0., 0.);
    dFe = hypot3d(GETXFE(sphere), GETYFE(sphere), GETZFE(sphere));
    Cal = log(dFe);

    viewer_->drawArrow(GETX(sphere), GETY(sphere), GETZ(sphere),
                       GETX(sphere) + Cal * GETXFE(sphere) / dFe,
                       GETY(sphere) + Cal * GETYFE(sphere) / dFe,
                       GETZ(sphere) + Cal * GETZFE(sphere), Cal / 10.);
  }

};

class LambdaSecond;

class LambdaFirst : public SiconosVisitor,
  public boost::enable_shared_from_this<LambdaFirst>
{

public:
  SP::LagrangianR relation;
  float x;
  float y;
  float z;
  float r;


  void visit(SP::DiskPlanR rel)
  {
    relation = rel;
  }

  void visit(SP::CircleCircleR rel)
  {
    relation = rel;
  }

  void visit(SP::DiskDiskR rel)
  {
    relation = rel;
  }

  void visit(SP::DiskMovingPlanR rel)
  {
    relation = rel;
  }

  void visit(SP::Circle circle)
  {
    x = GETX(circle);
    y = GETY(circle);
    r = circle->getRadius();
  }

  void visit(SP::Disk disk)
  {
    x = GETX(disk);
    y = GETY(disk);
    r = disk->getRadius();
  }

  void visit(SP::SphereLDS sphere)
  {
    x = GETX(sphere);
    y = GETY(sphere);
    z = GETZ(sphere);
    r = sphere->getRadius();
  }

  void accept(boost::shared_ptr<LambdaSecond>);

};



class LambdaSecond : public SiconosVisitor
{
private:
  float _w;
  float x1;
  float y1;
  float z1;
  float r1;
  SP::LagrangianR relation1;
  BodiesViewer* viewer_;

public:
  LambdaSecond(BodiesViewer& viewer, float w) : _w(w)
  {
    viewer_ = &viewer;
  };

  void visit(boost::shared_ptr<LambdaFirst> lambdastart)
  {
    x1 = lambdastart->x;
    y1 = lambdastart->y;
    z1 = lambdastart->z;
    r1 = lambdastart->r;
    relation1 = lambdastart->relation;
  };


  void visit(SP::Circle circle)
  {
    float x2 = GETX(circle);
    float y2 = GETY(circle);
    float r = circle->getRadius() + r1;
    float d = hypotf(x1 - x2, y1 - y2);

    glColor3f(.0f, .0f, .0f);
    viewer_->drawRec(x1, y1, x1 - (x2 - x1)*r / d, y1 - (y2 - y1)*r / d, _w);
  };

  void visit(SP::Disk disk)
  {

    float x2 = GETX(disk);
    float y2 = GETY(disk);
    float r = disk->getRadius() + r1;
    float d = hypotf(x1 - x2, y1 - y2);

    glColor3f(.0f, .0f, .0f);

    if (d > 0)
      viewer_->drawRec(x1, y1, x1 + (x2 - x1)*r / d, y1 + (y2 - y1)*r / d, _w);
    else
    {
      double jx = relation1->jachq()->getValue(0, 0);
      double jy = relation1->jachq()->getValue(0, 1);
      double dj = 2 * hypot(jx, jy);
      viewer_->drawRec(x1, y1, x1 - r * jx / dj, y1 - r * jy / dj, _w);
    }

  };

  void visit(SP::SphereLDS sphere)
  {

    float x2 = GETX(sphere);
    float y2 = GETY(sphere);
    float z2 = GETZ(sphere);
    float dx = x1 - x2;
    float dy = y1 - y2;
    float dz = z1 - z2;
    float r = sphere->getRadius() + r1;
    float d = sqrt(dx * dx + dy * dy + dz * dz);

    glColor3f(.0f, .0f, .0f);
    viewer_->drawPar(x1, y1, z1, x1 + (x2 - x1)*r / d, y1 + (y2 - y1)*r / d, z1 + (z2 - z1)*r / d, _w);

  };

};


#endif
