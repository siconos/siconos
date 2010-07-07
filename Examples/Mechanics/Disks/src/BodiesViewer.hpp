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
#include <qglviewer.h>
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
#define GETA4(C) C->q()->getValue(6)
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

#define ANSWER(T,CODE) \
  void visit(const T& ds)                       \
  {                                             \
    answer = ds . CODE;                         \
  }

#define ANSWER_V(T,CODE) \
  void visit(const T& ds)                       \
  {                                             \
    answer = CODE;                              \
  }

#define ANSWER_F(T,CODE) \
  void visit(const T& ds)                       \
  {                                             \
    answer = CODE(ds);                          \
  }


struct NeedNdof : public Question<unsigned int>
{
  ANSWER_V(Disk, 3);
  ANSWER_V(Circle, 3);
  ANSWER_V(SphereLDS, 6);
  ANSWER_V(SphereNEDS, 6);
};


struct NeedFExt : public Question<SP::SiconosVector>
{

  ANSWER(Disk, fExt());
  ANSWER(Circle, fExt());
  ANSWER(SphereLDS, fExt());
  ANSWER(SphereNEDS, fExt());
};

struct Needq : public Question<SP::SiconosVector>
{

  ANSWER(Disk, q());
  ANSWER(Circle, q());
  ANSWER(SphereLDS, q());
  ANSWER(SphereNEDS, q());

};


struct NeedRadius : public Question<double>
{

  ANSWER(Disk, getRadius());
  ANSWER(Circle, getRadius());
  ANSWER(SphereLDS, getRadius());
  ANSWER(SphereNEDS, getRadius());
};

struct NeedMassValue : public Question<double>
{
  ANSWER(Disk, mass()->getValue(0, 0));
  ANSWER(Circle, mass()->getValue(0, 0));
  ANSWER(SphereLDS, mass()->getValue(0, 0));
  ANSWER(SphereNEDS, massValue());
};

struct NeedJachq : public Question<SP::SiconosMatrix>
{
  ANSWER(LagrangianR, jachq());
  ANSWER(NewtonEulerR, jachq());
  ANSWER(DiskDiskR, jachq());
  ANSWER(DiskPlanR, jachq());
  ANSWER(DiskMovingPlanR, jachq());
  ANSWER(SphereLDSPlanR, jachq());
  ANSWER(SphereNEDSPlanR, jachq());
  ANSWER(SphereLDSSphereLDSR, jachq());
  ANSWER(SphereNEDSSphereNEDSR, jachq());
};


/* a drawing of a Siconos Lagrangian DS */
class Drawing
{

public:

  /* construction from a LagrangianDS */
  Drawing(SP::DynamicalSystem D)
  {
    DS_ = D;
    frame_.reset(new qglviewer::ManipulatedFrame());
    savedFExt_.reset(new SimpleVector(ask<NeedNdof>(*getDS())));
    selected_ = false;
  };

  ~Drawing() {};

  /* pointer to DS */
  SP::DynamicalSystem getDS()
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
  int getD()
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
    *savedFExt_ = *(ask<NeedFExt>(*getDS()));
  };
  void restoreFExt()
  {
    *(ask<NeedFExt>(*getDS())) = *savedFExt_;
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
  SP::DynamicalSystem DS_;
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
  void drawSphere(float x, float y, float z, float a, float b, float c, float d, float r, float *color);


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

  void visit(SP::SphereNEDS d)
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

  void visit(SP::SphereNEDS sphere)
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
                        GETA4(sphere),
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


  void visit(SP::SphereNEDS sphere)
  {
    float c[3], dFe, Cal;
    c[0] = .6 ;
    c[1] = .1;
    c[2] = .5;
    glLineWidth(3.);
    viewer_->drawSphere(GETX(sphere), GETY(sphere), GETZ(sphere),
                        GETA1(sphere), GETA2(sphere), GETA3(sphere), GETA4(sphere),
                        sphere->getRadius(), c);
    glColor3f(1., 0., 0.);
    dFe = hypot3d(GETXFE(sphere), GETYFE(sphere), GETZFE(sphere));
    Cal = log(dFe);

    //    viewer_->drawArrow(GETX(sphere), GETY(sphere), GETZ(sphere),
    //                       GETX(sphere)+Cal*GETXFE(sphere)/dFe,
    //                       GETY(sphere)+Cal*GETYFE(sphere)/dFe,
    //                       GETZ(sphere)+Cal*GETZFE(sphere), Cal/10.);
  }

};

#endif
