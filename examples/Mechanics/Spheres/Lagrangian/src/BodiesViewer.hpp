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

#ifndef BodiesViewer_hpp
#define BodiesViewer_hpp

/* Siconos */
#include <Disk.hpp>
#include <DiskDiskR.hpp>
#include <DiskPlanR.hpp>
#include <DiskMovingPlanR.hpp>
#include <Circle.hpp>
#include <CircleCircleR.hpp>
#include <SphereLDS.hpp>
#include <SphereLDSPlanR.hpp>
#include <SphereLDSSphereLDSR.hpp>
#include <SphereNEDS.hpp>
#include <SphereNEDSPlanR.hpp>
#include <SphereNEDSSphereNEDSR.hpp>
#include <Question.hpp>

#if WITH_BULLET
#include <BulletDS.hpp>
#include <BulletR.hpp>
#define IFBULLET(X) X
#else
#define IFBULLET(X)
#endif

#include "SiconosBodies.hpp"
#include "SiconosKernel.hpp"

/* QGLViewer */
#include <qglviewer.h>
#include <manipulatedCameraFrame.h>
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

/* convenient macros */


#define GETALLDS(M) M->simulation()->nonSmoothDynamicalSystem()->topology()->dSG(0)
#define GETNDS(M) GETALLDS(M)->size()

struct ForNdof : public Question<unsigned int>
{
  using SiconosVisitor::visit;

  ANSWER_V(Disk, 3);
  ANSWER_V(Circle, 3);
  ANSWER_V(SphereLDS, 6);
  ANSWER_V(SphereNEDS, 6);
  IFBULLET(ANSWER_V(BulletDS, 6));
};


struct ForFExt : public Question<SP::SiconosVector>
{
  using SiconosVisitor::visit;

  ANSWER(Disk, fExt());
  ANSWER(Circle, fExt());
  ANSWER(SphereLDS, fExt());
  ANSWER(SphereNEDS, fExt());
  IFBULLET(ANSWER(BulletDS, fExt()));
};


struct ForPosition : public Question<SP::SiconosVector>
{
  using SiconosVisitor::visit;

  ANSWER(Disk, q());
  ANSWER(Circle, q());
  ANSWER(SphereLDS, q());
  ANSWER(SphereNEDS, q());
  IFBULLET(ANSWER(BulletDS, q()));
};

struct ForRadius : public Question<double>
{
  using SiconosVisitor::visit;

  ANSWER(Disk, getRadius());
  ANSWER(Circle, getRadius());
  ANSWER(SphereLDS, getRadius());
  ANSWER(SphereNEDS, getRadius());
  IFBULLET(ANSWER_V(BulletDS, 0.)); // fix
};

struct ForMassValue : public Question<double>
{
  using SiconosVisitor::visit;

  ANSWER(Disk, mass()->getValue(0, 0));
  ANSWER(Circle, mass()->getValue(0, 0));
  ANSWER(SphereLDS, mass()->getValue(0, 0));
  ANSWER(SphereNEDS, scalarMass());
  IFBULLET(ANSWER(BulletDS, scalarMass()));
};

struct ForJachq : public Question<SP::SiconosMatrix>
{
  using SiconosVisitor::visit;

  ANSWER(LagrangianR, jachq());
  ANSWER(NewtonEulerR, jachq());
  ANSWER(DiskDiskR, jachq());
  ANSWER(DiskPlanR, jachq());
  ANSWER(DiskMovingPlanR, jachq());
  ANSWER(SphereLDSPlanR, jachq());
  ANSWER(SphereNEDSPlanR, jachq());
  ANSWER(SphereLDSSphereLDSR, jachq());
  ANSWER(SphereNEDSSphereNEDSR, jachq());
  IFBULLET(ANSWER(BulletR, jachq()));
};

struct ForContactForce : public Question<SP::SiconosVector>
{
  using SiconosVisitor::visit;

  IFBULLET(ANSWER(BulletR, contactForce()));
};




#define GETX(C) ask<ForPosition>(*C)->getValue(0)
#define GETY(C) ask<ForPosition>(*C)->getValue(1)
#define GETZ(C) ask<ForPosition>(*C)->getValue(2)
#define GETA1(C) ask<ForPosition>(*C)->getValue(3)
#define GETA2(C) ask<ForPosition>(*C)->getValue(4)
#define GETA3(C) ask<ForPosition>(*C)->getValue(5)
#define GETA4(C) ask<ForPosition>(*C)->getValue(6)
#define GETXFE(C) ask<ForFExt>(*C)->getValue(0)
#define GETYFE(C) ask<ForFExt>(*C)->getValue(1)
#define GETZFE(C) ask<ForFExt>(*C)->getValue(2)

#define GETVX(C) ask<ForVelocity>(*C)->getValue(0)
#define GETVY(C) ask<ForVelocity>(*C)->getValue(1)
#define GETVZ(C) ask<ForVelocity>(*C)->getValue(2)

#define GETVA1(C) ask<ForVelocity>(*C)->getValue(3)
#define GETVA2(C) ask<ForVelocity>(*C)->getValue(4)
#define GETVA3(C) ask<ForVelocity>(*C)->getValue(5)

#define GETRADIUS(C) ask<ForRadius>(*C)


enum SHAPE
{
  DISK,
  CIRCLE,
  SPHERE,
  IFBULLET(BULLET)
};

struct ForShape : public Question<SHAPE>
{
  using SiconosVisitor::visit;

  ANSWER_V(Disk, DISK);
  ANSWER_V(Circle, CIRCLE);
  ANSWER_V(SphereLDS, SPHERE);
  ANSWER_V(SphereNEDS, SPHERE);
  IFBULLET(ANSWER_V(BulletDS, BULLET));
};

/* dynamical system / figure association */
class QGLShape
{

public:

  /* construction from a LagrangianDS */
  QGLShape(SHAPE f, SP::DynamicalSystem D, const qglviewer::Frame* ref)
  {
    assert(D);

    figure_ = f;
    DS_ = D;
    frame_.reset(new qglviewer::ManipulatedFrame());
    frame_->setReferenceFrame(ref);
    savedFExt_ = ask<ForFExt>(*DS());
    selected_ = false;
    saved_ = true;
    type_ = Type::value(*DS());

  };

  ~QGLShape() {};

  /* pointer to DS */
  SP::DynamicalSystem DS() const
  {
    return DS_;
  };

  /* selection with mouse */
  bool selected()
  {
    return selected_ ;
  };
  void setSelection(bool v)
  {
    selected_ = v ;
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
    savedFExt_ = ask<ForFExt>(*DS());
  };
  void restoreFExt()
  {
    switch (Type::value(*DS()))
    {
    case Type::NewtonEulerDS :
    {
      std11::static_pointer_cast<NewtonEulerDS>(DS())
      ->setFExtPtr(std11::static_pointer_cast<SiconosVector>(savedFExt_));
      break;
    }
    case Type::LagrangianDS :
    {
      std11::static_pointer_cast<LagrangianDS>(DS())
      ->setFExtPtr(std11::static_pointer_cast<SiconosVector>(savedFExt_));
      break;
    };
    default:
    {};
    }
  };

  /* DS frame */
  qglviewer::ManipulatedFrame * frame()
  {
    return frame_.get();
  };

  SHAPE shape() const
  {
    return figure_;
  }

  Type::Siconos type() const
  {
    return type_;
  }

protected:
  int id_;
  bool selected_;
  bool saved_;
  SP::SiconosVector savedFExt_;
  boost::shared_ptr<qglviewer::ManipulatedFrame> frame_;

  SHAPE figure_;
  SP::DynamicalSystem DS_;
  int positionSize_;

  Type::Siconos type_;

};

TYPEDEF_SPTR(QGLShape);

/* QGLViewer main object */
class BodiesViewer : public QGLViewer
{

public:

#ifdef QT_INTERFACE
  Viewer(QWidget *parent);
#endif

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
  void drawPolyg(unsigned int n, double* coor, float *c);

  virtual void drawQGLShape(const QGLShape&);
  virtual void drawSelectedQGLShape(const QGLShape&);
  void insertQGLShape(SHAPE, SP::DynamicalSystem);


protected :

  void postSelection(const QPoint& point);
  virtual void init();
  virtual void animate();
  virtual void mousePressEvent(QMouseEvent *);
  virtual void mouseMoveEvent(QMouseEvent *);
  virtual void mouseReleaseEvent(QMouseEvent *);
  virtual void keyPressEvent(QKeyEvent *);

  virtual QString helpString() const;
  boost::shared_ptr<qglviewer::WorldConstraint> constraint_;

  SP::SiconosBodies Siconos_;

  qglviewer::Frame* referenceFrame_;

  std::vector<SP::QGLShape>  shapes_;

  void print(float x, float y, const char *s, int size);

  int NDS_;

  qglviewer::Vec selectedPoint_;

  int lastSelected_;

  bool myMouseBehavior_;

  bool stepSimulation_;
  bool stepNow_;

  long timeSiconos_;

  long timeGlob_;

  float _transparency;

};



#endif
