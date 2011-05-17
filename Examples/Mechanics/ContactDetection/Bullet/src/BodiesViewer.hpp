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

#include <BulletDS.hpp>
#include <BulletR.hpp>

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

/* convenient macros */


#define GETALLDS(M) M->model()->nonSmoothDynamicalSystem()->topology()->dSG(0)
#define GETNDS(M) GETALLDS(M)->size()

#ifndef WITH_LMGC
struct Lmgc2DPOLYG : public LagrangianDS
{
  SP::SimpleVector vertices() const
  {
    assert(0);
  };
  int vertices_number() const
  {
    assert(0);
  };
};
struct Lmgc2DDSK : public LagrangianDS
{
  double getRadius() const
  {
    assert(0);
  };
};

struct Lmgc2DR : public LagrangianR {};
#endif


struct ForVertices : public Question<SP::SimpleVector>
{
  ANSWER(Lmgc2DPOLYG, vertices());
};

struct ForVerticesNumber : public Question<int>
{
  ANSWER(Lmgc2DPOLYG, vertices_number());
};

struct ForNdof : public Question<unsigned int>
{
  ANSWER_V(Disk, 3);
  ANSWER_V(Circle, 3);
  ANSWER_V(SphereLDS, 6);
  ANSWER_V(SphereNEDS, 6);
  ANSWER_V(BulletDS, 6);
  ANSWER_V(Lmgc2DDSK, 3);
  ANSWER_V(Lmgc2DPOLYG, 3);
};


struct ForFExt : public Question<SP::SiconosVector>
{

  ANSWER(Disk, fExt());
  ANSWER(Circle, fExt());
  ANSWER(SphereLDS, fExt());
  ANSWER(SphereNEDS, fExt());
  ANSWER(BulletDS, fExt());
  ANSWER(Lmgc2DDSK, fExt());
  ANSWER(Lmgc2DPOLYG, fExt());
};

struct ForPosition : public Question<SP::SiconosVector>
{

  ANSWER(Disk, q());
  ANSWER(Circle, q());
  ANSWER(SphereLDS, q());
  ANSWER(SphereNEDS, q());
  ANSWER(BulletDS, q());
  ANSWER(Lmgc2DDSK, q());
  ANSWER(Lmgc2DPOLYG, q());

};


struct ForRadius : public Question<double>
{

  ANSWER(Disk, getRadius());
  ANSWER(Circle, getRadius());
  ANSWER(SphereLDS, getRadius());
  ANSWER(SphereNEDS, getRadius());
  ANSWER(Lmgc2DDSK, getRadius());
  ANSWER_V(Lmgc2DPOLYG, 0.);
  ANSWER_V(BulletDS, 0.); // fix
};

struct ForMassValue : public Question<double>
{
  ANSWER(Disk, mass()->getValue(0, 0));
  ANSWER(Circle, mass()->getValue(0, 0));
  ANSWER(SphereLDS, mass()->getValue(0, 0));
  ANSWER(SphereNEDS, massValue());
  ANSWER(BulletDS, massValue());
  ANSWER(Lmgc2DDSK, mass()->getValue(0, 0));;
  ANSWER(Lmgc2DPOLYG, mass()->getValue(0, 0));;
};

struct ForJachq : public Question<SP::SiconosMatrix>
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
  ANSWER(BulletR, jachq());
  ANSWER(Lmgc2DR, jachq());
};

struct ForContactForce : public Question<SP::SimpleVector>
{
  ANSWER(BulletR, contactForce());
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
  SPHERELDS,
  SPHERENEDS,
  POLYG,
  POLYH,
  BULLET
};

struct ForShape : public Question<SHAPE>
{
  ANSWER_V(Disk, DISK);
  ANSWER_V(Circle, CIRCLE);
  ANSWER_V(SphereLDS, SPHERE);
  ANSWER_V(SphereNEDS, SPHERE);
  ANSWER_V(Lmgc2DDSK, DISK);
  ANSWER_V(Lmgc2DPOLYG, POLYG);
  ANSWER_V(BulletDS, BULLET);
};

/* dynamical system / figure association */
class QGLShape
{

public:

  /* construction from a LagrangianDS */
  QGLShape(SHAPE f, SP::DynamicalSystem D)
  {
    assert(D);

    figure_ = f;
    DS_ = D;
    frame_.reset(new qglviewer::ManipulatedFrame());
    savedFExt_.reset(new SimpleVector(ask<ForNdof>(*DS())));
    selected_ = false;
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
    *savedFExt_ = *(ask<ForFExt>(*DS()));
  };
  void restoreFExt()
  {
    *(ask<ForFExt>(*DS())) = *savedFExt_;
  };

  /* DS frame */
  qglviewer::ManipulatedFrame * getFrame()
  {
    return frame_.get();
  };

  SHAPE kind() const
  {
    return figure_;
  }

protected:
  int id_;
  bool selected_;
  SP::SiconosVector savedFExt_;
  boost::shared_ptr<qglviewer::ManipulatedFrame> frame_;

  SHAPE figure_;
  SP::DynamicalSystem DS_;
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
  virtual void init() = 0;
  virtual void animate() = 0;
  virtual void mousePressEvent(QMouseEvent *) = 0;
  virtual void mouseMoveEvent(QMouseEvent *) = 0;
  virtual void mouseReleaseEvent(QMouseEvent *) = 0;

  virtual QString helpString() const = 0;
  boost::shared_ptr<qglviewer::WorldConstraint> constraint_;

  SP::SiconosBodies Siconos_;

  std::vector<SP::QGLShape>  shapes_;

  void print(float x, float y, const char *s, int size);

  int NDS_;

  qglviewer::Vec selectedPoint_;

  int lastSelected_;

  bool myMouseBehavior_;

  long timeSiconos_;

  long timeGlob_;

};



#endif
