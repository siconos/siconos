
#include "SiconosConfig.h"

#include "MechanicsIO.hpp"

#define DUMMY(X, Y) class X : public Y {}

#undef PROPOSED_CLASSES
#undef BULLET_CLASSES
#undef OCC_CLASSES
#undef MECHANISMS_CLASSES

#define PROPOSED_CLASSES() \
  REGISTER(BodyDS)

#include <BodyDS.hpp>

#define XBULLET_CLASSES() \
  REGISTER(BulletDS)     \
  REGISTER(BulletR) \
  REGISTER(BulletSpaceFilter)

#ifdef HAVE_BULLET
#include <BulletDS.hpp>
#include <BulletR.hpp>
#include <BulletSpaceFilter.hpp>
#else
#include <NewtonEulerDS.hpp>
#include <NewtonEulerFrom3DLocalFrameR.hpp>
#include <SpaceFilter.hpp>
DUMMY(BulletDS, NewtonEulerDS);
DUMMY(BulletR, NewtonEulerFrom3DLocalFrameR);
DUMMY(BulletSpaceFilter, SpaceFilter);
#endif

#define OCC_CLASSES() \
  REGISTER(OccBody) \
  REGISTER(OccR)
#ifdef HAVE_OCC
#include <OccBody.hpp>
#include <OccR.hpp>
#else
#include <NewtonEulerDS.hpp>
#include <NewtonEulerFrom3DLocalFrameR.hpp>
DUMMY(OccBody, NewtonEulerDS);
DUMMY(OccR, NewtonEulerFrom3DLocalFrameR);
#endif

#define MECHANISMS_CLASSES() \
  REGISTER(MBTB_FC3DContactRelation) \
  REGISTER(MBTB_ContactRelation)
#ifdef HAVE_MECHANISMS
#include <MBTB_FC3DContactRelation.hpp>
#include <MBTB_ContactRelation.hpp>
#else
#include <NewtonEulerFrom3DLocalFrameR.hpp>
#include <NewtonEulerFrom1DLocalFrameR.hpp>
DUMMY(MBTB_FC3DContactRelation, NewtonEulerFrom3DLocalFrameR);
DUMMY(MBTB_ContactRelation, NewtonEulerFrom1DLocalFrameR);
#endif

#define VISITOR_CLASSES() \
  REGISTER(DynamicalSystem)                     \
  REGISTER(LagrangianDS)                        \
  REGISTER(NewtonEulerDS)                       \
  REGISTER(LagrangianR)                         \
  REGISTER(Disk)                                \
  REGISTER(Circle)                              \
  REGISTER(NewtonEulerFrom1DLocalFrameR)        \
  REGISTER(NewtonEulerFrom3DLocalFrameR)        \
  REGISTER(PivotJointR)                         \
  REGISTER(KneeJointR)                          \
  REGISTER(PrismaticJointR)                     \
  PROPOSED_CLASSES()                            \
  OCC_CLASSES()                                 \
  XBULLET_CLASSES()                             \
  MECHANISMS_CLASSES()                          \
  REGISTER(NewtonEulerR)

#undef SICONOS_VISITABLES
#define SICONOS_VISITABLES() VISITOR_CLASSES()

#include <BlockVector.hpp>
#include <Question.hpp>

#include <LagrangianDS.hpp>
#include <NewtonEulerDS.hpp>

/* ... */
/* to be fixed: forward mess with mpl::is_base_of who needs fully
 * declared classes */
#include <SiconosKernel.hpp>

/* Mechanics visitables bodies */
#include "Circle.hpp"
#include "Disk.hpp"
#include "DiskDiskR.hpp"
#include "CircleCircleR.hpp"
#include "DiskPlanR.hpp"
#include "DiskMovingPlanR.hpp"
#include "SphereLDS.hpp"
#include "SphereLDSSphereLDSR.hpp"
#include "SphereNEDSSphereNEDSR.hpp"
#include "SphereLDSPlanR.hpp"
#include "SphereNEDS.hpp"
#include "SphereNEDSPlanR.hpp"
#include "ExternalBody.hpp"

#include <PivotJointR.hpp>
#include <KneeJointR.hpp>
#include <PrismaticJointR.hpp>

#include <VisitorMaker.hpp>

//#define DEBUG_MESSAGES 1
#include <debug.h>

using namespace Experimental;

struct GetPosition : public SiconosVisitor
{

  SP::SiconosVector result;

  template<typename T>
  void operator()(const T& ds)
  {
    result.reset(new SiconosVector(1+ds.q()->size()));
    result->setValue(0, ds.number());
    result->setBlock(1, *ds.q());
  }
};

struct GetVelocity : public SiconosVisitor
{

  SP::SiconosVector result;

  template<typename T>
  void operator()(const T& ds)
  {
    result.reset(new SiconosVector(1+ds.velocity()->size()));
    result->setValue(0, ds.number());
    result->setBlock(1, *ds.velocity());
  }
};

struct ForMu : public Question<double>
{
    ANSWER(NewtonImpactFrictionNSL, mu());
    ANSWER_V_NOUSING(NewtonImpactNSL, 0.);
};

/* template partial specilization is not possible inside struct, so we
 * need an helper function */
template<typename T>
void contactPointProcess(SiconosVector& answer,
                         const Interaction& inter,
                         const T& rel)
{

  answer.resize(14);
  const SiconosVector& posa = *rel.pc1();
  const SiconosVector& posb = *rel.pc2();
  const SiconosVector& nc = *rel.nc();
  const SimpleMatrix& jachqT = *rel.jachqT();
  double id = inter.number();
  double mu = ask<ForMu>(*inter.nonSmoothLaw());
  SiconosVector cf(jachqT.size(1));
  prod(*inter.lambda(1), jachqT, cf, true);
  answer.setValue(0, mu);

  DEBUG_PRINTF("posa(0)=%g\n", posa(0));
  DEBUG_PRINTF("posa(1)=%g\n", posa(1));
  DEBUG_PRINTF("posa(2)=%g\n", posa(2));


  answer.setValue(1, posa(0));
  answer.setValue(2, posa(1));
  answer.setValue(3, posa(2));
  answer.setValue(4, posb(0));
  answer.setValue(5, posb(1));
  answer.setValue(6, posb(2));
  answer.setValue(7, nc(0));
  answer.setValue(8, nc(1));
  answer.setValue(9, nc(2));
  answer.setValue(10, cf(0));
  answer.setValue(11, cf(1));
  answer.setValue(12, cf(2));
  answer.setValue(13, id);
};

template<>
void contactPointProcess<PivotJointR>(SiconosVector& answer,
                                      const Interaction& inter,
                                      const PivotJointR& rel)
{
};

template<>
void contactPointProcess<KneeJointR>(SiconosVector& answer,
                                     const Interaction& inter,
                                     const KneeJointR& rel)
{
};

template<>
void contactPointProcess<PrismaticJointR>(SiconosVector& answer,
                                          const Interaction& inter,
                                          const PrismaticJointR& rel)
{
};

struct ContactPointVisitor : public SiconosVisitor
{
  SP::Interaction inter;
  SiconosVector answer;

  template<typename T>
  void operator()(const T& rel)
  {
    contactPointProcess<T>(answer, *inter, rel);
  }

};

struct ContactPointDomainVisitor : public SiconosVisitor
{
  SP::Interaction inter;
  SiconosVector answer;

  template<typename T>
  void operator()(const T& rel)
  {
  }

};

template<>
void ContactPointDomainVisitor::operator()(const BulletR& rel)
{
  answer.resize(2);

  /*
   * TODO: contact point domain coloring (e.g. based on broadphase).
   * currently, domain = (x>0):1?0
   */
  answer.setValue(0, rel.pc1()->getValue(0) > 0);

  answer.setValue(1, inter->number());
}

template<typename T, typename G>
SP::SimpleMatrix MechanicsIO::visitAllVerticesForVector(const G& graph) const
{
  SP::SimpleMatrix result(new SimpleMatrix());
  typename G::VIterator vi, viend;
  unsigned int current_row;
  for(current_row=0,std11::tie(vi,viend)=graph.vertices();
      vi!=viend; ++vi, ++current_row)
  {
    T getter;
    graph.bundle(*vi)->accept(getter);
    const SiconosVector& data = *getter.result;
    result->resize(current_row+1, data.size());
    result->setRow(current_row, data);
  }
  return result;
}

template<typename T, typename G>
SP::SiconosVector MechanicsIO::visitAllVerticesForDouble(const G& graph) const
{
  SP::SiconosVector result(new SiconosVector(graph.vertices_number()));
  typename G::VIterator vi, viend;
  unsigned int current_row;
  for(current_row=0,std11::tie(vi,viend)=graph.vertices();
      vi!=viend; ++vi, ++current_row)
  {
    T getter;
    graph.bundle(*vi)->accept(getter);
    result->setValue(current_row, *getter.result);
  }
  return result;
}


SP::SimpleMatrix MechanicsIO::positions(const Model& model) const
{

  typedef
    Visitor < Classes < LagrangianDS, NewtonEulerDS >,
              GetPosition >::Make Getter;

  return visitAllVerticesForVector<Getter>
    (*(model.nonSmoothDynamicalSystem()->topology()->dSG(0)));
};


SP::SimpleMatrix MechanicsIO::velocities(const Model& model) const
{
  typedef
    Visitor < Classes < LagrangianDS, NewtonEulerDS >,
              GetVelocity>::Make Getter;

  return visitAllVerticesForVector<Getter>
    (*model.nonSmoothDynamicalSystem()->topology()->dSG(0));
}

SP::SimpleMatrix MechanicsIO::contactPoints(const Model& model,
                                            unsigned int index_set) const
{
  SP::SimpleMatrix result(new SimpleMatrix());
  InteractionsGraph::VIterator vi, viend;
  if (model.nonSmoothDynamicalSystem()->topology()->numberOfIndexSet() > 0)
  {
    InteractionsGraph& graph =
      *model.nonSmoothDynamicalSystem()->topology()->indexSet(index_set);
    unsigned int current_row;
    result->resize(graph.vertices_number(), 14);
    for(current_row=0, std11::tie(vi,viend) = graph.vertices();
        vi!=viend; ++vi, ++current_row)
    {
      DEBUG_PRINTF("process interaction : %p\n", &*graph.bundle(*vi));

      typedef Visitor < Classes <
                          NewtonEulerFrom1DLocalFrameR,
                          NewtonEulerFrom3DLocalFrameR,
                          PrismaticJointR,
                          KneeJointR,
                          PivotJointR>,
                        ContactPointVisitor>::Make ContactPointInspector;
      ContactPointInspector inspector;
      inspector.inter = graph.bundle(*vi);
      graph.bundle(*vi)->relation()->accept(inspector);
      const SiconosVector& data = inspector.answer;
      if (data.size() == 14) result->setRow(current_row, data);
    }
  }
  return result;
}

SP::SimpleMatrix MechanicsIO::domains(const Model& model) const
{
  SP::SimpleMatrix result(new SimpleMatrix());
  InteractionsGraph::VIterator vi, viend;
  if (model.nonSmoothDynamicalSystem()->topology()->numberOfIndexSet() > 0)
  {
    InteractionsGraph& graph =
      *model.nonSmoothDynamicalSystem()->topology()->indexSet(1);
    unsigned int current_row;
    result->resize(graph.vertices_number(), 2);
    for(current_row=0, std11::tie(vi,viend) = graph.vertices();
        vi!=viend; ++vi, ++current_row)
    {
      DEBUG_PRINTF("process interaction : %p\n", &*graph.bundle(*vi));

      typedef Visitor < Classes <
                          NewtonEulerFrom1DLocalFrameR,
                          NewtonEulerFrom3DLocalFrameR,
                          PrismaticJointR,
                          KneeJointR,
                          PivotJointR>,
                        ContactPointDomainVisitor>::Make DomainInspector;
      DomainInspector inspector;
      inspector.inter = graph.bundle(*vi);
      graph.bundle(*vi)->relation()->accept(inspector);
      const SiconosVector& data = inspector.answer;
      if (data.size() == 2) result->setRow(current_row, data);
    }
  }
  return result;
}
