#include "SiconosConfig.h"
#include "MechanicsIO.hpp"
#include "SiconosAlgebraProd.hpp"

#define DUMMY(X, Y) class X : public Y {}

#undef BULLET_CLASSES
#undef OCC_CLASSES
#undef MECHANISMS_CLASSES

#include <RigidBodyDS.hpp>
#include <RigidBody2dDS.hpp>

#define BULLET_CLASSES() \
  REGISTER(BulletR)\
  REGISTER(Bullet5DR)\
  REGISTER(Bullet2dR)\
  REGISTER(Bullet2d3DR)


#ifdef SICONOS_HAS_BULLET
#include <BulletR.hpp>
#include <Bullet5DR.hpp>
#include <Bullet2dR.hpp>
#include <Bullet2d3DR.hpp>
#else
#include <NewtonEuler3DR.hpp>
#include <NewtonEuler5DR.hpp>
#include <SpaceFilter.hpp>
DUMMY(BulletR, NewtonEuler3DR);
DUMMY(Bullet5DR, NewtonEuler5DR);
#include <Lagrangian2d2DR.cpp>
#include <Lagrangian2d3DR.cpp>
DUMMY(Bullet2dR, Lagrangian2d2DR);
DUMMY(Bullet2d3DR, Lagrangian2d3DR);
#endif

#define OCC_CLASSES() \
  REGISTER(OccBody) \
  REGISTER(OccR)
#ifdef SICONOS_HAS_OCE
#include <OccBody.hpp>
#include <OccR.hpp>
#else
#include <NewtonEulerDS.hpp>
#include <NewtonEuler3DR.hpp>
DUMMY(OccBody, NewtonEulerDS);
DUMMY(OccR, NewtonEuler3DR);
#endif

#define MECHANISMS_CLASSES() \
  REGISTER(MBTB_FC3DContactRelation) \
  REGISTER(MBTB_ContactRelation)

#ifdef HAVE_SICONOS_MECHANISMS
#include <MBTB_FC3DContactRelation.hpp>
#include <MBTB_ContactRelation.hpp>
#else
#include <NewtonEuler3DR.hpp>
#include <NewtonEuler1DR.hpp>
DUMMY(MBTB_FC3DContactRelation, NewtonEuler3DR);
DUMMY(MBTB_ContactRelation, NewtonEuler1DR);
#endif


/* all the classes that may be visited */
#define VISITOR_CLASSES()                       \
  REGISTER(DynamicalSystem)                     \
  REGISTER(LagrangianDS)                        \
  REGISTER(NewtonEulerDS)                       \
  REGISTER(LagrangianR)                         \
  REGISTER(Disk)                                \
  REGISTER(Circle)                              \
  REGISTER(Lagrangian2d2DR)                     \
  REGISTER(Lagrangian2d3DR)                     \
  REGISTER(NewtonEulerR)                        \
  REGISTER(NewtonEuler1DR)                      \
  REGISTER(NewtonEuler3DR)                      \
  REGISTER(NewtonEuler5DR)                      \
  REGISTER(PivotJointR)                         \
  REGISTER(KneeJointR)                          \
  REGISTER(PrismaticJointR)                     \
  REGISTER(RigidBodyDS)                         \
  REGISTER(RigidBody2dDS)                       \
  MECHANISMS_CLASSES()                          \
  OCC_CLASSES()                                 \
  BULLET_CLASSES()


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

/* all the visitable classes must have been included at this point */
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
  using SiconosVisitor::visit;
  void visit(const NewtonImpactFrictionNSL& nsl)
  {
    answer = nsl . mu();
  }
  void visit(const NewtonImpactRollingFrictionNSL& nsl)
  {
    answer = nsl . mu();
  }
  void visit(const NewtonImpactNSL& nsl)
  {
    answer = 0.;
  }
};


/* Get contact informations */
/* default: a visitor that do nothing */
struct ContactPointVisitor : public SiconosVisitor
{
  SP::Interaction inter;
  SiconosVector answer;

  template<typename T>
  void operator()(const T& rel)
  {
  }

};

/* then specializations : */
template<>
void ContactPointVisitor::operator()(const NewtonEuler3DR& rel)
{
  const SiconosVector& posa = *rel.pc1();
  const SiconosVector& posb = *rel.pc2();
  const SiconosVector& nc = *rel.nc();
  DEBUG_PRINTF("posa(0)=%g\n", posa(0));
  DEBUG_PRINTF("posa(1)=%g\n", posa(1));
  DEBUG_PRINTF("posa(2)=%g\n", posa(2));

  double id = inter->number();
  double mu = ask<ForMu>(*inter->nonSmoothLaw());
  const SimpleMatrix& jachqT = *rel.jachqT();
  SiconosVector cf(jachqT.size(1));
  prod(*inter->lambda(1), jachqT, cf, true);
  answer.resize(23);

  answer.setValue(0, mu);
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
  answer.setValue(13,inter->y(0)->getValue(0));
  answer.setValue(14,inter->y(0)->getValue(1));
  answer.setValue(15,inter->y(0)->getValue(2));
  answer.setValue(16,inter->y(1)->getValue(0));
  answer.setValue(17,inter->y(1)->getValue(1));
  answer.setValue(18,inter->y(1)->getValue(2));
  answer.setValue(19,inter->lambda(1)->getValue(0));
  answer.setValue(20,inter->lambda(1)->getValue(1));
  answer.setValue(21,inter->lambda(1)->getValue(2));
  answer.setValue(22, id);
}
/* then specializations : */
template<>
void ContactPointVisitor::operator()(const NewtonEuler5DR& rel)
{
  const SiconosVector& posa = *rel.pc1();
  const SiconosVector& posb = *rel.pc2();
  const SiconosVector& nc = *rel.nc();
  DEBUG_PRINTF("posa(0)=%g\n", posa(0));
  DEBUG_PRINTF("posa(1)=%g\n", posa(1));
  DEBUG_PRINTF("posa(2)=%g\n", posa(2));

  double id = inter->number();
  double mu = ask<ForMu>(*inter->nonSmoothLaw());
  const SimpleMatrix& jachqT = *rel.jachqT();
  SiconosVector cf(jachqT.size(1));
  prod(*inter->lambda(1), jachqT, cf, true);
  answer.resize(23);

  answer.setValue(0, mu);
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
  answer.setValue(13,inter->y(0)->getValue(0));
  answer.setValue(14,inter->y(0)->getValue(1));
  answer.setValue(15,inter->y(0)->getValue(2));
  answer.setValue(16,inter->y(1)->getValue(0));
  answer.setValue(17,inter->y(1)->getValue(1));
  answer.setValue(18,inter->y(1)->getValue(2));
  answer.setValue(19,inter->lambda(1)->getValue(0));
  answer.setValue(20,inter->lambda(1)->getValue(1));
  answer.setValue(21,inter->lambda(1)->getValue(2));
  answer.setValue(22, id);
}
template<>
void ContactPointVisitor::operator()(const Lagrangian2d2DR& rel)
{
  const SiconosVector& posa = *rel.pc1();
  const SiconosVector& posb = *rel.pc2();
  const SiconosVector& nc = *rel.nc();
  DEBUG_PRINTF("posa(0)=%g\n", posa(0));
  DEBUG_PRINTF("posa(1)=%g\n", posa(1));

  double id = inter->number();
  double mu = ask<ForMu>(*inter->nonSmoothLaw());
  const SimpleMatrix& jachq = *rel.jachq();
  SiconosVector cf(jachq.size(1));
  prod(*inter->lambda(1), jachq, cf, true);


  answer.resize(16);

  answer.setValue(0, mu);
  answer.setValue(1, posa(0));
  answer.setValue(2, posa(1));

  answer.setValue(3, posb(0));
  answer.setValue(4, posb(1));

  answer.setValue(5, nc(0));
  answer.setValue(6, nc(1));

  answer.setValue(7, cf(0));
  answer.setValue(8, cf(1));

  answer.setValue(9,inter->y(0)->getValue(0));
  answer.setValue(10,inter->y(0)->getValue(1));

  answer.setValue(11,inter->y(1)->getValue(0));
  answer.setValue(12,inter->y(1)->getValue(1));

  answer.setValue(13,inter->lambda(1)->getValue(0));
  answer.setValue(14,inter->lambda(1)->getValue(1));

  answer.setValue(15, id);
};

template<>
void ContactPointVisitor::operator()(const Lagrangian2d3DR& rel)
{
  const SiconosVector& posa = *rel.pc1();
  const SiconosVector& posb = *rel.pc2();
  const SiconosVector& nc = *rel.nc();
  DEBUG_PRINTF("posa(0)=%g\n", posa(0));
  DEBUG_PRINTF("posa(1)=%g\n", posa(1));

  double id = inter->number();
  double mu = ask<ForMu>(*inter->nonSmoothLaw());
  const SimpleMatrix& jachq = *rel.jachq();
  SiconosVector cf(jachq.size(1));
  prod(*inter->lambda(1), jachq, cf, true);


  answer.resize(16);

  answer.setValue(0, mu);
  answer.setValue(1, posa(0));
  answer.setValue(2, posa(1));

  answer.setValue(3, posb(0));
  answer.setValue(4, posb(1));

  answer.setValue(5, nc(0));
  answer.setValue(6, nc(1));

  answer.setValue(7, cf(0));
  answer.setValue(8, cf(1));

  answer.setValue(9,inter->y(0)->getValue(0));
  answer.setValue(10,inter->y(0)->getValue(1));

  answer.setValue(11,inter->y(1)->getValue(0));
  answer.setValue(12,inter->y(1)->getValue(1));

  answer.setValue(13,inter->lambda(1)->getValue(0));
  answer.setValue(14,inter->lambda(1)->getValue(1));

  answer.setValue(15, id);
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
void ContactPointDomainVisitor::operator()(const NewtonEuler3DR& rel)
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
  for(current_row=0,std::tie(vi,viend)=graph.vertices();
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
  for(current_row=0,std::tie(vi,viend)=graph.vertices();
      vi!=viend; ++vi, ++current_row)
  {
    T getter;
    graph.bundle(*vi)->accept(getter);
    result->setValue(current_row, *getter.result);
  }
  return result;
}


SP::SimpleMatrix MechanicsIO::positions(const NonSmoothDynamicalSystem& nsds) const
{

  typedef
  Visitor < Classes < LagrangianDS, NewtonEulerDS >,
          GetPosition >::Make Getter;

  return visitAllVerticesForVector<Getter>
         (*(nsds.topology()->dSG(0)));
};


SP::SimpleMatrix MechanicsIO::velocities(const NonSmoothDynamicalSystem& nsds) const
{
  typedef
  Visitor < Classes < LagrangianDS, NewtonEulerDS >,
          GetVelocity>::Make Getter;

  return visitAllVerticesForVector<Getter>
         (*nsds.topology()->dSG(0));
}

SP::SimpleMatrix MechanicsIO::contactPoints(const NonSmoothDynamicalSystem& nsds,
    unsigned int index_set) const
{
  SP::SimpleMatrix result(new SimpleMatrix());
  InteractionsGraph::VIterator vi, viend;
  if(nsds.topology()->numberOfIndexSet() > 0)
  {
    InteractionsGraph& graph =
      *nsds.topology()->indexSet(index_set);
    unsigned int current_row;
    result->resize(graph.vertices_number(), 25);

    int data_size =0;
    for(current_row=0, std::tie(vi,viend) = graph.vertices();
        vi!=viend; ++vi)
    {
      DEBUG_PRINTF("process interaction : %p\n", &*graph.bundle(*vi));

      /* create a visitor for specified classes */
      typedef Visitor < Classes <
      NewtonEuler1DR,
      NewtonEuler3DR,
      NewtonEuler5DR,
      Lagrangian2d2DR,
      Lagrangian2d3DR>,
      ContactPointVisitor>::Make ContactPointInspector;
      ContactPointInspector inspector;
      inspector.inter = graph.bundle(*vi);
      graph.bundle(*vi)->relation()->accept(inspector);
      SiconosVector& data = inspector.answer;
      data_size = data.size();

      if (data_size ==0)
      {
        // Nothing is done since the relation does not appear as a relation
        // related to a contact points (perhaps a joint)
      }
      else
      {
        // We add at the end the number of ds1 and ds2
        data.resize(data_size+2);
        DEBUG_EXPR(data.display(););
        DynamicalSystem& ds1 = *graph.properties(*vi).source;
        DynamicalSystem& ds2 = *graph.properties(*vi).target;

        data.setValue(data_size, ds1.number());
        data.setValue(data_size+1, ds2.number());
        DEBUG_EXPR(data.display(););
        if(result->size(1) != data.size())
        {
          result->resize(graph.vertices_number(), data.size());
        }
        result->setRow(current_row++, data);
        data_size +=2;
      }

    }
    result->resize(current_row, data_size);
    DEBUG_EXPR(result->display(););
  }

  return result;
}

SP::SimpleMatrix MechanicsIO::domains(const NonSmoothDynamicalSystem& nsds) const
{
  SP::SimpleMatrix result(new SimpleMatrix());
  InteractionsGraph::VIterator vi, viend;
  if(nsds.topology()->numberOfIndexSet() > 0)
  {
    InteractionsGraph& graph =
      *nsds.topology()->indexSet(1);
    unsigned int current_row;
    result->resize(graph.vertices_number(), 2);
    for(current_row=0, std::tie(vi,viend) = graph.vertices();
        vi!=viend; ++vi, ++current_row)
    {
      DEBUG_PRINTF("process interaction : %p\n", &*graph.bundle(*vi));

      typedef Visitor < Classes <
      NewtonEuler1DR,
      NewtonEuler3DR,
      PrismaticJointR,
      KneeJointR,
      PivotJointR>,
      ContactPointDomainVisitor>::Make DomainInspector;
      DomainInspector inspector;
      inspector.inter = graph.bundle(*vi);
      graph.bundle(*vi)->relation()->accept(inspector);
      const SiconosVector& data = inspector.answer;
      if(data.size() == 2) result->setRow(current_row, data);
    }
  }
  return result;
}



