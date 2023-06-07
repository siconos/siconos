#include "SiconosConfig.h"
#include "MechanicsIO.hpp"
#include "SiconosAlgebraProd.hpp"

#define DUMMY(X, Y) class X : public Y {}

#undef BULLET_CLASSES
#undef OCC_CLASSES
#undef MECHANISMS_CLASSES


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


// The following classes are common classes of mechanics/src
#include <RigidBodyDS.hpp>
#include <RigidBody2dDS.hpp>
#include <ContactR.hpp>
#include <Contact5DR.hpp>
#include <Contact2dR.hpp>
#include <Contact2d3DR.hpp>
#include <BodyShapeRecord.hpp>

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
  REGISTER(ContactR)                            \
  REGISTER(Contact5DR)                          \
  REGISTER(Contact2dR)                          \
  REGISTER(Contact2d3DR)                        \
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
#include "siconos_debug.h"

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

struct ForE : public Question<double>
{
  using SiconosVisitor::visit;
  void visit(const NewtonImpactFrictionNSL& nsl)
  {
    answer = nsl . en();
  }
  void visit(const NewtonImpactRollingFrictionNSL& nsl)
  {
    answer = nsl . en();
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
void ContactPointVisitor::operator()(const DiskDiskR& rel)
{
  VectorOfBlockVectors& DSlink = inter->linkToDSVariables();
  auto& q = *DSlink[LagrangianR::q0];

  double x1 = q(0); double y1 = q(1); const double r1 = rel.getRadius1();
  double x2 = q(3); double y2 = q(4); const double r2 = rel.getRadius2();
  double d = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));

  double ncx = d > 0 ? (x2-x1)/d : 0;
  double ncy = d > 0 ? (y2-y1)/d : 0;

  double cpax = x1 + ncx * r1;
  double cpay = y1 + ncy * r1;

  double cpbx = x2 - ncx * r2;
  double cpby = y2 - ncy * r2;

  double id = inter->number();
  double mu = ask<ForMu>(*inter->nonSmoothLaw());
  const SimpleMatrix& jachq = *rel.jachq();
  SiconosVector cf(jachq.size(1));
  prod(*inter->lambda(1), jachq, cf, true);

  answer.resize(16);

  answer.setValue(0, mu);
  answer.setValue(1, cpax);
  answer.setValue(2, cpay);

  answer.setValue(3, cpbx);
  answer.setValue(4, cpby);

  answer.setValue(5, ncx);
  answer.setValue(6, ncy);

  answer.setValue(7, cf(0));
  answer.setValue(8, cf(1));

  answer.setValue(9, inter->y(0)->getValue(0));
  answer.setValue(10, inter->y(0)->getValue(1));

  answer.setValue(11, inter->y(1)->getValue(0));
  answer.setValue(12, inter->y(1)->getValue(1));

  answer.setValue(13, inter->lambda(1)->getValue(0));
  answer.setValue(14, inter->lambda(1)->getValue(1));

  answer.setValue(15, id);

}

// CircleCircleR should be named DiskCircleR
template<>
void ContactPointVisitor::operator()(const CircleCircleR& rel)
{
  VectorOfBlockVectors& DSlink = inter->linkToDSVariables();
  auto& q = *DSlink[LagrangianR::q0];

  double x1 = q(0); double y1 = q(1); const double r1 = rel.getRadius1();
  double x2 = q(3); double y2 = q(4); const double r2 = rel.getRadius2();
  double d = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));

  double ncx = d > 0 ? (x2-x1)/d : 0;
  double ncy = d > 0 ? (y2-y1)/d : 0;

  double cpax, cpay, cpbx, cpby;
  if (r1 < r2) // disk1 inside circle2
  {
    cpax = x1 - ncx * r1;
    cpay = y1 - ncy * r1;

    cpbx = x2 - ncx * r2;
    cpby = y2 - ncy * r2;
  }
  else // disk2 inside circle1
  {
    cpbx = x2 + ncx * r2;
    cpby = y2 + ncx * r2;

    cpax = x1 + ncx * r1;
    cpay = y1 + ncy * r1;
  }
  double id = inter->number();
  double mu = ask<ForMu>(*inter->nonSmoothLaw());
  const SimpleMatrix& jachq = *rel.jachq();
  SiconosVector cf(jachq.size(1));
  prod(*inter->lambda(1), jachq, cf, true);

  answer.resize(16);

  answer.setValue(0, mu);
  answer.setValue(1, cpax);
  answer.setValue(2, cpay);

  answer.setValue(3, cpbx);
  answer.setValue(4, cpby);

  answer.setValue(5, ncx);
  answer.setValue(6, ncy);

  answer.setValue(7, cf(0));
  answer.setValue(8, cf(1));

  answer.setValue(9, inter->y(0)->getValue(0));
  answer.setValue(10, inter->y(0)->getValue(1));

  answer.setValue(11, inter->y(1)->getValue(0));
  answer.setValue(12, inter->y(1)->getValue(1));

  answer.setValue(13, inter->lambda(1)->getValue(0));
  answer.setValue(14, inter->lambda(1)->getValue(1));

  answer.setValue(15, id);

}

template<>
void ContactPointVisitor::operator()(const DiskPlanR& rel)
{
  auto& DSlink = inter->linkToDSVariables();
  const auto& q0 = *DSlink[LagrangianR::q0];

  auto x1 = q0(0); auto y1 = q0(1); auto r1 = rel.getRadius();
  auto A = rel.getA(); auto B = rel.getB(); auto C = rel.getC();


  double x2 = - (A*C - B*B * x1 + A*B * y1) / (A*A + B*B);
  double y2 = - (B*C - A*A * y1 + A*B * x1) / (A*A + B*B);

  double d = (*(inter->y()[1]))(0) + r1;

  double ncx = d > 0 ? (x2-x1)/d : 0;
  double ncy = d > 0 ? (y2-y1)/d : 0;

  double cpax = x1 + ncx * r1;
  double cpay = y1 + ncy * r1;

  double id = inter->number();
  double mu = ask<ForMu>(*inter->nonSmoothLaw());
  const SimpleMatrix& jachq = *rel.jachq();
  SiconosVector cf(jachq.size(1));
  prod(*inter->lambda(1), jachq, cf, true);

  answer.resize(16);

  answer.setValue(0, mu);
  answer.setValue(1, cpax);
  answer.setValue(2, cpay);

  answer.setValue(3, x2);
  answer.setValue(4, y2);

  answer.setValue(5, ncx);
  answer.setValue(6, ncy);

  answer.setValue(7, cf(0));
  answer.setValue(8, cf(1));

  answer.setValue(9, inter->y(0)->getValue(0));
  answer.setValue(10, inter->y(0)->getValue(1));

  answer.setValue(11, inter->y(1)->getValue(0));
  answer.setValue(12, inter->y(1)->getValue(1));

  answer.setValue(13, inter->lambda(1)->getValue(0));
  answer.setValue(14, inter->lambda(1)->getValue(1));

  answer.setValue(15, id);

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
        Lagrangian2d3DR,
        CircleCircleR,
        DiskDiskR,
        DiskPlanR>,
      ContactPointVisitor>::Make ContactPointInspector;
      ContactPointInspector inspector;
      inspector.inter = graph.bundle(*vi);
      graph.bundle(*vi)->relation()->accept(inspector);
      SiconosVector& data = inspector.answer;
      data_size = data.size();

      if(data_size ==0)
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

/* Get contact informations */
/* default: a visitor that do nothing */
struct ContactInfoVisitor : public SiconosVisitor
{
  SP::Interaction inter;
  // std::vector<int> answer; better with a vector of int
  SiconosVector answer;

  template<typename T>
  void operator()(const T& rel)
  {
  }
};

/* then specializations : */
template<>
void ContactInfoVisitor::operator()(const NewtonEuler3DR& rel)
{
  double id = inter->number();
  answer.resize(10);
  // answer[0]= id;
  // answer[1]= 0; // reserve for ds1.number
  // answer[2]= 0; // reserve for ds2.number
  answer.setValue(0,id);
  answer.setValue(1,0);
  answer.setValue(2,0);

}


template<>
void ContactInfoVisitor::operator()(const ContactR& rel)
{
  double id = inter->number();
  answer.resize(4);
  answer.setValue(0,id);
  answer.setValue(1,0);
  answer.setValue(2,0);
  if (rel.bodyShapeRecordB->staticBody)
  {
    answer.setValue(3, rel.bodyShapeRecordB->staticBody->number);
  }
  else
    answer.setValue(3, 0);
}

template<>
void ContactInfoVisitor::operator()(const Contact5DR& rel)
{
  double id = inter->number();
  answer.resize(4);
  answer.setValue(0,id);
  answer.setValue(1,0);
  answer.setValue(2,0);
  if (rel.bodyShapeRecordB->staticBody)
  {
    answer.setValue(3, rel.bodyShapeRecordB->staticBody->number);
  }
  else
    answer.setValue(3, 0);
}

template<>
void ContactInfoVisitor::operator()(const Contact2dR& rel)
{
  double id = inter->number();
  answer.resize(4);
  answer.setValue(0,id);
  answer.setValue(1,0);
  answer.setValue(2,0);
  if (rel.bodyShapeRecordB->staticBody)
  {
    answer.setValue(3, rel.bodyShapeRecordB->staticBody->number);
  }
  else
    answer.setValue(3, 0);
}

template<>
void ContactInfoVisitor::operator()(const Contact2d3DR& rel)
{
  double id = inter->number();
  answer.resize(4);
  answer.setValue(0,id);
  answer.setValue(1,0);
  answer.setValue(2,0);
  if (rel.bodyShapeRecordB->staticBody)
  {
    answer.setValue(3, rel.bodyShapeRecordB->staticBody->number);
  }
  else
    answer.setValue(3, 0);
}




SP::SimpleMatrix MechanicsIO::contactInfo(const NonSmoothDynamicalSystem& nsds,
    unsigned int index_set) const
{
  DEBUG_BEGIN("SP::SimpleMatrix MechanicsIO::contactInfo");
  SP::SimpleMatrix result(new SimpleMatrix());
  InteractionsGraph::VIterator vi, viend;
  if(nsds.topology()->numberOfIndexSet() > 0)
  {
    InteractionsGraph& graph =
      *nsds.topology()->indexSet(index_set);
    unsigned int current_row;
    result->resize(graph.vertices_number(), 4);

    int data_size =0;
    for(current_row=0, std::tie(vi,viend) = graph.vertices();
        vi!=viend; ++vi)
    {
      DEBUG_PRINTF("process interaction : %p\n", &*graph.bundle(*vi));

      /* create a visitor for specified classes */
      typedef Visitor < Classes < NewtonEuler3DR,
                                  ContactR,
                                  Contact5DR,
                                  Contact2dR,
                                  Contact2d3DR>,
                        ContactInfoVisitor>::Make ContactInfoInspector;
      ContactInfoInspector inspector;
      inspector.inter = graph.bundle(*vi);
      graph.bundle(*vi)->relation()->accept(inspector);
      SiconosVector& data = inspector.answer;
      data_size = data.size();

      if(data_size ==0)
      {
        // Nothing is done since the relation does not appear as a relation
        // related to a contact points (perhaps a joint)
      }
      else
      {

        // We add at the end the number of ds1 and ds2
        DEBUG_EXPR(data.display(););
        DynamicalSystem& ds1 = *graph.properties(*vi).source;
        DynamicalSystem& ds2 = *graph.properties(*vi).target;
        data.setValue(1, ds1.number());
        data.setValue(2, ds2.number());
      }
      if(result->size(1) != data.size())
      {
        result->resize(graph.vertices_number(), data.size());
      }
      result->setRow(current_row++, data);

    }
    result->resize(current_row, data_size);
    DEBUG_EXPR(result->display(););

  }
  DEBUG_END("SP::SimpleMatrix MechanicsIO::contactInfo");

  return result;
}

/* Get contact work information */
/* default: a visitor that do nothing */

struct ContactContactWorkVisitor : public SiconosVisitor
{
  SP::Interaction inter;
  // std::vector<int> answer; better with a vector of int
  SiconosVector answer;
  double tol;
  template<typename T>
  void operator()(const T& rel)
  {
  }
};

/* then specializations : */
template<>
void ContactContactWorkVisitor::operator()(const NewtonEuler3DR& rel)
{
  double id = inter->number();
  answer.resize(6);
}

static void compute_contact_work_and_status(SP::Interaction inter, double tol, SiconosVector& answer) {
  double mu = ask<ForMu>(*inter->nonSmoothLaw());
  double e = ask<ForE>(*inter->nonSmoothLaw());
  // Compute normal contact work
  double vn_minus =  inter->y_k(1).getValue(0);
  double vn_plus = inter->y(1)->getValue(0);
  double pn  = inter->lambda(1)->getValue(0);

  double normal_contact_work = 0.5 * ( vn_minus + vn_plus)*pn;
  answer.setValue(1,normal_contact_work);

  // Compute tangent contact work of impulse

  double vt_1_minus =  inter->y_k(1).getValue(1);
  double vt_2_minus =  inter->y_k(1).getValue(2);
  double vt_1_plus = inter->y(1)->getValue(1);
  double vt_2_plus = inter->y(1)->getValue(2);

  double pt_1  = inter->lambda(1)->getValue(1);
  double pt_2  = inter->lambda(1)->getValue(2);

  double tangent_contact_work = 0.5 * ( vt_1_minus + vt_1_plus)*pt_1 + 0.5 * (vt_2_minus + vt_2_plus)*pt_2;
  answer.setValue(2,tangent_contact_work);

  // Compute directly work dissipated by friction impulse
  double theta=1/2.;
  double norm_vt_plus =  sqrt(vt_1_plus*vt_1_plus+vt_2_plus*vt_2_plus);
  double norm_vt_minus = sqrt(vt_1_minus*vt_1_minus+vt_2_minus*vt_2_minus);

  double friction_dissipation = mu* (theta * norm_vt_plus  + (1-theta)*norm_vt_minus) * pn;
  answer.setValue(3,friction_dissipation);
  // compute contact status

  double norm_pt = sqrt(pt_1*pt_1 + pt_2*pt_2);

  if ( (pn < tol ) and (vn_plus + e * vn_minus > tol) )
    answer.setValue(4,0);// take-off = 0
  else if ( (pn > tol ) and (vn_plus + e * vn_minus  < tol) )
    {
      if (     (norm_pt  - mu * pn >  tol))
	{
	  //std::cout << "WARNING: the impulse is outside the Coulomb cone  " << std::endl;
	  answer.setValue(4,-3);// outside the cone = -3
	}
      else if ( (norm_pt  - mu*pn < - tol))
	{
	  //std::cout << "the impulse is in the *interior* of  the Coulomb cone  " << std::endl;
	  //std::cout << "norm_vt_plus  " << norm_vt_plus << std::endl;
	  if (norm_vt_plus > tol)
	    {
	      //std::cout << "WARNING: but the norm of vt is not zero  " << std::endl;
	      answer.setValue(4,-2);// sticking with a non zero slifing velocity = -2
	    }
	  answer.setValue(4,1);// sticking = 1
	}
      else
	{
	  //std::cout << "the impulse is on the *boundary* of the Coulomb cone  " << std::endl;
	  //std::cout << "norm_vt_plus  " << norm_vt_plus << std::endl;
	  answer.setValue(4,2); // sliding = 2
	}
    }
  else
    answer.setValue(4,-1);// undetermined = -1

  if ( (pn > tol ) and (vn_minus  > tol) )
    {
      std::cout << "WARNING: we apply the impact law of positive velocity " << std::endl;
      std::cout << "pn " << pn << " vn minus " << vn_minus << " vn plus " << vn_plus
		<< " normal_contact_work " << normal_contact_work
		<< " -e * vn_minus   " << -e*vn_minus 
		<< std::endl;
      answer.setValue(5, normal_contact_work);
    }
  // double id = inter->number();
  // std::cout << "\nid "<< id << std::endl;
  // std::cout << " e "<< e  << " mu "<< mu << std::endl;
  // std::cout << " tol "<< tol<< std::endl;
  // std::cout << "vn_plus "<< vn_plus << std::endl;
  // std::cout << "vn_minus "<< vn_minus << std::endl;
  // std::cout << "pn "<< pn << std::endl;
  // std::cout << "normal_contact_work  "<< normal_contact_work  << std::endl;

  // std::cout << "vt_plus "<< vt_1_plus << " " << vt_2_plus <<  std::endl;
  // std::cout << "vt_minus "<< vt_1_minus << " " << vt_2_minus <<  std::endl;
  // std::cout << "pt "<< pt_1 << " "  << pt_2 << std::endl;
  // std::cout << "tangent_contact_work  "<< tangent_contact_work  << std::endl;

  // std::cout << "friction_dissipation  "<< friction_dissipation << std::endl;

  // std::cout << "norm_pt  "<< norm_pt  << std::endl;
  // std::cout << "norm_pt - mu* pn  "<< norm_pt -mu *pn   << std::endl;
  // std::cout << "vn_plus + e * vn_minus  " << vn_plus + e * vn_minus   << std::endl;
  // std::cout << "status   "<<   answer.getValue(4) << std::endl;
}


template<>
void ContactContactWorkVisitor::operator()(const ContactR& rel)
{
  double id = inter->number();
  answer.resize(6);
  answer.setValue(0,id);
  compute_contact_work_and_status(inter,  tol, answer);
}

template<>
void ContactContactWorkVisitor::operator()(const Contact5DR& rel)
{
  double id = inter->number();
  answer.resize(6);
  answer.setValue(0,id);


  compute_contact_work_and_status(inter,  tol, answer);

}

template<>
void ContactContactWorkVisitor::operator()(const Contact2dR& rel)
{
  double id = inter->number();
  answer.resize(6);
  answer.setValue(0,id);


  compute_contact_work_and_status(inter,  tol, answer);

}

template<>
void ContactContactWorkVisitor::operator()(const Contact2d3DR& rel)
{
  double id = inter->number();
  answer.resize(6);
  answer.setValue(0,id);


  compute_contact_work_and_status(inter,  tol, answer);

}


SP::SimpleMatrix MechanicsIO::contactContactWork(const NonSmoothDynamicalSystem& nsds,
						 unsigned int index_set, double tol) const
{
  DEBUG_BEGIN("SP::SimpleMatrix MechanicsIO::contactContactWork");
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
      typedef Visitor < Classes < NewtonEuler3DR,
                                  ContactR,
                                  Contact5DR,
                                  Contact2dR,
                                  Contact2d3DR>,
                        ContactContactWorkVisitor>::Make ContactContactWorkInspector;
      ContactContactWorkInspector inspector;
      inspector.inter = graph.bundle(*vi);
      inspector.tol = tol;
      graph.bundle(*vi)->relation()->accept(inspector);
      SiconosVector& data = inspector.answer;
      data_size = data.size();

      if(data_size ==0)
      {
        // Nothing is done since the relation does not appear as a relation
        // related to a contact points (perhaps a joint)
      }
      else
      {

      }
      if(result->size(1) != data.size())
      {
        result->resize(graph.vertices_number(), data.size());
      }
      result->setRow(current_row++, data);

    }
    result->resize(current_row, data_size);
    DEBUG_EXPR(result->display(););

  }
  DEBUG_END("SP::SimpleMatrix MechanicsIO::contactContactWork");

  //result->display();
  return result;
}
