
#include <IOConfig.h>

#define HAVE_SICONOS_MECHANICS
#include <VisitorMaker.hpp>

#include <SpaceFilter.hpp>
#include <BlockVector.hpp>
#include <Interaction.hpp>

#include <Question.hpp>

#ifdef HAVE_BULLET
#include <BulletDS.hpp>
#include <BulletR.hpp>
#include <BulletSpaceFilter.hpp>
#include <btBulletCollisionCommon.h>
#endif


#include "MechanicsIO.hpp"

#include <SiconosGraph.hpp>
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


using namespace Alternative;

struct GetPosition : public SiconosVisitor
{

  SP::SiconosVector result;

  template<typename T>
  void operator()(const T& ds)
  {
    result = ds.q();
  }
};

struct GetVelocity : public SiconosVisitor
{

  SP::SiconosVector result;

  template<typename T>
  void operator()(const T& ds)
  {
    result = ds.velocity();
  }
};

struct ForMu : public Question<double>
{
  ANSWER(NewtonImpactFrictionNSL, mu());
};

struct ContactPointVisitor : public SiconosVisitor
{
  const Interaction& inter;
  const VectorOfVectors& workVectors;
  SiconosVector answer;

  ContactPointVisitor(const Interaction& inter,
                      const VectorOfVectors& workVectors) :
    inter(inter), workVectors(workVectors) {};

#ifdef HAVE_BULLET
  void visit(const BulletR& rel)
  {
    answer.resize(14);
    btManifoldPoint& cp = *rel.contactPoint();
    const btVector3& posa = cp.getPositionWorldOnA();
    const btVector3& posb = cp.getPositionWorldOnB();
    const SiconosVector& nc = *rel.nc();
    const SimpleMatrix& jachqT = *rel.jachqT();
    double id = (size_t) &*rel.contactPoint();
    double mu = ask<ForMu>(*inter.nslaw());
    SiconosVector cf(workVectors[NewtonEulerR::p1]->size());
    prod(*inter.lambda(1), jachqT, cf, true);
    answer.setValue(0, mu);
    answer.setValue(1, posa[0]);
    answer.setValue(2, posa[1]);
    answer.setValue(3, posa[2]);
    answer.setValue(4, posb[0]);
    answer.setValue(5, posb[1]);
    answer.setValue(6, posb[2]);
    answer.setValue(7, nc(0));
    answer.setValue(8, nc(1));
    answer.setValue(9, nc(2));
    answer.setValue(10, cf(0));
    answer.setValue(11, cf(1));
    answer.setValue(12, cf(2));
    answer.setValue(13, id);

  }
#endif
};

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

SP::SimpleMatrix MechanicsIO::contactPoints(const Model& model) const
{
  SP::SimpleMatrix result(new SimpleMatrix());
  InteractionsGraph::VIterator vi, viend;
  if (model.nonSmoothDynamicalSystem()->topology()->numberOfIndexSet() > 0)
  {
    InteractionsGraph& graph =
      *model.nonSmoothDynamicalSystem()->topology()->indexSet(1);
    unsigned int current_row;
    result->resize(graph.vertices_number(), 14);
    for(current_row=1, std11::tie(vi,viend) = graph.vertices();
        vi!=viend; ++vi, ++current_row)
    {
      const Interaction& inter = *graph.bundle(*vi);
      const VectorOfVectors& workVectors = *graph.properties(*vi).workVectors;
      ContactPointVisitor visitor(inter, workVectors);
      inter.relation()->accept(visitor);
      const SiconosVector& data = visitor.answer;
      result->setRow(current_row, data);
    }
  }
  return result;
}
