#define FROM_IMPL
#include <IOConfig.h>

#define HAVE_SICONOS_MECHANICS
#include <VisitorMaker.hpp>

#include <SpaceFilter.hpp>
#include <BlockVector.hpp>

#ifdef HAVE_BULLET
#include <BulletDS.hpp>
#include <BulletR.hpp>
#include <BulletSpaceFilter.hpp>
#include <btBulletCollisionCommon.h>
#endif


#include "MechanicsIO.hpp"

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



struct ContactPointVisitor : public SiconosVisitor
{
  const Interaction& inter;
  SiconosVector answer;

  ContactPointVisitor(Interaction& inter) : inter(inter) {};

#ifdef HAVE_BULLET
  void visit(const BulletR& rel)
  {
    answer.resize(9);
    btManifoldPoint& cp = *rel.contactPoint();
    const btVector3& posa = cp.getPositionWorldOnA();
    const SiconosVector& nc = *rel.nc();
    const SimpleMatrix& jachqT = *rel.jachqT();
    SiconosVector cf(7);
    prod(*inter.lambda(1), jachqT, cf, true);
    answer.setValue(0, posa[0]);
    answer.setValue(1, posa[1]);
    answer.setValue(2, posa[2]);
    answer.setValue(3, nc(0));
    answer.setValue(4, nc(1));
    answer.setValue(5, nc(2));
    answer.setValue(6, cf(0));
    answer.setValue(7, cf(1));
    answer.setValue(8, cf(2));
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
  DynamicalSystemsGraph::EIterator ei, eiend;
  const DynamicalSystemsGraph& graph = 
    *model.nonSmoothDynamicalSystem()->topology()->dSG(0);
  unsigned int current_row;
  for(current_row=1,std11::tie(ei,eiend)=graph.edges();
      ei!=eiend; ++ei, ++current_row)
  {
    Interaction& inter = *graph.bundle(*ei);
    ContactPointVisitor visitor(inter);
    inter.relation()->accept(visitor);
    const SiconosVector& data = visitor.answer;
    result->resize(current_row, data.size());
    result->setRow(current_row, data);
  }
  return result;
}

