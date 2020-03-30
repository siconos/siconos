#ifndef ADDONS_HPP
#define ADDONS_HPP
#include "SiconosPointers.hpp"
#include "SimpleMatrix.hpp"
#include "DynamicalSystem.hpp"

// => we need swig iterators for bgl iterators

std::vector<SP::DynamicalSystem> dynamicalSystems(SP::DynamicalSystemsGraph dsg)
{
  std::vector<SP::DynamicalSystem> r = std::vector<SP::DynamicalSystem>();
  DynamicalSystemsGraph::VIterator vi, viend;
  for (boost::tie(vi, viend) = dsg->vertices(); vi != viend; ++vi)
  {
    r.push_back(dsg->bundle(*vi));
  };
  return r;
};

std::vector<SP::Interaction> interactions(SP::InteractionsGraph dsg)
{
  std::vector<SP::Interaction> r = std::vector<SP::Interaction>();
  InteractionsGraph::VIterator vi, viend;
  for (boost::tie(vi, viend) = dsg->vertices(); vi != viend; ++vi)
  {
    r.push_back(dsg->bundle(*vi));
  };
  return r;
};

std::vector<SP::DynamicalSystem> dynamicalSystemsVector()
{
  return std::vector<SP::DynamicalSystem>();
}

std::vector<std::pair<SP::DynamicalSystem, SP::DynamicalSystem> >
graphLayout(SP::DynamicalSystemsGraph dsg)
{

  std::vector<std::pair<SP::DynamicalSystem, SP::DynamicalSystem> > r =
    std::vector<std::pair<SP::DynamicalSystem, SP::DynamicalSystem> >();

  DynamicalSystemsGraph::EIterator ei, eiend;

  for (boost::tie(ei, eiend) = dsg->edges(); ei != eiend; ++ei)
  {
    std::pair<SP::DynamicalSystem, SP::DynamicalSystem>
    p(dsg->bundle(dsg->source(*ei)),
      dsg->bundle(dsg->target(*ei)));
    r.push_back(p);
  };
  return r;
};

std::vector<std::pair<unsigned int, unsigned int> >
graphLayoutInt(SP::DynamicalSystemsGraph dsg)
{

  std::vector<std::pair<unsigned int, unsigned int> > r =
    std::vector<std::pair<unsigned int, unsigned int> >();

  DynamicalSystemsGraph::EIterator ei, eiend;

  for (boost::tie(ei, eiend) = dsg->edges(); ei != eiend; ++ei)
  {
    std::pair<unsigned int, unsigned int>
    p(dsg->bundle(dsg->source(*ei))->number(),
      dsg->bundle(dsg->target(*ei))->number());
    r.push_back(p);
  };
  return r;
};

/* access to edges & vertices from python */
struct graphAccess  
{
  std::shared_ptr<_InteractionsGraph> graph;
  std::vector<SP::Interaction> vertices;
  std::vector<SP::DynamicalSystem> edges;
  
  void update()
  {
    vertices.clear();
    edges.clear();
    InteractionsGraph::VIterator ui, uiend;
    for (boost::tie(ui,uiend) = graph->vertices(); ui != uiend; ++ui)
    {
      vertices.push_back(graph->bundle(*ui));
    };
    
    InteractionsGraph::EIterator ei, eiend;
    for (boost::tie(ei,eiend) = graph->edges(); ei != eiend; ++ei)
    {
      edges.push_back(graph->bundle(*ei));
    };
    
    
  };
  
  graphAccess(std::shared_ptr<_InteractionsGraph> ig) : graph(ig)
  {
    update();
  };
  
  ~graphAccess()
  {
    std::vector<SP::Interaction>().swap(vertices);
    std::vector<SP::DynamicalSystem>().swap(edges);
  }
  
};

#endif
