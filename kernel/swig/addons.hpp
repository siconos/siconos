#ifndef ADDONS_HPP
#define ADDONS_HPP
#include "SiconosPointers.hpp"
#include "SiconosMatrix.hpp"
#include "DynamicalSystem.hpp"

// => we need swig iterators for bgl iterators

std::vector<std::shared_ptr<siconos::modeling::DynamicalSystem>> dynamicalSystems(std::shared_ptr<siconos::graphs::DynamicalSystemGraph> dsg)
{
  std::vector<std::shared_ptr<siconos::modeling::DynamicalSystem>> r = std::vector<std::shared_ptr<siconos::modeling::DynamicalSystem>>();
  DynamicalSystemsGraph::VIterator vi, viend;
  for (boost::tie(vi, viend) = dsg->vertices(); vi != viend; ++vi)
  {
    r.push_back(dsg->bundle(*vi));
  };
  return r;
};

std::vector<std::shared_ptr<siconos::modeling::Interaction>> interactions(std::shared_ptr<siconos::graphs::InteractionsGraph> dsg)
{
  std::vector<std::shared_ptr<siconos::modeling::Interaction>> r = std::vector<std::shared_ptr<siconos::modeling::Interaction>>();
  InteractionsGraph::VIterator vi, viend;
  for (boost::tie(vi, viend) = dsg->vertices(); vi != viend; ++vi)
  {
    r.push_back(dsg->bundle(*vi));
  };
  return r;
};

std::vector<std::shared_ptr<siconos::modeling::DynamicalSystem>> dynamicalSystemsVector()
{
  return std::vector<std::shared_ptr<siconos::modeling::DynamicalSystem>>();
}

std::vector<std::pair<std::shared_ptr<siconos::modeling::DynamicalSystem>, std::shared_ptr<siconos::modeling::DynamicalSystem>> >
graphLayout(std::shared_ptr<siconos::graphs::DynamicalSystemGraph> dsg)
{

  std::vector<std::pair<std::shared_ptr<siconos::modeling::DynamicalSystem>, std::shared_ptr<siconos::modeling::DynamicalSystem>> > r =
    std::vector<std::pair<std::shared_ptr<siconos::modeling::DynamicalSystem>, std::shared_ptr<siconos::modeling::DynamicalSystem>> >();

  DynamicalSystemsGraph::EIterator ei, eiend;

  for (boost::tie(ei, eiend) = dsg->edges(); ei != eiend; ++ei)
  {
    std::pair<std::shared_ptr<siconos::modeling::DynamicalSystem>, std::shared_ptr<siconos::modeling::DynamicalSystem>>
    p(dsg->bundle(dsg->source(*ei)),
      dsg->bundle(dsg->target(*ei)));
    r.push_back(p);
  };
  return r;
};

std::vector<std::pair<unsigned int, unsigned int> >
graphLayoutInt(std::shared_ptr<siconos::graphs::DynamicalSystemGraph> dsg)
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
  std::vector<std::shared_ptr<siconos::modeling::Interaction>> vertices;
  std::vector<std::shared_ptr<siconos::modeling::DynamicalSystem>> edges;
  
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
    std::vector<std::shared_ptr<siconos::modeling::Interaction>>().swap(vertices);
    std::vector<std::shared_ptr<siconos::modeling::DynamicalSystem>>().swap(edges);
  }
  
};

#endif
