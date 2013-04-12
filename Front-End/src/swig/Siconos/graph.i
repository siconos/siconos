// -*- c++ -*-

%{
#include <SiconosGraph.hpp>
%}

%extend InteractionsGraph
{
  const std::vector<SP::Interaction> vertices()
  {
    std::vector<SP::Interaction> r;
    InteractionsGraph::VIterator ui, uiend;
    for (boost::tie(ui,uiend) = $self->vertices(); ui != uiend; ++ui)
    {
      r.push_back($self->bundle(*ui));
    };
    return r;
  };


  const std::vector<SP::DynamicalSystem> edges()
  {
    std::vector<SP::DynamicalSystem> r;
    InteractionsGraph::EIterator ui, uiend;
    for (boost::tie(ui,uiend) = $self->edges(); ui != uiend; ++ui)
    {
      r.push_back($self->bundle(*ui));
    };
    return r;
  };

}


%extend DynamicalSystemsGraph
{
  const std::vector<SP::DynamicalSystem> vertices()
  {
    std::vector<SP::DynamicalSystem> r;
    DynamicalSystemsGraph::VIterator ui, uiend;
    for (boost::tie(ui,uiend) = $self->vertices(); ui != uiend; ++ui)
    {
      r.push_back($self->bundle(*ui));
    };
    return r;
  };


  const std::vector<SP::Interaction> edges()
  {
    std::vector<SP::Interaction> r;
    DynamicalSystemsGraph::EIterator ui, uiend;
    for (boost::tie(ui,uiend) = $self->edges(); ui != uiend; ++ui)
    {
      r.push_back($self->bundle(*ui));
    };
    return r;
  };

}



 /* swig has difficulties with this macro in SiconosProperties */
#undef INSTALL_GRAPH_PROPERTIES
#define INSTALL_GRAPH_PROPERTIES(X,Y)

%include "SiconosGraph.hpp"

TYPEDEF_SPTR(_DynamicalSystemsGraph);
%feature("director") _DynamicalSystemsGraph;
%shared_ptr( SiconosGraph<std11::shared_ptr<DynamicalSystem>, 
                          std11::shared_ptr<Interaction>, 
                          SystemProperties, InteractionProperties, 
                          GraphProperties >);

TYPEDEF_SPTR(_InteractionsGraph);
%feature("director") _InteractionsGraph;
%shared_ptr( SiconosGraph<std11::shared_ptr<Interaction>, 
                          std11::shared_ptr<DynamicalSystem>, 
                          InteractionProperties, SystemProperties, 
                          GraphProperties >);

TYPEDEF_SPTR(DynamicalSystemsGraph);
%feature("director") DynamicalSystemsGraph;
%shared_ptr(DynamicalSystemsGraph);

TYPEDEF_SPTR(InteractionsGraph);
%feature("director") InteractionsGraph;
%shared_ptr(InteractionsGraph);

// must be specified after %shared_ptr, if ever needed
%template(_DynamicalSystemsGraph) SiconosGraph<
  std11::shared_ptr<DynamicalSystem>, 
  std11::shared_ptr<Interaction>, 
  SystemProperties, InteractionProperties, 
  GraphProperties >;

%template(SP_DynamicalSystemsGraph) std11::shared_ptr<
  SiconosGraph<std11::shared_ptr<DynamicalSystem>, 
               std11::shared_ptr<Interaction>, 
               SystemProperties, InteractionProperties, 
               GraphProperties > >;

%template(_InteractionsGraph) SiconosGraph<
  std11::shared_ptr<Interaction>, 
  std11::shared_ptr<DynamicalSystem>, 
  InteractionProperties, SystemProperties, 
  GraphProperties >;

%template(SP_InteractionsGraph) std11::shared_ptr<
  SiconosGraph<std11::shared_ptr<Interaction>, 
               std11::shared_ptr<DynamicalSystem>, 
               InteractionProperties, SystemProperties, 
               GraphProperties > >;

%ignore DynamicalSystemsGraph::vertices;
%ignore DynamicalSystemsGraph::edges;
%ignore InteractionsGraph::vertices;
%ignore InteractionsGraph::edges;

%template (ig_vertices) std::vector<std11::shared_ptr<Interaction> >;
%template (ig_edges) std::vector<std11::shared_ptr<DynamicalSystem> >;
