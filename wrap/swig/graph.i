// -*- c++ -*-

%{
#include <SimulationGraphs.hpp>

  struct no_property{};
%}

//%include <SimulationGraphs.hpp>

%extend Siconos::Properties
{
  Siconos::Properties::reference __getitem__(const Siconos::Properties::key_type& v)
  {
    return (*self[v]);
  };
};

%extend Siconos::SubProperties
{
  Siconos::SubProperties::reference __getitem__(const Siconos::SubProperties::key_type& v)
  {
    return (*self[v]);
  };
};

struct InteractionsGraph{};
%extend InteractionsGraph
{

  // const std::vector<InteractionsGraph::VDescriptor> vertices()
  // {
  //   std::vector<InteractionsGraph::VDescriptor> r;
  //   InteractionsGraph::VIterator ui, uiend;
  //   for (boost::tie(ui,uiend) = $self->vertices(); ui != uiend; ++ui)
  //   {
  //     r.push_back(*ui);
  //   };
  //   return r;
  // };

  // const std::vector<InteractionsGraph::EDescriptor> edges()
  // {
  //   std::vector<InteractionsGraph::EDescriptor> r;
  //   InteractionsGraph::EIterator ui, uiend;
  //   for (boost::tie(ui,uiend) = $self->edges(); ui != uiend; ++ui)
  //   {
  //     r.push_back(*ui);
  //   };
  //   return r;
  // };

  const std::vector<SP::Interaction> interactions()
  {
    std::vector<SP::Interaction> r;
    InteractionsGraph::VIterator ui, uiend;
    for (boost::tie(ui,uiend) = $self->vertices(); ui != uiend; ++ui)
    {
      r.push_back($self->bundle(*ui));
    };
    return r;
  };


  const std::vector<SP::DynamicalSystem> dynamicalSystems()
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

struct DynamicalSystemsGraph{};
%extend DynamicalSystemsGraph
{

  // const std::vector<DynamicalSystemsGraph::VDescriptor> vertices()
  // {
  //   std::vector<DynamicalSystemsGraph::VDescriptor> r;
  //   DynamicalSystemsGraph::VIterator ui, uiend;
  //   for (boost::tie(ui,uiend) = $self->vertices(); ui != uiend; ++ui)
  //   {
  //     r.push_back(*ui);
  //   };
  //   return r;
  // };


  // const std::vector<DynamicalSystemsGraph::EDescriptor> edges()
  // {
  //   std::vector<DynamicalSystemsGraph::EDescriptor> r;
  //   DynamicalSystemsGraph::EIterator ui, uiend;
  //   for (boost::tie(ui,uiend) = $self->edges(); ui != uiend; ++ui)
  //   {
  //     r.push_back(*ui);
  //   };
  //   return r;
  // };

  const std::vector<SP::DynamicalSystem> dynamicalSystems()
  {
    std::vector<SP::DynamicalSystem> r;
    DynamicalSystemsGraph::VIterator ui, uiend;
    for (boost::tie(ui,uiend) = $self->vertices(); ui != uiend; ++ui)
    {
      r.push_back($self->bundle(*ui));
    };
    return r;
  };

  const std::vector<SP::Interaction> interactions()
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
// %import <boost/config.hpp>
// %import <boost/version.hpp>
// %import <boost/graph/graph_utility.hpp>
// %import <boost/graph/adjacency_list.hpp>
// #if (BOOST_VERSION >= 104000)
// %import <boost/property_map/property_map.hpp>
// #else
// %import <boost/property_map.hpp>
// #endif

 //%import <boost/static_assert.hpp>

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

%feature("director") DynamicalSystemsGraph;
%shared_ptr(DynamicalSystemsGraph);

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

%template (ig_interactions) std::vector<std11::shared_ptr<Interaction> >;
%template (ig_dynamicalSystems) std::vector<std11::shared_ptr<DynamicalSystem> >;


/* missing in generated file, why ? */
// %{
//   namespace swig {
//     template<>
//     struct traits<void>
//     {
//       typedef value_category category;
//       static const char* type_name() { return "void"; }
//     };
//   }
// %}


// /* this 2 should be sufficients */
// %template (dsg_edescr) std::vector<
//   SiconosGraph<std11::shared_ptr<DynamicalSystem>,
//                std11::shared_ptr<Interaction>,
//                SystemProperties, InteractionProperties,
//                GraphProperties >::EDescriptor >;

// %template (ig_vdescr) std::vector<
//   SiconosGraph<std11::shared_ptr<Interaction>,
//                std11::shared_ptr<DynamicalSystem>,
//                InteractionProperties, SystemProperties,
//                GraphProperties >::VDescriptor >;



