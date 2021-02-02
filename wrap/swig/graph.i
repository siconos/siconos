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

// Must be before declaration of struct InteractionsGraph
// otherwise swig doesn't associate shared_ptr type!
%shared_ptr(InteractionsGraph);
%feature("director") InteractionsGraph;

struct InteractionsGraph{};
%extend InteractionsGraph
{

  PyObject* interactions()
  {
    PyObject* py_tuple = PyTuple_New($self->size());
    InteractionsGraph::VIterator ui, uiend;
    PyObject* resultobj;
    size_t i = 0;
    for (boost::tie(ui,uiend) = $self->vertices(); ui != uiend; ++ui, ++i)
    {
      SP::Interaction * nptr = new SP::Interaction($self->bundle(*ui));
      resultobj = SWIG_NewPointerObj(%as_voidptr(nptr), $descriptor(SP::Interaction *), SWIG_POINTER_OWN);
      PyTuple_SetItem(py_tuple, i, resultobj);
    };
    return py_tuple;
  };


  PyObject* dynamicalSystems()
  {
    PyObject* py_tuple = PyTuple_New($self->edges_number());
    InteractionsGraph::EIterator ui, uiend;
    PyObject* resultobj;
    size_t i = 0;
    for (boost::tie(ui,uiend) = $self->edges(); ui != uiend; ++ui, ++i)
    {
      %convert_sp_ds($self->bundle(*ui), resultobj);
      PyTuple_SetItem(py_tuple, i, resultobj);
    };
    return py_tuple;
  };

}

// Must be before declaration of struct DynamicalSystemsGraph
// otherwise swig doesn't associate shared_ptr type!
%feature("director") DynamicalSystemsGraph;
%shared_ptr(DynamicalSystemsGraph);

struct DynamicalSystemsGraph{};
%extend DynamicalSystemsGraph
{

  PyObject* interactions()
  {
    PyObject* py_tuple = PyTuple_New($self->edges_number());
    DynamicalSystemsGraph::EIterator ui, uiend;
    PyObject* resultobj;
    size_t i = 0;
    for (boost::tie(ui,uiend) = $self->edges(); ui != uiend; ++ui, ++i)
    {
      SP::Interaction * nptr = new SP::Interaction($self->bundle(*ui));
      resultobj = SWIG_NewPointerObj(%as_voidptr(nptr), $descriptor(SP::Interaction *), SWIG_POINTER_OWN);
      PyTuple_SetItem(py_tuple, i, resultobj);
    };
    return py_tuple;
  };


  PyObject* dynamicalSystems()
  {
    PyObject* py_tuple = PyTuple_New($self->size());
    InteractionsGraph::VIterator ui, uiend;
    PyObject* resultobj;
    size_t i = 0;
    for (boost::tie(ui,uiend) = $self->vertices(); ui != uiend; ++ui, ++i)
    {
      %convert_sp_ds($self->bundle(*ui), resultobj);
      PyTuple_SetItem(py_tuple, i, resultobj);
    };
    return py_tuple;
  };

}

%typemap(out) (std::vector<SP::DynamicalSystem>)
{
  PyObject* py_tuple = PyTuple_New($1.size());
  if (!py_tuple) SWIG_fail;
  PyObject* tmpobj;

  for (size_t i = 0; i < $1.size(); ++i)
  {
    %convert_sp_ds($1.at(i), tmpobj);
    PyTuple_SetItem(py_tuple, i, tmpobj);
  }

  $result = py_tuple;
}

%typemap(out) (std::vector<SP::Interaction>)
{
  PyObject* py_tuple = PyTuple_New($1.size());
  if (!py_tuple) SWIG_fail;
  PyObject* tmpobj;

  for (size_t i = 0; i < $1.size(); ++i)
  {
    SP::Interaction * nptr = new SP::Interaction($1.at(i));
    tmpobj = SWIG_NewPointerObj(%as_voidptr(nptr), $descriptor(SP::Interaction *), SWIG_POINTER_OWN);
    PyTuple_SetItem(py_tuple, i, tmpobj);
  }

  $result = py_tuple;
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
%shared_ptr( SiconosGraph<std::shared_ptr<DynamicalSystem>,
                          std::shared_ptr<Interaction>,
                          DynamicalSystemProperties, InteractionProperties,
                          GraphProperties >);

TYPEDEF_SPTR(_InteractionsGraph);
%feature("director") _InteractionsGraph;
%shared_ptr( SiconosGraph<std::shared_ptr<Interaction>,
                          std::shared_ptr<DynamicalSystem>,
                          InteractionProperties, DynamicalSystemProperties,
                          GraphProperties >);

// must be specified after %shared_ptr, if ever needed
%template(_DynamicalSystemsGraph) SiconosGraph<
  std::shared_ptr<DynamicalSystem>,
  std::shared_ptr<Interaction>,
  DynamicalSystemProperties, InteractionProperties,
  GraphProperties >;

%template(SP_DynamicalSystemsGraph) std::shared_ptr<
  SiconosGraph<std::shared_ptr<DynamicalSystem>,
               std::shared_ptr<Interaction>,
               DynamicalSystemProperties, InteractionProperties,
               GraphProperties > >;

%template(_InteractionsGraph) SiconosGraph<
  std::shared_ptr<Interaction>,
  std::shared_ptr<DynamicalSystem>,
  InteractionProperties, DynamicalSystemProperties,
  GraphProperties >;

%template(SP_InteractionsGraph) std::shared_ptr<
  SiconosGraph<std::shared_ptr<Interaction>,
               std::shared_ptr<DynamicalSystem>,
               InteractionProperties, DynamicalSystemProperties,
               GraphProperties > >;

%ignore DynamicalSystemsGraph::vertices;
%ignore DynamicalSystemsGraph::edges;
%ignore InteractionsGraph::vertices;
%ignore InteractionsGraph::edges;


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
//   SiconosGraph<std::shared_ptr<DynamicalSystem>,
//                std::shared_ptr<Interaction>,
//                DynamicalSystemProperties, InteractionProperties,
//                GraphProperties >::EDescriptor >;

// %template (ig_vdescr) std::vector<
//   SiconosGraph<std::shared_ptr<Interaction>,
//                std::shared_ptr<DynamicalSystem>,
//                InteractionProperties, DynamicalSystemProperties,
//                GraphProperties >::VDescriptor >;



