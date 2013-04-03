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
