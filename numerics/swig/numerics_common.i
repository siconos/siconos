%include "SolverOptions.h"
%include "relay_cst.h"
%include "AVI_cst.h"
 //%include "SOCLCP_cst.h"
%include "Friction_cst.h"
%include "lcp_cst.h"
%include "MCP_cst.h"
%include "NCP_cst.h"
%include "mlcp_cst.h"
%include "VI_cst.h"
%include "GenericMechanical_cst.h"
%include "Newton_methods.h"
%include "projectionOnCone.h"
%include "projectionOnRollingCone.h"

%extend SolverOptions
{
  SolverOptions(int id)
  {
    SolverOptions *SO;
    SO = solver_options_create(id);
    return SO;
  }
  
  SolverOptions()
  {
    SolverOptions *SO;
    SO = (SolverOptions *) malloc(sizeof(SolverOptions));
    return SO;
  }

  ~SolverOptions()
  {
    solver_options_clear(&$self);
  }

};
