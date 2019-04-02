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

%include "projectionOnCone.h"
%include "projectionOnRollingCone.h"

%extend SolverOptions
{
  SolverOptions(enum FRICTION_SOLVER id)
  {
    SolverOptions *SO;
    SO = (SolverOptions *) malloc(sizeof(SolverOptions));

    /* cf Friction_cst.h */
    if(id >= 400 && id < 500)
    {
      fc2d_setDefaultSolverOptions(SO, id);
    }
    else if (id >= 500 && id < 600)
    {
      fc3d_setDefaultSolverOptions(SO, id);
    }
    else if (id >= 600 && id < 700)
    {
      gfc3d_setDefaultSolverOptions(SO, id);
    }
    else
    {
      SWIG_Error(SWIG_RuntimeError, "Unknown friction contact problem solver");
      free(SO);
      return NULL;
    }


    return SO;
  }
  
  SolverOptions(SecondOrderConeLinearComplementarityProblem* soclcp, enum SOCLCP_SOLVER id)
  {
    SolverOptions *SO;
    SO = (SolverOptions *) malloc(sizeof(SolverOptions));

    if (id >= 1100 && id < 1200)
    {
      soclcp_setDefaultSolverOptions(SO, id);
    }
    else
    {
      SWIG_Error(SWIG_RuntimeError, "Unknown SOCLCP solver");
      free(SO);
      return NULL;
    }


    return SO;
  }

  SolverOptions()
  {
    SolverOptions *SO;
    SO = (SolverOptions *) malloc(sizeof(SolverOptions));
    return SO;
  }

  SolverOptions(LinearComplementarityProblem* lcp, enum LCP_SOLVER id)
  {
    SolverOptions *SO;
    SO = (SolverOptions *) malloc(sizeof(SolverOptions));
    solver_options_set(SO, id);
    return SO;
  }

  SolverOptions(MixedLinearComplementarityProblem* mlcp, enum MLCP_SOLVER id)
  {
    SolverOptions *SO;
    SO = (SolverOptions *) malloc(sizeof(SolverOptions));
    SO->solverId=id;
    mixedLinearComplementarity_setDefaultSolverOptions(mlcp, SO);
    return SO;
  }

  SolverOptions(MixedComplementarityProblem* mlcp, enum MCP_SOLVER id)
  {
    SolverOptions *SO;
    SO = (SolverOptions *) malloc(sizeof(SolverOptions));
    SO->solverId=id;
    mixedComplementarity_setDefaultSolverOptions(mlcp, SO);
    return SO;
  }

  SolverOptions(MixedComplementarityProblem2* mcp, enum MCP_SOLVER id)
  {
    SolverOptions *SO;
    SO = (SolverOptions *) malloc(sizeof(SolverOptions));
    solver_options_set(SO, id);
    return SO;
  }

  SolverOptions(NonlinearComplementarityProblem* ncp, enum NCP_SOLVER id)
  {
    SolverOptions *SO;
    SO = (SolverOptions *) malloc(sizeof(SolverOptions));
    solver_options_set(SO, id);
    return SO;
  }

  SolverOptions(VariationalInequality* vi, enum VI_SOLVER id)
  {
    SolverOptions *SO;
    SO = (SolverOptions *) malloc(sizeof(SolverOptions));
    solver_options_set(SO, id);
    return SO;
  }

  SolverOptions(AffineVariationalInequalities* avi, enum AVI_SOLVER id)
  {
    SolverOptions *SO;
    SO = (SolverOptions *) malloc(sizeof(SolverOptions));
    SO->solverId=id;
    solver_options_set(SO, id);
    return SO;
  }

  ~SolverOptions()
  {
    solver_options_delete($self);
    free($self);
  }

};
