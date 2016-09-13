%include "SolverOptions.h"

%extend SolverOptions
{
  SolverOptions()
  {
    SolverOptions *SO;
    SO = (SolverOptions *) malloc(sizeof(SolverOptions));
    return SO;
  }

  ~SolverOptions()
  {
    solver_options_delete($self);
    free($self);
  }

};


