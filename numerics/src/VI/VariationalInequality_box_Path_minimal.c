/* Test file */
#include "SolverOptions.h"
#include "VariationalInequality.h"
#include "VI_cst.h"
#include "solver_registry.h"
#include "numerics_errors.h"
#include "utils/numerics_errors.h"

void vi_box_path(VariationalInequality* problem, double *z, double *F, int *info, SolverOptions* options);

static int vi_box_path_solve_wrap(void* problem, double* z, double* F, SolverOptions* options) {
    int info = NUMERICS_OK;
    vi_box_path((VariationalInequality*)problem, z, F, &info, options);
    return info;
}

static void vi_box_path_set_default(SolverOptions* options) {
    (void)options;
    /* No special defaults needed - uses standard SolverOptions */
}

REGISTER_SOLVER(SICONOS_VI_BOX_PATH, "VI_BOX_PATH",
                             "Box VI solver based on PATH solver",
                             NULL,   /* init_wrap */
                             vi_box_path_solve_wrap,
                             NULL,   /* free_wrap */
                             NULL,   /* err_fn */
                             vi_box_path_set_default,
                             1000,   /* default_max_iter */
                             1e-4,   /* default_tol */
                             0)      /* is_local_solver */
