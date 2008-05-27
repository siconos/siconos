#ifndef STANDALONE_H
#define STANDALONE_H

/*! \page PathFerrisInterface Interface to Path Solver

  For details on Path Solver, see:
  "The Path Solver: A Non-Monotone Stabilization Scheme for Mixed Complementarity Problems",\n
  S.P. Dirkse, M.C. Ferris, Sept. 93 \n
  or
  "Algorithms and Environments for Complementarity", Todd S. Munson


  To call Path Solver, you must: \n\n
   - provide the following two functions:\n\n
    - funcEval - evaluate the function at z, placing the result in f.
      Return the number of domain violations. \n\n
    - jacEval  - evaluate the Jacobian at z, placing the result in
      col_start, col_len, row, and data.  Return the number
      of domain violations.\n\n
   - call pathMain() with the right arguments

 Note: all indices in the Jacobian begin with one.  The Jacobian is stored
 by columns.  Finally, if the structure of the Jacobian remains constant,
 col_start, col_len, and row need only be filled once.


*/

/*!\file Standalone_Path.h
  \brief general interface to call Path Solver

 */

/** pointer to function used to call funcEval */
typedef int (*FuncEvalPtr)(int, double*, double*);

/** pointer to function used to call jacEval */
typedef int (*JacEvalPtr)(int, int, double*, int*, int*, int*, double*);

#ifdef __cplusplus
extern "C" {
#endif
  /** Routine used to call Path solver -
      The main path routine takes the following as arguments:
      \param[in] n      - the number of variables in the problem
      \param[in] nnz    - the number of nonzeros in the jacobian
      \param[in,out] status - final status of the problem, one of the following:\n
      1 - solved \n
      2 - stationary point found\n
      3 - major iteration limit\n
      4 - cumulative minor iteration limit\n
      5 - time limit\n
      6 - user interrupt\n
      7 - bound error (lb is not less than ub)\n
      8 - domain error (could not find a starting point\n)
      9 - internal error\n
      \param[in,out] z  IN: a starting point for the problem OUT: final point
      \param[in,out] f - final function evaluation
      \param[in] lb     - lower bounds on the variables
      \param[in] ub     - upper bounds on the variables\n
      On input, z, f, lb, and ub must be vectors with at least n elements.
      For infinity, use 1e20.
  */
  void pathMain(int n, int nnz, int *status,
                double *z, double *f, double *lb, double *ub);

  /** Set the function used to compute f(z)
      \param fPtr the pointed function
   */
  void setFuncEval(FuncEvalPtr);

  /** Set the function used to compute \f$ \nabla_zf(z) \f$
      \param jfPtr the pointed function
  */
  void setJacEval(JacEvalPtr);


#ifdef __cplusplus
}
#endif
#endif
