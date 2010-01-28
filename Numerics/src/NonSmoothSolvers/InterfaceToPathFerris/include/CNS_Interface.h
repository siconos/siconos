/*****************************************************************************/
/* CNS_Interface.h                                                           */
/*                                                                           */
/* DESCRIPTION                                                               */
/*   The constrained nonlinear system is to find an x_star such that:        */
/*     F(x_star) = 0                                                         */
/*     lb <= x_start <= ub                                                   */
/*   where F is a given function from R^n to R^n, and lb and ub are          */
/*   prescribed lower and upper bounds.                                      */
/*                                                                           */
/*   The function, jacobian, lower, and upper bounds are communicated        */
/*   through a CNS_Interface, which must be set in an allocated CNS          */
/*   structure using the CNS_SetInterface( ) function before calling any     */
/*   algorithms or conversion routines.                                      */
/*****************************************************************************/

#ifndef CNS_INTERFACE_H
#define CNS_INTERFACE_H

#include "Types.h"
#include "Presolve_Interface.h"

struct _CNS;
typedef struct _CNS CNS;

/*****************************************************************************/
/* Allocation and deallocation functions.                                    */
/*****************************************************************************/
/*                                                                           */
/* - CNS_Create -  allocate a CNS structure with the given estimate of the   */
/*                 problem size and number of nonzeros in the Jacobian.      */
/*                 If the estimated size passed into the creation routine    */
/*                 as smaller than the actual problem size, the code will    */
/*            reallocate memory.                                        */
/*                                                                           */
/* - CNS_Destroy - deallocate the given CNS structure.  Further use of the   */
/*                 structure is prohibited.                                  */
/*                                                                           */
/* - CNS_Size    - reallocate the CNS structure with the new problem size    */
/*                 and number of nonzeros given.  The reallocation only      */
/*                 occurs if the new size is larger than the current size of */
/*                 the given CNS structure.                                  */
/*                                                                           */
/*****************************************************************************/

FUN_DECL(CNS *)CNS_Create(Int maxModSize, Int maxModNNZ);
FUN_DECL(Void) CNS_Destroy(CNS *c);
FUN_DECL(Void) CNS_Size(CNS *c, Int maxModSize, Int maxModNNZ);

/*****************************************************************************/
/* Access routines - used to obtain information about the problem.  Do NOT   */
/*                   directly modify any of the data provided!  This can     */
/*                   cause really bad problems.                              */
/*****************************************************************************/
/*                                                                           */
/* - CNS_GetX - obtain the current x vector                                  */
/* - CNS_GetL - obtain the current lower bound vector                        */
/* - CNS_GetU - obtain the current upper bound vector                        */
/* - CNS_GetB - obtain the current basis indication vector                   */
/*                                                                           */
/* - CNS_GetF - obtain the current function evaluation vector                */
/* - CNS_GetJ - obtain the current Jacobian evaluation matrix                */
/*                                                                           */
/*****************************************************************************/

FUN_DECL(Double *) CNS_GetX(CNS *c);
FUN_DECL(Double *) CNS_GetL(CNS *c);
FUN_DECL(Double *) CNS_GetU(CNS *c);
FUN_DECL(Int *)    CNS_GetB(CNS *c);

FUN_DECL(Double *) CNS_GetF(CNS *c);
FUN_DECL(Void)     CNS_GetJ(CNS *c, Int *nz, Int **col, Int **len,
                            Int **row, Double **data);

/*****************************************************************************/
/* CNS_Interface declaration.                                                */
/*****************************************************************************/
/*                                                                           */
/* - interface_data is a user defined piece of data passed as the first      */
/*   argument to all of the interface functions.                             */
/*                                                                           */
/* - problem_size should fill in the size of the problem.  The number of     */
/*   nonzeros estimate must be an upper bound on the actual number of        */
/*   nonzeros that will be encountered.                                      */
/*                                                                           */
/* - bounds should fill in an initial point, and the lower and upper         */
/*   bounds on the problem.  The input lower and upper bounds are            */
/*   initialized to minus and plus infinity, respectively.  Only modify      */
/*   them to indicate the finite bounds on your problem.                     */
/*                                                                           */
/* - function_evaluations should perform a function evaluation at x,         */
/*   placing the result in f and returning the number of domain violations   */
/*   (i.e. divisions by zero)                                                */
/*                                                                           */
/* - jacobian_evaluation should:                                             */
/*       1)  Perform a function evaluation if wantf is true.                 */
/*       2)  Fill in the jacobian of F at x in compressed sparse column      */
/*           format, filling in nnz.  col[i] is the column start in          */
/*           [1, ..., nnz] in the row/data vector, len[i] is the number of   */
/*           nonzeros in the column, row[j] is the row index in [1, ..., n]  */
/*           data[j] is the value of the element.  NOTE: this uses FORTRAN   */
/*           style indices in col and row!  The function returns the number  */
/*           of domain violations encountered.                               */
/*                                                                           */
/* - start and finish are not required, but will be called at the start and  */
/*   end of the path call respectively.                                      */
/*                                                                           */
/* - variable_name and constraint_name are not required but can be used to   */
/*   tell the names of the variable/constraint.  The index given is between  */
/*   1 and n.  I.e. it is a FORTRAN style index.                             */
/*                                                                           */
/* - basis is not required, but is a function used to communicate an initial */
/*   basis to the code.                                                      */
/*                                                                           */
/*****************************************************************************/

typedef struct
{
  Void *interface_data;

  Void(CB_FPTR problem_size)(Void *id, Int *size, Int *nnz);
  Void(CB_FPTR bounds)(Void *id, Int size,
                       Double *x, Double *l, Double *u);

  Int(CB_FPTR function_evaluation)(Void *id, Int n, Double *x, Double *f);
  Int(CB_FPTR jacobian_evaluation)(Void *id, Int n, Double *x,
                                   Int wantf, Double *f,
                                   Int *nnz, Int *col, Int *len,
                                   Int *row, Double *data);

  /***************************************************************************/
  /* The following functions are not required.  If they are not provided,    */
  /* simply fill in NULL for the value.  The code reacts appropriately in    */
  /* such circumstances.                                                     */
  /***************************************************************************/

  Void(CB_FPTR start)(Void *id);
  Void(CB_FPTR finish)(Void *id, Double *x);
  Void(CB_FPTR variable_name)(Void *id, Int variable,
                              Char *buffer, Int buf_size);
  Void(CB_FPTR constraint_name)(Void *id, Int constr,
                                Char *buffer, Int buf_size);
  Void(CB_FPTR basis)(Void *id, Int size, Int *basX);
} CNS_Interface;

/*****************************************************************************/
/* Interface functions.                                                      */
/*****************************************************************************/
/*                                                                           */
/* CNS_SetInterface - set the CNS_Interface for an allocated CNS.  This      */
/*                    function must be called before using the CNS in an     */
/*                    algorithm or conversion routine.                       */
/*                                                                           */
/* CNS_GetInterface - obtain the CNS_Interface for the provided CNS.         */
/*                                                                           */
/* CNS_SetPresolveInterface - set the Presolve_Interface for an allocated    */
/*                            CNS.  If preprocessing is going to be used,    */
/*                            this function must be called before using the  */
/*                            CNS in an algorithm or conversion routine.     */
/*                                                                           */
/* CNS_GetPresolveInterface - obtain the Presolve_Interface for the provided */
/*                            CNS.                                           */
/*                                                                           */
/*****************************************************************************/

FUN_DECL(Void) CNS_SetInterface(CNS *c, CNS_Interface *i);
FUN_DECL(Void) CNS_GetInterface(CNS *c, CNS_Interface *i);

FUN_DECL(Void) CNS_SetPresolveInterface(CNS *c, Presolve_Interface *i);
FUN_DECL(Void) CNS_GetPresolveInterface(CNS *c, Presolve_Interface *i);

/*****************************************************************************/
/* Jacobian flag functions - used to provide information about the structure */
/*                           of the Jacobian provided by the interface.      */
/*                           Make sure that if you set a flag, your jacobian */
/*                           evaluation routine obeys the convention.        */
/*                           Failure to do so can lead to problems.          */
/*                                                                           */
/* NOTE: in order for the presolve to work, you MUST have a constant         */
/*       Jacobian structure.                                                 */
/*****************************************************************************/
/*                                                                           */
/* CNS_Jacobian_Structure_Constant - if set to True, the structure of the    */
/*                                   Jacobian does not change when a new     */
/*                                   evaluation is performed; only the data  */
/*                                   changes.  A constant Jacobian structure */
/*                                   is a prerequisite for preprocessing.    */
/*                                                                           */
/* CNS_Jacobian_Data_Contiguous    - if set to True, the Jacobian data is    */
/*                                   stored continuously from [1..nnz] in    */
/*                                   data and row arrays.                    */
/*                                                                           */
/* CNS_Jacobian_Diagonal           - if set to True, each column in the      */
/*                                   Jacobian has an element on the          */
/*                                   diagonal.                               */
/*                                                                           */
/* CNS_Jacobian_Diagonal_First     - if set to True, for each column in the  */
/*                                   Jacobian having a diagonal element, the */
/*                                   diagonal element is first in the column */
/*                                   listing.                                */
/*                                                                           */
/*****************************************************************************/

FUN_DECL(Void) CNS_Jacobian_Structure_Constant(CNS *c, Boolean b);
FUN_DECL(Void) CNS_Jacobian_Data_Contiguous(CNS *c, Boolean b);
FUN_DECL(Void) CNS_Jacobian_Diagonal(CNS *c, Boolean b);
FUN_DECL(Void) CNS_Jacobian_Diagonal_First(CNS *c, Boolean b);

#endif

