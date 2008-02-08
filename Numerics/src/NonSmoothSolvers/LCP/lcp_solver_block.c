/* Siconos-Numerics version 2.1.1, Copyright INRIA 2005-2007.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifndef MEXFLAG
#include "NonSmoothDrivers.h"
#endif
#include "LA.h"

/** lcp_solver_block is a generic interface allowing the call of one of the block LCP solvers.\n
 *
 * \param blmat    On enter, the sparse block matrix M of the LCP
 * \param q        On enter, the vector of doubles of the LCP
 * \param pt       On enter, a union which contains a LCP structure
 *
 * \param z        On enter/return, a vector of doubles which contains the initial iterate for the LCP(q,M) and returns the solution of the problem.
 * \param w        On return, a vector of doubles which returns the solution of the problem.
 * \param it_end   On return, an integer which returns the final number of block iterations
 * \param itt_end  On return, an integer which returns the total number of local LCP iterations or pivots
 * \param res      On return, a double which returns the final value of error criteria.
 *
 * \return info    Integer identifiant for the solver result\n
 *                 0 : convergence\n
 *                 >0 : no convergence (see solver for specific info value)\n
 *
 * lcp_solver_block is a generic interface which consider LCP with a sparse block structure. The global LCP is solved
 * as a succession of local LCP solved via lcp_solver.\n
 *
 * list Keywords to call solvers:
 *
 *   - Lemke    for lcp_lexicolemke
 *   - PGS      for lcp_pgs
 *   - RPGS     for lcp_rpgs
 *   - CPG      for lcp_cpg
 *   - QP       for lcp_qp
 *   - NSQP     for lcp_nsqp
 *   - Latin    for lcp_latin
 *   - Newton   for lcp_newton_min
 *
 * \author Mathieu Renouf & Pascal Denoyelle
 */
int lcp_GS_SBM(LinearComplementarity_Problem* problem, double *z, double *w, Solver_Options* options, Numerics_Options* global_options, int *it_end, int *itt_end , double *res)
{
  /*   /\* Notes:  */
  /*      - we suppose that the trivial solution case has been checked before, */
  /*      since this function is supposed to be called from lcp_driver(). */
  /*      - Input matrix M of the problem is supposed to be sparse-block with no null row (ie no rows with all blocks equal to null) */
  /*   *\/ */

  /*   if(problem->M->storageType!=1) */
  /*     { */
  /*       fprintf(stderr,"lcp_GS_SBS error: wrong storage type for input matrix M of the LCP.\n"); */
  /*       exit(EXIT_FAILURE); */
  /*     } */
  /*   /\* Number of non-null blocks in blmat *\/ */
  /*   int nbOfNonNullBlocks = blmat->nbblocks; */
  /*   if (nbOfNonNullBlocks < 1)  */
  /*     { */
  /*       fprintf(stderr,"Numerics lcp_solver_SBM error: empty M matrix (all blocks = NULL).\n"); */
  /*       exit(EXIT_FAILURE); */
  /*     } */

  /*   /\* Solver parameters*\/ */
  /*   int itermax = options->iparam[0]; */
  /*   double tol = options->dparam[0]; */

  /*   /\* Matrix M/vector q of the LCP *\/ */
  /*   SparseBlockStructuredMatrix* blmat = problem->M->matrix1;  */
  /*   double * q = problem->q; */

  /*   /\* Local problem *\/ */
  /*   LinearComplementarity_Problem * local_problem = malloc(sizeof(*local_problem)); */
  /*   local_problem->storageType = 0; // dense storage */
  /*   local_problem->M->matrix0 = NULL; */
  /*   local_problem->q = NULL; */


  /*   /\* Loop over the rows of blocks in blmat *\/ */
  /*   for( i = 0, i< blmat->size; ++i) */
  /*     { */
  /*       /\* row/col positions of the current block *\/ */
  /*       rowcurbl = blmat->RowIndex[i]; */
  /*       colcurbl = blmat->ColumnIndex[i]; */

  /*       /\* Gets diagonal block *\/ */
  /*       if (rowcurbl == colcurbl)  */
  /*  adrbldiag = blmat->block[i]; */

  /*       /\* extra diagonal blocks *\/ */
  /*       else  */
  /*  { */

  /*  } */


  /*     } */


























  /*   double *adrbldiag; */

  /*   int rowprecbl,rowcurbl,indicrow,rowsize; */
  /*   int colprecbl,colcurbl,indiccol,colsize; */
  /*   int i,j,k; */

  /*   /\* Initialize output *\/ */
  /*   *it_end   = 0; */
  /*   *itt_end  = 0; */
  /*   *res      = 0.0; */
  /*   int totaliter = 0; */

  /*   /\* Number of blocks in a row *\/ */
  /*   int nbBlocksInARow = blmat->size; */

  /*   /\* Dim. of the Matrix M (number of scalar elements (not blocks!) in a row *\/ */
  /*   int n = blmat->blocksize[nbBlocksInARow-1]; */

  /*   /\* Size of the biggest block - Used to allocate memory for rhs *\/ */
  /*   int blsizemax = 0;  */
  /*   for( i = 0 ; i < nbBlocksInARow ; i++ )  */
  /*     { */
  /*       k = blmat->blocksize[i]; */
  /*       if (k > blsizemax) blsizemax = k; */
  /*     } */


  /*   /\* Unused ... *\/ */
  /*   /\*   double qs = DNRM2( n , q , incx ); *\/ */
  /*   /\*   double den = 1.0/qs; *\/ */
  /*   /\* Initialization of z and w *\/ */
  /*   //   for( i = 0 ; i < n ; i++ ) z[i] = 0.; */
  /*   //   dcopy_( (integer *)&n , q , (integer *)&incx , w , (integer *)&incy ); */


  /*   /\* Memory allocation for RHS *\/ */
  /*   double *rhs = ( double* )malloc( blsizemax * sizeof( double ) ); */
  /*   { */
  /*     fprintf(stderr,"lcp_GS_SBS error: allocation failed for rhs.\n"); */
  /*     exit(EXIT_FAILURE); */
  /*   } */

  /*   /\* Termination value *\/ */
  /*   int info= 1; */
  /*   int info1 = 1; */
  /*   /\* Current iteration number *\/ */
  /*   int iter = 0; */
  /*   /\* resulting error *\/ */
  /*   double err  = 1.; */

  /*   int incx = 1,incy = 1; */

  /*   /\* ========== Iterations loop ==========*\/ */
  /*   while( ( iter < itermax ) && ( info ) ) */
  /*     { */
  /*       iter++; */

  /*       /\* dcopy_( (integer *)&n , w , (integer *)&incx , ww , (integer *)&incy );*\/ */

  /*       rowprecbl = -1; */
  /*       rowsize = 0; */
  /*       indicrow = 0; */

  /*       /\* ==== Iterations over non null blocks ==== *\/ */
  /*       for( i = 0 ; i < nbOfNonNullBlocks ; i++ ) */
  /*  { */
  /*    /\* row/col positions of the current block *\/ */
  /*    rowcurbl = blmat->RowIndex[i]; */
  /*    colcurbl = blmat->ColumnIndex[i]; */

  /*    /\* == First case: new row of blocks == *\/ */
  /*    if  (rowcurbl != rowprecbl)  */
  /*      { */
  /*        /\* If not the first row of blocks ... *\/ */
  /*        if (rowprecbl != -1) */
  /*    { */
  /*      if (adrbldiag != NULL) */
  /*        { */
  /*          /\* Local LCP resolution  *\/ */

  /*          options->iparam[1] = 0; */
  /*          local_problem->M->matrix0 = adrbldiag; */
  /*          local_problem->M->size0 = rowsize; */
  /*          local_problem->M->size1 = rowsize; */
  /*          local_problem->q = rhs; */

  /*          info1 = lcp_driver(local_problem, , &z[indicrow] , &w[indicrow], options, global_options); */
  /*          totaliter += options->iparam[1]; */
  /*          /\* If solver failed ...*\/ */
  /*          if (info1 > 0)  */
  /*      { */
  /*        free(rhs); */
  /*        fprintf(stderr,"lcp_GS_SBS error: local LCP solver failed at global iteration %d.\n for block(%d,%d).n", iter,rowprecbl,rowprecbl); */
  /*        exit(EXIT_FAILURE); */
  /*      } */
  /*        } */
  /*      else // if the diagonal block is null  */
  /*        { */
  /*          free(rhs); */
  /*          fprintf(stderr,"lcp_GS_SBS error: null diagonal block for row number %d.\n", rowcurbl); */
  /*          exit(EXIT_FAILURE); */
  /*        } */
  /*    } */

  /*        adrbldiag = NULL; */
  /*        indiccol = blmat->blocksize[colcurbl-1]; */
  /*        colprecbl = colcurbl; */

  /*        indicrow += rowsize; */
  /*        indicrow += blmat->blocksize[rowcurbl-1]; */
  /*        indicrow -= blmat->blocksize[rowprecbl]; */

  /*        rowprecbl = rowcurbl; */
  /*        rowsize = blmat->blocksize[rowcurbl] -  blmat->blocksize[rowcurbl-1]; */
  /*        DCOPY(rowsize, &q[indicrow], incx, rhs, incy); */
  /*      } */

  /*    /\* Gets diagonal block *\/ */
  /*    if (rowcurbl == colcurbl)  */
  /*      adrbldiag = blmat->block[i]; */

  /*    /\* extra diagonal blocks *\/ */
  /*    else  */
  /*      {    */
  /*        indiccol += blmat->blocksize[colcurbl-1]; */
  /*        indiccol -=  blmat->blocksize[colprecbl]; */
  /*        colprecbl = colcurbl; */

  /*        /\* Dim of the current block *\/ */
  /*        rowsize = blmat->blocksize[rowcurbl] - blmat->blocksize[rowcurbl-1]; */
  /*        colsize = blmat->blocksize[colcurbl] - blmat->blocksize[colcurbl-1]; */

  /*        DGEMV( LA_NOTRANS , rowsize , colsize , 1.0 , blmat->block[i], rowsize , &z[indiccol], incx , 1.0 , rhs , incy ); */
  /*      } */
  /*  } */
  /*       /\* ==== End of iterations over non null blocks ==== *\/ */

  /*       if (adrbldiag != NULL) */
  /*  { */
  /*    /\* Local LCP resolution  *\/ */
  /*    options->iparam[1] = 0; */
  /*    /\* Local_problem */
  /*       M = adrbldiag */
  /*       q = rhs  */
  /*       size = rowsize */
  /*       options = pt */
  /*    *\/ */
  /*    info1 = lcp_driver(local_problem, , &z[indicrow] , &w[indicrow], options, global_options); */

  /*    totaliter += options->iparam[1]; */

  /*    if (info1 > 0)  */
  /*      { */
  /*        if (verbose > 0)  */
  /*    { */
  /*      printf(" Sub LCP solver failed at global iteration %d in lcp_solver_block",iter); */
  /*      printf(" for block(%d,%d)\n",rowprecbl,rowprecbl); */
  /*    } */
  /*        free(rhs); */
  /*        return 1; */
  /*      } */
  /*  } */
  /*       else  */
  /*  { */
  /*    printf("NULL diagonal block in LCP !!!\n"); */
  /*    free(rhs); */
  /*    return 1; */
  /*  } */

  /*       /\* **** Criterium convergence **** *\/ */

  /*       info = filter_result_LCP_block(blmat,q,z,tol,verbose,w); */
  /*       /\*    qs   = -1; */

  /*       daxpy_( (integer *)&n , &qs , w , (integer *)&incx , ww , (integer *)&incy ); */

  /*       num = dnrm2_( (integer *)&n, ww , (integer *)&incx ); */
  /*       err = num*den;*\/ */
  /*       /\* **** ********************* **** *\/ */
  /*     } */
  /*   /\* ========== End of Iterations loop ==========*\/ */

  /*   *it_end  = iter; */
  /*   *itt_end = totaliter; */
  /*   *res     = err; */

  /*   /\*  free(ww);*\/ */
  /*   free(rhs); */

  /*   if(verbose > 0 ) */
  /*     { */
  /*       /\*    if( err > tol ){*\/ */
  /*       /\*      printf( " The residue is : %g \n", err );*\/ */
  /*       /\*      return 1;*\/ */
  /*       if (info) */
  /*  printf( " No convergence of lcp_solver_block after %d iterations\n" , iter ); */

  /*       else */
  /*  printf(" Convergence of lcp_solver_block after %d iterations\n" , iter ); */
  /*       /\*      printf( " The residue is : %g \n", err );*\/ */
  /*     } */
  /*   /\*  else{ */
  /*       if( err > tol ) return 1; */
  /*       }*\/ */

  /*   /\*  info = filter_result_LCP_block(blmat,q,z,tol,pt->lcp.chat,w);*\/ */

  /*   return info; */
  return 0;
}
