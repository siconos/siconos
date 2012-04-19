/* Siconos-Numerics, Copyright INRIA 2005-2011.
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
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */

/*!\file lcp_pred.h
  \brief Header for lcp solvers with extract-predict mechanism (Unstable, under devel.)

*/
#ifndef LCP_pred_H
#define LCP_pred_H

#include "LCP_Solvers.h"

/** structure used to define the lcp problem - Deprecated.
 */
typedef struct
{
} method_lcp;


#ifdef __cplusplus
extern "C"
{
#endif

  /**
     lcp_solver_pred is a generic interface allowing the call of one of the LCP solvers.
     - At first, the signs of q elements are checked to detect the trivial case of positive q.\n
     - Then it tries to find z by assuming that the indices of the non-zero elements
     are the same as the previous solution (previous solutions with trivial case of q positive excepted of course).\n
     - If q is not positive and prediction failed, the regular LCP solver is called.\n

     \param[in] vec          On enter, a (\f$n \times n\f$)-vector of doubles which contains the components of the LCP matrix with a Fortran storage.
     \param[in] q            On enter, a n-vector of doubles which contains the components of the constant right hand side vector.
     \param[in] n            On enter, an integer which represents the dimension of the LCP problem.
     \param[in] ptvec          On enter, a union containing the LCP structure.
     \n \n
     \param[in,out] z        On enter, an initial guess for iterative LCP solvers.\n
     On return, a n-vector of doubles which contains the solution of the problem.
     \param[out] w           On return, a n-vector of doubles which contains the complementary solution of the problem.
     \n
     \param[in]     firsttime      At 1, forces the regular LCP solver to be used (for initialization purpose).
     \param[out]    soltype        On return, indicates how the solution was found (0 : no sol,1 : q positive,2 : prediction,3 : LCP solver)
     \param[in,out] indic          The set of indices of non-zero z values in ascending order
     \param[in,out] indicop        The complementary set of indices of "indic".
     \param[in,out] submatlcp      The submatrix of M defined by "indic".
     \param[in,out] submatlcpop    The submatrix of M defined by "indicop".
     \param[in,out] ipiv           Pivot indices in LU factorization of "submatlcp".
     \param[in,out] sizesublcp     "submatlcp" size.
     \param[in,out] sizesublcpop   "submatlcpop" size.
     \param subq
     \param bufz
     \param newz
     \param workspace
     \return integer
     - 0 : successful\n
     - >0 : otherwise (see specific solvers for more information about the log info)

     \author Nineb Sheherazade & Mathieu Renouf & Pascal Denoyelle
  */
  int lcp_solver_pred(double *vec, double *q , int *n , method_lcp * ptvec ,
                      double *z , double *w ,
                      int firsttime, int *soltype , int *indic , int *indicop ,
                      double *submatlcp ,
                      double *submatlcpop,
                      int *ipiv, int *sizesublcp, int *sizesublcpop,
                      double *subq, double *bufz, double *newz, double *workspace);
  /**
   *  lcp_solver_block_pred_vec solves a LCP with a matrix stored block by block.\n
   *  - It iterates by solving successively diagonal sub-LCPs (Gauss-Seidel block by block).\n
   *  - "pred" means that for each diagonal sub-LCP (Mi,qi), the position of non-zero terms zi_k~ of its previous solution \n
   *  in block iteration are stored as well as the corresponding submatrix Mi~ of this diagonal block. \n
   *  Thus Mi~ * zi_k~ = -qi_k~. The solver tries first to solve the new sub-LCP (Mi,qi_k+1) \n
   *  by solving the linear system Mi~ * zi_k+1~ = -qi_k+1~. If this does not work, the regular LCP algorithm is called.\n
   *  - "vec" means that a set of LCP solver methods is provided, and applied individually to the corresponding diagonal sub-LCP.\n
   *  If there are less solvers than diagonal blocks, the last solver is applied to all last diagonal blocks.\n
   *
   * \param blmat    On enter, the sparse block matrix M of the LCP
   * \param blmatpred   On enter, the position of non-zero z values block by block and all corresponding submatrices
   * \param nbmethod On enter, the number of available solving methods provided
   * \param q        On enter, the vector of doubles of the LCP
   * \param ptvec       On enter, a pointer on the list of solving methods
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
   * \param maxiterglob
   * \param tolglob
   * lcp_solver_block_pred_vec is a generic interface which considers LCP with a sparse block structure. The global LCP is solved
   * as a succession of local LCP solved via lcp_solver_pred.\n
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
   * Data file example:\n
   * If we consider the matrix M and the right-hand-side q defined as
   *
   * \f$
   * M=\left[\begin{array}{cccc|cc|cc}
   *          1 & 2 & 0 & 4   & 3 &-1   & 0 & 0\\
   *          2 & 1 & 0 & 0   & 4 & 1   & 0 & 0\\
   *          0 & 0 & 1 &-1   & 0 & 0   & 0 & 0\\
   *          5 & 0 &-1 & 6   & 0 & 6   & 0 & 0\\
   *          \hline
   *          0 & 0 & 0 & 0   & 1 & 0   & 0 & 5\\
   *          0 & 0 & 0 & 0   & 0 & 2   & 0 & 2\\
   *          \hline
   *          0 & 0 & 2 & 1   & 0 & 0   & 2 & 2\\
   *          0 & 0 & 2 & 2   & 0 & 0   & -1 & 2\\
   *        \end{array}\right] \quad, q=\left[\begin{array}{c}-1\\-1\\0\\-1\\\hline 1\\0\\\hline -1\\2\end{array}\right].
   * \f$
   *
   * then
   * - the number of non null blocks is 6 (blmat.nbblocks) and the number of diagonal blocks is 3 (blmat.size)
   * - the vector blmat.blocksize of diagonal blocks sizes is equal to [4,2,2]
   * - the vectors blmat.ColumnIndex and blmat.RowIndex of non null blocks indices are equal to [0,1,1,2,0,2] and [0,0,1,1,2,2]
   * - the blmat.block contains all non null block matrices stored in Fortran order (column by column) as\n
   *   blmat.block[0] = {1,2,0,5,2,1,0,0,0,0,1,-1,4,0,-1,6}\n
   *   blmat.block[1] = {3,4,0,0,-1,1,0,6}\n
   *   ...\n
   *   blmat.block[5] = {2,-1,2,2}
   * \author Pascal Denoyelle
   */
  int lcp_solver_block_pred_vec(SparseBlockStructuredMatrix *blmat,
                                SparseBlockStructuredMatrixPred *blmatpred,
                                int nbmethod,
                                int maxiterglob, double tolglob,
                                double *q, method_lcp ** ptvec,
                                double *z, double *w,
                                int *it_end, int *itt_end , double *res);

#ifdef __cplusplus
}
#endif

#endif
