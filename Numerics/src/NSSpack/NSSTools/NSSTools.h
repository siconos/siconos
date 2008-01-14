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

#ifndef NSSTOOLS_H
#define NSSTOOLS_H

/*!\file NSSTools.h
  Header to collect basic tools, structures definition or any usefull things for NSSpack

*/

/*!\struct SparseBlockStructuredMatrix

\brief To store sparse block matrices with square diagonal blocks

\param nbblocks         : the total number of non null blocks
\param **block          : *block contains the double values of one block in Fortran storage (column by column)
**block is the list of non null blocks
\param size             : the number of blocks along a row (or column)
\param *blocksize       : the list of the sizes of diagonal (square) blocks
\param *RowIndex        : the list of *block row indices (first row = 0)
\param *ColumnIndex     : the list of *block column indices (first column = 0)
*/

typedef struct
{
  int nbblocks;
  double **block;
  int size;
  int *blocksize;
  int *RowIndex;
  int *ColumnIndex;
} SparseBlockStructuredMatrix;

#ifdef __cplusplus
extern "C" {
#endif

  /** Destructor for SparseBlockStructuredMatrix objects
      \param SparseBlockStructuredMatrix, the matrix to be destroyed.
   */
  void freeSpBlMat(SparseBlockStructuredMatrix *blmat);

  /** Search for the max. element of a vector
      \param[in] x, the vector
      \param[in-out] solution, value of the greatest element of x
      \param[in] n, size of x
  */
  void max_part(double*, double*, int);

  /** compare two double a and b, and return the max.
   *  \param a, double*
   *  \param b, double*
   *  \param c, double*, the max
   */
  void maxf(double*, double*, double*);

  /** Search for the min. element of a vector
      \param[in] x, the vector
      \param[in-out] solution, value of the smallest element of x
      \param[in] n, size of x
  */
  void min_part(double*, double*, int);

  /** compare two double a and b, and return the min.
   *  \param a, double*
   *  \param b, double*
   *  \param c, double*, the min
   */
  void minf(double*, double*, double*);

  /** Positive part values of the components of a vector
      \param[in] x, the vector
      \param[in-out] solution, vector of positive part values of x components
      \param[in] n, size of x
  */
  void pos_part(double*, double*, int);

  /** Absolute values of the components of a vector
      \param[in] x, the vector
      \param[in-out] solution, vector of absolute values of x components
      \param[in] n, size of x
  */
  void abs_part(double*, double*, int);

  /**
      Input na, a, nb, b
      Output nc, c
      a and b: interger vectors in increasing order
      c : vector of integers of a that are not in b.
      \author Nineb Sheherazade & Dureisseix David.
  */
  void diffns(int *na, int *a, int *nb, int * b, int *nc, int *c);

  /** */
  void sortsn_(int*ddl_i, int *sort, int *n);

#ifdef __cplusplus
}
#endif

#endif /* NSSPACK_H */
