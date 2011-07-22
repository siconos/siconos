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



#ifndef NUMERICSMATRIX_TEST_FUNCTION_H
#define NUMERICSMATRIX_TEST_FUNCTION_H

#ifdef __cplusplus
extern "C"
{
#endif
  int test_BuildNumericsMatrix(NumericsMatrix** MM);
  int test_prodNumericsMatrix(NumericsMatrix** MM);
  int test_prodNumericsMatrixNumericsMatrix(NumericsMatrix** MM);
  int test_subRowprod(NumericsMatrix* M1, NumericsMatrix* M2);
  int test_rowProdNoDiag(NumericsMatrix* M1, NumericsMatrix* M2);
  int test_subRowprodNonSquare(NumericsMatrix* M3, NumericsMatrix* M4);
  int test_rowProdNoDiagNonSquare(NumericsMatrix* M3, NumericsMatrix* M4);
  int test_SBMRowToDense(SparseBlockStructuredMatrix *M);
  int test_RowPermutationSBM(SparseBlockStructuredMatrix *M);
  int test_ColPermutationSBM(SparseBlockStructuredMatrix *M);
#ifdef __cplusplus
}
#endif

#endif


