/* Siconos-Kernel, Copyright INRIA 2005-2011.
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
#include<iostream>
//using namespace std;

// ===== Linear DS  ====

// Plugins for A, B (matrices), u and f. See LinearDS.h

extern "C" void computeA(unsigned int *sizeOfA, double* APtr, const double *time)
{
  /* input parameter : time, sizeOfA
   * output parameter : APtr (pointer to SiconosMatrix)
   */
  unsigned int n = *sizeOfA;
  for (unsigned int i = 0; i < n; i++)
  {
    for (unsigned int j = 0; j < n; j++)
      APtr[i + j * n] = (*time);
  }


}

extern "C" void computeb(unsigned int *sizeOfB, double* b, const double *time)
{
  /* input parameter : time
   * output parameter : sizeOfB (size of the vector b); b (pointer to b vector)
   */
  unsigned int n = *sizeOfB;
  for (unsigned int i = 0; i < n; i++)
    b[i] = (*time);
}

extern "C" void computeU(unsigned int *sizeOfU, double* u, const double *time)
{
  /* input parameter : time
   * output parameter: sizeOfU (size of the vector u); u (pointer to u vector)
   */
  unsigned int n = *sizeOfU;
  for (unsigned int i = 0; i < n; i++)
    u[i] = (*time);
}

extern "C" void computeE(unsigned int* rowsOfE, unsigned int* colOfE, double* EPtr, const double* time)
{
  /* input parameter : time
   * output parameter : rowsOfE and colOfE (number of lines and columns in E); E (pointer to SiconosMatrix)
   */
  unsigned int n1 = *rowsOfE;
  unsigned int n2 = *colOfE;
  for (unsigned int i = 0; i < n1; i++)
  {
    for (unsigned int j = 0; j < n2; j++)
      EPtr[i + j * n1] = (*time);
  }


}
