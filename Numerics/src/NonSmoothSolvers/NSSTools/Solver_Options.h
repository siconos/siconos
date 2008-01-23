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

#ifndef Solver_Options_H
#define Solver_Options_H

/*!\file Solver_Options.h
  \brief Structure used to send options (name, parameters ...) to a specific solver-driver (mainly from Kernel to Numerics).
  \author Franck Perignon
*/

/** Structure used to send options (name, parameters ...) to a specific solver-driver (mainly from Kernel to Numerics).
    \param notSet, int equal to false(0) if the parameters below have not been set (ie need to read default values) else true(1)
    \param name of the solver
    \param a list of int parameters (depends on each solver, see solver doc.)
    \param a list of double parameters (depends on each solver, see solver doc.)
    \param int to check storage type (0: double*, 1: SparseBlockStructuredMatrix)
*/
typedef struct
{
  int notSet;
  char solverName[64];
  int * iparam;
  double * dparam;
  int storageType;
} Solver_Options;

#ifdef __cplusplus
extern "C" {
#endif

  /** Read default parameters values for a solver and save them in a Solver_Options structure
      \param[in] driverName, type of the considered problem \n
      0: LCP\n
      1: MLCP\n
      2: FrictionContact2D\n
      3: FrictionContact3D\n
      4: Relay\n
      5: QP\n
      6: NCP
      \param[out] options, structure used to save the parameters
   */
  void readSolverOptions(int, Solver_Options*);

#ifdef __cplusplus
}
#endif



#endif
