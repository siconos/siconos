/* Siconos-Kernel version 1.1.3, Copyright INRIA 2005-2006.
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
#ifndef FrictionContact3D_H
#define FrictionContact3D_H

#include "FrictionContact.h"

/** \class FrictionContact3D
 *  \brief This class is devoted to the formalization and the resolution of the
 * Linear Complementarity Problem (FrictionContact3D)
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 1.1.3.
 *  \date (Creation) Dec 15, 2005
 *
 *
 * This class is derived from FrictionContact one.
 * It provides specific function for 3D friction problems, and corresponds to pfc3D Numerics function.
 *
 */
class FrictionContact3D : public FrictionContact
{

public:

  /** \fn FrictionContact3D()
   *  \brief default constructor
   */
  FrictionContact3D();

  /** \fn FrictionContact3D(OneStepNSProblemXML*, Strategy*=NULL)
   *  \brief xml constructor
   *  \param OneStepNSProblemXML* : the XML linked-object
   *  \param Strategy *: the strategy that owns the problem (optional)
   */
  FrictionContact3D(OneStepNSProblemXML*, Strategy* = NULL);

  /** \fn FrictionContact3D(Strategy * , const std::string& =DEFAULT_SOLVER, const unsigned int& = DEFAULT_ITER, const double& = DEFAULT_TOL,
   *                        const std::string & = DEFAULT_NORMTYPE, const double & = DEFAULT_SEARCHDIR);
   *  \brief constructor from data
   *  \param Strategy *: the strategy that owns this problem
   *  \param string: solver name (optional)
   *  \param int : MaxIter (optional) required if a solver is given
   *  \param double : Tolerance (optional) -> for NLGS, Gcp, Latin
   *  \param string : NormType (optional) -> never used at the time
   *  \param double : SearchDirection (optional) -> for Latin
   */
  FrictionContact3D(Strategy * , const std::string& = DEFAULT_SOLVER, const unsigned int& = DEFAULT_ITER, const double& = DEFAULT_TOL,
                    const std::string & = DEFAULT_NORMTYPE, const double & = DEFAULT_SEARCHDIR);

  /** \fn FrictionContact3D(Strategy*, Sovler*)
   *  \brief constructor from data
   *  \param Strategy *: the strategy that owns this problem
   *  \param Solver* : pointer to object that contains solver algorithm and formulation
   */
  FrictionContact3D(Strategy * , Solver*);

  // --- Destructror ---
  ~FrictionContact3D();

  // --- Others functions ---

  /** \fn void computeQ (void)
   *  \brief compute vector q
   *  \param double : current time
   */
  void computeQ(const double& time);

  /** \fn FrictionContact3D* convert (OneStepNSProblem* ds)
   *  \brief encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param OneStepNSProblem* : the one step problem which must be converted
   * \return a pointer on the problem if it is of the right type, NULL otherwise
   */
  static FrictionContact3D* convert(OneStepNSProblem* osnsp);
};

#endif // FrictionContact3D_H
