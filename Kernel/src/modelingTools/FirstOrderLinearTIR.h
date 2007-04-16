/* Siconos-Kernel version 2.0.1, Copyright INRIA 2005-2006.
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
/*! \file FirstOrderLinearTIR.h

*/
#ifndef FirstOrderLinearTIR_H
#define FirstOrderLinearTIR_H

#include "FirstOrderLinearR.h"

/** Linear Time Invariant Relation, derived from class FirstOrderR
 *
 * \author SICONOS Development Team - copyright INRIA
 *  \version 2.0.1.
 *  \date Apr 27, 2004
 *
 *  Linear Relation for First Order Dynamical Systems, with:
 *
 * \f[
 * y = C x + Fz + D \lambda + e \\
 *
 * R = B \lambda
 * \f]
 *
 */
class FirstOrderLinearTIR : public FirstOrderLinearR
{

private:

  /** Default (private) constructor
   */
  FirstOrderLinearTIR();

public:

  /** xml constructor
   *  \param RelationXML* : the XML object corresponding
   */
  FirstOrderLinearTIR(RelationXML*);

  /** create the Relation from a set of data
   *  \param SiconosMatrix : the matrix C
   *  \param SiconosMatrix : the matrix B
   *  \exception RuntimeException
   */
  FirstOrderLinearTIR(const SiconosMatrix& , const SiconosMatrix&);

  /** create the Relation from a set of data
   *  \param SiconosMatrix : C
   *  \param SiconosMatrix : D
   *  \param SiconosMatrix : F
   *  \param SimpleVector  : e
   *  \param SiconosMatrix : B
   *  \exception RuntimeException
   */
  FirstOrderLinearTIR(const SiconosMatrix& , const SiconosMatrix& ,
                      const SiconosMatrix& , const SimpleVector& ,
                      const SiconosMatrix&);

  /** create the Relation from a set of data
   *  \param pointer to SiconosMatrix : the matrix C
   *  \param pointer to SiconosMatrix : the matrix B
   *  \exception RuntimeException
   */
  FirstOrderLinearTIR(SiconosMatrix* , SiconosMatrix*);

  /** create the Relation from a set of data
   *  \param pointer to SiconosMatrix : C
   *  \param pointer to SiconosMatrix : D
   *  \param pointer to SiconosMatrix : F
   *  \param pointer to SimpleVector  : e
   *  \param pointer to SiconosMatrix : B
   *  \exception RuntimeException
   */
  FirstOrderLinearTIR(SiconosMatrix* , SiconosMatrix* ,
                      SiconosMatrix* , SimpleVector* ,
                      SiconosMatrix*);

  /** destructor
   */
  ~FirstOrderLinearTIR();

  // --- OTHER FUNCTIONS ---

  /** Computes y
   *  \param double : current time
   *  \param unsigned int: number of the derivative to compute, optional, default = 0.
   */
  void computeOutput(double, unsigned int = 0);

  /** Computes yFree AND save it into y
   *  \param double : current time
   *  \param unsigned int: number of the derivative to compute, optional, default = 0.
   */
  void computeFreeOutput(double = 0, unsigned int = 0);

  /** Computes lambda
   *  \param double : current time
   *  \param unsigned int: "derivative" order of lambda used to compute input
   */
  void computeInput(double, unsigned int);

  /** encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param Relation * : the relation which must be converted
   * \return a pointer on the relation if it is of the right type, NULL otherwise
   */
  static FirstOrderLinearTIR* convert(Relation *r);
};

#endif
