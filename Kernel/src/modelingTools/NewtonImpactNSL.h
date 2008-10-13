/* Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
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
/*! \file NewtonImpactNSL.h

*/
#ifndef NEWTONIMPACTNSL_H
#define NEWTONIMPACTNSL_H

#include "NonSmoothLaw.h"
class NonSmoothLaw;

/** Newton impact Non Smooth Law
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) June 29, 2004
 *
 * This class formalizes the Newton Impact law together with a complementarity condition. i.e.
 * \f[
 * \left\{\begin{array}{l}
 * y \geq 0, \lambda \geq 0, y^{T} \lambda=0\\
 *  if y \leq 0 \quad \mbox{then} \quad \dot y(t^{+}) - e \dot y(t^{-}) \geq 0, \quad  \lambda \geq 0, (\dot y(t^{+}) - e \dot y(t^{-}))^{T} \lambda=0
 * \end{array}\right.
 * \f]
 *
 * nsLawSize is equal to 1.
 *
 */
class NewtonImpactNSL : public NonSmoothLaw
{

private:
  /** The Newton normal coefficient of restitution  */
  double e;

public:

  /** default constructor
   */
  NewtonImpactNSL();

  /** constructor with XML object of the NewtonImpactNSL
  *  \param NonSmoothLawXML* : the XML object corresponding
  */
  NewtonImpactNSL(SP::NonSmoothLawXML);

  /** constructor with the value of the NewtonImpactNSL attributes
  *  \param a double value e
  */
  NewtonImpactNSL(double);

  /** destructor
   */
  ~NewtonImpactNSL();

  /** check the ns law to see if it is verified
  *  \return a boolean value whioch determines if the NS Law is verified
  */
  bool isVerified() const;

  /** getter of e
  *  \return the value of e
  */
  inline const double getE() const
  {
    return e;
  };

  /** setter of e
  *  \param a double to set e
  */
  inline void setE(double newVal)
  {
    e = newVal;
  };

  /** copy the data of the NonSmoothLaw in the XML tree
  */
  void saveNonSmoothLawToXML();

  /** print the data to the screen
  */
  void display() const;

  /** encapsulates an operation of dynamic casting. Needed by Python interface.
  *  \param NonSmoothLaw* : the law which must be converted
  * \return a pointer on the law if it is of the right type, NULL otherwise
  */
  static NewtonImpactNSL* convert(NonSmoothLaw* nsl);
};

#endif // NewtonImpactNSL_H
