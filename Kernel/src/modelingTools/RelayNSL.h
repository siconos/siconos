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
/*! \file RelayNSL.h

*/
#ifndef RELAYNSLAW_H
#define RELAYNSLAW_H

#include "NonSmoothLaw.h"

class NonSmoothLaw;

/** Relay Non Smooth Law
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) Apr 27, 2004
 *
 *
 */
class RelayNSL : public NonSmoothLaw
{

private:

  /** default constructor
   */
  RelayNSL();

public:

  /** basic constructor
      \param law size
   */
  RelayNSL(unsigned int);

  /** constructor with XML object of the RelayNSL
  *  \param NonSmoothLawXML* : the XML object corresponding
  */
  RelayNSL(SP::NonSmoothLawXML);

  /** constructor with the value of the RelayNSL attributes
  *  \param a double value c
  *  \param a double value d
  *  \param unsigned int: size of the ns law
  */
  RelayNSL(double, double, unsigned int);

  ~RelayNSL();

  /** check the ns law to see if it is verified
  *  \return a boolean value whioch determines if the NS Law is verified
  */
  bool isVerified() const;

  /** to get lb
  *  \return the value of lb
  */
  inline const double getLb() const
  {
    return _lb;
  };

  /** to set c
  *  \param a double
  */
  inline void setLb(double newVal)
  {
    _lb = newVal;
  };

  /** to get ub
  *  \return the value of ub
  */
  inline const double getUb() const
  {
    return _ub;
  };


  /** to set ub
  *  \param a double
  */
  inline void setUb(double newVal)
  {
    _ub = newVal;
  };

  /** copy the data of the NonSmoothLaw to the XML tree
  */
  void saveNonSmoothLawToXML();

  /** print the data to the screen
  */
  void display() const;

  /** encapsulates an operation of dynamic casting. Needed by Python interface.
  *  \param NonSmoothLaw* : the law which must be converted
  * \return a pointer on the law if it is of the right type, NULL otherwise
  */
  static RelayNSL* convert(NonSmoothLaw* nsl);

  /** visitors hook
   */
  ACCEPT_STD_VISITORS();

private:
  /** represent the lower bound of the Relay */
  double _lb;

  /** represent the upper bound of the Relay*/
  double _ub;

};

#endif // RELAYNSLAW_H
