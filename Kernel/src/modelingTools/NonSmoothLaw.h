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

/*! \file NonSmoothLaw.h

*/

#ifndef NSLAW_H
#define NSLAW_H

#include "SiconosConst.h"
#include "RuntimeException.h"

#include "SiconosPointers.hpp"
#include "SiconosVisitor.hpp"

class Interaction;

/** Non Smooth Laws Base Class
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) May 05, 2004
 *
 *
 *
 */
class NonSmoothLaw
{
protected:

  /** type of the NonSmoothLaw */
  std::string nsLawType;

  /** "size" of the NonSmoothLaw */
  unsigned int size;

  /** the XML pbject linked to the NonSmoothLaw to read XML data */
  SP::NonSmoothLawXML nslawxml;

  /** default constructor
   */
  NonSmoothLaw() {};

  /** copy constructor (private=> no copy nor pass-by value allowed)
   */
  NonSmoothLaw(const NonSmoothLaw&) {};

public:

  /** basic constructor
  * \param unsigned int, the non smooth law size
  */
  NonSmoothLaw(unsigned int);

  /** constructor with XML object of the NonSmoothLaw
  *  \param NonSmoothLawXML* : the XML object corresponding
  */
  NonSmoothLaw(SP::NonSmoothLawXML);

  /** destructor
  */
  virtual ~NonSmoothLaw();

  /** check if the NS law is verified
  *  \return a boolean value whioch determines if the NS Law is verified
  */
  virtual bool isVerified() const = 0;

  /** get the NonSmoothLawXML of the NonSmoothLaw
  *  \return the pointer on the NonSmoothLawXML of the NonSmoothLaw
  */
  inline SP::NonSmoothLawXML getNonSmoothLawXML()
  {
    return nslawxml;
  }

  /** set the NonSmoothLawXML of the NonSmoothLaw
  *  \param NonSmoothLawXML* : the pointer to set nslawxml
  */
  inline void setNonSmoothLawXML(SP::NonSmoothLawXML newNslawxml)
  {
    nslawxml = newNslawxml;
  }

  /** to get the size
  *  \return an unsigned int
  */
  inline const unsigned int getNsLawSize() const
  {
    return size;
  }

  /** set the size of the nsLaw
  *  \param an unsigned int
  */
  inline void setNsLawSize(unsigned int newVal)
  {
    size = newVal;
  }

  /** copy the data of the NonSmoothLaw to the XML tree
  *  \exception RuntimeException
  */
  virtual void saveNonSmoothLawToXML() = 0;

  /** display the data of the NonSmoothLaw on the standard output
  *  \exception RuntimeException
  */
  virtual void display() const = 0;

  /** visitors hook
   */
  VIRTUAL_ACCEPT_VISITORS(NonSmoothLaw);

};

DEFINE_SPTR(NonSmoothLaw);

#endif
