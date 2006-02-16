/* Siconos-Kernel version 1.1.1, Copyright INRIA 2005-2006.
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

/** \class NonSmoothLaw
 *  \brief this class contains non-smooth law
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 1.1.1.
 *  \date (Creation) May 05, 2004
 *
 *
 *
 */

#ifndef NSLAW_H
#define NSLAW_H

#include "NonSmoothLawXML.h"
#include "Interaction.h"
#include "SiconosConst.h"

const std::string COMPLEMENTARITYCONDITIONNSLAW = "ComplementarityNSL";
const std::string NEWTONIMPACTNSLAW = "NewtonImpactNSL";
const std::string RELAYNSLAW = "RelayNSL";
const std::string NEWTONIMPACTFRICTIONNSLAW = "NewtonImpactFrictionNSL";
class Interaction;
class NonSmoothLawXML;

class NonSmoothLaw
{
public:

  /** \fn NonSmoothLaw()
   *  \brief default constructor
   */
  NonSmoothLaw();

  /** \fn NonSmoothLaw(NonSmoothLawXML*)
   *  \brief constructor with XML object of the NonSmoothLaw
   *  \param NonSmoothLawXML* : the XML object corresponding
   */
  NonSmoothLaw(NonSmoothLawXML*);

  virtual ~NonSmoothLaw();

  /** \fn bool isVerified()
   *  \brief check if the NS law is verified
   *  \return a boolean value whioch determines if the NS Law is verified
   */
  virtual bool isVerified() const = 0;

  /** \fn inline NonSmoothLawXML* getNonSmoothLawXML()
   *  \brief get the NonSmoothLawXML of the NonSmoothLaw
   *  \return the pointer on the NonSmoothLawXML of the NonSmoothLaw
   */
  inline NonSmoothLawXML* getNonSmoothLawXML()
  {
    return nslawxml;
  }

  /** \fn inline void setNonSmoothLawXML( NonSmoothLawXML* nslawxml )
   *  \brief set the NonSmoothLawXML of the NonSmoothLaw
   *  \param NonSmoothLawXML* : the pointer to set nslawxml
   */
  inline void setNonSmoothLawXML(NonSmoothLawXML* newNslawxml)
  {
    nslawxml = newNslawxml;
  }

  /** \fn inline string getType()
   *  \brief allows to get the type of the NonSmoothLaw
   *  \return string : the type of the NonSmoothLaw
   */
  inline const std::string getType() const
  {
    return nsLawType;
  }

  /** \fn inline string setType(const string &)
   *  \brief set the type of the NonSmoothLaw
   *  \param: string, the type of the NonSmoothLaw
   */
  inline void setType(const std::string& newType)
  {
    nsLawType = newType;
  }

  /** \fn void saveNonSmoothLawToXML()
   *  \brief copy the data of the NonSmoothLaw to the XML tree
   *  \exception RuntimeException
   */
  virtual void saveNonSmoothLawToXML() = 0;

  /** \fn void display()
   *  \brief display the data of the NonSmoothLaw on the standard output
   *  \exception RuntimeException
   */
  virtual void display() const = 0;

protected:

  /** type of the NonSmoothLaw */
  std::string nsLawType;

  /** the XML pbject linked to the NonSmoothLaw to read XML data */
  NonSmoothLawXML *nslawxml;

};

#endif
