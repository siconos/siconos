/* Siconos-Kernel, Copyright INRIA 2005-2012.
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
/*! \file
*/

#ifndef __LsodarOSIXML__
#define __LsodarOSIXML__

#include "SiconosPointers.hpp"
#include "OneStepIntegratorXML.hpp"

/** XML management for LsodarOSI
 *
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 3.0.0.
 *   \date 05/17/2004
 *
 *
 *
 * LsodarOSIXML allows to manage data of a LsodarOSI DOM tree.
 */
class LsodarOSIXML : public OneStepIntegratorXML
{
public:

  LsodarOSIXML();

  /** Build a LsodarOSIXML object from a DOM tree describing LsodarOSI OneStepIntegrator
  *   \param LsodarOSINode : the LsodarOSI DOM tree
  *   \param map<int, bool> definedDSNumbers : to know if DS numbers are not used by another OneStepIntegrator
  */
  LsodarOSIXML(xmlNode * LsodarOSINode);
private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(LsodarOSIXML);

};

DEFINE_SPTR(LsodarOSIXML)
#endif
