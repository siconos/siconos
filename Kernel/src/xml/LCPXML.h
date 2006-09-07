/* Siconos-Kernel version 1.3.0, Copyright INRIA 2005-2006.
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
/** \class LCPXML
 *   \brief This class manages Lagrangian LCP data
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 1.3.0.
 *   \date 05/18/2004
 *
 *
 *
 * LCPXML XML data management for LCP NSProblem
 *
 */

#ifndef __LCPXML__
#define __LCPXML__

#include "OneStepNSProblemXML.h"

class OneStepNSProblem;
class LCPXML : public OneStepNSProblemXML
{
private:

  //Nodes
  xmlNode * MNode;
  xmlNode * qNode;

public:
  /** \fn LCPXML()
   *   \brief Default constructor
   */
  LCPXML();

  /** \fn LCPXML(xmlNode * LCPNode)
   *   \brief Build a LCPXML object from a DOM tree describing a LCP
   *   \param LCPNode : the LCP DOM tree
   *   \exception XMLException : if a property of the LCP lacks in the DOM tree
   */
  LCPXML(xmlNode * LCPNode);

  /** \fn LCPXML()
   *   \brief Destructor
   */
  ~LCPXML();

  /** \fn SimpleMatrix getM()
   *   \brief Return M
   *   \return The M SimpleMatrix of the LCP
   */
  inline SimpleMatrix getM() const
  {
    return  SiconosDOMTreeTools::getSiconosMatrixValue(MNode);
  }

  /** \fn SimpleVector getQ()
   *   \brief Return vector q
   *   \return SimpleVector : q vector of the LCP
   */
  inline SimpleVector getQ() const
  {
    return SiconosDOMTreeTools::getSiconosVectorValue(qNode);
  }

  /** \fn void setM(const SiconosMatrix &m)
   *   \brief set matrix M
   *   \param The M SiconosMatrix to set
   */
  void setM(const SiconosMatrix &);

  /** \fn void setQ(const SiconosVector &q)
   *   \brief set vector q
   *   \param The q SiconosVector to set
   */
  void setQ(const SiconosVector&);

  /** \fn bool hasM()
   *  \brief returns true if MNode is defined
   *  \return true if MNode is defined
   */
  inline bool hasM() const
  {
    return (MNode != NULL);
  }

  /** \fn bool hasQ()
   *  \brief returns true if qNode is defined
   *  \return true if qNode is defined
   */
  inline bool hasQ() const
  {
    return (qNode != NULL);
  }
};


#endif
