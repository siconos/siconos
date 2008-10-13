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

/*! \file LagrangianLinearRXML.h
\brief xml management for Lagrangian linear relations.
*/
#ifndef __LLRelationXML__
#define __LLRelationXML__

#include "LagrangianRXML.h"
#include "SimpleMatrix.h"
#include "SimpleVector.h"

class SimpleMatrix;
class SimpleVector;

/** XML management for LagrangianRXML
 *
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 3.0.0.
 *   \date 05/25/2004
 *
 *
 *
 * LagrangianLinearRXML allows to manage data of a LLRelation DOM tree.
 */
class LagrangianLinearRXML : public LagrangianRXML
{

private:

  /**xml node for H matrix. */
  xmlNodePtr HNode;
  /**xml node for b vector. */
  xmlNodePtr bNode;
  /**xml node for D matrix. */
  xmlNodePtr DNode;
  /**xml node for F matrix. */
  xmlNodePtr FNode;

  /** Default constructor */
  LagrangianLinearRXML();

public:

  /** Basic constructor.
   *   \param xml pointer to the LagrangianLinearR data.
   */
  LagrangianLinearRXML(xmlNode * LLRelationNode);

  /** Destructor */
  ~LagrangianLinearRXML();

  /** return true if H is given in xmlfile
   */
  inline bool hasH() const
  {
    return HNode;
  }

  /** Return the H of the LLRelationXML
   *   \return The H SimpleMatrix of the LLRelationXML
   */
  inline SimpleMatrix getH()
  {
    return  SiconosDOMTreeTools::getSiconosMatrixValue(HNode);
  }

  /** Change the H matrix value (in xml file or external data file switch his origin position)
   *   \param SiconosMatrix matrix : the new value for H matrix
   */
  void setH(const SiconosMatrix&);

  /** return true if b is given in xmlfile
   */
  inline bool hasB() const
  {
    return bNode;
  }

  /** Return b vector of the LLRelationXML
   *   \return SimpleVector : b vector of the LLRelationXML
   */
  inline SimpleVector getB()
  {
    return  SiconosDOMTreeTools::getSiconosVectorValue(bNode);
  }

  /** Change the b vector value (in xml file or external data file switch his origin position)
   *   \param SiconosVector vector : the new value for b vector
   */
  void setB(const SiconosVector&);

  /** return true if D is given in xmlfile
   */
  inline bool hasD() const
  {
    return DNode;
  }

  /** Return the D of the LLRelationXML
   *   \return The D SimpleMatrix of the LLRelationXML
   */
  inline SimpleMatrix getD()
  {
    return  SiconosDOMTreeTools::getSiconosMatrixValue(DNode);
  }

  /** Change the D matrix value (in xml file or external data file switch his origin position)
   *   \param SiconosMatrix matrix : the new value for D matrix
   */
  void setD(const SiconosMatrix&);

  /** return true if F is given in xmlfile
   */
  inline bool hasF() const
  {
    return FNode;
  }

  /** Return the F of the LLRelationXML
   *   \return The F SimpleMatrix of the LLRelationXML
   */
  inline SimpleMatrix getF()
  {
    return  SiconosDOMTreeTools::getSiconosMatrixValue(FNode);
  }

  /** Change the F matrix value (in xml file or external data file switch his origin position)
   *   \param SiconosMatrix matrix : the new value for F matrix
   */
  void setF(const SiconosMatrix&);
};

#endif
