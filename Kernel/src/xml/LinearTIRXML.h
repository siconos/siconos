/* Siconos-Kernel version 2.0.0, Copyright INRIA 2005-2006.
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

/*! \file LinearTIRXML.h

*/

#ifndef __LTIRelationXML__
#define __LTIRelationXML__

#include "RelationXML.h"
#include "SimpleMatrix.h"
#include "SimpleVector.h"

const std::string  LTIR_C = "C";
const std::string  LTIR_D = "D";
const std::string  LTIR_F = "F";
const std::string  LTIR_E = "e";
const std::string  LTIR_B = "B";

class SimpleMatrix;
class SimpleVector;

/** XML management for LinearTIR
 *
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 2.0.0.
 *   \date 05/13/2004
 *
 *
 *
 * LinearTIRXML allows to manage data of a LTIRelation DOM tree.
 */
class LinearTIRXML : public RelationXML
{
public:
  LinearTIRXML();

  /** Build a LinearTIRXML object from a DOM tree describing a Relation with LTI type
  *   \param LinearTIRXML : the LinearTIR DOM tree
  *   \exception XMLException : if a property of the LinearTI Relation lacks in the DOM tree
  */
  LinearTIRXML(xmlNode * LTIRelationNode);

  ~LinearTIRXML();

  /** Return the C of the LTIRelationXML
  *   \return The C SimpleMatrix of the LTIRelationXML
  */
  inline SimpleMatrix getC()
  {
    return SiconosDOMTreeTools::getSiconosMatrixValue(CNode);
  }

  /** Return the D of the LTIRelationXML
  *   \return The D SimpleMatrix of the LTIRelationXML
  */
  inline SimpleMatrix getD()
  {
    return SiconosDOMTreeTools::getSiconosMatrixValue(DNode);
  }

  /** Return the F of the LTIRelationXML
  *   \return The F SimpleMatrix of the LTIRelationXML
  */
  inline SimpleMatrix getF()
  {
    return SiconosDOMTreeTools::getSiconosMatrixValue(FNode);
  }

  /** Return e of the LTIRelationXML
   *   \return SimpleVector : e of LTIRelationXML
   */
  inline SimpleVector getE()
  {
    return SiconosDOMTreeTools::getSiconosVectorValue(eNode);
  }

  /** Return the B of the LTIRelationXML
  *   \return The B SimpleMatrix of the LTIRelationXML
  */
  inline SimpleMatrix getB()
  {
    return SiconosDOMTreeTools::getSiconosMatrixValue(BNode);
  }

  /** Change the C matrix values (in xml file or external data file switch his origin position)
  *   \param SiconosMatrix matrix : the new value for C matrix
  */
  void setC(const SiconosMatrix&);

  /** Change the D matrix values (in xml file or external data file switch his origin position)
  *   \param SiconosMatrix matrix : the new value for D matrix
  */
  void setD(const SiconosMatrix&);

  /** Change the F matrix values (in xml file or external data file switch his origin position)
  *   \param SiconosMatrix matrix : the new value for F matrix
  */
  void setF(const SiconosMatrix&);

  /** Change the e Vector values (in xml file or external data file switch his origin position)
  *   \param SiconosVector *vector : new value of e
  */
  void setE(const SiconosVector&);

  /** Change the B matrix values (in xml file or external data file switch his origin position)
  *   \param SiconosMatrix matrix : the new value for B matrix
  */
  void setB(const SiconosMatrix&);

  /** return true if XXnode exists */
  inline bool hasC() const
  {
    return (CNode != NULL);
  }
  inline bool hasD() const
  {
    return (DNode != NULL);
  }
  inline bool hasF() const
  {
    return (FNode != NULL);
  }
  inline bool hasE() const
  {
    return (eNode != NULL);
  }
  inline bool hasB() const
  {
    return (BNode != NULL);
  }

private:

  //Nodes
  xmlNode * CNode;
  xmlNode * DNode;
  xmlNode * FNode;
  xmlNode * eNode;
  xmlNode * BNode;
};


#endif
