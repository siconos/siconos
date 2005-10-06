/* Siconos version 1.0, Copyright INRIA 2005.
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

/** \class LinearTIRXML
 *   \brief This class manages LTIR Relation data
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 1.0
 *   \date 05/13/2004
 *
 *
 *
 * LinearTIRXML allows to manage data of a LTIRelation DOM tree.
 */


#ifndef __LTIRelationXML__
#define __LTIRelationXML__

#include "RelationXML.h"

const std::string  LTIR_C = "C";
const std::string  LTIR_D = "D";
const std::string  LTIR_F = "F";
const std::string  LTIR_E = "e";
const std::string  LTIR_B = "B";
const std::string  LTIR_A = "a";

class LinearTIRXML : public RelationXML
{
public:
  LinearTIRXML();

  /** \fn LinearTIRXML(xmlNode * LTIRelationNode)
   *   \brief Build a LinearTIRXML object from a DOM tree describing a Relation with LTI type
   *   \param LinearTIRXML : the LinearTIR DOM tree
   *   \exception XMLException : if a property of the LinearTI Relation lacks in the DOM tree
   */
  LinearTIRXML(xmlNode * LTIRelationNode);

  ~LinearTIRXML();

  /** \fn SiconosMatrix getC()
   *   \brief Return the C of the LTIRelationXML
   *   \return The C SiconosMatrix of the LTIRelationXML
   */
  inline SiconosMatrix getC()
  {
    return SiconosDOMTreeTools::getSiconosMatrixValue(CNode);
  }

  /** \fn SiconosMatrix getD()
   *   \brief Return the D of the LTIRelationXML
   *   \return The D SiconosMatrix of the LTIRelationXML
   */
  inline SiconosMatrix getD()
  {
    return SiconosDOMTreeTools::getSiconosMatrixValue(DNode);
  }

  /** \fn SiconosMatrix getF()
   *   \brief Return the F of the LTIRelationXML
   *   \return The F SiconosMatrix of the LTIRelationXML
   */
  inline SiconosMatrix getF()
  {
    return SiconosDOMTreeTools::getSiconosMatrixValue(FNode);
  }

  /** \fn SimpleVector getE()
   *   \brief Return e of the LTIRelationXML
   *   \return SimpleVector : e of LTIRelationXML
   */
  inline SimpleVector getE()
  {
    return SiconosDOMTreeTools::getSiconosVectorValue(eNode);
  }

  /** \fn SiconosMatrix getB()
   *   \brief Return the B of the LTIRelationXML
   *   \return The B SiconosMatrix of the LTIRelationXML
   */
  inline SiconosMatrix getB()
  {
    return SiconosDOMTreeTools::getSiconosMatrixValue(BNode);
  }

  /** \fn SimpleVector getA()
   *   \brief Return a of the LTIRelationXML
   *   \return SimpleVector : a of LTIRelationXML
   */
  inline SimpleVector getA()
  {
    return SiconosDOMTreeTools::getSiconosVectorValue(aNode);
  }

  /** \fn void setC(SiconosMatrix )
   *   \brief Change the C matrix values (in xml file or external data file switch his origin position)
   *   \param SiconosMatrix matrix : the new value for C matrix
   */
  void setC(const SiconosMatrix&);

  /** \fn void setD(SiconosMatrix )
   *   \brief Change the D matrix values (in xml file or external data file switch his origin position)
   *   \param SiconosMatrix matrix : the new value for D matrix
   */
  void setD(const SiconosMatrix&);

  /** \fn void setF(SiconosMatrix )
   *   \brief Change the F matrix values (in xml file or external data file switch his origin position)
   *   \param SiconosMatrix matrix : the new value for F matrix
   */
  void setF(const SiconosMatrix&);

  /** \fn void setE(SiconosVector )
   *   \brief Change the e Vector values (in xml file or external data file switch his origin position)
   *   \param SiconosVector *vector : new value of e
   */
  void setE(const SiconosVector&);

  /** \fn void setB(SiconosMatrix )
    *   \brief Change the B matrix values (in xml file or external data file switch his origin position)
    *   \param SiconosMatrix matrix : the new value for B matrix
    */
  void setB(const SiconosMatrix&);

  /** \fn void setA(SiconosVector )
   *   \brief Change the a Vector values (in xml file or external data file switch his origin position)
   *   \param SiconosVector *vector : new value of a
   */
  void setA(const SiconosVector&);

  /** \fn bool hasXX()
   * \brief return true if XXnode exists */
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
  inline bool hasA() const
  {
    return (aNode != NULL);
  }

private:

  //Nodes
  xmlNode * CNode;
  xmlNode * DNode;
  xmlNode * FNode;
  xmlNode * eNode;
  xmlNode * BNode;
  xmlNode * aNode;

};


#endif
