/* Siconos-Kernel, Copyright INRIA 2005-2010.
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

/*! \file LinearRXML.h
  \brief xml management for first order linear relations (LinearR).
*/

#ifndef LinearRXML_H
#define LinearRXML_H

#include "RelationXML.hpp"
#include "SimpleMatrix.hpp"
#include "SimpleVector.hpp"

/** XML management for Linear Relations (FirstOrder or Lagrangian)

    \author SICONOS Development Team - copyright INRIA
    \version 3.0.0.
    \date 05/13/2004

    For a First Order Linear Relation:
    \code
    <FirstOrderRelation type="LinearR">
    <C matrixRowSize='1' matrixColSize='3'>
    <row> 1.0  0.0  0.0</row>
    </C>
    </FirstOrderRelation>
    \endcode
    For a Lagrangian Linear Relation:
    \code
    <LagrangianRelation type="LinearR">
    <C matrixRowSize='1' matrixColSize='3'>
    <row> 1.0  0.0  0.0</row>
    </C>
    </LagrangianRelation>
    \endcode
    *
    */
class LinearRXML : public RelationXML
{

private:

  /** C matrix node */
  xmlNodePtr CNode;
  /** D matrix node */
  xmlNodePtr DNode;
  /** F matrix node */
  xmlNodePtr FNode;
  /** e vector node */
  xmlNodePtr eNode;
  /** B matrix node */
  xmlNodePtr BNode;

  /** Default constructor */
  LinearRXML();

public:

  /** Build a LinearRXML object from a DOM tree describing a Relation with LTI type
   *   \param xml pointer to relation data.
   */
  LinearRXML(xmlNodePtr);

  /** Destructor*/
  ~LinearRXML();

  /** return true if CNode exists */
  inline bool hasC() const
  {
    return (CNode);
  }
  /** return true if DNode exists */
  inline bool hasD() const
  {
    return (DNode);
  }
  /** return true if FNode exists */
  inline bool hasF() const
  {
    return (FNode);
  }
  /** return true if eNode exists */
  inline bool hasE() const
  {
    return (eNode);
  }
  /** return true if BNode exists */
  inline bool hasB() const
  {
    return (BNode);
  }

  /** Return true if C is calculated from a plugin */
  inline bool isCPlugin() const
  {
    return xmlHasProp(CNode, (xmlChar *)"matrixPlugin");
  }
  /** Return true if C is calculated from a plugin */
  inline bool isDPlugin() const
  {
    return xmlHasProp(DNode, (xmlChar *)"matrixPlugin");
  }
  /** Return true if C is calculated from a plugin */
  inline bool isFPlugin() const
  {
    return xmlHasProp(FNode, (xmlChar *)"matrixPlugin");
  }
  /** Return true if e is calculated from a plugin */
  inline bool isEPlugin() const
  {
    return xmlHasProp(eNode, (xmlChar *)"vectorPlugin");
  }
  /** Return true if C is calculated from a plugin */
  inline bool isBPlugin() const
  {
    return xmlHasProp(BNode, (xmlChar *)"matrixPlugin");
  }

  /** Return matrix C
   *   \return a SimpleMatrix
   */
  inline SimpleMatrix getC()
  {
    if (isCPlugin())
      XMLException::selfThrow("LinearRXML - getC : no matrix input, C is loaded from a plugin");
    return SiconosDOMTreeTools::getSiconosMatrixValue(CNode);
  }

  /** Return matrix D
   *   \return a SimpleMatrix
   */
  inline SimpleMatrix getD()
  {
    if (isDPlugin())
      XMLException::selfThrow("LinearRXML - getD : no matrix input, D is loaded from a plugin");
    return SiconosDOMTreeTools::getSiconosMatrixValue(DNode);
  }

  /** Return matrix F
   *   \return a SimpleMatrix
   */
  inline SimpleMatrix getF()
  {
    if (isFPlugin())
      XMLException::selfThrow("LinearRXML - getF : no matrix input, F is loaded from a plugin");
    return SiconosDOMTreeTools::getSiconosMatrixValue(FNode);
  }

  /** Return vector e.
   *  \return a SimpleVector
   */
  inline SimpleVector getE()
  {
    if (isEPlugin())
      XMLException::selfThrow("LinearRXML - getE : no vector input, e is loaded from a plugin");
    return SiconosDOMTreeTools::getSiconosVectorValue(eNode);
  }

  /** Return matrix B
   *   \return a SimpleMatrix
   */
  inline SimpleMatrix getB()
  {
    if (isBPlugin())
      XMLException::selfThrow("LinearRXML - getB : no matrix input, B is loaded from a plugin");
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

  /** return the name of the plug-in used for C
   * \return a string
   */
  inline const std::string getCPlugin() const
  {
    if (!isCPlugin())
      XMLException::selfThrow("LinearRXML - getCPlugin : C is not loaded from a plugin");
    return  SiconosDOMTreeTools::getStringAttributeValue(CNode, "matrixPlugin");
  }

  /** return the name of the plug-in used for D
   * \return a string
   */
  inline const std::string getDPlugin() const
  {
    if (!isDPlugin())
      XMLException::selfThrow("LinearRXML - getDPlugin : D is not loaded from a plugin");
    return  SiconosDOMTreeTools::getStringAttributeValue(DNode, "matrixPlugin");
  }

  /** return the name of the plug-in used for F
   * \return a string
   */
  inline const std::string getFPlugin() const
  {
    if (!isFPlugin())
      XMLException::selfThrow("LinearRXML - getFPlugin : F is not loaded from a plugin");
    return  SiconosDOMTreeTools::getStringAttributeValue(FNode, "matrixPlugin");
  }

  /** return the name of the plug-in used for e
   *   \return a string
   */
  inline const std::string getEPlugin() const
  {
    if (!isEPlugin())
      XMLException::selfThrow("LinearRXML - getEPlugin : e is not loaded from a plugin.");
    return  SiconosDOMTreeTools::getStringAttributeValue(eNode, "vectorPlugin");
  }

  /** return the name of the plug-in used for B
   * \return a string
   */
  inline const std::string getBPlugin() const
  {
    if (!isBPlugin())
      XMLException::selfThrow("LinearRXML - getBPlugin : B is not loaded from a plugin");
    return  SiconosDOMTreeTools::getStringAttributeValue(BNode, "matrixPlugin");
  }

  /** to save the plug-in for C.
   *   \param a string (name of the plug-in)
   */
  void setCPlugin(const std::string&);

  /** to save the plug-in for D.
   *   \param a string (name of the plug-in)
   */
  void setDPlugin(const std::string&);

  /** to save the plug-in for F.
   *   \param a string (name of the plug-in)
   */
  void setFPlugin(const std::string&);

  /** to save the plug-in for e.
   *   \param a string (name of the plug-in)
   */
  void setEPlugin(const std::string&);

  /** to save the plug-in for B.
   *   \param a string (name of the plug-in)
   */
  void setBPlugin(const std::string&);

};

TYPEDEF_SPTR(LinearRXML);

#endif
