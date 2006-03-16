/* Siconos-Kernel version 1.1.3, Copyright INRIA 2005-2006.
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

/** \class LinearDSIOXML
*   \brief This class manages LinearDSIO DSInputOutput data
*  \author SICONOS Development Team - copyright INRIA
*   \version 1.1.3.
*   \date 17/01/2005
*
*
*
* LinearTIRXML allows to manage data of a LinearDSIO in the DOM tree.
*/


#ifndef __LinearDSIOXML__
#define __LinearDSIOXML__

#include "DSInputOutputXML.h"

class SimpleMatrix;


const std::string LINEARDSIO_A = "A";
const std::string LINEARDSIO_B = "B";


class LinearDSIOXML : public DSInputOutputXML
{
public:
  LinearDSIOXML();

  /** \fn LinearDSIOXML(xmlNode * , vector<int> )
  *   \brief Build a DSInputOutputXML object from a DOM tree describing a DSInputOutput
  *   \param xmlNode* : the DSInputOutput DOM tree
  //    *   \param vector<int>  : vector of DSXML numbers to verify DS concerned by the DSInputOutput (identified by number) exists
  */
  LinearDSIOXML(xmlNode * dsioNode/*, vector<int> definedDSNumbers */);
  ~LinearDSIOXML();

  /** \fn SimpleMatrix getA()
  *   \brief Return the A of the LinearDSIOXML
  *   \return The A SimpleMatrix of the LinearDSIOXML
  */
  inline SimpleMatrix getA()
  {
    return SiconosDOMTreeTools::getSiconosMatrixValue(ANode);
  }

  /** \fn SimpleMatrix getB()
  *   \brief Return the B of the LinearDSIOXML
  *   \return The B SimpleMatrix of the LinearDSIOXML
  */
  inline SimpleMatrix getB()
  {
    return SiconosDOMTreeTools::getSiconosMatrixValue(BNode);
  }

  //    /** \fn SimpleMatrix getE()
  //    *   \brief Return the E of the LinearDSIOXML
  //    *   \return The E SimpleMatrix of the LinearDSIOXML
  //    */
  //    inline SimpleMatrix getE()
  //    {
  //      return SiconosDOMTreeTools::getSiconosMatrixValue(this->ENode);
  //    }
  //
  //    /** \fn SimpleVector getA()
  //    *   \brief Return a of the LinearDSIOXML
  //    *   \return SimpleVector : a of LinearDSIOXML
  //    */
  //    inline /*SiconosVector*/SimpleVector getA()
  //    {
  //      return SiconosDOMTreeTools::getSiconosVectorValue(this->aNode);
  //    }
  //
  //    /** \fn void setC(SiconosMatrix *matrix)
  //    *   \brief Change the C matrix values (in xml file or external data file switch his origin position)
  //    *   \param SiconosMatrix matrix : the new value for C matrix
  //    */
  //    void setC(SiconosMatrix *matrix);
  //
  //    /** \fn void setD(SiconosMatrix *matrix)
  //    *   \brief Change the D matrix values (in xml file or external data file switch his origin position)
  //    *   \param SiconosMatrix matrix : the new value for D matrix
  //    */
  //    void setD(SiconosMatrix *matrix);


  /** \fn void setA(SiconosMatrix *matrix)
  *   \brief Change the A matrix values (in xml file or external data file switch his origin position)
  *   \param SiconosMatrix matrix : the new value for A matrix
  */
  void setA(SiconosMatrix *matrix);

  /** \fn void setB(SiconosMatrix *matrix)
  *   \brief Change the B matrix values (in xml file or external data file switch his origin position)
  *   \param SiconosMatrix matrix : the new value for B matrix
  */
  void setB(SiconosMatrix *matrix);


private:
  //Nodes
  xmlNode * ANode;
  xmlNode * BNode;
  //    xmlNode * ENode;
  //    xmlNode * aNode;
};


#endif
