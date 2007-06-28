/* Siconos-Kernel version 2.1.0, Copyright INRIA 2005-2006.
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
/*! \file
*/

#ifndef __TIMEDISCRETISATIONXML__
#define __TIMEDISCRETISATIONXML__


#include "SiconosDOMTreeTools.h"
#include "SimpleVector.h"

class TimeDiscretisation;

const std::string  TD_H = "h";
const std::string  TD_N = "N";
const std::string  TD_TK = "tk";
const std::string  TD_HMIN = "hMin";
const std::string  TD_HMAX = "hMax";
const std::string  TD_ISCONSTANT = "isConstant";


/** XML management for Time Discretisation
 *
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 2.1.0.
 *   \date 05/17/2004
 *
 */
class TimeDiscretisationXML
{
public:
  TimeDiscretisationXML();

  /** Build a TimeDiscretisationXML object from a DOM tree describing a TimeDiscretisation
  *   \param timeDiscretisationNode : the TimeDiscretisation DOM tree
  *   \exception XMLException : if a property of the TimeDiscretisation lacks in the DOM tree
  */
  TimeDiscretisationXML(xmlNode * TimeDiscretisationNode);


  /** Gets the rootNode of the TimeDiscretisationXML
  *   \return xmlNode* : the rootNode
  */
  inline xmlNode* getRootNode()
  {
    return rootNode;
  }

  /** sets the rootNode of the TimeDiscretisationXML
  *   \param xmlNode* : the rootNode to set
  */
  inline void setRootNode(xmlNode*node)
  {
    rootNode = node;
  }

  /** Return the h value of the TimeDiscretisation
  *   \return The h double value of the TimeDiscretisation
  */
  inline double getH() const
  {
    return  SiconosDOMTreeTools::getContentValue<double>(hNode);
  }

  /** allows to set the h value of the TimeDiscretisation
  *   \param The h double value of the TimeDiscretisation
  */
  inline void setH(const double& d)
  {
    if (!hasH())
    {
      hNode = SiconosDOMTreeTools::createDoubleNode(rootNode, TD_H, d);
    }
    else SiconosDOMTreeTools::setDoubleContentValue(hNode, d);
  }

  /** return the N value of the TimeDiscretisation
  */
  inline int getN() const
  {
    return SiconosDOMTreeTools::getContentValue<int>(NNode);
  }

  /** allows to set the N value of the TimeDiscretisation
  *   \param The N int value of the TimeDiscretisation
  */
  inline void setN(const int& i)
  {
    if (!hasN())
    {
      NNode = SiconosDOMTreeTools::createIntegerNode(rootNode, TD_N, i);
    }
    else SiconosDOMTreeTools::setIntegerContentValue(NNode, i);
  }

  /** Return the tk values of the TimeDiscretisation
  *   \return SimpleVector : tk vector of the TimeDiscretisation
  */
  inline SimpleVector getTk() const
  {
    return SiconosDOMTreeTools::getSiconosVectorValue(tkNode);
  }

  /** set the tk Node value of the TimeDiscretisation
  *   \param The tk SiconosVector of the TimeDiscretisation
  */
  inline void setTkNode(const SimpleVector& v)
  {
    if (!hasTk())
    {
      tkNode = SiconosDOMTreeTools::createVectorNode(rootNode, TD_TK, v);
    }
    else SiconosDOMTreeTools::setSiconosVectorNodeValue(tkNode, v);
  }

  /** Return the hMin value of the TimeDiscretisation / -1.0 if not defined
  *   \return The hMin double value of the TimeDiscretisation
  */
  inline double getHMin()const
  {
    if (hasHMin())
      return SiconosDOMTreeTools::getContentValue<double>(hMinNode);
    else return -1.0;
  }

  /** allows to set the hMin value of the TimeDiscretisation
  *   \param The hMin double value of the TimeDiscretisation
  */
  inline void setHMin(const double& d)
  {
    if (!hasHMin())
    {
      hMinNode = SiconosDOMTreeTools::createDoubleNode(rootNode, TD_HMIN, d);
    }
    else SiconosDOMTreeTools::setDoubleContentValue(hMinNode, d);
  }

  /** Return the hMax value of the TimeDiscretisation / -1.0 if not defined
  *   \return The hMax double value of the TimeDiscretisation
  */
  inline double getHMax() const
  {
    if (hasHMax())
      return SiconosDOMTreeTools::getContentValue<double>(hMaxNode);
    else return -1.0;
  }

  /** allows to set the hMax value of the TimeDiscretisation
  *   \param The hMax double value of the TimeDiscretisation
  */
  inline void setHMax(const double& d)
  {
    if (!hasHMax())
    {
      hMaxNode = SiconosDOMTreeTools::createDoubleNode(rootNode, TD_HMAX, d);
    }
    else SiconosDOMTreeTools::setDoubleContentValue(hMaxNode, d);
  }


  /** defines if the TimeDiscretisation is constant or not
  *   \param bool : true if TimeDiscretisation is constant, false otherwise
  */
  inline void setConstant(const bool& b)
  {
    if (SiconosDOMTreeTools::hasAttributeValue(rootNode, TD_ISCONSTANT))
    {
      SiconosDOMTreeTools::setBooleanAttributeValue(rootNode, TD_ISCONSTANT, b);
    }
    else
    {
      SiconosDOMTreeTools::createBooleanAttribute(rootNode, TD_ISCONSTANT, b);
    }
  }

  /** Return true if the TimeDiscretisation is constant
  *   \return A boolean value : true if TimeDiscretisation is constant, false otherwise
  */
  inline bool isConstant() const
  {
    return SiconosDOMTreeTools::getAttributeValue<bool>(rootNode, TD_ISCONSTANT);
  }

  /** returns true if hNode is defined
  *  \return true if hNode is defined
  */
  inline bool hasH() const
  {
    return (hNode != NULL);
  }

  /** returns true if NNode is defined
  *  \return true if NNode is defined
  */
  inline bool hasN() const
  {
    return (NNode != NULL);
  }

  /** returns true if tkNode is defined
  *  \return true if tkNode is defined
  */
  inline bool hasTk() const
  {
    return (tkNode != NULL);
  }

  /** returns true if hMinNode is defined
  *  \return true if hMinNode is defined
  */
  inline bool hasHMin() const
  {
    return (hMinNode != NULL);
  }

  /** returns true if hMaxNode is defined
  *  \return true if hMaxNode is defined
  */
  inline bool hasHMax() const
  {
    return (hMaxNode != NULL);
  }

  /** makes the operations to create the TimeDiscretisation of the SimulationXML
  *   \param xmlNode* : the root node for the TimeDiscretisationXML
  *   \param NSDS* : the TimeDiscretisation of this TimeDiscretisationXML
  */
  void updateTimeDiscretisationXML(xmlNode*, TimeDiscretisation*);

private:

  //Nodes
  xmlNode * rootNode;
  xmlNode * hNode;
  xmlNode * NNode;
  xmlNode * tkNode;
  xmlNode * hMinNode;
  xmlNode * hMaxNode;

};


#endif
