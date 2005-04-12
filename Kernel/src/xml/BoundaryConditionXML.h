#ifndef BOUNDARYCONDITIONXML_H
#define BOUNDARYCONDITIONXML_H


#include <libxml/tree.h>
#include "SiconosDOMTreeTools.h"


//const string BC_TYPE = "type";
//const string BC_NLINEAR = "NLinear";
//const string BC_LINEAR = "Linear";
//const string BC_PERIODIC = "Periodic";
#include "XMLTagsName.h"

/** \class BoundaryConditionXML
 *  \brief describes the boundary conditions for a BVP NSDS
*  \author SICONOS Development Team - copyright INRIA
 *  \version 1.0
 *  \date May 24, 2004
 *
 */
class BoundaryConditionXML
{
public:

  BoundaryConditionXML();
  BoundaryConditionXML(xmlNode *);
  ~BoundaryConditionXML();

  /** \fn int getType()
   *  \brief Return the type of the BoundaryConditionXML
   *  \return The string type of the BoundaryConditionXML
   */
  inline string getType()
  {
    string type((char*)rootBCNode->name);
    return type;
    //return SiconosDOMTreeTools::getStringAttributeValue(this->rootBCNode, BC_TYPE);
  }

protected:
  /* node of the DOM tree containing the attributes of the BouindaryCondition
   * The BoundaryCondition is formed like that : <BoundaryCondition> <[type]> <[type]/> <BoundaryCondition/>
   * the node here corrrespond to the tag <[type]>
   */
  xmlNode *rootBCNode;
};

#endif // BOUNDARYCONDITIONXML_H
