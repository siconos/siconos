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
//$Log: BoundaryConditionXML.h,v $
//Revision 1.12  2005/02/24 15:50:20  jbarbier
//- LCP prepared to changes needed for several interactions
//
//- new function for the SiconosMatrices to copy a block matrix into another matrix
//
//- tag names of BoundaryConditionXML, DSInputOutputXML, DSXML, InteractionXML, LagrangianLinearRXML, LagrangianDSXML put in XMLTagNames.h
//
//Revision 1.11  2005/01/18 10:35:17  jbarbier
//- attribute "r" no longer used for Moreau integrator
//
//- modificatoin in the tests for Moreau integrator
//
//- file XMLTagsName.h for further use to regroup all xml tags name...
//
//Revision 1.10  2005/01/11 17:08:30  jbarbier
//- last modification about the BoundaryCondition
//<BoundaryCondition>
//  <[type]>
//  <[type]/>
//<BoundaryCondition/>
//
//- modification of the xml files for this modification
//
//- version 1.2 of the xml schema
//
//Revision 1.9  2004/09/10 08:04:48  jbarbier
//- XML save available for BoundaryCondition and Interaction
//
//Revision 1.8  2004/07/29 14:25:41  jbarbier
