//$Id: PeriodicBCXML.h,v 1.6 2004/09/10 08:05:24 jbarbier Exp $

/** \class PeriodicBCXML
*   \brief This class manages Periodic BC data part
*   \author J. Blanc-Tranchant
*   \version 1.0
*   \date 05/25/2004
*
*
* $Date: 2004/09/10 08:05:24 $
* $Revision: 1.6 $
* $Author: jbarbier $
* $Source: /CVS/Siconos/SICONOS/src/xml/PeriodicBCXML.h,v $
*
* PeriodicBCXML allows to manage data of a PeriodicBC DOM tree.
*/




#ifndef __PERIODICBCXML__
#define __PERIODICBCXML__

#include <libxml/tree.h>

#include "BoundaryConditionXML.h"

class PeriodicBCXML : public BoundaryConditionXML
{
public:

  PeriodicBCXML();

  /** \fn PeriodicBCXML(xmlNode * PeriodicBCNode)
  *   \brief Build a PeriodicBCXML object from a DOM tree describing a PeriodicBC
  *   \param xmlNode * PeriodicBCNode : the PeriodicBC DOM tree
  */
  PeriodicBCXML(xmlNode * PeriodicBCNode);

  ~PeriodicBCXML();


  /** \fn void updateBoundaryConditionXML( xmlNode* node)//, BoundaryCondition* bc )
   *  \brief makes the operations to add a BoundaryCondition to the DynamicalSystemXML
   *  \param xmlNode* : the root node of this BoundaryCondition
  //     *  \param BoundaryCondition* : the BoundaryCondition of the DS
   */
  void updateBoundaryConditionXML(xmlNode* node/*, BoundaryCondition* bc */);
};


#endif
//$Log: PeriodicBCXML.h,v $
//Revision 1.6  2004/09/10 08:05:24  jbarbier
//- XML save available for BoundaryCondition and Interaction
//
//Revision 1.5  2004/07/29 14:25:45  jbarbier
//- $Log: PeriodicBCXML.h,v $
//- Revision 1.6  2004/09/10 08:05:24  jbarbier
//- - XML save available for BoundaryCondition and Interaction
//- and $Id: PeriodicBCXML.h,v 1.6 2004/09/10 08:05:24 jbarbier Exp $ added
//
