//$Id: NLinearBCXML.h,v 1.6 2004/09/10 08:04:51 jbarbier Exp $

/** \class NLinearBCXML
*   \brief This class manages None Linear BC data part
*   \author J. Blanc-Tranchant
*   \version 1.0
*   \date 05/25/2004
*
*
* $Date: 2004/09/10 08:04:51 $
* $Revision: 1.6 $
* $Author: jbarbier $
* $Source: /CVS/Siconos/SICONOS/src/xml/NLinearBCXML.h,v $
*
* NLinearBCXML allows to manage data of a NLinearBC DOM tree.
*/




#ifndef __NLINEARBCXML__
#define __NLINEARBCXML__

#include <libxml/tree.h>

#include "BoundaryConditionXML.h"

class NLinearBCXML : public BoundaryConditionXML
{
public:

  NLinearBCXML();

  /** \fn NLinearBCXML(xmlNode * NLinearBCNode)
  *   \brief Build a NLinearBCXML object from a DOM tree describing a NLinearBC
  *   \param xmlNode * NLinearBCNode : the NLinearBC DOM tree
  */
  NLinearBCXML(xmlNode * NLinearBCNode);

  ~NLinearBCXML();


  /** \fn void updateBoundaryConditionXML( xmlNode* node)// BoundaryCondition* bc )
   *  \brief makes the operations to add a BoundaryCondition to the DynamicalSystemXML
   *  \param xmlNode* : the root node of this BoundaryCondition
  //     *  \param BoundaryCondition* : the BoundaryCondition of the DS
   */
  void updateBoundaryConditionXML(xmlNode* node/*, BoundaryCondition* bc */);
};


#endif
//$Log: NLinearBCXML.h,v $
//Revision 1.6  2004/09/10 08:04:51  jbarbier
//- XML save available for BoundaryCondition and Interaction
//
//Revision 1.5  2004/07/29 14:25:43  jbarbier
//- $Log: NLinearBCXML.h,v $
//- Revision 1.6  2004/09/10 08:04:51  jbarbier
//- - XML save available for BoundaryCondition and Interaction
//- and $Id: NLinearBCXML.h,v 1.6 2004/09/10 08:04:51 jbarbier Exp $ added
//
