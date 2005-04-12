
/** \class NLinearBCXML
*   \brief This class manages None Linear BC data part
*  \author SICONOS Development Team - copyright INRIA
*   \version 1.0
*   \date 05/25/2004
*
*
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
