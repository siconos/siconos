//$Id: LagrangianECXML.h,v 1.2 2005/01/26 13:50:40 jbarbier Exp $
#ifndef LAGRANGIANECXML_H
#define LAGRANGIANECXML_H

#include "EqualityConstraintXML.h"

/** \class LagrangianECXML
 *  \brief object to manage XML data of an Lagrangian EqualityConstraint
 *  \author Jean-Michel Barbier
 *  \version 0.1
 *  \date 17/01/2005
 *
 * $Date: 2005/01/26 13:50:40 $
 * $Revision: 1.2 $
 * $Author: jbarbier $
 * $Source: /CVS/Siconos/SICONOS/src/xml/LagrangianECXML.h,v $
 *
 */

class LagrangianECXML: public EqualityConstraintXML
{
public:

  LagrangianECXML();

  /** \fn LagrangianECXML(xmlNode * , vector<int> )
  *   \brief Build a EqualityConstraintXML object from a DOM tree describing a EqualityConstraint
  *   \param xmlNode* : the EqualityConstraint DOM tree
  *   \param vector<int>  : vector of DSXML numbers to verify DS concerned by the EqualityConstraint (identified by number) exists
  */
  LagrangianECXML(xmlNode*, vector<int>);
  ~LagrangianECXML();
};

#endif // LAGRANGIANECXML_H

//$Log£