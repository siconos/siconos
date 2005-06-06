#ifndef LAGRANGIANECXML_H
#define LAGRANGIANECXML_H

#include "EqualityConstraintXML.h"

/** \class LagrangianECXML
 *  \brief object to manage XML data of an Lagrangian EqualityConstraint
*  \author SICONOS Development Team - copyright INRIA
 *  \version 0.1
 *  \date 17/01/2005
 *
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
  LagrangianECXML(xmlNode*, std::vector<int>);
  ~LagrangianECXML();
};

#endif // LAGRANGIANECXML_H

