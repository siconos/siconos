#ifndef LAGRANGIANLINEARECXML_H
#define LAGRANGIANLINEARECXML_H

#include "EqualityConstraintXML.h"

/** \class LagrangianLinearECXML
 *  \brief object to manage XML data of an Lagrangian EqualityConstraint
*  \author SICONOS Development Team - copyright INRIA
 *  \version 0.1
 *  \date 17/01/2005
 *
 *
 */

class LagrangianLinearECXML: public EqualityConstraintXML
{
public:

  LagrangianLinearECXML();

  /** \fn LagrangianLinearECXML(xmlNode * , vector<int> )
  *   \brief Build a EqualityConstraintXML object from a DOM tree describing a EqualityConstraint
  *   \param xmlNode* : the EqualityConstraint DOM tree
  *   \param vector<int>  : vector of DSXML numbers to verify DS concerned by the EqualityConstraint (identified by number) exists
  */
  LagrangianLinearECXML(xmlNode*, vector<int>);
  ~LagrangianLinearECXML();
};

#endif // LAGRANGIANLINEARECXML_H

