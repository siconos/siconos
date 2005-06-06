#ifndef LINEARECXML_H
#define LINEARECXML_H

#include "EqualityConstraintXML.h"

/** \class LinearECXML
 *  \brief object to manage XML data of an Linear EqualityConstraint
*  \author SICONOS Development Team - copyright INRIA
 *  \version 0.1
 *  \date 17/01/2005
 *
 *
 */

class LinearECXML : public EqualityConstraintXML
{
public:

  LinearECXML();

  /** \fn LinearECXML(xmlNode * , vector<int> )
  *   \brief Build a EqualityConstraintXML object from a DOM tree describing a EqualityConstraint
  *   \param xmlNode* : the EqualityConstraint DOM tree
  *   \param vector<int>  : vector of DSXML numbers to verify DS concerned by the EqualityConstraint (identified by number) exists
  */
  LinearECXML(xmlNode*, std::vector<int>);
  virtual ~LinearECXML();
};

#endif // LINEARECXML_H

