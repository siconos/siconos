#ifndef LINEARTIECXML_H
#define LINEARTIECXML_H

#include "LinearECXML.h"

/** \class LinearTIECXML
 *  \brief object to manage XML data of an Linear Time Invariant EqualityConstraint
*  \author SICONOS Development Team - copyright INRIA
 *  \version 0.1
 *  \date 17/01/2005
 *
 *
 */

class LinearTIECXML : public LinearECXML
{
public:

  LinearTIECXML();

  /** \fn LinearTIECXML(xmlNode * , vector<int> )
  *   \brief Build a EqualityConstraintXML object from a DOM tree describing a EqualityConstraint
  *   \param xmlNode* : the EqualityConstraint DOM tree
  *   \param vector<int>  : vector of DSXML numbers to verify DS concerned by the EqualityConstraint (identified by number) exists
  */
  LinearTIECXML(xmlNode*, vector<int>);
  ~LinearTIECXML();
};

#endif // LINEARTIECXML_H

