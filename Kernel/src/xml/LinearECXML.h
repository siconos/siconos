//$Id: LinearECXML.h,v 1.2 2005/01/26 13:50:40 jbarbier Exp $
#ifndef LINEARECXML_H
#define LINEARECXML_H

#include "EqualityConstraintXML.h"

/** \class LinearECXML
 *  \brief object to manage XML data of an Linear EqualityConstraint
 *  \author Jean-Michel Barbier
 *  \version 0.1
 *  \date 17/01/2005
 *
 * $Date: 2005/01/26 13:50:40 $
 * $Revision: 1.2 $
 * $Author: jbarbier $
 * $Source: /CVS/Siconos/SICONOS/src/xml/LinearECXML.h,v $
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
  LinearECXML(xmlNode*, vector<int>);
  ~LinearECXML();
};

#endif // LINEARECXML_H

//$Log: LinearECXML.h,v $
//Revision 1.2  2005/01/26 13:50:40  jbarbier
//
//- loading of an XML input file now loads EqualityConstraints and DSInputOutputs
//
//Revision 1.1  2005/01/17 10:56:27  jbarbier
//- classes EqualityConstraint and DSInputOutput added with inherited classes
//
//- classes EqualityConstraintXML and DSInputOutputXML added with inherited classes
//