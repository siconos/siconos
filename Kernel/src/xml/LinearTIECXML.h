//$Id: LinearTIECXML.h,v 1.3 2005/01/26 13:50:40 jbarbier Exp $
#ifndef LINEARTIECXML_H
#define LINEARTIECXML_H

#include "LinearECXML.h"

/** \class LinearTIECXML
 *  \brief object to manage XML data of an Linear Time Invariant EqualityConstraint
 *  \author Jean-Michel Barbier
 *  \version 0.1
 *  \date 17/01/2005
 *
 * $Date: 2005/01/26 13:50:40 $
 * $Revision: 1.3 $
 * $Author: jbarbier $
 * $Source: /CVS/Siconos/SICONOS/src/xml/LinearTIECXML.h,v $
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

//$Log: LinearTIECXML.h,v $
//Revision 1.3  2005/01/26 13:50:40  jbarbier
//
//- loading of an XML input file now loads EqualityConstraints and DSInputOutputs
//
//Revision 1.2  2005/01/17 14:09:34  jbarbier
//- LagrangianECXML class added
//
//Revision 1.1  2005/01/17 10:56:27  jbarbier
//- classes EqualityConstraint and DSInputOutput added with inherited classes
//
//- classes EqualityConstraintXML and DSInputOutputXML added with inherited classes
//