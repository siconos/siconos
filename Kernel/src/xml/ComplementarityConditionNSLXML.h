//$Id: ComplementarityConditionNSLXML.h,v 1.4 2004/07/29 14:25:42 jbarbier Exp $

/** \class ComplementarityConditionNSLXML
*   \brief This class manages ComplementarityConditionNSLXML data part
*   \author J. Blanc-Tranchant
*   \version 1.0
*   \date 05/14/2004
*
*
* $Date: 2004/07/29 14:25:42 $
* $Revision: 1.4 $
* $Author: jbarbier $
* $Source: /CVS/Siconos/SICONOS/src/xml/ComplementarityConditionNSLXML.h,v $
*
* ComplementarityConditionNSLXML allows to manage data of a CCNSLaw DOM tree.
*/


#ifndef _ComplementarityConditionNSLXML_
#define _ComplementarityConditionNSLXML_

#include "NonSmoothLawXML.h"


class ComplementarityConditionNSLXML : public NonSmoothLawXML
{
public:
  ComplementarityConditionNSLXML();
  ComplementarityConditionNSLXML(xmlNode*);
  ~ComplementarityConditionNSLXML();

};


#endif
//$Log: ComplementarityConditionNSLXML.h,v $
//Revision 1.4  2004/07/29 14:25:42  jbarbier
//- $Log$ and $Id$ added
//
