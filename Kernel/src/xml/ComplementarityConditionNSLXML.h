
/** \class ComplementarityConditionNSLXML
*   \brief This class manages ComplementarityConditionNSLXML data part
*  \author SICONOS Development Team - copyright INRIA
*   \version 1.0
*   \date 05/14/2004
*
*
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
