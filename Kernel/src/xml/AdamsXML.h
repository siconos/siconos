
/** \class AdamsXML
*   \brief This class manages Adams data part
*  \author SICONOS Development Team - copyright INRIA
*   \version 1.0
*   \date 05/17/2004
*
*
* AdamsXML allows to manage data of a Adams DOM tree.
*/


#ifndef __ADAMSXML__
#define __ADAMSXML__

#include <libxml/tree.h>
#include "OneStepIntegratorXML.h"


//using namespace std;

//const string ADAMS_R = "r";


class AdamsXML : public OneStepIntegratorXML
{
public:

  AdamsXML();

  /** \fn AdamsXML(xmlNode * AdamsNode)
  *   \brief Build a AdamsXML object from a DOM tree describing Adams OneStepIntegrator
  *   \param AdamsNode : the Adams DOM tree
  *   \param map<int, bool> definedDSNumbersMap : to know if DS numbers are not used by another OneStepIntegrator
  */
  AdamsXML(xmlNode * AdamsNode, map<int, bool> definedDSNumbersMap);


private:

  //Nodes
};


#endif
