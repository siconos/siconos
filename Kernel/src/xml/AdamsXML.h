
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


using namespace std;

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
//$Log: AdamsXML.h,v $
//Revision 1.10  2004/09/23 14:09:24  jbarbier
//- modification of the integrators, the attribute r is always optional.
//
//- modification of the LagrangianNonLinearR. computeInput and computeOutput are
//required.
//
//Revision 1.9  2004/09/14 13:49:55  jbarbier
//- files added in sample/ to run run the main_siconos test program
//
//- all the platform can now be saved in an XML file when it is created manually
//
//Revision 1.8  2004/08/09 15:00:54  jbarbier
//- changes in the cardinality of some attributes of the DynamicalSystem,
//OneStepIntegrator
//
//- modifications in classes Moreau, Lsodar, Adams for these new cardinalities
//
//- corrections in the test xml files
//
//Revision 1.7  2004/07/29 14:25:41  jbarbier
