
/** \class NSLawXML
*   \brief This class manages NSLaw data part
*  \author SICONOS Development Team - copyright INRIA
*   \version 1.0
*   \date 05/13/2004
*
*
*
* NSLawXML allows to manage data of a NSLaw DOM tree.
*/


#ifndef __NSLawXML__
#define __NSLawXML__

#include <string>
#include <libxml/tree.h>
#include "SiconosDOMTreeTools.h"

#include "NonSmoothLaw.h"


using namespace std;

class NonSmoothLaw;

//const string NSLAW_TYPE = "type";
#include "XMLTagsName.h"


class NonSmoothLawXML
{
public:
  NonSmoothLawXML();
  NonSmoothLawXML(xmlNode*);
  ~NonSmoothLawXML();

  /** \fn int getType()
  *   \brief Return the type of the NSLawXML
  *   \return The string type of the NSLawXML
  */
  inline string getType()
  {
    //return SiconosDOMTreeTools::getStringAttributeValue(this->rootNSLawXMLNode, NSLAW_TYPE);
    string type((char*)this->rootNSLawXMLNode->name);
    return type;
  }

  /** \fn xmlNode* getNode()
  *   \brief Return the node of the NonSmoothLawXML
  *   \return xmlNode* : the node of the NonSmoothLawXML in the DOM tree
  */
  inline xmlNode* getNode() const
  {
    return this->rootNSLawXMLNode;
  }

  /** \fn void updateNonSmoothLawXML( xmlNode* node, NonSmoothLaw* nsl )
  *   \brief makes the operations to create the NonSmoothLaw of the Interaction
  *   \param xmlNode* : the root node of the NonSmoothLawXML
  *   \param Relation* : the NonSmoothLaw of this NonSmoothLawXML
  */
  void updateNonSmoothLawXML(xmlNode* node, NonSmoothLaw* nsl);

protected:
  xmlNode * rootNSLawXMLNode;
};


#endif
//$Log: NonSmoothLawXML.h,v $
//Revision 1.7  2005/03/08 12:41:39  jbarbier
//- constant variables files modified :
//Some constants added in SiconosConst
//
//all global tag of the modeling tools are in XMLTagsName, other tags are specific to an XML class
//
//Revision 1.6  2004/12/08 12:49:39  jbarbier
//- changes in the XML Schema, respect of the recommandations of the W3C
//version 1.1
//
//- changes in all balises DS, Relation, NSLaw, OneStepIntegrator, OneStepNSProblem
//in the XML files into specific names like LagrangianNLDS, LinearSystemDS, ...
//for the DS
//
//Revision 1.5  2004/09/14 13:49:58  jbarbier
//- files added in sample/ to run run the main_siconos test program
//
//- all the platform can now be saved in an XML file when it is created manually
//
//Revision 1.4  2004/07/29 14:25:44  jbarbier
