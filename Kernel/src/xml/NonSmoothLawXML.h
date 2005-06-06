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

#include "SiconosDOMTreeTools.h"
#include "NonSmoothLaw.h"

class NonSmoothLaw;

//const string NSLAW_TYPE = "type";
class NonSmoothLawXML
{
public:
  NonSmoothLawXML();
  NonSmoothLawXML(xmlNode*);
  virtual ~NonSmoothLawXML();

  /** \fn int getType()
  *   \brief Return the type of the NSLawXML
  *   \return The string type of the NSLawXML
  */
  inline std::string getType()
  {
    //return SiconosDOMTreeTools::getStringAttributeValue(this->rootNSLawXMLNode, NSLAW_TYPE);
    std::string type((char*) rootNSLawXMLNode->name);
    return type;
  }

  /** \fn xmlNode* getNode()
  *   \brief Return the node of the NonSmoothLawXML
  *   \return xmlNode* : the node of the NonSmoothLawXML in the DOM tree
  */
  inline xmlNode* getNode() const
  {
    return rootNSLawXMLNode;
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
