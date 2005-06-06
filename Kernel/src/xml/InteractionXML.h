/** \class InteractionXML
 *   \brief This class manages Interaction data part
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 1.0
 *   \date (Creation) 04/12/2004
 *

 *
 * InteractionXML allows to manage data of a Interaction DOM tree.
 */


#ifndef __INTERACTIONXML__
#define __INTERACTIONXML__

#include "SiconosDOMTreeTools.h"
#include "Interaction.h"

#include "NonSmoothLawXML.h"
#include "DynamicalSystem.h"
#include "RelationXML.h"

class DynamicalSystem;
class Interaction;
class NonSmoothLawXML;
class RelationXML;

const std::string INTERACTION_STATUS = "Status";
const std::string INTERACTION_Y = "y";
const std::string INTERACTION_LAMBDA = "lambda";
const std::string INTERACTION_NINTER = "nInter";
const std::string INTERACTION_DS_CONCERNED = "DS_Concerned";
const std::string INTERACTION_INTERACTWITHDS_NUMBER = "interactsWithDS_Number";

class InteractionXML
{
public:
  InteractionXML();

  /** \fn InteractionXML(xmlNode * rootNodeInteraction, int number)
   *   \brief Build a InteractionXML object from a DOM tree describing a Interaction
   *   \param xmlNode * interactionNode : the Interaction DOM tree
   *   \param vector<int> definedDSNumbers : vector of DSXML numbers to verify DS concerned by the interaction (identified by number) exists
   */
  InteractionXML(xmlNode * interactionNode, std::vector<int> definedDSNumbers);

  ~InteractionXML();


  /** \fn xmlNode* getInteractionXMLNode()
   *   \brief Return the root node of the InteractionXML
   *   \return xmlNode* : the root node
   */
  inline xmlNode* getInteractionXMLNode()
  {
    return rootInteractionXMLNode;
  }

  /** \fn int getNumber()
   *   \brief Return the number of the InteractionXML
   *   \return The number of the InteractionXML
   */
  inline int getNumber()
  {
    return SiconosDOMTreeTools::getIntegerAttributeValue(this->rootInteractionXMLNode, NUMBER_ATTRIBUTE);
  }

  /** \fn void setNumber(int i)
   *   \brief allows to save the number of the InteractionXML
   *   \return The number to save
   */
  inline void setNumber(int i)
  {
    SiconosDOMTreeTools::setIntegerAttributeValue(this->rootInteractionXMLNode, NUMBER_ATTRIBUTE, i);
  }

  /** \fn bool hasId()
   *  \brief return true if idNode is defined
   *  \return true if idNode is defined
   */
  inline bool hasId()
  {
    return (this->idNode != NULL);
  }

  /** \fn string getId()
   *   \brief Return the id of the InteractionXML
   *   \return The string id of the InteractionXML
   */
  inline std::string getId()
  {
    return SiconosDOMTreeTools::getStringContentValue(this->idNode);
  }

  /** \fn void setId(string s)
   *   \brief allows to save the id of the InteractionXML
   *   \return The string id to save
   */
  inline void setId(std::string s)
  {
    if (hasId() == false)
    {
      idNode = SiconosDOMTreeTools::createStringNode(this->rootInteractionXMLNode, ID_ATTRIBUTE, s);
    }
    else SiconosDOMTreeTools::setStringContentValue(this->idNode, s);
  }

  /** \fn int getNumber()
   *   \brief Return the number of the InteractionXML
   *   \return The number of the InteractionXML
   */
  inline int getNInter()
  {
    return SiconosDOMTreeTools::getIntegerContentValue(this->nInterNode);
  }

  /** \fn void setNumber(int i)
   *   \brief allows to save the number of the InteractionXML
   *   \return The number to save
   */
  inline void setNInter(int nInter)
  {
    if (this->nInterNode == false)
    {
      this->nInterNode = SiconosDOMTreeTools::createIntegerNode(this->rootInteractionXMLNode, INTERACTION_NINTER, nInter);
    }
    else SiconosDOMTreeTools::setIntegerContentValue(this->nInterNode, nInter);
  }

  /** \fn vector<int> getStatus()
   *   \brief Return the status of the InteractionXML
   *   \return vector<int> status : the status of the InteractionXML
   */
  inline std::vector<int> getStatus()
  {
    return SiconosDOMTreeTools::getVectorIntContentValue(this->statusNode);
  }

  /** \fn void setStatus(vector<int> i)
   *   \brief allows to save the status of the InteractionXML
   *   \return vector<int> status : the status to save
   */
  inline void setStatus(std::vector<int> i)
  {
    if (this->statusNode == false)
    {
      this->statusNode = SiconosDOMTreeTools::createVectorIntNode(this->rootInteractionXMLNode, INTERACTION_STATUS, i);
    }
    else SiconosDOMTreeTools::setVectorIntContentValue(this->statusNode, i);
  }

  /** \fn bool hasY()
   *  \brief return true if yNode is defined
   *  \return true if yNode is defined
   */
  inline bool hasY()
  {
    return (this->yNode != NULL);
  }

  /** \fn SimpleVector getY()
   *   \brief Return y vector of the InteractionXML
   *   \return SimpleVector : the y of the InteractionXML
   */
  inline /*SiconosVector*/SimpleVector getY()
  {
    return SiconosDOMTreeTools::getSiconosVectorValue(this->yNode);
  }

  /** \fn void setY(const SiconosVector *v)
   *   \brief allows to save the y of the InteractionXML
   *   \param SiconosVector* : the y to save
   */
  inline void setY(const SiconosVector& v)
  {
    if (hasY() == false)
    {
      yNode = SiconosDOMTreeTools::createVectorNode(rootInteractionXMLNode, INTERACTION_Y, v);
    }
    else SiconosDOMTreeTools::setSiconosVectorNodeValue(this->yNode, v);
  }

  /** \fn bool hasLambda()
   *  \brief return true if lambdaNode is defined
   *  \return true if lambdaNode is defined
   */
  inline bool hasLambda()
  {
    return (this->lambdaNode != NULL);
  }

  /** \fn SimpleVector getLambda()
   *   \brief Return the lambda vector of the InteractionXML
   *   \return SimpleVector : the lambda of the InteractionXML
   */
  inline /*SiconosVector*/ SimpleVector getLambda()
  {
    return SiconosDOMTreeTools::getSiconosVectorValue(this->lambdaNode);
  }

  /** \fn void setLambda(const SiconosVector& v)
   *   \brief allows to save the lambda of the InteractionXML
   *   \return SiconosVector* : the lambda to save
   */
  inline void setLambda(const SiconosVector& v)
  {
    if (hasLambda() == false)
      lambdaNode = SiconosDOMTreeTools::createVectorNode(this->rootInteractionXMLNode, INTERACTION_LAMBDA, v);
    else SiconosDOMTreeTools::setSiconosVectorNodeValue(this->lambdaNode, v);
  }


  /** \fn vector< vector<int> > getDSConcerned()
   *   \brief Return the DSs concerned by the InteractionXML
   *   \return the 2D integer vector who contains the DSs concerned coulpes by the InteractionXML
   */
  inline std::vector< std::vector<int> > getDSConcerned()
  {
    return this->DSCouples;
  }

  /** \fn void setDSConcerned( vector<DynamicalSystem*> )
   *   \brief allows to set the dynamical systems which are interacting together with this interaction
   *   \param vector<DynamicalSystem*> : the dynamical systems in interaction
   */
  void setDSConcerned(std::vector<DynamicalSystem*>);


  /** \fn RelationXML* getRelationXML()
   *   \brief Return the relationXML of the InteractionXML
   *   \return The relationXML of the InteractionXML
   */
  inline RelationXML* getRelationXML()
  {
    return relationXML;
  }


  /** \fn NonSmoothLawXML* getNonSmoothLawXML()
   *   \brief Return the NonSmoothLawXML of the InteractionXML
   *   \return The NonSmoothLawXML of the InteractionXML
   */
  inline NonSmoothLawXML* getNonSmoothLawXML()
  {
    return nSLawXML;
  }

  /** \fn void updateInteractionXML( xmlNode* node, Interaction* inter );
   *   \brief makes the operations to add an Interaction to the NSDS
   *   \param xmlNode* : the root node of the InteractionXML
   *   \param Interaction* : the Interaction of this InteractionXML
   */
  void updateInteractionXML(xmlNode* node, Interaction* inter);

  /** \fn void loadInteraction( Interaction* )
   *   \brief loads the data of the Interaction into the InteractionXML (in the DOM tree)
   *   \param NSDS* : the Interaction of this InteractionXML
   */
  void loadInteraction(Interaction*);

  /** \fn bool hasAll()
   *  \brief All is an attribute of the DS_Concerned tag
   *  \return bool : true if attribute all is defined
   */
  inline bool hasAll()
  {
    if (SiconosDOMTreeTools::hasAttributeValue(this->dsConcernedNode, ALL_ATTRIBUTE))
      return SiconosDOMTreeTools::getBooleanAttributeValue(this->dsConcernedNode, ALL_ATTRIBUTE);
    else return false;
  }

  /** \fn void setAll(bool all)
   *   \brief Allows to modify the attribute "all" of the DS_concerned tag
   *   \param bool : the value to assign to the attribute
   */
  inline void setAll(bool all)
  {
    if (this->hasAll() == false)
    {
      if (all == true)
        xmlNewProp(this->dsConcernedNode, (xmlChar*)ALL_ATTRIBUTE.c_str(), (xmlChar*)"true");
    }
    else
    {
      if (all == false)
        xmlRemoveProp(xmlHasProp(this->dsConcernedNode, (xmlChar*)INTERACTION_DS_CONCERNED.c_str()));
    }
  }


private:

  //Nodes
  xmlNode * rootInteractionXMLNode;

  xmlNode * idNode;
  xmlNode * nInterNode;
  xmlNode * statusNode;
  xmlNode * yNode;
  xmlNode * lambdaNode;
  xmlNode * isActiveNode;
  xmlNode * dsConcernedNode;

  //Couples of DSs (DS numbers)
  std::vector< std::vector<int> > DSCouples;

  //Relation
  RelationXML *relationXML;

  //NSLAW
  NonSmoothLawXML *nSLawXML;

  //Methods

  /** \fn loadInteractionProperties(xmlNode * , vector<int>)
   *   \brief load the different properties of a Interaction
   *   \param xmlNode * nodeInteraction : the DOM tree node of the concern Interaction
   *   \param vector<int> definedDSNumbers : vector of DSXML numbers to verify DS concerned by the interaction (identified by number) exists
   *   \exception XMLException : if a property of the Interaction lacks in the DOM tree
   */
  void loadInteractionProperties(xmlNode * interactionNode, std::vector<int> definedDSNumbers);

  /** \fn loadInteractionConcernedDS(xmlNode * , vector<int>)
   *   \brief load the DSs concerned by this interaction
   *   \param xmlNode * DSConcernedNode : the DOM tree node of DS concerned by the interaction
   *   \param vector<int> definedDSNumbers : vector of DSXML numbers to verify DS concerned by the interaction (identified by number) exists
   *   \exception XMLException : if a DS number not exists
   */
  void loadInteractionConcernedDS(xmlNode * DSConcernedNode, std::vector<int> definedDSNumbers);


};

#endif
