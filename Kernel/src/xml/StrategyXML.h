
/** \class StrategyXML
*   \brief This class manages Strategy data part
*  \author SICONOS Development Team - copyright INRIA
*   \version 1.0
*   \date 05/17/2004
*
*
*
* StrategyXML allows to get OneStepIntegratorXMLs and OneStepNSProblemXMLs from a DOM tree.
*/


#ifndef __STRATEGYXML__
#define __STRATEGYXML__

#include <string>
#include <vector>
#include <map>
#include <libxml/tree.h>


#include "OneStepIntegratorXML.h"
#include "OneStepNSProblemXML.h"

#include "TimeDiscretisationXML.h"
#include "SiconosDOMTreeTools.h"

#include "Strategy.h"
#include "XMLTagsName.h"



//using namespace std;

class Strategy;
class TimeDiscretisation;
class OneStepIntegratorXML;
class OneStepNSProblemXML;
class TimeDiscretisationXML;


class StrategyXML
{
public:
  StrategyXML();
  ~StrategyXML();

  /** \fn StrategyXML(xmlNode * rootStrategyNode)
  *   \brief Build a StrategyXML object from a DOM tree describing a Strategy - verifying DS indicated by numbers in OneStepIntegrators exist
  *   \param rootNSDSNode : the NSDS DOM tree
  *   \param definedNumberDS : the numbers of the defined DS in the Model
  *   \param definedNumberInteraction : the numbers of the defined Interaction in the Model
  *   \exception XMLException : if a property of the LagrangianLinearR lacks in the DOM tree
  */
  StrategyXML(xmlNode * rootStrategyNode, vector<int> definedNumberDS, vector<int> definedNumberInteraction);


  /** \fn vector<OneStepIntegratorXML*> getOneStepIntegratorXML()
  *   \brief Ables to have the OneStepIntegratorXMLs of a StrategyXML
  *   \return the a vector containing OneStepIntegratorXMLs
  */
  inline vector<OneStepIntegratorXML*> getOneStepIntegratorXML()
  {
    return this->oneStepIntegratorXMLVector;
  }

  /** \fn bool hasOneStepIntegratorXML()
  *   \brief Allows to know if the StrategyXML owns one or more OneStepIntegratorXML
  *   \return bool : false if no OneStepIntegratorXML exists in this StratgeyXML
  *                  true otherwise
  */
  inline bool hasOneStepIntegratorXML()
  {
    return (this->oneStepIntegratorXMLVector.size() >= 1);
  }


  /** \fn * OneStepNSProblemXML getOneStepNSProblemXML()
  *   \brief Return the OneStepNSProblemXML pointer of the StrategyXML
  *   \return the OneStepNSProblemXML pointer of the StrategyXML ; NULL if StrategyXML does not have
  */
  inline OneStepNSProblemXML * getOneStepNSProblemXML()
  {
    return this->oneStepNSProblemXML;
  }

  /** \fn TimeDiscretisationXML* getTimeDiscretisationXML()
  *   \brief Return the TimeDiscretisationXML of the StrategyXML
  *   \return the TimeDiscretisationXML of the StrategyXML
  */
  inline TimeDiscretisationXML* getTimeDiscretisationXML()
  {
    return this->timeDiscretisationXML;
  }

  /** \fn inline string getStrategyXMLType()
  *   \brief Return the type of the Strategy
  *   \return the type of the StrategyXML
  */
  inline string getStrategyXMLType()
  {
    return SiconosDOMTreeTools::getStringAttributeValue(this->strategyNode, TYPE_ATTRIBUTE);
  }

  /** \fn inline vector<int> getDefinedNumberInteractionVector()
   *  \brief allows to get the number of the Interactions which are already defined
   *  \return vector<int> : the vector of the number of the Interactions already defined
   */
  inline vector<int> getDefinedNumberInteractionVector()
  {
    return this->definedNumberInteractionVector;
  }

  /** \fn inline xmlNode* getNode()
   *  \brief allows to get the root node of the Strategy
   *  \return xmlNode : the rootNode of the Strategy
   */
  inline xmlNode* getNode()
  {
    return (xmlNode*)this->strategyNode;
  }


  /** \fn bool hasOneStepNSProblemXML()
  *   \brief determines if the Strategy has a OneStepNSProblemXML
  *   \return bool :  false if the oneStepNSProblemXML* is NULL
  */
  inline bool hasOneStepNSProblemXML()
  {
    return (this->oneStepNSProblemXML != NULL);
  }

  /** \fn void updateStrategyXML( xmlNode* node, Strategy* str )
  *   \brief makes the operations to create a StrategyXML to the SiconosModelXML
  *   \param xmlNode* : the root node of the StrategyXML
  *   \param Strategy* : the Strategy of this StrategyXML
  */
  void updateStrategyXML(xmlNode* node, Strategy* str);

  /** \fn void loadStrategy( Strategy* )
  *   \brief loads the data of the Strategy into the StrategyXML (in the DOM tree)
  *   \param Strategy* : the Strategy of this StrategyXML
  */
  void loadStrategy(Strategy*);


private:
  xmlNode *strategyNode;


  /* vector of OneStepIntegratorXML* */
  vector<OneStepIntegratorXML*> oneStepIntegratorXMLVector;

  /* OneStepNSProblemXML - maybe not defined */
  OneStepNSProblemXML *oneStepNSProblemXML;


  /* TimeDiscretisationXML */
  TimeDiscretisationXML *timeDiscretisationXML;

  //---
  /* Map of availables DS : DS defined in the model / for OneStepIntegrator : to know if a DS is already used by another OneStepIntegrator or if well defined */
  map<int, bool> DSAvailabilityMap;

  /* Vector of defined Interaction in the model / for OneStepNSProblem : to know if an interactiion is well defined*/
  vector<int> definedNumberInteractionVector;
  //--


  /** \fn void loadStrategy(xmlNode * rootStrategyNode)
  *   \brief Load the StrategyXML : OneStepIntegratorXMLs, OneStepNSProblemXML and TimeDiscretisationXML components
  *   \exception XMLException : if a property of the Strategy lacks in the DOM tree
  */
  void loadStrategyXML();


  /** \fn void loadOneStepIntegratorXML(xmlNode * rootOneStepIntegratorNode)
  *   \brief Build OneStepIntegratorXML objects from a DOM tree describing OneStepIntegrators
  *   \param rootOneStepIntegratorNode : the OneStepIntegrators DOM tree
  *   \exception XMLException : if a property of OneStepIntegrator lacks in the DOM tree
  */
  void loadOneStepIntegratorXML(xmlNode *rootOneStepIntegratorNode);


  /** \fn void loadOneStepNSProblemXML(xmlNode * rootOneStepNSProblemNode)
  *   \brief Build OneStepNSProblemXML object from a DOM tree describing OneStepNSProblem
  *   \param rootOneStepNSProblemXMLrNode : the OneStepNSProblem DOM tree
  *   \exception XMLException : if a property of OneStepNSProblem lacks in the DOM tree
  */
  void loadOneStepNSProblemXML(xmlNode * rootOneStepNSProblemNode);

  /** \fn void loadTimeDiscretisationXML(xmlNode * rootTimeDiscretisationNode)
  *   \brief Build TimeDiscretisationXML object from a DOM tree describing TimeDiscretisation
  *   \param rootTimeDiscretisationNode : the TimeDiscretisation DOM tree
  */
  void loadTimeDiscretisationXML(xmlNode * rootTimeDiscretisationNode);
};



#endif
