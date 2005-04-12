
/** \class NSDSXML
*   \brief This class manages NSDS data part
*  \author SICONOS Development Team - copyright INRIA
*   \version 1.0
*   \date 04/04/2004
*
*
*
* NSDSXML allows to get DSXMLs and InteractionXMLs from a DOM tree.
*/


#ifndef __NSDSXML__
#define __NSDSXML__

#include <string>
#include <vector>
#include <map>
#include <libxml/tree.h>


#include "NonSmoothDynamicalSystem.h"

#include "DSXML.h"
#include "InteractionXML.h"
#include "EqualityConstraintXML.h"
#include "DSInputOutputXML.h"
#include "SiconosDOMTreeTools.h"

#include "XMLTagsName.h"


//using namespace std;


class NonSmoothDynamicalSystem;
class SiconosModelXML;
class DSXML;
class InteractionXML;
class EqualityConstraintXML;
class DSInputOutputXML;


const string NSDS_BVP = "bvp";


class NSDSXML
{
public:
  NSDSXML();

  /** \fn NSDSXML(xmlNode * rootNSDSNode)
  *   \brief Build an NSDSXML object from a DOM tree describing an NSDS
  *   \param rootNSDSNode : the NSDS DOM tree
  */
  NSDSXML(xmlNode * rootNSDSNode);

  ~NSDSXML();

  ///* return a DSs Map */
  //inline map<int, DSXML> getDSMap() const;

  /** \fn xmlNode* getNSDSXMLNode()
  *   \brief Return the root node of the NSDSXML
  *   \return xmlNode* : the root node
  */
  inline xmlNode* getNSDSXMLNode()
  {
    return this->NSDSNode;
  }

  /** \fn DSXML* getDSXML(int number)
  *   \brief Return the DSXML with id number
  *   \param int number : the number of the DSXML to return
  *   \return the DSXML of number number NULL if doesn't exist
  */
  DSXML* getDSXML(int number);

  /** \fn InteractionXML* getInteractionXML(int number)
  *   \brief Return the InteracionXML with id number
  *   \param number : int number : the number of the InteractionXML to return
  *   \exception XMLException
  *   \return the InteractionXML of number number, NULL if doesn't exist
  */
  InteractionXML* getInteractionXML(int number);

  /** \fn EqualityConstraintXML* getEqualityConstraintXML(int number)
  *   \brief Return the EqualityConstraintXML with id number
  *   \param number : int number : the number of the EqualityConstraintXML to return
  *   \exception XMLException
  *   \return the EqualityConstraintXML of number number, NULL if doesn't exist
  */
  EqualityConstraintXML* getEqualityConstraintXML(int number);


  /** \fn inline vector<int> getDSNumbers();
  *   \brief Allows to know the defined DS
  *   \exception XMLException
  *   \return vector DS numbers
  */
  inline vector<int> getDSNumbers()
  {
    return this->definedDSNumbers;
  }

  /** \fn inline vector<int> getInteractionNumbers();
  *   \brief Allows to know the defined interactions
  *   \return vector Interactions integer numbers
  */
  inline vector<int> getInteractionNumbers()
  {
    return this->definedInteractionNumbers;
  }

  /** \fn inline vector<int> getEqualityConstraintNumbers();
  *   \brief Allows to know the defined EqualityConstraints
  *   \return vector EqualityConstraints integer numbers
  */
  inline vector<int> getEqualityConstraintNumbers()
  {
    return this->definedEqualityConstraintNumbers;
  }


  /** \fn bool isBVP()
  *   \brief Allows to know if the NSDS is BVP or not
  *   \return True if the NSDS is BVP false otherwise
  */
  inline bool isBVP()
  {
    if (SiconosDOMTreeTools::hasAttributeValue(this->NSDSNode, NSDS_BVP))
      return SiconosDOMTreeTools::getBooleanAttributeValue(this->NSDSNode, NSDS_BVP);
    else return false;
  }

  /** \fn void setBVP(bool)
  *   \brief Allows to define if the NSDS is BVP
  *   \param True if the NSDS is BVP false otherwise
  */
  inline void setBVP(bool b)
  {
    if (!(SiconosDOMTreeTools::hasAttributeValue(this->NSDSNode, NSDS_BVP)))
    {
      //if( b ) SiconosDOMTreeTools::setBooleanAttributeValue(this->NSDSNode, NSDS_BVP, true);
      if (b == true) xmlNewProp(this->NSDSNode, (xmlChar*)NSDS_BVP.c_str(), (xmlChar*)"true");
    }
    else
    {
      if (b == false) xmlRemoveProp(xmlHasProp(this->NSDSNode, (xmlChar*)NSDS_BVP.c_str()));
    }
  }

  /** \fn void updateNSDSXML(xmlNode*, NonSmoothDynamicalSystem*)
  *   \brief makes the operations to add a NSDS to the SiconosModelXML
  *   \param xmlNode* : the root node for the NSDSXML
  *   \param NonSmoothDynamicalSystem* : the NonSmoothDynamicalSystem of this NSDSXML
  */
  void updateNSDSXML(xmlNode*, NonSmoothDynamicalSystem*);

  /** \fn void loadNonSmoothDynamicalSystem( NonSmoothDynamicalSystem* )
  *   \brief loads the data of the NSDS into the NSDSXML (in the DOM tree)
  *   \param NonSmoothDynamicalSystem* : the NSDS of this NSDSXML
  */
  void loadNonSmoothDynamicalSystem(NonSmoothDynamicalSystem*);


private:
  xmlNode *NSDSNode;

  /* Map of DSs */
  map<int, DSXML*> DSXMLMap;

  /* Map of interactions */
  map<int, InteractionXML*> interactionXMLMap;

  /* Map of EqualityConstraints */
  map<int, EqualityConstraintXML*> equalityConstraintXMLMap;

  /* Map of DSInputOutputs */
  map<int, DSInputOutputXML*> dsInputOutputXMLMap;


  /* vector of DSInputOutput numbers*/
  vector<int> definedDSInputOutputNumbers;

  /* vector of DS numbers*/
  vector<int> definedDSNumbers;

  /* vector of Interaction numbers*/
  vector<int> definedInteractionNumbers;

  /* vector of EqualityConstraint numbers*/
  vector<int> definedEqualityConstraintNumbers;


  /** \fn loadNonSmoothDynamicalSystem()
  *   \brief Load the NonSmoothDynamicalSystem : Interactions and DSs
  *   \exception XMLException
  */
  void loadNonSmoothDynamicalSystem();

  /** \fn loadDSXML(xmlNode * rootDSNode)
  *   \brief Builds DSXML objects from a DOM tree describing DSs
  *   \param xmlNode* : the DSs DOM tree
  *   \exception XMLException : if a property of the NSDS lacks in the DOM tree
  */
  void loadDSXML(xmlNode * rootDSNode);

  /** \fn loadInteractionXML(xmlNode * rootInteractionNode)
  *   \brief Builds InteractionXML objects from a DOM tree describing Interactions
  *   \param xmlNode* : the Interactions DOM tree
  *   \exception XMLException : if a number relating to an Interaction declares in the NSDS is already used
  */
  void loadInteractionXML(xmlNode * rootInteractionNode);

  /** \fn void loadEqualityConstraintXML(xmlNode * rootECNode)
  *   \brief Builds EqualityConstraintXML objects from a DOM tree describing EqualityConstraints
  *   \param xmlNode* : the EqualityConstraints DOM tree
  *   \exception XMLException : if a number relating to an EqualityConstraint declares in the NSDS is already used
  */
  void loadEqualityConstraintXML(xmlNode * rootECNode);

  /** \fn void loadDSInputOutputXML(xmlNode * )
  *   \brief Builds DSInputOutputXML objects from a DOM tree describing DSInputOutputs
  *   \param xmlNode* : the DSInputOutputs DOM tree
  *   \exception XMLException : if a number relating to an DSInputOutput declares in the NSDS is already used
  */
  void loadDSInputOutputXML(xmlNode * rootdsioNode);

  /** \fn map<int, DSInputOutputXML*> getDSInputOutputXMLRelatingToDS( int number )
  *   \brief selects the DSInputOutputXML objects relating to a specific DynamicalSystem
  *   \param map<int, DSInputOutputXML*> : the map containing the DSInputOutputXML for a specific DynamicalSystem
  */
  map<int, DSInputOutputXML*> getDSInputOutputXMLRelatingToDS(int number);

};



#endif
