
/** \class OneStepNSProblemXML
*   \brief This class manages OneStepNSProblem data part
*  \author SICONOS Development Team - copyright INRIA
*   \version 1.0
*   \date 05/17/2004
*
*
*
* OneStepNSProblemXML allows to manage data of a OneStepNSProblem DOM tree.
*/


#ifndef __OneStepNSProblemXML__
#define __OneStepNSProblemXML__


#include <vector>
#include <string>
#include <map>
#include <libxml/tree.h>

//#include "SiconosVector.h"
#include "NewSiconosVector.h"
#include "SiconosMatrix.h"
#include "SiconosDOMTreeTools.h"

//#include "KernelDefaultConfig.h"

using namespace std;

//Tags
const string OSNSP_N = "n";
const string OSNSP_INTERACTION_CONCERNED = "Interaction_Concerned";
const string OSNSP_SOLVER = "Solver";

const string OSNSP_TOLERANCE = "tolerance";
const string OSNSP_MAXITER = "maxIter";
const string OSNSP_NORMTYPE = "normType";
const string OSNSP_SEARCHDIRECTION = "searchDirection";

const string OSNSP_LCPSOLVING = "LcpSolving";
const string OSNSP_RPSOLVING = "RelayPrimalSolving";
const string OSNSP_RDSOLVING = "RelayDualSolving";
const string OSNSP_CFPSOLVING = "ContactFrictionPrimalSolving";
const string OSNSP_CFDSOLVING = "ContactFrictionDualSolving";
const string OSNSP_LEMKE = "Lemke";
const string OSNSP_GSNL = "Gsnl";
const string OSNSP_GCP = "Gcp";
const string OSNSP_LATIN = "Latin";
//const string OSNSP_lemke = "lemke";
//const string OSNSP_gsnl = "gsnl";
//const string OSNSP_gcp = "gcp";
//const string OSNSP_latin = "latin";

#include "XMLTagsName.h"


extern string   DefaultSolver;
extern string   DefaultAlgoName;
extern string   DefaultAlgoNormType;
extern double   DefaultAlgoTolerance;
extern int    DefaultAlgoMaxIter;
extern double   DefaultAlgoSearchDirection;


class OneStepNSProblemXML
{
public:

  OneStepNSProblemXML();

  /** \fn OneStepNSProblemXML(xmlNode * OneStepNSProblemNode)
  *   \brief Build a OneStepNSProblemXML object from a DOM tree describing a OneStepNSProblem
  *   \param OneStepNSProblemNode : the OneStepNSProblem DOM tree
  *   \param vector<int> definedInteractionNumbers : the Interaction numbers effectivly defined in the model
  *   \exception XMLException : if a property of the OneStepNSProblemXML lacks in the DOM tree
  */
  OneStepNSProblemXML(xmlNode * oneStepNSProblemXMLNode, vector<int> definedInteractionNumbers);


  /** \fn int getN()
  *   \brief Return the n value of the OneStepNSProblemXML
  *   \return The n integer of the OneStepNSProblemXML
  */
  inline int getN()
  {
    return SiconosDOMTreeTools::getIntegerContentValue(this->nNode);
  }

  /** \fn void setN(int n)
  *   \brief allows to save the n value of the OneStepNSProblemXML
  *   \param The n integer of the OneStepNSProblemXML
  */
  inline void setN(int n)
  {
    if (this->hasN() == false)
    {
      this->nNode = SiconosDOMTreeTools::createIntegerNode(this->rootNSProblemXMLNode, OSNSP_N, n);
    }
    else SiconosDOMTreeTools::setIntegerContentValue(this->nNode, n);
  }

  /** \fn bool hasN()
   *  \brief returns true if nNode is defined
   *  \return true if nNode is defined
   */
  inline bool hasN()
  {
    return (this->nNode != NULL);
  }


  /** \fn vector<int> getInteractionConcerned()
  *   \brief Return the Interaction numbers of the OneStepNSProblem
  *   \return The Interaction numbers vector of the OneStepNSProblem
  */
  inline vector<int> getInteractionConcerned()
  {
    return this->interactionNumbersVector;
  }

  /** \fn void setInteractionConcerned(vector<int>, bool all=false)
  *   \brief allows to save the Interaction numbers concerned by this OneStepNSProblem
  * This function will manage operations in the XML DOM tree
  * if the interactionConcernedNode is NULL all is created
  * else, the interactions concerned are defined in the XML input file and it only need to check if more interactions are defined
  *   \param vectir<int> : The Interaction numbers vector to save
  *   \param bool : defines if all the interaction of the nsds are concerned by this OneStepNSProblem. Default value is false.
  */
  void setInteractionConcerned(vector<int>, bool all = false);

  /** \fn int getType()
  *   \brief Return the type of the OneStepNSProblemXML
  *   \return The string type of the OneStepNSProblemXML
  */
  inline string getType()
  {
    //return SiconosDOMTreeTools::getStringAttributeValue(this->rootNSProblemXMLNode, OSNSP_TYPE);
    string type((char*)this->rootNSProblemXMLNode->name);
    return type;
  }

  /** \fn bool hasAll()
   *  \brief All is an attribute of the DS_Concerned tag
   *  \return bool : true if attribute all is defined
   */
  inline bool hasAll()
  {
    if (SiconosDOMTreeTools::hasAttributeValue(this->interactionConcernedNode, ALL_ATTRIBUTE))
      return SiconosDOMTreeTools::getBooleanAttributeValue(this->interactionConcernedNode, ALL_ATTRIBUTE);
    else return false;
  }

  /** \fn void setAll(bool all)
  *   \brief Allows to modify the attribute "all" of the Interaction_concerned tag
  *   \param bool : the value to assign to the attribute
  */
  inline void setAll(bool all)
  {
    if (this->hasAll() == false)
    {
      if (all == true)
        xmlNewProp(this->interactionConcernedNode, (xmlChar*)ALL_ATTRIBUTE.c_str(), (xmlChar*)"true");
    }
    else
    {
      if (all == false)
        xmlRemoveProp(xmlHasProp(this->interactionConcernedNode, (xmlChar*)ALL_ATTRIBUTE.c_str()));
    }
  }


  /*============================================================
   *      Solver tags and attributes of the OneStepNSProblem
   *               ( according to SICONOS/Numerics )
   * ===========================================================/
  /** \fn string getSolver()
  *   \brief Return the kind of solver of the OneStepNSProblem (LcpSolving, RelayDualSolving, ...)
  *   \return string : the type of solver
  */
  inline string getSolver()
  {
    if (this->solverNode != NULL)
    {
      xmlNode *node = SiconosDOMTreeTools::findNodeChild(this->solverNode);
      string type((char*)node->name);
      return type;
    }
    /*
     * for the moment, Solver tag is not required, so "default" is returned when solverNode == NULL
     */
    return DefaultSolver;
    //XMLException::selfThrow("OneStepNSProblemXML - getSolver : solverNode == NULL");
  }

  /** \fn string getSolverAlgorithmName()
  *   \brief Return the kind of algorithm of the solver used
  *   \return string : the type of algorithm, "default" is returned if no algorithm is defined
  */
  inline string getSolverAlgorithmName()
  {
    if (this->solverAlgorithmNode != NULL)
    {
      string type((char*)this->solverAlgorithmNode->name);
      return type;
    }
    else
    {
      cout << "Warning : No algorithm defined for the Solver." << endl;
      return DefaultAlgoName;
    }
  }

  /** \fn double getSolverAlgorithmTolerance()
  *   \brief Return the tolerance of the algorithm used
  *   \return double : the tolerance of algorithm, -1.0 is returned if no algorithm is defined
  */
  inline double getSolverAlgorithmTolerance()
  {
    double res = DefaultAlgoTolerance;
    if ((this->solverAlgorithmNode != NULL)
        && xmlHasProp(this->solverAlgorithmNode, (xmlChar*)OSNSP_TOLERANCE.c_str()))
    {
      string type = SiconosDOMTreeTools::getStringAttributeValue(this->solverAlgorithmNode, OSNSP_TOLERANCE);
      res = atof(type.c_str());
    }
    else
      cout << "Warning : No tolerance defined for the Solver." << endl;

    return res;
  }

  /** \fn int getSolverAlgorithmMaxIter()
  *   \brief Return the maximum number of iteration the algorithm can do
  *   \return double : the maxIter number of algorithm, -1 is returned if no algorithm is defined
  */
  inline int getSolverAlgorithmMaxIter()
  {
    int res = DefaultAlgoMaxIter;
    if ((this->solverAlgorithmNode != NULL)
        && xmlHasProp(this->solverAlgorithmNode, (xmlChar*)OSNSP_MAXITER.c_str()))
      res = SiconosDOMTreeTools::getIntegerAttributeValue(this->solverAlgorithmNode, OSNSP_MAXITER);
    else
      cout << "Warning : No maxIter defined for the Solver." << endl;

    return res;
  }

  /** \fn double getSolverAlgorithmSearchDirection()
  *   \brief Return the search direction of the latin algorithm used
  *   \return double : the search direction of the algorithm, -1.0 is returned if no algorithm is defined
  */
  inline double getSolverAlgorithmSearchDirection()
  {
    double res = DefaultAlgoSearchDirection;
    if ((this->solverAlgorithmNode != NULL)
        && xmlHasProp(this->solverAlgorithmNode, (xmlChar*)OSNSP_SEARCHDIRECTION.c_str()))
      res = SiconosDOMTreeTools::getDoubleAttributeValue(this->solverAlgorithmNode, OSNSP_SEARCHDIRECTION);
    else
      cout << "Warning : No search direction defined for the Solver." << endl;

    return res;
  }

  /** \fn string getSolverAlgorithmNormType()
  *   \brief Return the snorm type of the latin algorithm used
  *   \return double : the norm type of the algorithm, "default" is returned if no algorithm is defined
  */
  inline string getSolverAlgorithmNormType()
  {
    string res = DefaultAlgoNormType;
    if ((this->solverAlgorithmNode != NULL)
        && xmlHasProp(this->solverAlgorithmNode, (xmlChar*)OSNSP_NORMTYPE.c_str()))
      res = SiconosDOMTreeTools::getStringAttributeValue(this->solverAlgorithmNode, OSNSP_NORMTYPE);
    else
      cout << "Warning : No norm type defined for the Solver." << endl;

    return res;
  }

  /** \fn bool hasSolver()
   *  \brief checks if tag Solver exists
   *  \return bool : true if tag Solver exists
   */
  inline bool hasSolver()
  {
    return (this->solverNode != NULL);
  }

  /** \fn bool hasSolverAlgorithm()
   *  \brief checks if tag Solver contains a solving method
   *  \return bool : true if tag Solver contains a solving method
   */
  inline bool hasSolverAlgorithm()
  {
    return (this->solverAlgorithmNode != NULL);
  }


  /** \fn void setSolver(string name)
  *   \brief set the kind of solver of the OneStepNSProblem (LcpSolving, RelayDualSolving, ...)
  *   \param string : the type of solver
  *   \param string : the name of the method used for by the solver
  *   \param string : the norm type used by the solver
  *   \param double : the tolerance parameter used by the solver
  *   \param int : the maximum iteration parameter used by the solver
  *   \param double : the search direction parameter used by the solver
  */
  void setSolver(string name, string methodName, string normType,
                 double tolerance, int maxIter, double searchDirection);


  //  /** \fn void updateOneStepNSProblemXML( xmlNode* node, OneStepNSProblemXML* str )
  //  *   \brief makes the operations to create a OneStepNSProblemXML to the StrategyXML
  //  *   \param xmlNode* : the root node of the OneStepNSProblemXML
  //  *   \param OneStepNSProblem* : the OneStepNSProblem of this OneStepNSProblemXML
  //  */
  //  virtual void updateOneStepNSProblemXML( xmlNode* node, OneStepNSProblem* osnspb );


  //  void loadOneStepNSProblem( /*OneStepNSProblem *osnsp*/ );



protected:
  /** root of the model problem (contains n, M, q, interaction concerned, ...) */
  xmlNode* rootNSProblemXMLNode;

  /** root node name "OneStepNSProblem" */
  xmlNode* rootNode;

  xmlNode* interactionConcernedNode;
  xmlNode* solverNode;
  xmlNode* solverAlgorithmNode;

private:

  //Nodes
  xmlNode * nNode;

  //Interactions (Interaction numbers)
  vector<int> interactionNumbersVector;

  //Methods

  /** \fn loadOneStepNSProblemConcernedInteraction(xmlNode * interactionConcernedNode,map<int, bool> definedInteractionNumbers)
  *   \brief load the Interaction numbers for  of the OneStepNSProblem
  *   \param xmlNode * InteractionConcernedNode : the DOM tree node of the concerned Interaction
  *   \param vector<int> definedInteractionNumbers : the Interaction numbers effectivly defined in the model
  *   \exception XMLException
  */
  void loadOneStepNSProblemConcernedInteraction(xmlNode * interactionConcernedNode, vector<int> definedInteractionNumbers);
};


#endif
//$Log: OneStepNSProblemXML.h,v $
//Revision 1.23  2005/03/08 14:23:45  jbarbier
//- modification of constant variables :
//in the XML module, main tags of the XML objects of the strategy are in XMLTagsName.h
//
//in simualtion tools, some constants have been moved to SiconosConst.h
//
//Revision 1.22  2005/01/25 10:33:15  jbarbier
//- modifications for test purpose
//
//Revision 1.21  2005/01/25 09:27:18  jbarbier
//- save of Solver tag in the OneStepNSProblem tag available when saving without XML input file and with partial XML input file
//
//Revision 1.20  2005/01/24 14:33:03  jbarbier
//- OneStepNSProblem > Solver tag is available and managed in the XML part
//
//- tests added on OneStepNSProblem > Solver tag
//
//Revision 1.19  2005/01/14 11:51:25  jbarbier
//- attribute "all" of the OneStepNSProblem terminated
//in OneStepIntegrator and Interaction the attribute is available
//
//Revision 1.18  2004/12/08 12:49:39  jbarbier
//- changes in the XML Schema, respect of the recommandations of the W3C
//version 1.1
//
//- changes in all balises DS, Relation, NSLaw, OneStepIntegrator, OneStepNSProblem
//in the XML files into specific names like LagrangianDS, LinearSystemDS, ...
//for the DS
//
//Revision 1.17  2004/09/27 13:27:14  jbarbier
//
//- Siconos schema renamed : SiconosModelSchema-V1.0.xsd
//
//- new required tags of the model : title, author, description, date, xmlSchema.
//They replace previous attributes author, description and date of the Model.
//
//Revision 1.16  2004/09/23 14:09:24  jbarbier
//- modification of the integrators, the attribute r is always optional.
//
//- modification of the LagrangianNonLinearR. computeInput and computeOutput are
//required.
//
//Revision 1.15  2004/09/15 13:23:14  jbarbier
//- corrections in the OneStepNSProblem, for the XML save. The list of interaction
//linked to the onestepnsproblem is now saved correctly. It is updated before
//during the creation process.
//
//Revision 1.14  2004/09/14 13:49:59  jbarbier
//- files added in sample/ to run run the main_siconos test program
//
//- all the platform can now be saved in an XML file when it is created manually
//
//Revision 1.13  2004/09/10 11:26:28  charlety
//
//_ Integration of the new version of the SiconosVector in the platform. the class simpleVector is used mostly to replace old SiconosVector. When a vector can be composite or simple, like the state of a dynamical system, a pointer on SiconosVector is used, and the vector is initialized simple or composite when the system is initialized.
//
//_ All the tests which worked with the previous version of the vector are OK with the new version.
//
//_ Example SICONOS and bouncingBall are OK
//
//_ some comments have still to be adapted to NewSiconosVector .
//
//_ This version of NewSiconosVector could be called 0.9. some details have to be fixed, it will be done before the end of September.
//
//Revision 1.12  2004/07/29 14:25:44  jbarbier
