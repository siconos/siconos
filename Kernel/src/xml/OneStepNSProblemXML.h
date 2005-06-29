
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

#include "SiconosDOMTreeTools.h"
#include "OneStepNSProblem.h"

//Tags
const std::string OSNSP_N = "n";
const std::string  OSNSP_INTERACTION_CONCERNED = "Interaction_Concerned";
const std::string  OSNSP_SOLVER = "Solver";

const std::string  OSNSP_TOLERANCE = "tolerance";
const std::string  OSNSP_MAXITER = "maxIter";
const std::string  OSNSP_NORMTYPE = "normType";
const std::string  OSNSP_SEARCHDIRECTION = "searchDirection";

const std::string  OSNSP_LCPSOLVING = "LcpSolving";
const std::string  OSNSP_RPSOLVING = "RelayPrimalSolving";
const std::string  OSNSP_RDSOLVING = "RelayDualSolving";
const std::string  OSNSP_CFPSOLVING = "ContactFrictionPrimalSolving";
const std::string  OSNSP_CFDSOLVING = "ContactFrictionDualSolving";
const std::string  OSNSP_LEMKE = "Lemke";
const std::string  OSNSP_GSNL = "Gsnl";
const std::string  OSNSP_GCP = "Gcp";
const std::string  OSNSP_LATIN = "Latin";
//const std::string  OSNSP_lemke = "lemke";
//const std::string  OSNSP_gsnl = "gsnl";
//const std::string  OSNSP_gcp = "gcp";
//const std::string  OSNSP_latin = "latin";

#include "XMLTagsName.h"


extern std::string    DefaultSolver;
extern std::string    DefaultAlgoName;
extern std::string    DefaultAlgoNormType;
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
  OneStepNSProblemXML(xmlNode * oneStepNSProblemXMLNode, std::vector<int> definedInteractionNumbers);

  virtual ~OneStepNSProblemXML();

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
  inline std::vector<int> getInteractionConcerned()
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
  void setInteractionConcerned(std::vector<int>, bool = false);

  /** \fn int getType()
  *   \brief Return the type of the OneStepNSProblemXML
  *   \return The string type of the OneStepNSProblemXML
  */
  inline std::string  getType()
  {
    //return SiconosDOMTreeTools::getStringAttributeValue(this->rootNSProblemXMLNode, OSNSP_TYPE);
    std::string  type((char*)this->rootNSProblemXMLNode->name);
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


  //============================================================
  //      Solver tags and attributes of the OneStepNSProblem
  //               ( according to SICONOS/Numerics )
  //============================================================
  /** \fn string getSolver()
  *   \brief Return the kind of solver of the OneStepNSProblem (LcpSolving, RelayDualSolving, ...)
  *   \return string : the type of solver
  */
  inline std::string  getSolver()
  {
    if (this->solverNode != NULL)
    {
      xmlNode *node = SiconosDOMTreeTools::findNodeChild(this->solverNode);
      std::string  type((char*)node->name);
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
  inline std::string  getSolverAlgorithmName()
  {
    if (this->solverAlgorithmNode != NULL)
    {
      std::string type((char*)this->solverAlgorithmNode->name);
      return type;
    }
    else
    {
      std::cout << "Warning : No algorithm defined for the Solver." << std::endl;
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
      std::string  type = SiconosDOMTreeTools::getStringAttributeValue(this->solverAlgorithmNode, OSNSP_TOLERANCE);
      res = atof(type.c_str());
    }
    else
      std::cout << "Warning : No tolerance defined for the Solver." << std::endl;

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
      std::cout << "Warning : No maxIter defined for the Solver." << std::endl;

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
      std::cout << "Warning : No search direction defined for the Solver." << std::endl;

    return res;
  }

  /** \fn string getSolverAlgorithmNormType()
  *   \brief Return the snorm type of the latin algorithm used
  *   \return double : the norm type of the algorithm, "default" is returned if no algorithm is defined
  */
  inline std::string  getSolverAlgorithmNormType()
  {
    std::string  res = DefaultAlgoNormType;
    if ((this->solverAlgorithmNode != NULL)
        && xmlHasProp(this->solverAlgorithmNode, (xmlChar*)OSNSP_NORMTYPE.c_str()))
      res = SiconosDOMTreeTools::getStringAttributeValue(this->solverAlgorithmNode, OSNSP_NORMTYPE);
    else
      std::cout << "Warning : No norm type defined for the Solver." << std::endl;

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
  void setSolver(std::string  name, std::string  methodName, std::string  normType,
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
  std::vector<int> interactionNumbersVector;

  //Methods

  /** \fn loadOneStepNSProblemConcernedInteraction(xmlNode * interactionConcernedNode,map<int, bool> definedInteractionNumbers)
  *   \brief load the Interaction numbers for  of the OneStepNSProblem
  *   \param xmlNode * InteractionConcernedNode : the DOM tree node of the concerned Interaction
  *   \param vector<int> definedInteractionNumbers : the Interaction numbers effectivly defined in the model
  *   \exception XMLException
  */
  void loadOneStepNSProblemConcernedInteraction(xmlNode * interactionConcernedNode, std::vector<int> definedInteractionNumbers);
};


#endif
