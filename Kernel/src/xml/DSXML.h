// $Id: DSXML.h,v 1.48 2005/03/14 16:05:27 jbarbier Exp $

/** \class DSXML
*   \brief This class manages DS data part
*   \author J. Blanc-Tranchant
*   \version 1.0
*   \date 04/04/2004
*
* $Date: 2005/03/14 16:05:27 $
* $Revision: 1.48 $
* $Author: jbarbier $
* $Source: /CVS/Siconos/SICONOS/src/xml/DSXML.h,v $
*
* DSXML allows to manage data of a DS DOM tree.
*/


#ifndef __DSXML__
#define __DSXML__


#include <vector>
#include <map>
#include <string>
#include <libxml/tree.h>

//#include "SiconosVector.h"
#include "NewSiconosVector.h"
#include "SimpleVector.h"
#include "SiconosMatrix.h"
#include "SiconosMemory.h"
#include "SiconosMemoryXML.h"
#include "SiconosDOMTreeTools.h"
#include "BoundaryConditionXML.h"
#include "NLinearBCXML.h"
#include "LinearBCXML.h"
#include "PeriodicBCXML.h"
#include "DSInputOutputXML.h"

#include "DynamicalSystem.h"
#include "BoundaryCondition.h"

#include "XMLException.h"

using namespace std;

class DynamicalSystem;
class BoundaryCondition;
class DSInputOutputXML;
//class LagrangianNLDS;
//class LagrangianTIDS;
//class LinearSystemDS;

//Tags
const string DS_N = "n";
const string DS_X0 = "x0";
const string DS_X = "x";
const string DS_XDOT = "xDot";
const string DS_R = "R";
const string DS_XMEMORY = "xMemory";
const string DS_XDOTMEMORY = "xDotMemory";
const string DS_RMEMORY = "rMemory";
const string DS_STEPSINMEMORY = "StepsInMemory";
const string DS_VECTORFIELD = "vectorField";
const string DS_COMPUTEJACOBIANX = "computeJacobianX";
#include "XMLTagsName.h"

class DSXML
{
public:

  DSXML();

  ~DSXML();

  /** \fn DSXML(xmlNode * DSNode)
  *   \brief Build a DSXML object from a DOM tree describing a DS
  *   \param xmlNode * DSNode : the DS DOM tree
  *   \param bool isBVP : if NSDS is IBP DS have boundary condition
  */
  DSXML(xmlNode * DSNode, bool isBVP);

  /** \fn int getNumber()
  *   \brief Return the number of the DSXML
  *   \return The integer number of the DSXML
  */
  inline int getNumber()
  {
    return SiconosDOMTreeTools::getIntegerAttributeValue(this->rootDSXMLNode, NUMBER_ATTRIBUTE);
  }


  /** \fn int getType()
  *   \brief Return the type of the DSXML
  *   \return The string type of the DSXML
  */
  inline string getType()
  {
    //return SiconosDOMTreeTools::getStringAttributeValue(this->rootDSXMLNode, TYPE_ATTRIBUTE);
    string res((char*)this->rootDSXMLNode->name);
    return res;
  }


  /** \fn string getId()
  *   \brief Return the id of the DSXML
  *   \return The string id of the DSXML
  */
  inline string getId()
  {
    return  SiconosDOMTreeTools::getStringContentValue(this->idNode);
  }

  /** \fn void setId(string s)
  *   \brief allows to save the id of the DSXML
  *   \param Integer : The string s of the DSXML
  */
  inline void setId(string s)
  {
    if (this->hasId() == false)
    {
      this->idNode = SiconosDOMTreeTools::createStringNode(this->rootDSXMLNode, ID_ATTRIBUTE, s);
    }
    else SiconosDOMTreeTools::setStringContentValue(this->idNode, s);
  }

  /** \fn int getN()
  *   \brief Return the n of the DSXML
  *   \return The integer n of the DSXML
  */
  inline int getN()
  {
    return  SiconosDOMTreeTools::getIntegerContentValue(this->nNode);
  }

  /** \fn void setN(int n)
  *   \brief allows to save the n of the DSXML
  *   \param Integer : The integer n of the DSXML
  */
  inline void setN(int n)
  {
    if (this->hasN() == false)
    {
      this->nNode = SiconosDOMTreeTools::createIntegerNode(this->rootDSXMLNode, DS_N, n);
    }
    else SiconosDOMTreeTools::setIntegerContentValue(this->nNode, n);
  }


  /** \fn SimpleVector getX0()
  *   \brief Returns the X0 vector of the DSXML
  *   \return SimpleVector : X0 vector of the DSXML
  */
  inline /*SiconosVector*/SimpleVector getX0()
  {
    return  SiconosDOMTreeTools::getSiconosVectorValue(this->x0Node);
  }

  /** \fn void setX0(SiconosVector *v)
  *   \brief allows to set the X0 of the DSXML
  *   \param The X0 SiconosVector to save
  */
  inline void setX0(SiconosVector *v)
  {
    if (this->hasX0() == false)
    {
      this->x0Node = SiconosDOMTreeTools::createVectorNode(this->rootDSXMLNode, DS_X0, v);
    }
    else SiconosDOMTreeTools::setSiconosVectorValue(this->x0Node, v);
  }

  /** \fn SimpleVector getX()
  *   \brief Returns the X vector of the DSXML
  *   \return SimpleVector : X vector of the DSXML
  */
  inline /*SiconosVector*/SimpleVector getX()
  {
    OUT("inline /*SiconosVector*/SimpleVector getX()");
    return  SiconosDOMTreeTools::getSiconosVectorValue(this->xNode);
  }

  /** \fn void setX(SiconosVector *v)
  *   \brief allows to save the X of the DSXML
  *   \return The X SiconosVector to save
  */
  inline void setX(SiconosVector *v)
  {
    if (this->hasX() == false)
    {
      this->xNode = SiconosDOMTreeTools::createVectorNode(this->rootDSXMLNode, DS_X, v);
    }
    else SiconosDOMTreeTools::setSiconosVectorValue(this->xNode, v);
  }

  /** \fn SimpleVector getXDot()
  *   \brief Returns the XDot vector of the DSXML
  *   \return SimpleVector : XDot vector of the DSXML
  */
  inline /*SiconosVector*/SimpleVector getXDot()
  {
    return  SiconosDOMTreeTools::getSiconosVectorValue(this->xDotNode);
  }

  /** \fn void setXDot(SiconosVector *v)
  *   \brief allows to save the XDot of the DSXML
  *   \param The XDot SiconosVector to save
  */
  inline void setXDot(SiconosVector *v)
  {
    if (this->hasXDot() == false)
    {
      this->xDotNode = SiconosDOMTreeTools::createVectorNode(this->rootDSXMLNode, DS_XDOT, v);
    }
    else SiconosDOMTreeTools::setSiconosVectorValue(this->xDotNode, v);
  }


  /** \fn SiconosMemoryXML* getXMemoryXML()
  *   \brief Returns the xMemoryXML* of the DSXML
  *   \return SiconosMemoryXML*
  */
  inline SiconosMemoryXML* getXMemoryXML()
  {
    return this->xMemoryXML;
  }

  /** \fn void setXMemory(SiconosMemory* smem)
  *   \brief allows to save the XMemory of the DSXML
  *   \param SiconosMemory* smem : SiconosMemory to save
  */
  inline void setXMemory(SiconosMemory* smem)
  {
    if (this->hasXMemory() == false)
    {
      this->xMemoryXML = new SiconosMemoryXML(NULL, this->rootDSXMLNode, DS_XMEMORY);
      this->xMemoryNode = this->xMemoryXML->getSiconosMemoryXMLNode();

      this->xMemoryXML->setSiconosMemorySize(smem->getMemorySize());
      this->xMemoryXML->setSiconosMemoryVector(smem->getVectorMemory());
    }
    else
    {
      this->xMemoryXML->setSiconosMemorySize(smem->getMemorySize());
      this->xMemoryXML->setSiconosMemoryVector(smem->getVectorMemory());
    }
  }


  /** \fn SiconosMemoryXML* getXDotMemoryXML()
  *   \brief Returns the xDotMemoryXML* of the DSXML
  *   \return SiconosMemoryXML*
  */
  inline SiconosMemoryXML* getXDotMemoryXML()
  {
    return this->xDotMemoryXML;
  }

  /** \fn void setXDotMemory(SiconosMemory* smem)
  *   \brief allows to save the xDotMemory of the DSXML
  *   \param SiconosMemory* smem : SiconosMemory to save
  */
  inline void setXDotMemory(SiconosMemory* smem)
  {

    if (this->hasXDotMemory() == false)
    {
      this->xDotMemoryXML = new SiconosMemoryXML(NULL, this->rootDSXMLNode, DS_XDOTMEMORY);
      this->xDotMemoryNode = this->xDotMemoryXML->getSiconosMemoryXMLNode();

      this->xDotMemoryXML->setSiconosMemorySize(smem->getMemorySize());
      this->xDotMemoryXML->setSiconosMemoryVector(smem->getVectorMemory());
    }
    else
    {
      this->xDotMemoryXML->setSiconosMemorySize(smem->getMemorySize());
      this->xDotMemoryXML->setSiconosMemoryVector(smem->getVectorMemory());
    }
  }

  /** \fn void setR(SiconosVector *r)
  *   \brief allows to save the R of the DSXML
  *   \param SiconosVector R of the DSXML
  */
  inline void setR(SiconosVector *r)
  {
    if (this->hasR() == false)
    {
      this->rNode = SiconosDOMTreeTools::createVectorNode(this->rootDSXMLNode, DS_R, r);
    }
    else SiconosDOMTreeTools::setSiconosVectorValue(this->rNode, r);
  }

  /** \fn SiconosMemoryXML* getRMemoryXML()
  *   \brief Returns the rMemoryXML* of the DSXML
  *   \return SiconosMemoryXML*
  */
  inline SiconosMemoryXML* getRMemoryXML()
  {
    return this->rMemoryXML;
  }

  /** \fn void setRMemory(SiconosMemory* smem)
  *   \brief allows to save the rMemory of the DSXML
  *   \param SiconosMemory* smem : SiconosMemory to save
  */
  inline void setRMemory(SiconosMemory* smem)
  {
    if (this->hasRMemory() == false)
    {
      this->rMemoryXML = new SiconosMemoryXML(NULL, this->rootDSXMLNode, DS_RMEMORY);
      this->rMemoryNode = this->rMemoryXML->getSiconosMemoryXMLNode();

      this->rMemoryXML->setSiconosMemorySize(smem->getMemorySize());
      this->rMemoryXML->setSiconosMemoryVector(smem->getVectorMemory());
    }
    else
    {
      this->rMemoryXML->setSiconosMemorySize(smem->getMemorySize());
      this->rMemoryXML->setSiconosMemoryVector(smem->getVectorMemory());
    }
  }

  /** \fn int getStepsInMemory()
  *   \brief Returns the steps in memory for the DSXML
  *   \return The integer number of steps in memory for the DSXML
  */
  inline int getStepsInMemory()
  {
    return  SiconosDOMTreeTools::getIntegerContentValue(this->stepsInMemoryNode);
  }

  /** \fn inline void setStepsInMemory(int nb)
  *   \brief allows to save the steps in memory for the DSXML
  *   \param The integer number of steps in memory to save
  */
  inline void setStepsInMemory(int nb)
  {
    if (this->hasStepsInMemory() == false)
    {
      this->stepsInMemoryNode = SiconosDOMTreeTools::createIntegerNode(this->rootDSXMLNode, DS_STEPSINMEMORY, nb);
    }
    else SiconosDOMTreeTools::setIntegerContentValue(this->stepsInMemoryNode, nb);
  }

  /** \fn BoundaryConditionXML * getBoundaryConditionXML()
  *   \brief Returns the BoundaryConditionXML pointer of the DSXML
  *   \return the BoundaryConditionXML pointer of the DSXML ; NULL if DSXML does not have
  */
  inline BoundaryConditionXML * getBoundaryConditionXML()
  {
    return this->boundaryConditionXML;
  }


  /** \fn int getVectorFieldPlugin()
  *   \brief Returns the plugin for the DSXML
  *   \return string which defines the plugin for the DSXML
  */
  inline string getVectorFieldPlugin()
  {
    return  SiconosDOMTreeTools::getStringAttributeValue(this->vectorFieldNode, PLUGIN_ATTRIBUTE);
  }

  /** \fn inline void setVectorFieldPlugin(string plugin)
  *   \brief allows to save the the vectorFieldPlugin for the DSXML
  *   \param The string corresponding to the plugin to save
  */
  inline void setVectorFieldPlugin(string plugin)
  {
    if (this->hasVectorFieldPlugin() == false)
    {
      this->vectorFieldNode = SiconosDOMTreeTools::createSingleNode(this->rootDSXMLNode, DS_VECTORFIELD);
      xmlNewProp(this->vectorFieldNode, (xmlChar*)(PLUGIN_ATTRIBUTE.c_str()), (xmlChar*)plugin.c_str());
    }
    else SiconosDOMTreeTools::setStringAttributeValue(this->vectorFieldNode, PLUGIN_ATTRIBUTE, plugin);
  }

  /** \fn int getComputeJacobianXPlugin()
  *   \brief Returns the plugin for the DSXML
  *   \return string which defines the plugin for the DSXML
  */
  inline string getComputeJacobianXPlugin()
  {
    return  SiconosDOMTreeTools::getStringAttributeValue(this->computeJacobianXNode, PLUGIN_ATTRIBUTE);
  }

  /** \fn inline void setComputeJacobianXPlugin(string plugin)
  *   \brief allows to save the the vectorFieldPlugin for the DSXML
  *   \param The string corresponding to the plugin to save
  */
  inline void setComputeJacobianXPlugin(string plugin)
  {
    if (this->hasComputeJacobianXPlugin() == false)
    {
      this->computeJacobianXNode = SiconosDOMTreeTools::createSingleNode(this->rootDSXMLNode, DS_COMPUTEJACOBIANX);
      xmlNewProp(this->computeJacobianXNode, (xmlChar*)(PLUGIN_ATTRIBUTE.c_str()), (xmlChar*)plugin.c_str());
    }
    else SiconosDOMTreeTools::setStringAttributeValue(this->computeJacobianXNode, PLUGIN_ATTRIBUTE, plugin);
  }

  /** \fn bool hasN()
   *  \brief returns true if nNode is defined
   *  \return true if nNode is defined
   */
  inline bool hasN()
  {
    return (this->nNode != NULL);
  }

  /** \fn bool hasId()
   *  \brief returns true if idNode is defined
   *  \return true if idNode is defined
   */
  inline bool hasId()
  {
    return (this->idNode != NULL);
  }

  /** \fn bool hasX()
   *  \brief returns true if xNode is defined
   *  \return true if xNode is defined
   */
  inline bool hasX()
  {
    return (this->xNode != NULL);
  }

  /** \fn bool hasXDot()
   *  \brief returns true if xDotNode is defined
   *  \return true if xDotNode is defined
   */
  inline bool hasXDot()
  {
    return (this->xDotNode != NULL);
  }

  /** \fn bool hasXMemory()
   *  \brief returns true if xMemoryNode is defined
   *  \return true if xMemoryNode is defined
   */
  inline bool hasXMemory()
  {
    return (this->xMemoryNode != NULL);
  }

  /** \fn bool hasXDotMemory()
   *  \brief returns true if xDotMemoryNode is defined
   *  \return true if xDotMemoryNode is defined
   */
  inline bool hasXDotMemory()
  {
    return (this->xDotMemoryNode != NULL);
  }

  /** \fn bool hasX0()
   *  \brief returns true if x0Node is defined
   *  \return true if x0Node is defined
   */
  inline bool hasX0()
  {
    return (this->x0Node != NULL);
  }

  /** \fn bool hasStepsInMemory()
   *  \brief returns true if stepsInMemoryNode is defined
   *  \return true if stepsInMemoryNode is defined
   */
  inline bool hasStepsInMemory()
  {
    return (this->stepsInMemoryNode != NULL);
  }

  /** \fn bool hasR()
   *  \brief returns true if R is defined
   *  \return true if R is defined
   */
  inline bool hasR()
  {
    return (this->rNode != NULL);
  }

  /** \fn bool hasRMemory()
   *  \brief returns true if RMemory is defined
   *  \return true if RMemory is defined
   */
  inline bool hasRMemory()
  {
    return (this->rMemoryNode != NULL);
  }

  /** \fn bool hasBoundaryCondition()
   *  \brief returns true if boundaryConditionNode is defined
   *  \return true if boundaryConditionNode is defined
   */
  inline bool hasBoundaryCondition()
  {
    return (this->boundaryConditionNode != NULL);
  }

  /** \fn bool hasComputeJacobianXPlugin()
   *  \brief returns true if computeJacobianXNode is defined
   *  \return true if computeJacobianXNode is defined
   */
  inline bool hasComputeJacobianXPlugin()
  {
    return (this->computeJacobianXNode != NULL);
  }

  /** \fn bool hasVectorFieldPlugin()
   *  \brief returns true if vectorFieldNode is defined
   *  \return true if vectorFieldNode is defined
   */
  inline bool hasVectorFieldPlugin()
  {
    return (this->vectorFieldNode != NULL);
  }

  /** \fn void updateDynamicalSystemXML(xmlNode*, DynamicalSystem*, BoundaryCondition*)
  *   \brief makes the operations to add a DynamicalSystem to the NSDSXML
  *   \param xmlNode* : the root node of this DynamicalSystem
  *   \param DynamicalSystem* : the DynamicalSystem of this DSXML
  *   \param BoundaryCondition* : the BoundaryCondition of the DS if the NSDS is BVP (optional)
  */
  void updateDynamicalSystemXML(xmlNode*, DynamicalSystem*, BoundaryCondition* bc = NULL);

  /** \fn void loadDS( DynamicalSystem* )
  *   \brief loads the depending data of the DynamicalSystem into the DSXML (the BoundaryCondition if exists)
  *   \param DynamicalSystem* : the DynamicalSystem of this DSXML
  */
  void loadDS(DynamicalSystem*);


  /** \fn DSInputOutputXML* getDSInputOutputXML(int number)
  *   \brief Return the DSInputOutputXML with id number
  *   \param number : int number : the number of the DSInputOutputXML to return
  *   \exception XMLException
  *   \return the DSInputOutputXML of number number, NULL if doesn't exist
  */
  DSInputOutputXML* getDSInputOutputXML(int number);

  /** \fn inline vector<int> getDSInputOutputNumbers();
  *   \brief Allows to know the defined DSInputOutputs
  *   \return vector DSInputOutputs integer numbers
  */
  inline vector<int> getDSInputOutputNumbers()
  {
    return this->definedDSInputOutputNumbers;
  }

  /** \fn void loadDSInputOutputXML(xmlNode * )
  *   \brief Builds DSInputOutputXML objects from a DOM tree describing DSInputOutputs
  *   \param xmlNode* : the DSInputOutputs DOM tree
  *   \exception XMLException : if a number relating to an DSInputOutput declares in the NSDS is already used
  */
  void loadDSInputOutputXML(xmlNode * rootdsioNode);


protected:
  xmlNode * rootDSXMLNode;
  xmlNode * parentNode;

  //Object
  BoundaryConditionXML * boundaryConditionXML;  //Maybe not defined (if not BVP NSDS)

  /** \fn loadDSProperties(xmlNode * DSnode)
  *   \brief load the different properties of a DS
  *   \param bool isBVP : if NSDS is BVP DS have boundary condition
  *   \exception XMLException : if a property of the DS lacks in the DOM tree
  */
  void loadDSProperties(bool isBVP);


  /** \fn void loadBoundaryConditionXML(xmlNode * rootBoundaryConditionNode)
  *   \brief Build BoundaryConditionXML object from a DOM tree describing BoundaryCondition
  *   \param rootBoundaryConditionXMLNode : the BoundaryCondition DOM tree
  *   \exception XMLException : if the type of the BoundaryCondition given in the DOM tree  does not exist
  */
  void loadBoundaryConditionXML(xmlNode * rootBoundaryConditionNode);

  SiconosMemoryXML * xMemoryXML;
  SiconosMemoryXML * xDotMemoryXML;
  SiconosMemoryXML * rMemoryXML;

  /* Map of DSInputOutputs */
  map<int, DSInputOutputXML*> dsInputOutputXMLMap;

  /* vector of DSInputOutput numbers*/
  vector<int> definedDSInputOutputNumbers;

private:
  //Nodes

  xmlNode * idNode;
  xmlNode * nNode;
  xmlNode * x0Node;
  xmlNode * xNode;
  xmlNode * xDotNode;
  xmlNode * xMemoryNode;
  xmlNode * xDotMemoryNode;
  xmlNode * stepsInMemoryNode;
  xmlNode * vectorFieldNode;
  xmlNode * computeJacobianXNode;
  xmlNode * boundaryConditionNode;
  xmlNode * dsInputOutputNode;
  xmlNode * rNode;
  xmlNode * rMemoryNode;

  //    //Object
  //    BoundaryConditionXML * boundaryConditionXML;  //Maybe not defined (if not BVP NSDS)

  //Methods

  //    /** \fn loadDSProperties(xmlNode * DSnode)
  //    *   \brief load the different properties of a DS
  //    *   \param bool isBVP : if NSDS is BVP DS have boundary condition
  //    *   \exception XMLException : if a property of the DS lacks in the DOM tree
  //    */
  //    void loadDSProperties(bool isBVP);
  //
  //
  //    /** \fn void loadBoundaryConditionXML(xmlNode * rootBoundaryConditionNode)
  //    *   \brief Build BoundaryConditionXML object from a DOM tree describing BoundaryCondition
  //    *   \param rootBoundaryConditionXMLNode : the BoundaryCondition DOM tree
  //    *   \exception XMLException : if the type of the BoundaryCondition given in the DOM tree  does not exist
  //    */
  //    void loadBoundaryConditionXML(xmlNode * rootBoundaryConditionNode);

};

#endif
//$Log: DSXML.h,v $
//Revision 1.48  2005/03/14 16:05:27  jbarbier
//- manual creation of DSInputOutput saving OK
//
//- in progress for EqualityConstraint
//
//Revision 1.47  2005/03/09 15:30:36  jbarbier
//- add of LagrangianEC class
//
//- in progress : implementation of the EqualityConstraint and DSInputOutput - create methods
//
//Revision 1.46  2005/03/08 12:41:38  jbarbier
//- constant variables files modified :
//Some constants added in SiconosConst
//
//all global tag of the modeling tools are in XMLTagsName, other tags are specific to an XML class
//
//Revision 1.45  2005/03/04 15:35:26  jbarbier
//- README files added for some samples
//
//- beginning of the refactoring of XML module constants
//
//Revision 1.44  2005/02/24 15:50:20  jbarbier
//- LCP prepared to changes needed for several interactions
//
//- new function for the SiconosMatrices to copy a block matrix into another matrix
//
//- tag names of BoundaryConditionXML, DSInputOutputXML, DSXML, InteractionXML, LagrangianLinearRXML, LagrangianNLDSXML put in XMLTagNames.h
//
//Revision 1.43  2005/01/18 17:07:44  charlety
//
//_ added autotools makefiles for sample directory
//
//Revision 1.42  2005/01/18 10:35:17  jbarbier
//- attribute "r" no longer used for Moreau integrator
//
//- modificatoin in the tests for Moreau integrator
//
//- file XMLTagsName.h for further use to regroup all xml tags name...
//
//Revision 1.41  2004/12/08 12:49:38  jbarbier
//- changes in the XML Schema, respect of the recommandations of the W3C
//version 1.1
//
//- changes in all balises DS, Relation, NSLaw, OneStepIntegrator, OneStepNSProblem
//in the XML files into specific names like LagrangianNLDS, LinearSystemDS, ...
//for the DS
//
//Revision 1.40  2004/09/27 08:24:26  charlety
//
//_ Modifications in doxygen comments.
//
//Revision 1.39  2004/09/10 11:26:26  charlety
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
//Revision 1.38  2004/09/10 08:04:48  jbarbier
//- XML save available for BoundaryCondition and Interaction
//
//Revision 1.37  2004/08/23 14:30:02  jbarbier
//- All the dynamical systems can be created in a comand program and added to a
//NSDS. The save is OK, but the creation of the boundary conditions is not yet
//finished.
//
//Revision 1.36  2004/08/20 15:26:45  jbarbier
//- creation of a Model and save in the XML is ok
//- creation of a NSDS and save in the XML is ok
//- creation of a NonLinearSystemDS and save in the XML is OK
//
//Revision 1.35  2004/08/10 12:04:30  jbarbier
//- save of the plugin's name for fInt
//
//Revision 1.34  2004/08/04 11:03:23  jbarbier
//- about the SiconosMemory : when a SiconosMemory has a maxSize greater than the
//number of steps in memory required by an integrator, the oldest SiconosVector
//are deleted
//
//- the way to initialize the SiconosMemory by the integrator has been updated to
//match with these changes
//
//Revision 1.33  2004/08/03 12:07:12  jbarbier
//- all test on th eModel are successfull
//
//- new tests on the Model with the opening of XML file
//
//- link TimeDiscretisation -> Strategy
//
//- attribute T of the Model is now optional
//
//Revision 1.32  2004/08/02 09:26:26  jbarbier
//- xml save for SiconosMemory corrected
//- temporary operation in Moreau::integrate because of the current version of
//SiconosVector
//
//Revision 1.31  2004/07/30 14:37:15  jbarbier
//- saving methods for DynamicalSystemXML and LagrangianNLDSXML
//
//Revision 1.30  2004/07/29 14:25:42  jbarbier
//- $Log: DSXML.h,v $
//- Revision 1.48  2005/03/14 16:05:27  jbarbier
//- - manual creation of DSInputOutput saving OK
//-
//- - in progress for EqualityConstraint
//-
//- Revision 1.47  2005/03/09 15:30:36  jbarbier
//- - add of LagrangianEC class
//-
//- - in progress : implementation of the EqualityConstraint and DSInputOutput - create methods
//-
//- Revision 1.46  2005/03/08 12:41:38  jbarbier
//- - constant variables files modified :
//- Some constants added in SiconosConst
//-
//- all global tag of the modeling tools are in XMLTagsName, other tags are specific to an XML class
//-
//- Revision 1.45  2005/03/04 15:35:26  jbarbier
//- - README files added for some samples
//-
//- - beginning of the refactoring of XML module constants
//-
//- Revision 1.44  2005/02/24 15:50:20  jbarbier
//- - LCP prepared to changes needed for several interactions
//-
//- - new function for the SiconosMatrices to copy a block matrix into another matrix
//-
//- - tag names of BoundaryConditionXML, DSInputOutputXML, DSXML, InteractionXML, LagrangianLinearRXML, LagrangianNLDSXML put in XMLTagNames.h
//-
//- Revision 1.43  2005/01/18 17:07:44  charlety
//-
//- _ added autotools makefiles for sample directory
//-
//- Revision 1.42  2005/01/18 10:35:17  jbarbier
//- - attribute "r" no longer used for Moreau integrator
//-
//- - modificatoin in the tests for Moreau integrator
//-
//- - file XMLTagsName.h for further use to regroup all xml tags name...
//-
//- Revision 1.41  2004/12/08 12:49:38  jbarbier
//- - changes in the XML Schema, respect of the recommandations of the W3C
//- version 1.1
//-
//- - changes in all balises DS, Relation, NSLaw, OneStepIntegrator, OneStepNSProblem
//- in the XML files into specific names like LagrangianNLDS, LinearSystemDS, ...
//- for the DS
//-
//- Revision 1.40  2004/09/27 08:24:26  charlety
//-
//- _ Modifications in doxygen comments.
//-
//- Revision 1.39  2004/09/10 11:26:26  charlety
//-
//- _ Integration of the new version of the SiconosVector in the platform. the class simpleVector is used mostly to replace old SiconosVector. When a vector can be composite or simple, like the state of a dynamical system, a pointer on SiconosVector is used, and the vector is initialized simple or composite when the system is initialized.
//-
//- _ All the tests which worked with the previous version of the vector are OK with the new version.
//-
//- _ Example SICONOS and bouncingBall are OK
//-
//- _ some comments have still to be adapted to NewSiconosVector .
//-
//- _ This version of NewSiconosVector could be called 0.9. some details have to be fixed, it will be done before the end of September.
//-
//- Revision 1.38  2004/09/10 08:04:48  jbarbier
//- - XML save available for BoundaryCondition and Interaction
//-
//- Revision 1.37  2004/08/23 14:30:02  jbarbier
//- - All the dynamical systems can be created in a comand program and added to a
//- NSDS. The save is OK, but the creation of the boundary conditions is not yet
//- finished.
//-
//- Revision 1.36  2004/08/20 15:26:45  jbarbier
//- - creation of a Model and save in the XML is ok
//- - creation of a NSDS and save in the XML is ok
//- - creation of a NonLinearSystemDS and save in the XML is OK
//-
//- Revision 1.35  2004/08/10 12:04:30  jbarbier
//- - save of the plugin's name for fInt
//-
//- Revision 1.34  2004/08/04 11:03:23  jbarbier
//- - about the SiconosMemory : when a SiconosMemory has a maxSize greater than the
//- number of steps in memory required by an integrator, the oldest SiconosVector
//- are deleted
//-
//- - the way to initialize the SiconosMemory by the integrator has been updated to
//- match with these changes
//-
//- Revision 1.33  2004/08/03 12:07:12  jbarbier
//- - all test on th eModel are successfull
//-
//- - new tests on the Model with the opening of XML file
//-
//- - link TimeDiscretisation -> Strategy
//-
//- - attribute T of the Model is now optional
//-
//- Revision 1.32  2004/08/02 09:26:26  jbarbier
//- - xml save for SiconosMemory corrected
//- - temporary operation in Moreau::integrate because of the current version of
//- SiconosVector
//-
//- Revision 1.31  2004/07/30 14:37:15  jbarbier
//- - saving methods for DynamicalSystemXML and LagrangianNLDSXML
//- and $Id: DSXML.h,v 1.48 2005/03/14 16:05:27 jbarbier Exp $ added
//
