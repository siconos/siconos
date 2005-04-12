//$IId$

/** \class TimeDiscretisationXML
*   \brief This class manages Time Discretisation data
*  \author SICONOS Development Team - copyright INRIA
*   \version 1.0
*   \date 05/17/2004
*
*
*
* TimeDiscretisationXML allows to manage data of a TimeDiscretisation DOM tree.
*/


#ifndef __TIMEDISCRETISATIONXML__
#define __TIMEDISCRETISATIONXML__


#include <libxml/tree.h>
//#include "SiconosVector.h"
#include "NewSiconosVector.h"
#include "SimpleVector.h"
#include "SiconosDOMTreeTools.h"

#include "TimeDiscretisation.h"



//using namespace std;

class TimeDiscretisation;


const string TD_H = "h";
const string TD_N = "N";
const string TD_TK = "tk";
const string TD_HMIN = "hMin";
const string TD_HMAX = "hMax";
const string TD_ISCONSTANT = "isConstant";


class TimeDiscretisationXML
{
public:
  TimeDiscretisationXML();

  /** \fn TimeDiscretisationXML(xmlNode * TimeDiscretisationNode)
  *   \brief Build a TimeDiscretisationXML object from a DOM tree describing a TimeDiscretisation
  *   \param timeDiscretisationNode : the TimeDiscretisation DOM tree
  *   \exception XMLException : if a property of the TimeDiscretisation lacks in the DOM tree
  */
  TimeDiscretisationXML(xmlNode * TimeDiscretisationNode);


  /** \fn xmlNode* getRootNode()
  *   \brief Gets the rootNode of the TimeDiscretisationXML
  *   \return xmlNode* : the rootNode
  */
  inline xmlNode* getRootNode()
  {
    return this->rootNode;
  }

  /** \fn void setRootNode(xmlNode*)
  *   \brief sets the rootNode of the TimeDiscretisationXML
  *   \param xmlNode* : the rootNode to set
  */
  inline void setRootNode(xmlNode*node)
  {
    this->rootNode = node;
  }

  /** \fn inline double getH()
  *   \brief Return the h value of the TimeDiscretisation
  *   \return The h double value of the TimeDiscretisation
  */
  inline double getH()
  {
    return  SiconosDOMTreeTools::getDoubleContentValue(this->hNode);
  }

  /** \fn void setH(double d)
  *   \brief allows to set the h value of the TimeDiscretisation
  *   \param The h double value of the TimeDiscretisation
  */
  inline void setH(double d)
  {
    if (this->hasH() == false)
    {
      this->hNode = SiconosDOMTreeTools::createDoubleNode(this->rootNode, TD_H, d);
    }
    else SiconosDOMTreeTools::setDoubleContentValue(this->hNode, d);
  }

  /** \fn inline void getN(int i)
  *   \brief return the N value of the TimeDiscretisation
  *   \param The N int value of the TimeDiscretisation
  */
  inline int getN()
  {
    return SiconosDOMTreeTools::getIntegerContentValue(this->NNode);
  }

  /** \fn inline void setN(int i)
  *   \brief allows to set the N value of the TimeDiscretisation
  *   \param The N int value of the TimeDiscretisation
  */
  inline void setN(int i)
  {
    if (this->hasN() == false)
    {
      this->NNode = SiconosDOMTreeTools::createIntegerNode(this->rootNode, TD_N, i);
    }
    else SiconosDOMTreeTools::setIntegerContentValue(this->NNode, i);
  }

  /** \fn inline SimpleVector getTk()
  *   \brief Return the tk values of the TimeDiscretisation
  *   \return SimpleVector : tk vector of the TimeDiscretisation
  */
  inline SimpleVector* getTk()
  {
    /*SiconosVector*/SimpleVector *v = new
    SimpleVector(SiconosDOMTreeTools::getSiconosVectorValue(this->tkNode));
    return  v; //SiconosDOMTreeTools::getSiconosVectorValue(this->tkNode);
  }

  /** \fn void setTk(SiconosVector *v)
  *   \brief allows to set the tk values of the TimeDiscretisation
  *   \param The tk SiconosVector of the TimeDiscretisation
  */
  inline void setTk(SiconosVector *v)
  {
    //SiconosDOMTreeTools::setSiconosVectorValue(this->tkNode, *v);
    if (this->hasTk() == false)
    {
      this->tkNode = SiconosDOMTreeTools::createVectorNode(this->rootNode, TD_TK, v);
    }
    else SiconosDOMTreeTools::setSiconosVectorValue(this->tkNode, v);
  }

  /** \fn inline double getHMin()
  *   \brief Return the hMin value of the TimeDiscretisation / -1.0 if not defined
  *   \return The hMin double value of the TimeDiscretisation
  */
  inline double getHMin()
  {
    if (hMin)
      return SiconosDOMTreeTools::getDoubleContentValue(this->hMinNode);
    return -1.0;
  }

  /** \fn void setHMin(double d)
  *   \brief allows to set the hMin value of the TimeDiscretisation
  *   \param The hMin double value of the TimeDiscretisation
  */
  inline void setHMin(double d)
  {
    if (this->hasHMin() == false)
    {
      this->hMinNode = SiconosDOMTreeTools::createDoubleNode(this->rootNode, TD_HMIN, d);
    }
    else SiconosDOMTreeTools::setDoubleContentValue(this->hMinNode, d);
  }

  /** \fn inline double getHMax()
  *   \brief Return the hMax value of the TimeDiscretisation / -1.0 if not defined
  *   \return The hMax double value of the TimeDiscretisation
  */
  inline double getHMax()
  {
    if (hMax)
      return SiconosDOMTreeTools::getDoubleContentValue(this->hMaxNode);
    return -1.0;
  }

  /** \fn void setHMax(double d)
  *   \brief allows to set the hMax value of the TimeDiscretisation
  *   \param The hMax double value of the TimeDiscretisation
  */
  inline void setHMax(double d)
  {
    if (this->hasHMax() == false)
    {
      this->hMaxNode = SiconosDOMTreeTools::createDoubleNode(this->rootNode, TD_HMAX, d);
    }
    else SiconosDOMTreeTools::setDoubleContentValue(this->hMaxNode, d);
  }


  /** \fn inline void setConstant(bool)
  *   \brief defines if the TimeDiscretisation is constant or not
  *   \param bool : true if TimeDiscretisation is constant, false otherwise
  */
  inline void setConstant(bool b)
  {
    if (SiconosDOMTreeTools::hasAttributeValue(this->rootNode, TD_ISCONSTANT))
    {
      SiconosDOMTreeTools::setBooleanAttributeValue(this->rootNode, TD_ISCONSTANT, b);
    }
    else
    {
      SiconosDOMTreeTools::createBooleanAttribute(this->rootNode, TD_ISCONSTANT, b);
    }
  }

  /** \fn inline bool isConstant()
  *   \brief Return true if the TimeDiscretisation is constant
  *   \return A boolean value : true if TimeDiscretisation is constant, false otherwise
  */
  inline bool isConstant()
  {
    return SiconosDOMTreeTools::getBooleanAttributeValue(this->rootNode, TD_ISCONSTANT);
  }

  //    /** \fn void setHConstant(bool b)
  //    *   \brief allows to set the value of constant
  //    *   \return A boolean value : true if h is constant, false otherwise
  //    */
  //    inline void setHConstant(bool b)
  //    {
  //      SiconosDOMTreeTools::setBooleanAttributeValue(this->hNode, TD_ISCONSTANT, b);
  //    }

  /** \fn bool hasH()
   *  \brief returns true if hNode is defined
   *  \return true if hNode is defined
   */
  inline bool hasH()
  {
    return (this->hNode != NULL);
  }

  /** \fn bool hasN()
   *  \brief returns true if NNode is defined
   *  \return true if NNode is defined
   */
  inline bool hasN()
  {
    return (this->NNode != NULL);
  }

  /** \fn bool hasTk()
   *  \brief returns true if tkNode is defined
   *  \return true if tkNode is defined
   */
  inline bool hasTk()
  {
    return (this->tkNode != NULL);
  }

  /** \fn bool hasHMin()
   *  \brief returns true if hMinNode is defined
   *  \return true if hMinNode is defined
   */
  inline bool hasHMin()
  {
    return (this->hMinNode != NULL);
  }

  /** \fn bool hasHMax()
   *  \brief returns true if hMaxNode is defined
   *  \return true if hMaxNode is defined
   */
  inline bool hasHMax()
  {
    return (this->hMaxNode != NULL);
  }

  /** \fn void updateTimeDiscretisationXML(xmlNode*, TimeDiscretisation*)
  *   \brief makes the operations to create the TimeDiscretisation of the StrategyXML
  *   \param xmlNode* : the root node for the TimeDiscretisationXML
  *   \param NSDS* : the TimeDiscretisation of this TimeDiscretisationXML
  */
  void updateTimeDiscretisationXML(xmlNode*, TimeDiscretisation*);

private:

  //Nodes
  xmlNode * rootNode;
  xmlNode * hNode;
  xmlNode * NNode;
  xmlNode * tkNode;
  xmlNode * hMinNode;
  xmlNode * hMaxNode;

  bool hMin;
  bool hMax;

  //Methods

  /** \fn void loadTimeDiscretisationProperties(xmlNode * timeDiscretiationNode)
  *   \brief load the different properties of a TimeDiscretisation
  *   \param xmlNode * timeDiscretisationNode : the DOM tree node of the concern TimeDiscretisation
  *   \exception XMLException
  */
  void loadTimeDiscretisationProperties(xmlNode * timeDiscretisationNode);
};


#endif
