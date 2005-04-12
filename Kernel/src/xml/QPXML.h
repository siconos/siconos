
/** \class QPXML
*   \brief This class manages Lagrangian QP data
*  \author SICONOS Development Team - copyright INRIA
*   \version 1.0
*   \date 05/18/2004
*
*
*
* QPXML allows to manage data of a QP DOM tree.
*/


#ifndef __QPXMLDEF__
#define __QPXMLDEF__



#include <string>
#include <libxml/tree.h>

#include "OneStepNSProblemXML.h"
#include "OneStepNSProblem.h"

#include "SiconosDOMTreeTools.h"
//#include "SiconosVector.h"
#include "NewSiconosVector.h"
#include "SimpleVector.h"
#include "SiconosMatrix.h"

//using namespace std;

class OneStepNSProblem;

const string QP_Q = "Q";
const string QP_P = "p";


class QPXML : public OneStepNSProblemXML
{
public:
  QPXML();

  /** \fn QPXML(xmlNode * QPNode)
  *   \brief Build a QPXML object from a DOM tree describing a QP
  *   \param QPNode : the QP DOM tree
  *   \param vector<int> definedInteractionNumbers : the Interaction numbers effectivly defined in the model
  *   \exception XMLException : if a property of the QP lacks in the DOM tree
  */
  QPXML(xmlNode * QPNode, vector<int> definedInteractionNumbers);

  /** \fn SiconosMatrix getQ()
  *   \brief Return Q
  *   \return The Q SiconosMatrix of the QP
  */
  inline SiconosMatrix getQ()
  {
    return  SiconosDOMTreeTools::getSiconosMatrixValue(this->QNode);
  }

  /** \fn SimpleVector getP()
  *   \brief Return p
  *   \return SimpleVector :  vector p of the QP
  */
  inline /*SiconosVector*/SimpleVector getP()
  {
    return SiconosDOMTreeTools::getSiconosVectorValue(this->pNode);
  }

  /** \fn void setQ(SiconosMatrix *m)
  *   \brief allows  to save Q
  *   \param The Q SiconosMatrix to save
  */
  inline void setQ(SiconosMatrix *m)
  {
    if (this->hasQ() == false)
    {
      this->QNode = SiconosDOMTreeTools::createMatrixNode(this->rootNSProblemXMLNode, QP_Q, m);
    }
    else SiconosDOMTreeTools::setSiconosMatrixValue(this->QNode, m);
  }

  /** \fn void setP(SiconosVector *v)
  *   \brief allows to save p
  *   \param SimpleVector* : vector p to save
  */
  inline void setP(SiconosVector *v)
  {
    if (this->hasP() == false)
    {
      this->pNode = SiconosDOMTreeTools::createVectorNode(this->rootNSProblemXMLNode, QP_P, v);
    }
    else SiconosDOMTreeTools::setSiconosVectorValue(this->pNode, v);
  }

  /** \fn bool hasP()
   *  \brief returns true if pNode is defined
   *  \return true if pNode is defined
   */
  inline bool hasP()
  {
    return (this->pNode != NULL);
  }

  /** \fn bool hasQ()
   *  \brief returns true if QNode is defined
   *  \return true if QNode is defined
   */
  inline bool hasQ()
  {
    return (this->QNode != NULL);
  }

  /** \fn void updateOneStepNSProblemXML( xmlNode* node, OneStepNSProblemXML* str )
  *   \brief makes the operations to create a OneStepNSProblemXML to the StrategyXML
  *   \param xmlNode* : the root node of the OneStepNSProblemXML
  *   \param OneStepNSProblem* : the OneStepNSProblem of this OneStepNSProblemXML
  */
  void updateOneStepNSProblemXML(xmlNode* node, OneStepNSProblem* osnspb);


private:
  //Nodes
  xmlNode * QNode;
  xmlNode * pNode;

};

#endif
