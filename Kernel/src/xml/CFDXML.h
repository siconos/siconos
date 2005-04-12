/** \class CFDXML
*   \brief This class manages Lagrangian CFD data
*  \author SICONOS Development Team - copyright INRIA
*   \version 1.0
*   \date 12/04/2005
*
*
*
* CFDXML allows to manage data of a CFD DOM tree.
*/


#ifndef CFDXML_H
#define CFDXML_H



#include <string>
#include <libxml/tree.h>

#include "OneStepNSProblemXML.h"
#include "OneStepNSProblem.h"

#include "SiconosDOMTreeTools.h"
#include "NewSiconosVector.h"
#include "SimpleVector.h"
#include "SiconosMatrix.h"


class OneStepNSProblem;

const string CFD_M = "M";
const string CFD_Q = "q";


class CFDXML : public OneStepNSProblemXML
{
public:
  CFDXML();

  /** \fn CFDXML(xmlNode * CFDNode)
  *   \brief Build a CFDXML object from a DOM tree describing a CFD
  *   \param CFDNode : the CFD DOM tree
  *   \param vector<int> definedInteractionNumbers : the Interaction numbers effectivly defined in the model
  *   \exception XMLException : if a property of the CFD lacks in the DOM tree
  */
  CFDXML(xmlNode * CFDNode, vector<int> definedInteractionNumbers);

  /** \fn SiconosMatrix getM()
  *   \brief Return M
  *   \return The M SiconosMatrix of the CFD
  */
  inline SiconosMatrix getM()
  {
    return  SiconosDOMTreeTools::getSiconosMatrixValue(this->MNode);
  }

  /** \fn SimpleVector getQ()
  *   \brief Return vector q
  *   \return SimpleVector : q vector of the CFD
  */
  inline /*SiconosVector*/SimpleVector getQ()
  {
    return SiconosDOMTreeTools::getSiconosVectorValue(this->qNode);
  }

  /** \fn void setM(SiconosMatrix *m)
  *   \brief allows to save M
  *   \return The M SiconosMatrix to save
  */
  inline void setM(SiconosMatrix *m)
  {
    if (this->hasM() == false)
    {
      this->MNode = SiconosDOMTreeTools::createMatrixNode(this->rootNSProblemXMLNode, CFD_M, m);
    }
    else SiconosDOMTreeTools::setSiconosMatrixValue(this->MNode, m);
  }

  /** \fn void setQ(SiconosVector *v)
  *   \brief allows to save q
  *   \return The q SiconosVector to save
  */
  inline void setQ(SiconosVector *v)
  {
    if (this->hasQ() == false)
    {
      this->qNode = SiconosDOMTreeTools::createVectorNode(this->rootNSProblemXMLNode, CFD_Q, v);
    }
    else SiconosDOMTreeTools::setSiconosVectorValue(this->qNode, v);
  }

  /** \fn bool hasM()
   *  \brief returns true if MNode is defined
   *  \return true if MNode is defined
   */
  inline bool hasM()
  {
    return (this->MNode != NULL);
  }

  /** \fn bool hasQ()
   *  \brief returns true if qNode is defined
   *  \return true if qNode is defined
   */
  inline bool hasQ()
  {
    return (this->qNode != NULL);
  }

  /** \fn void updateOneStepNSProblemXML( xmlNode* node, OneStepNSProblemXML* str )
  *   \brief makes the operations to create a OneStepNSProblemXML to the StrategyXML
  *   \param xmlNode* : the root node of the OneStepNSProblemXML
  *   \param OneStepNSProblem* : the OneStepNSProblem of this OneStepNSProblemXML
  */
  void updateOneStepNSProblemXML(xmlNode* node, OneStepNSProblem* osnspb);


private:

  //Nodes
  xmlNode * MNode;
  xmlNode * qNode;
};


#endif
