/** \class DFC_2DXML
*   \brief This class manages Lagrangian DFC_2D data
*  \author SICONOS Development Team - copyright INRIA
*   \version 1.0
*   \date 12/04/2005
*
*
*
* DFC_2DXML allows to manage data of a DFC_2D DOM tree.
*/


#ifndef DFC_2DXML_H
#define DFC_2DXML_H

#include "OneStepNSProblemXML.h"

class OneStepNSProblem;

const std::string DFC_2D_M = "M";
const std::string DFC_2D_Q = "q";


class DFC_2DXML : public OneStepNSProblemXML
{
public:
  DFC_2DXML();

  /** \fn DFC_2DXML(xmlNode * DFC_2DNode)
  *   \brief Build a DFC_2DXML object from a DOM tree describing a DFC_2D
  *   \param DFC_2DNode : the DFC_2D DOM tree
  *   \param vector<int> definedInteractionNumbers : the Interaction numbers effectivly defined in the model
  *   \exception XMLException : if a property of the DFC_2D lacks in the DOM tree
  */
  DFC_2DXML(xmlNode * DFC_2DNode, std::vector<int> definedInteractionNumbers);

  /** \fn SiconosMatrix getM()
  *   \brief Return M
  *   \return The M SiconosMatrix of the DFC_2D
  */
  inline SiconosMatrix getM()
  {
    return  SiconosDOMTreeTools::getSiconosMatrixValue(this->MNode);
  }

  /** \fn SimpleVector getQ()
  *   \brief Return vector q
  *   \return SimpleVector : q vector of the DFC_2D
  */
  inline /*SiconosVector*/SimpleVector getQ()
  {
    return SiconosDOMTreeTools::getSiconosVectorValue(this->qNode);
  }

  /** \fn void setM(const SiconosMatrix &m)
  *   \brief save M
  *   \return The M SiconosMatrix to save
  */
  inline void setM(const SiconosMatrix& m)
  {
    if (this->hasM() == false)
    {
      this->MNode = SiconosDOMTreeTools::createMatrixNode(this->rootNSProblemXMLNode, DFC_2D_M, m);
    }
    else SiconosDOMTreeTools::setSiconosMatrixNodeValue(this->MNode, m);
  }

  /** \fn void setQ(const SiconosVector& v)
  *   \brief allows to save q
  *   \return The q SiconosVector to save
  */
  inline void setQ(const SiconosVector&v)
  {
    if (this->hasQ() == false)
    {
      this->qNode = SiconosDOMTreeTools::createVectorNode(this->rootNSProblemXMLNode, DFC_2D_Q, v);
    }
    else SiconosDOMTreeTools::setSiconosVectorNodeValue(this->qNode, v);
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
