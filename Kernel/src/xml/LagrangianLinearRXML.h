
/** \class LagrangianLinearRXML
 *   \brief This class manages LagrangianLinear Relation data
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 1.0
 *   \date 05/25/2004
 *
 *
 *
 * LagrangianLinearRXML allows to manage data of a LLRelation DOM tree.
 */

#ifndef __LLRelationXML__
#define __LLRelationXML__

#include "RelationXML.h"

const std::string  LLR_H = "H";
const std::string  LLR_B = "b";

class SiconosMatrix;
class SiconosVector;
class SimpleVector;

class LagrangianLinearRXML : public RelationXML
{
public:

  LagrangianLinearRXML();

  /** \fn LagrangianLinearRXML(xmlNode * LLRelationNode)
   *   \brief Build a LagrangianLinearRXML object from a DOM tree describing a Relation with LL type
   *   \param LagrangianLinearRXML : the LagrangianLinearR DOM tree
   *   \exception XMLException : if a property of the LagrangianLinear Relation lacks in the DOM tree
   */
  LagrangianLinearRXML(xmlNode * LLRelationNode);

  ~LagrangianLinearRXML();

  /** \fn SiconosMatrix getH()
   *   \brief Return the H of the LLRelationXML
   *   \return The H SiconosMatrix of the LLRelationXML
   */
  inline SiconosMatrix getH()
  {
    return  SiconosDOMTreeTools::getSiconosMatrixValue(HNode);
  }


  /** \fn SimpleVector getB()
   *   \brief Return b vector of the LLRelationXML
   *   \return SimpleVector : b vector of the LLRelationXML
   */
  inline /*SiconosVector*/SimpleVector getB()
  {
    return  SiconosDOMTreeTools::getSiconosVectorValue(bNode);
  }

  /** \fn void setH(SiconosMatrix *matrix)
   *   \brief Change the H matrix value (in xml file or external data file switch his origin position)
   *   \param SiconosMatrix matrix : the new value for H matrix
   */
  void setH(const SiconosMatrix&);

  /** \fn void setB(SiconosVector *vector)
   *   \brief Change the b vector value (in xml file or external data file switch his origin position)
   *   \param SiconosVector vector : the new value for b vector
   */
  void setB(const SiconosVector&);


  /** \fn bool hasB() const
   *   \brief return true if b is given in xmlfile
   */
  inline bool hasB() const
  {
    return (!(bNode == NULL));
  }

private:


  //Nodes
  xmlNode * HNode;
  xmlNode * bNode;

};


#endif
