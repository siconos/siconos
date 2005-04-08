
/** \class LagrangianLinearECXML
*   \brief This class manages LagrangianLinear Relation data
*  \author SICONOS Development Team - copyright INRIA
*   \version 1.0
*   \date 05/25/2004
*
*
*
* LagrangianLinearECXML allows to manage data of a LagrangianLinearEC DOM tree.
*/


#ifndef _LagrangianLinearECXML_
#define _LagrangianLinearECXML_


#include <libxml/tree.h>

#include "EqualityConstraintXML.h"
#include "LagrangianLinearECXML.h"
#include "SiconosMatrix.h"
#include "NewSiconosVector.h"
#include "SimpleVector.h"
#include "SiconosDOMTreeTools.h"


using namespace std;


const string LLEC_H = "H";
const string LLEC_B = "b";
//#include "XMLTagsName.h"


class LagrangianLinearECXML : public EqualityConstraintXML
{
public:

  LagrangianLinearECXML();

  /** \fn LagrangianLinearECXML(xmlNode * , vector<int> )
  *   \brief Build a EqualityConstraintXML object from a DOM tree describing a EqualityConstraint
  *   \param xmlNode* : the EqualityConstraint DOM tree
  *   \param vector<int>  : vector of DSXML numbers to verify DS concerned by the EqualityConstraint (identified by number) exists
  */
  LagrangianLinearECXML(xmlNode*, vector<int>);

  ~LagrangianLinearECXML();

  /** \fn SiconosMatrix getH()
  *   \brief Return the H of the LagrangianLinearECXML
  *   \return The H SiconosMatrix of the LagrangianLinearECXML
  */
  inline SiconosMatrix getH()
  {
    return  SiconosDOMTreeTools::getSiconosMatrixValue(this->HNode);
  }


  /** \fn SimpleVector getB()
  *   \brief Return b vector of the LLRelationXML
  *   \return SimpleVector : b vector of the LLRelationXML
  */
  inline /*SiconosVector*/SimpleVector getB()
  {
    return  SiconosDOMTreeTools::getSiconosVectorValue(this->bNode);
  }

  /** \fn void setH(SiconosMatrix *matrix)
  *   \brief Change the H matrix value (in xml file or external data file switch his origin position)
  *   \param SiconosMatrix matrix : the new value for H matrix
  */
  void setH(SiconosMatrix *matrix);

  /** \fn void setB(SiconosVector *vector)
  *   \brief Change the b vector value (in xml file or external data file switch his origin position)
  *   \param SiconosVector vector : the new value for b vector
  */
  void setB(SiconosVector *vector);


private:

  //Nodes
  xmlNode * HNode;
  xmlNode * bNode;

};

#endif
