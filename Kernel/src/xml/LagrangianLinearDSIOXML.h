#ifndef LAGRANGIANLINEARDSIOXML_H
#define LAGRANGIANLINEARDSIOXML_H

#include "DSInputOutputXML.h"
#include <libxml/tree.h>
#include "SiconosMatrix.h"
#include "NewSiconosVector.h"
#include "SimpleVector.h"
#include "SiconosDOMTreeTools.h"

using namespace std;

const string LLDSIO_H = "H";
const string LLDSIO_B = "b";

class LagrangianLinearDSIOXML : public DSInputOutputXML
{
public:

  LagrangianLinearDSIOXML();
  virtual ~LagrangianLinearDSIOXML();

  /** \fn LagrangianLinearDSIOXML(xmlNode * dsioNode)
  *   \brief Build a LagrangianLinearDSIOXML object from a DOM tree describing a DSIO with LagrangianLinear type
  *   \param LagrangianLinear : the LagrangianLinearDSIO DOM tree
  *   \exception XMLException : if a property of the LagrangianLinear DSIO lacks in the DOM tree
  */
  LagrangianLinearDSIOXML(xmlNode * dsioNode);

  //////////////////////////////////////////////

  /** \fn SiconosMatrix getH()
  *   \brief Return the H of the LagrangianLinearDSIOXML
  *   \return The H SiconosMatrix of the LagrangianLinearDSIOXML
  */
  inline SiconosMatrix getH()
  {
    return  SiconosDOMTreeTools::getSiconosMatrixValue(this->HNode);
  }


  /** \fn SimpleVector getB()
  *   \brief Return b vector of the LagrangianLinearDSIOXML
  *   \return SimpleVector : b vector of the LagrangianLinearDSIOXML
  */
  inline SimpleVector getB()
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

#endif // LAGRANGIANLINEARDSIOXML_H
