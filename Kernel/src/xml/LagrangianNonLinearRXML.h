
/** \class LagrangianLinearRXML
 *   \brief This class manages LagrangianLinear Relation data
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 1.0
 *   \date 08/12/2004
 *
 *
 *
 * LagrangianNonLinearRXML allows to manage data of a LNLRelation DOM tree.
 */

#ifndef __LNLRelationXML__
#define __LNLRelationXML__

#include "RelationXML.h"

class LagrangianNonLinearRXML : public RelationXML
{
public:

  LagrangianNonLinearRXML();

  /** \fn LagrangianNonLinearRXML(xmlNode * LNLRelationNode)
   *   \brief Build a LagrangianNonLinearRXML object from a DOM tree describing a Relation with LNL type
   *   \param LagrangianNonLinearRXML : the LagrangianNonLinearR DOM tree
   *   \exception XMLException : if a property of the LagrangianLinear Relation lacks in the DOM tree
   */
  LagrangianNonLinearRXML(xmlNode * LNLRelationNode);

  ~LagrangianNonLinearRXML();

  /** \fn int getComputeInputPlugin()
   *   \brief Returns the computeInput plugin of the Relation
   *   \return string which defines the plugin
   */
  std::string  getComputeInputPlugin() const;

  /** \fn int getComputeOutputPlugin()
   *   \brief Returns the computeOutput plugin of the Relation
   *   \return string which defines the plugin
   */
  std::string  getComputeOutputPlugin() const;

private:
  //Nodes
};


#endif
