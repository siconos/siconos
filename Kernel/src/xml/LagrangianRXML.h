
/** \class LagrangianRXML
 *   \brief This class manages Lagrangian Relation data
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 1.0
 *   \date 08/12/2004
 *
 *
 *
 * LagrangianRXML allows to manage data of a LNLRelation DOM tree.
 */

#ifndef __LAGRANGIANRelationXML__
#define __LAGRANGIANRelationXML__

#include "RelationXML.h"

class LagrangianRXML : public RelationXML
{
public:

  LagrangianRXML();

  /** \fn LagrangianRXML(xmlNode * LNLRelationNode)
   *   \brief Build a LagrangianRXML object from a DOM tree describing a Relation with LNL type
   *   \param LagrangianRXML : the LagrangianR DOM tree
   *   \exception XMLException : if a property of the LagrangianLinear Relation lacks in the DOM tree
   */
  LagrangianRXML(xmlNode * LNLRelationNode);

  ~LagrangianRXML();

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
