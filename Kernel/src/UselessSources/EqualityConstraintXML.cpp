/* Siconos-Kernel version 2.1.1, Copyright INRIA 2005-2006.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
*/
#include "EqualityConstraintXML.h"
using namespace std;

EqualityConstraintXML::EqualityConstraintXML()
{
  this->GNode = NULL;
  this->computeInputNode = NULL;
  this->computeOutputNode = NULL;
}

EqualityConstraintXML::EqualityConstraintXML(xmlNodePtr ecNode, vector<int> definedDSNumbers)
{
  xmlNodePtr node;
  //  string type ( (char*)ecNode->name );
  this->rootNode = ecNode;

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootNode, EQUALITYCONSTRAINT_G)) != NULL)
  {
    this->GNode = node;
  }
  else
  {
    XMLException::selfThrow("EqualityConstraintXML - EqualityConstraintXML(xmlNodePtr ecNode) error : tag " + EQUALITYCONSTRAINT_G + " not found.");
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootNode, COMPUTE_INPUT_TAG)) != NULL)
    this->computeInputNode = node;
  else
  {
    //      //if( SiconosDOMTreeTools::getStringAttributeValue(this->rootRelationXMLNode, RELATION_TYPE) == RELATION_LNL )
    //      if( LAGRANGIAN_NON_LINEAR_RELATION_TAG == type )
    //      {
    //      XMLException::selfThrow("RelationXML - RelationXML::RelationXML(xmlNodePtr relationNode) error : tag " + COMPUTE_INPUT_TAG + " not found.");
    //      }
    /*else*/ cout << "EqualityConstraintXML - EqualityConstraintXML::EqualityConstraintXML - Warning : tag " << COMPUTE_INPUT_TAG << " not found. This is attribute is optional for this relation" << endl;
    this->computeInputNode = NULL;
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootNode, COMPUTE_OUTPUT_TAG)) != NULL)
  {
    this->computeOutputNode = node;
  }
  else
  {
    //      //if( SiconosDOMTreeTools::getStringAttributeValue(this->rootRelationXMLNode, RELATION_TYPE) == RELATION_LNL )
    //      if( LAGRANGIAN_NON_LINEAR_RELATION_TAG == type )
    //      {
    //      XMLException::selfThrow("RelationXML - RelationXML::RelationXML(xmlNodePtr relationNode) error : tag " + COMPUTE_OUTPUT_TAG + " not found.");
    //      }
    /*else*/ cout << "EqualityConstraintXML - EqualityConstraintXML::EqualityConstraintXML - Warning : tag " << COMPUTE_OUTPUT_TAG << " not found. This is attribute is optional for this relation" << endl;
    this->computeOutputNode = NULL;
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootNode, EQUALITYCONSTRAINT_DSIO_CONCERNED)) != NULL)
  {
    this->dsioConcernedNode = node;
    loadECConcernedDSIO(node/*, definedDSNumbers*/);
  }
  else
  {
    XMLException::selfThrow("EqualityConstraintXML - EqualityConstraintXML(xmlNodePtr ecNode) error : tag " + EQUALITYCONSTRAINT_DSIO_CONCERNED + " not found.");
  }
}

EqualityConstraintXML::~EqualityConstraintXML()
{}


void EqualityConstraintXML::loadECConcernedDSIO(xmlNodePtr  DSIOConcernedNode/*, vector<int> definedDSNumbers*/)
{
  xmlNodePtr DSIOnode;
  int number;
  int size = 0;
  int i = 0;

  if ((DSIOnode = SiconosDOMTreeTools::findNodeChild((const xmlNodePtr)DSIOConcernedNode, DSINPUTOUTPUT_TAG)) == NULL)
  {
    XMLException::selfThrow("EqualityConstraintXML - loadECConcernedDSIO error : at least one couple of " + DSINPUTOUTPUT_TAG + " must be declared in " + EQUALITYCONSTRAINT_DSIO_CONCERNED + " tag.");
  }

  size = SiconosDOMTreeTools::getNodeChildrenNumber(DSIOConcernedNode);
  while ((DSIOnode != NULL) && (i < size))
  {
    number = SiconosDOMTreeTools::getAttributeValue<int>(DSIOnode, NUMBER_ATTRIBUTE);

    // \todo : verifying that the DS are defined before
    //    j = 0;
    //    //Verifying that the DS number exists
    //    while ((j<definedDSNumbers.size()) && (definedDSNumbers[j]!=number)) {  j++; }
    //
    //    if (j==definedDSNumbers.size())
    //    {
    //      char errorMsg[1024];
    //      sprintf(errorMsg, "DSInputOutputXML - loadDSInputOutputConcernedDS error : in a tag %s you define couple of DS with a DS number who doesn't exist : %d.", INTERACTION_DS_CONCERNED.c_str(), number1);
    //      XMLException::selfThrow(errorMsg);
    //    }

    this->definedDSIONumbers.push_back(number);

    DSIOnode = SiconosDOMTreeTools::findFollowNode(DSIOnode, DSINPUTOUTPUT_TAG);

    i++;
  }
}

void EqualityConstraintXML::setDSIOConcerned(vector<int> dsioConcerned)
{
  if (this->dsioConcernedNode == NULL)
  {
    xmlNodePtr  node;

    //    // conversion of the vector<DynamicalSystem*> in vector< vector<int> >
    //    this->dsVector.clear();
    //    for( i = 0; i<dsConcerned.size(); i=i+2 )
    //    {
    //      number = dsConcerned[i]->getNumber();
    //      this->dsVector.push_back( number );
    //    }

    //save in the DOM tree
    char num[32];

    /*
     * creation of the DS_Concerned node
     */
    this->dsioConcernedNode = xmlNewChild(this->rootNode, NULL, (xmlChar*)EQUALITYCONSTRAINT_DSIO_CONCERNED.c_str(), NULL);

    for (unsigned int i = 0; i < this->definedDSIONumbers.size(); i++)
    {
      node = xmlNewChild(this->dsioConcernedNode, NULL, (xmlChar*)DSINPUTOUTPUT_TAG.c_str(), NULL);
      sprintf(num, "%i", this->definedDSIONumbers[i]);
      xmlNewProp(node, (xmlChar*)NUMBER_ATTRIBUTE.c_str(), (xmlChar*)num);
    }
  }
  else
  {
    /* \todo : when EqualityConstraint have been given in the XML input file and that the user has added new ones */
  }
}

void EqualityConstraintXML::updateEqualityConstraintXML(xmlNodePtr  node, EqualityConstraint* ec)
{
  rootNode = node;
}
