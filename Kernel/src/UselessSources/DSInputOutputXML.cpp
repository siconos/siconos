/* Siconos-Kernel version 2.1.0, Copyright INRIA 2005-2006.
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
#include "DSInputOutputXML.h"
using namespace std;

DSInputOutputXML::DSInputOutputXML()
{
  this->computeInputNode = NULL;
  this->computeOutputNode = NULL;
  //  this->HNode = NULL;
}

DSInputOutputXML::DSInputOutputXML(xmlNodePtr dsioNode/*, vector<int> definedDSNumbers*/)
{
  xmlNodePtr node;
  string type((char*)dsioNode->name);
  this->rootDSIOXMLNode = dsioNode;

  if (type == NON_LINEAR_DSIO_TAG || type == LAGRANGIAN_DSIO_TAG)
  {
    if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSIOXMLNode, COMPUTE_INPUT_TAG)) != NULL)
      this->computeInputNode = node;
    else
      XMLException::selfThrow("DSInputOutputXML - DSInputOutputXML(xmlNodePtr dsioNode) error : tag " + COMPUTE_INPUT_TAG + " not found.");

    if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSIOXMLNode, COMPUTE_OUTPUT_TAG)) != NULL)
      this->computeOutputNode = node;
    else
      XMLException::selfThrow("DSInputOutputXML - DSInputOutputXML(xmlNodePtr dsioNode) error : tag " + COMPUTE_OUTPUT_TAG + " not found.");
  }
  else
  {
    this->computeInputNode = NULL;
    this->computeOutputNode = NULL;
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootDSIOXMLNode, DS_CONCERNED)) != NULL)
  {
    this->dsConcernedNode = node;
    loadDSIOConcernedDS(node/*, definedDSNumbers*/);
  }
  else
  {
    XMLException::selfThrow("DSInputOutputXML - DSInputOutputXML(xmlNodePtr dsioNode) error : tag " + DS_CONCERNED + " not found.");
  }
}

DSInputOutputXML::~DSInputOutputXML()
{}


void DSInputOutputXML::loadDSIOConcernedDS(xmlNodePtr  DSConcernedNode/*, vector<int> definedDSNumbers*/)
{
  xmlNodePtr DSnode;
  int number;
  int size = 0;
  int i = 0;

  if ((DSnode = SiconosDOMTreeTools::findNodeChild((const xmlNodePtr)DSConcernedNode, DYNAMICAL_SYSTEM_TAG)) == NULL)
  {
    XMLException::selfThrow("DSInputOutputXML - loadDSIOConcernedDS error : at least one couple of " + DYNAMICAL_SYSTEM_TAG + " must be declared in " + DS_CONCERNED + " tag.");
  }

  size = SiconosDOMTreeTools::getNodeChildrenNumber(DSConcernedNode);
  while ((DSnode != NULL) && (i < size))
  {
    number = SiconosDOMTreeTools::getAttributeValue<int>(DSnode, NUMBER_ATTRIBUTE);

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

    this->definedDSNumbers.push_back(number);

    DSnode = SiconosDOMTreeTools::findFollowNode(DSnode, DYNAMICAL_SYSTEM_TAG);

    i++;
  }
}

void DSInputOutputXML::updateDSInputOutputXML(xmlNodePtr  node, DSInputOutput* dsio)
{
  this->rootDSIOXMLNode = node;
}

void DSInputOutputXML::setDSConcerned(vector<int> dsConcerned)
{
  if (this->dsConcernedNode == NULL)
  {
    xmlNodePtr node;

    //save in the DOM tree
    char num[32];

    /*
     * creation of the DS_Concerned node
     */
    this->dsConcernedNode = xmlNewChild(this->rootDSIOXMLNode, NULL, (xmlChar*)DS_CONCERNED.c_str(), NULL);

    for (unsigned int i = 0; i < this->definedDSNumbers.size(); i++)
    {
      node = xmlNewChild(this->dsConcernedNode, NULL, (xmlChar*)DYNAMICAL_SYSTEM_TAG.c_str(), NULL);
      sprintf(num, "%i", this->definedDSNumbers[i]);
      xmlNewProp(node, (xmlChar*)NUMBER_ATTRIBUTE.c_str(), (xmlChar*)num);
    }
  }
  else
  {
    /* \todo : when DSIO have been given in the XML input file and that the user has added new ones */
  }
}

