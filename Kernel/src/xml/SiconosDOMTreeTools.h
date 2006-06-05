/* Siconos-Kernel version 1.1.4, Copyright INRIA 2005-2006.
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
/** \class SiconosDOMTreeTools
 *   \brief This class is a sort of tools box to get and set elements from nodes DOM tree
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 1.1.4.
 *   \date 04/06/2004
 *
 *
 * SiconosDOMTreeTools allows to manage data of Siconos nodes DOM tree, like set or get double, vector etc.
 */

#ifndef __SICONOSDOMTREETOOLS__
#define __SICONOSDOMTREETOOLS__

#include "SiconosConst.h"
#include "SimpleVector.h"
#include "SimpleMatrix.h"

#include "XMLException.h"
#include "XMLTagsName.h"
#include "RuntimeException.h"
#include "check.h"
#include <libxml/tree.h>
#include <vector>
#include <map>
#include <sstream>
#include <iostream>

const std::string FILE_STORAGE = "ascii";
const std::string SDTT_VECTOR = "Vector";
const std::string SDTT_MATRIX = "Matrix";
const std::string SDTT_VECTORSIZE = "vectorSize";
const std::string SDTT_MATRIXSIZE = "matrixSize";
const std::string SDTT_MATRIXCOLSIZE = "matrixColSize";
const std::string SDTT_MATRIXROWSIZE = "matrixRowSize";
const std::string SDTT_VECTORFILE = "vectorFile";
const std::string SDTT_VECTORPLUGIN = "vectorPlugin";
const std::string SDTT_MATRIXFILE = "matrixFile";
const std::string SDTT_ROW = "row";
const unsigned int MATRIX_MAX_SIZE = 10;
const unsigned int VECTOR_MAX_SIZE = 10;

class SiconosDOMTreeTools
{
public:

  /** \fn static SimpleVector getSiconosVectorValue(xmlNodePtr  vectorNode)
   *   \brief Return a SimpleVector read from a vector-type node
   *   \param xmlNodePtr : the vector node to be read
   *   \return A simpleVector
   */
  static SimpleVector getSiconosVectorValue(const xmlNodePtr);

  /** \fn static template<class T> void getVector(const xmlNodePtr vectorNode, vector<T>& outputVector)
   *   \brief read a vector-type node and save value of type T into a vector<T>
   *   \param xmlNodePtr: the vector node to be read
   *   \param vector<T> in-out parameter
   *   \return a vector<int>
   */
  template<class T> static void getVector(const xmlNodePtr vectorNode, std::vector<T>& outputVector)
  {
    if (vectorNode == NULL)
      XMLException::selfThrow("SiconosDOMTreeTools - getVector, node == NULL ");

    // vector loading, 3 options:
    //  - size of the vector + list of values
    //  - read from a file
    //  - read from a plug-in

    //  vector<int> outputVector; // the vector to be read ...
    unsigned int size;        // ... and its size

    // if attributes "vectorFile" is present
    if (xmlHasProp((xmlNodePtr)vectorNode, (xmlChar *)SDTT_VECTORFILE.c_str()))
    {
      // Case 1: read from an ascii file
      std::string fileName = SiconosDOMTreeTools::getStringAttributeValue(vectorNode, SDTT_VECTORFILE);
      std::ifstream in(fileName.c_str());
      in >> size;
      outputVector.reserve(size);
      T tmp;
      for (unsigned int i = 0; i < size; i++)
      {
        in >> tmp ;
        outputVector.push_back(tmp);
      }
    } //else if attribute vectorPlugin is present
    else if (xmlHasProp((xmlNodePtr)vectorNode, (xmlChar *)SDTT_VECTORPLUGIN.c_str()))
    {
      // Case 2: read using a plug-in
      XMLException::selfThrow("SiconosDOMTreeTools - getIntVector using a plug-in, not yet implemented");
    }
    else // if attribute vectorSize is present
    {
      xmlChar* tmp = xmlNodeGetContent((xmlNodePtr)vectorNode);
      std::string vOut = (char*)tmp;
      xmlFree(tmp);
      string2Vector(vOut, outputVector);
    }

    // return outputVector;
  };

  /** \fn static SimpleMatrix getSiconosMatrixValue(xmlNodePtr  siconosMatrixNode)
   *   \brief Return a SimpleMatrix computed from a siconos matrix node
   *   \param siconosMatrixNode : the matrix node you want to get in SimpleMatrix type
   *   \return A Sicono   else cout<<"!!! /!\ the number of dynamical system of the interaction is uneven, a couple of dynamical system in interaction is not complete."<<endl;
   sMatrix
  */
  static SimpleMatrix getSiconosMatrixValue(const xmlNodePtr  siconosMatrixNode);


  //    /** \fn static vector<SiconosVector*> getVectorMemoryValue(const xmlNodePtr  memoryNode)
  //    *   \brief Return a vector of SiconosVector computed from a memory node
  //    *   \param memoryNode : the memory node you want to get in a vector of SiconosVector type
  //    *   \return A  vector of SiconosVector
  //    */
  //    static vector<SiconosVector*> getVectorMemoryValue(const xmlNodePtr  memoryNode);
  //
  //    /** \fn static SiconosMemory getMemoryValue(const xmlNodePtr  memoryNode);
  //    *   \brief Return a SiconosMemory of SiconosVector computed from a memory node
  //    *   \param memoryNode : the memory node you want to get in a vector of SiconosVector type
  //    *   \return SiconosMemory
  //    */
  //    static SiconosMemory getMemoryValue(const xmlNodePtr  memoryNode);



  //-----------------

  /** \fn static const bool hasAttributeValue(const xmlNodePtr  node, const std::string attributeName) const
   *   \brief Returns the boolean value which shows if an attribut is defined or not
   *   \param node : the node who contents the attribute you want to if it exists
   *   \param attributeName : the attribute of the node you want to know if it exists
   *   \return true if the attribute exists, otherwise false
   */
  static const bool hasAttributeValue(const xmlNodePtr, const std::string&);


  /** \fn static string getStringAttributeValue(xmlNodePtr node, string attributeName)
   *   \brief Return the int value of the attribute attributeName of the node node
   *   \param node : the node who contents the attribute you want
   *   \param attributeName : the attribute of the node you want to have the int value
   *   \return The string value of the attribute attributeName contents in the node node
   */
  static std::string getStringAttributeValue(const xmlNodePtr node, const std::string& attributeName);

  /** \fn   template<typename T> static T getAttributeValue(const xmlNodePtr node, const std::string& attributeName)
   *   \brief Return the T-type value of the attribute attributeName of the node node
   *   \param node : the node who contents the attribute you want
   *   \param attributeName : attribute name
   *   \return type T object
   */
  template<typename T> static T getAttributeValue(const xmlNodePtr node, const std::string& attributeName)
  {
    if (!xmlHasProp((xmlNodePtr)node, (xmlChar *)attributeName.c_str()))
      XMLException::selfThrow("SiconosDOMTreeTools - getAttributeValue : the attribute " + attributeName + " doesn't exist in tag " + (char*)node->name);
    xmlChar* tmp = xmlGetProp((xmlNodePtr)node, (xmlChar *)(attributeName.c_str()));
    std::string vOut = (char*)tmp;
    xmlFree(tmp);
    // \todo Pb with bool and isstringstream -> find solution ...
    if (vOut == "true") vOut = "1";
    if (vOut == "false") vOut = "0";
    std::istringstream iss(vOut);
    T value;
    iss >> value;
    return value;
  };

  /** \fn static string getStringContentValue(xmlNodePtr  node)
   *   \brief Return the string content of the node node
   *   \param node : the node you want the string content
   *   \return The string value of the content of the node node
   */
  static std::string getStringContentValue(const xmlNodePtr  node);

  /** \fn   template<typename T> static T getContentValue(const xmlNodePtr node)
   *   \brief Return the type T content of a node
   *   \param node : the node you want the content
   *   \return a T
   */
  template<typename T> static T getContentValue(const xmlNodePtr node)
  {
    xmlChar* tmp = xmlNodeGetContent(node);
    std::string vOut = (char*)tmp;
    xmlFree(tmp);
    // \todo Pb with bool and isstringstream -> find solution ...
    if (vOut == "true") vOut = "1";
    if (vOut == "false") vOut = "0";
    std::istringstream iss(vOut);
    T value;
    iss >> value;
    return value;
  };

  /** \fn static void setSiconosVectorNodeValue(xmlNodePtr  siconosVectorNode, const SiconosVector& v)
   *   \brief Change values of a siconosVectorNode from a SiconosVector
   * if the vector is greater than VectorMaxSize, it will be saved in an external file
   * else, if the user already used an external file, this file will still be used
   *       else the vector is stored in the xml input/output file
   *   \param siconosVectorNode : the vector node you want to set
   *   \param v : the vector you want to copy the value in the siconosVectorNode
   *   \exception XMLException
   */
  static void setSiconosVectorNodeValue(const xmlNodePtr  siconosVectorNode, const SiconosVector &v);

  /** \fn static void setSiconosMatrixValue(xmlNodePtr  siconosMatrixNode, const SiconosMatrix& m)
   *   \brief Change values of a siconosMatrixNode from a SiconosMatrix
   *   \param siconosMatrixNode : the matrix node you want to set
   *   \param m : the matrix you want to copy the value in the siconosVectorNode
   *   \exception XMLException
   */
  static void setSiconosMatrixNodeValue(const xmlNodePtr  siconosMatrixNode, const SiconosMatrix& m);

  //    /** \fn static void setVectorMemoryValue(xmlNodePtr  memoryNode, vector<SiconosVector*> memory)
  //    *   \brief Change values of a memoryNode from a vector<SiconosVector>
  //    *   \param memoryNode : the memory node you want to set
  //    *   \param memory : the memory you want to copy the value in the memoryNode
  //    *   \exception XMLException
  //    */
  //    static void setVectorMemoryValue(const xmlNodePtr  memoryNode, const vector<SiconosVector*> memory);
  //
  //    /** \fn static void setMemoryValue(const xmlNodePtr  memoryNode, const SiconosMemory & memory)
  //    *   \brief Change values of a memoryNode from a SiconosMemory
  //    *   \param memoryNode : the memory node you want to set
  //    *   \param memory : the memory you want to copy the value in the memoryNode
  //    *   \exception XMLException
  //    */
  //    static void setMemoryValue( xmlNodePtr  memoryNode, const SiconosMemory & memory);

  /** \fn static void setStringAttributeValue(xmlNodePtr  node, char * attributeName, string value)
   *   \brief Set a string value at a node attribute
   *   \param node : the concern node
   *   \param attributeName : the concern attribute you want to set a string value
   *   \param value : the string value to set
   */
  static void setStringAttributeValue(const xmlNodePtr  node, const std::string attributeName, const std::string value);

  /** \fn static void setIntegerAttributeValue(xmlNodePtr  node, char * attributeName, int value)
   *   \brief Set a integer value at a node attribute
   *   \param node : the concern node
   *   \param attributeName : the concern attribute you want to set a integer value
   *   \param value : the integer value to set
   */
  static void setIntegerAttributeValue(const xmlNodePtr  node, const std::string attributeName, const int value);

  /** \fn static void setDoublAttributeeValue(xmlNodePtr  node, char * attributeName, double value)
   *   \brief Set a double value at a node attribute
   *   \param node : the concern node
   *   \param attributeName : the concern attribute you want to set a double value
   *   \param value : the double value to set
   */
  static void setDoubleAttributeValue(const xmlNodePtr  node, const std::string attributeName, const double value);

  /** \fn static void setBooleanAttributeeValue(xmlNodePtr  node, char * attributeName, bool value)
   *   \brief Set a boolean value at a node attribute
   *   \param node : the concern node
   *   \param attributeName : the concern attribute you want to set a boolean value
   *   \param value : the boolean value to set
   */
  static void setBooleanAttributeValue(const xmlNodePtr  node, const std::string attributeName, const bool value);

  /** \fn static void setIntegerContentValue(xmlNodePtr  node, int value)
   *   \brief Set an integer content at a node
   *   \param node : the node you want to set the integer content
   *   \param value : the integer value to set
   */
  static void setIntegerContentValue(const xmlNodePtr  node, const int value);

  /** \fn static void setDoubleContentValue(xmlNodePtr  node, double value)
   *   \brief Set a double content at a node
   *   \param node : the node you want to set the double content
   *   \param value : the double value to set
   */
  static void setDoubleContentValue(const xmlNodePtr  node, const double value);

  /** \fn static void setStringContentValue(xmlNodePtr  node, string value)
   *   \brief Set a string content at a node
   *   \param node : the node you want to set the string content
   *   \param value : the string value to set
   */
  static void setStringContentValue(const xmlNodePtr  node, const std::string value);

  /** \fn void setVectorIntContentValue(const xmlNodePtr  vectorNode, const vector<int> v)
   *   \brief set a vector<int> at a node
   *   \param vectorNode : the vector node you want to set in vector<int> type
   *   \param vector<int> : the vector of int to set
   */
  static void setVectorIntContentValue(const xmlNodePtr  vectorNode, const std::vector<int> v);

  //-----------------

  /** \fn static xmlNodePtr  createMatrixNode(xmlNodePtr , const string&, const SiconosMatrix& )
   *  \brief creates a new node in the DOM tree to save a SiconosMatrix
   *  \param xmlNode : the root node of the XML object calling this function
   *  \param string : the name of the balise of this SiconosMatrix
   *  \param SiconosMatrix : the SiconosMatrix to save in the XML file
   *  \return xmlNodePtr  : the node created
   */
  static xmlNodePtr  createMatrixNode(xmlNodePtr , const std::string&, const SiconosMatrix&);

  /** \fn static xmlNodePtr  createVectorNode(xmlNodePtr , const string&, SiconosVector& )
   *  \brief creates a new node in the DOM tree to save a SiconosVector
   *  \param xmlNode : the root node of the XML object calling this function
   *  \param string : the name of the balise of this SiconosVector
   *  \param SiconosVector& : the SiconosVector to save in the XML file
   *  \return xmlNodePtr  : the node created
   */
  static xmlNodePtr  createVectorNode(xmlNodePtr , const std::string&, const  SiconosVector&);

  /** \fn static xmlNodePtr  createVectorIntNode(xmlNodePtr , const string, vector<int> )
   *  \brief creates a new node in the DOM tree to save a SiconosVector
   *  \param xmlNode : the root node of the XML object calling this function
   *  \param string : the name of the balise of this SiconosVector
   *  \param vector<int> : the vector of int to save in the XML file
   *  \return xmlNodePtr  : the node created
   */
  static xmlNodePtr  createVectorIntNode(xmlNodePtr , const std::string, std::vector<int> v);

  //    /** \fn static xmlNodePtr  createVectorNode(xmlNodePtr , const string, vector<SiconosVector*> )
  //     *  \brief creates a new node in the DOM tree to save a SiconosVector
  //     *  \param xmlNode : the root node of the XML object calling this function
  //     *  \param string : the name of the balise of this vector
  //     *  \param vector<SiconosVector*> : the vector to save in the XML file
  //     *  \return xmlNodePtr  : the node created
  //     */
  //    static xmlNodePtr  createVectorMemoryNode(xmlNodePtr , const string, vector<SiconosVector*> );

  //    /** \fn static xmlNodePtr  createSiconosMemoryNode(xmlNodePtr , const string, SiconosMemory* )
  //     *  \brief creates a new node in the DOM tree to save a SiconosMemory
  //     *  \param xmlNode : the root node of the XML object calling this function
  //     *  \param string : the name of the balise of this SiconosMemory
  //     *  \param SiconosMemory* : the SiconosMemory to save in the XML file
  //     *  \return xmlNodePtr  : the node created
  //     */
  //    static xmlNodePtr  createSiconosMemoryNode(xmlNodePtr , const string, SiconosMemory* );

  /** \fn static xmlNodePtr  createDoubleNode(xmlNodePtr , const string, const double)
   *  \brief creates a new node in the DOM tree to save a double value
   *  \param xmlNode : the root node of the XML object calling this function
   *  \param string : the name of the balise of this double value
   *  \param double : the double value to save in the XML file
   *  \return xmlNodePtr  : the node created
   */
  static xmlNodePtr  createDoubleNode(xmlNodePtr , const std::string, const double);

  /** \fn static xmlNodePtr  createIntegerNode(xmlNodePtr , const string, const int)
   *  \brief creates a new node in the DOM tree to save an integer value
   *  \param xmlNode : the root node of the XML object calling this function
   *  \param string : the name of the balise of this integer value
   *  \param integer : the integer value to save in the XML file
   *  \return xmlNodePtr  : the node created
   */
  static xmlNodePtr  createIntegerNode(xmlNodePtr , const std::string, const int);

  /** \fn static xmlNodePtr  createBooleanNode(xmlNodePtr , const string, const bool)
   *  \brief creates a new node in the DOM tree to save a boolean value
   *  \param xmlNode : the root node of the XML object calling this function
   *  \param string : the name of the balise of this boolean value
   *  \param bool : the boolean value to save in the XML file
   *  \return xmlNodePtr  : the node created
   */
  static xmlNodePtr  createBooleanNode(xmlNodePtr , const std::string, const bool);

  /** \fn static xmlNodePtr  createStringNode(xmlNodePtr , const string, const string)
   *  \brief creates a new node in the DOM tree to save a string
   *  \param xmlNode : the root node of the XML object calling this function
   *  \param string : the name of the balise of this string
   *  \param string : the string to save in the XML file
   *  \return xmlNodePtr  : the node created
   */
  static xmlNodePtr  createStringNode(xmlNodePtr , const std::string, const std::string);

  /** \fn static xmlNodePtr  createSingleNode(xmlNodePtr , const string)
   *  \brief creates a new node in the DOM tree to save a string
   *  \param xmlNode : the root node of the XML object calling this function
   *  \param string : the name of the balise of this string
   *  \return xmlNodePtr  : the node created
   */
  static xmlNodePtr  createSingleNode(xmlNodePtr , const std::string);

  /** \fn static xmlNodePtr  createStringAttribute(xmlNodePtr , const string, const string)
   *  \brief creates a new attribute in the DOM tree to save a string
   *  \param xmlNode : the root node of the XML object calling this function
   *  \param string : the name of the balise of this string
   *  \param string : the string to save in the XML file
   *  \return xmlNodePtr  : the node created
   */
  static void createStringAttribute(xmlNodePtr , const std::string, const std::string);

  /** \fn static void createBooleanAttribute(xmlNodePtr , const string, const bool)
   *  \brief adds an attribute to a node in the DOM tree to save a boolean value
   *  \param xmlNode : the root node of the XML object calling this function
   *  \param string : the name of the attribute of this boolean value
   *  \param bool : the boolean value to save in the XML file
   *  \return xmlNodePtr  : the node created
   */
  static void createBooleanAttribute(xmlNodePtr , const std::string, const bool);

  //-----------------

  /** \fn static xmlNodePtr  findNodeChild(xmlNodePtr  node, string childNodeName)
   *   \brief Find the child node childNodeName of the node 'node'
   *   \param node : the node you want to search a element
   *   \param childNodeName : the name of the node you search
   *   \return the xmlNode or NULL if not found
   */
  static xmlNodePtr  findNodeChild(const xmlNodePtr  node, const std::string& childNodeName);

  /** \fn static xmlNodePtr  findNodeChild(xmlNodePtr  node)
   *   \brief Find the child node childNodeName of the node 'node'
   *   \param node : the node you want to search a element
   *   \return the xmlNode or NULL if not found
   */
  static xmlNodePtr  findNodeChild(const xmlNodePtr  node);

  /** \fn static xmlNodePtr  findFollowNode(xmlNodePtr  startNode, string followNodeName)
   *   \brief Find the first node with name followNodeName since the startNode node
   *   \param const xmlNodePtr  node : the node you want to begin search
   *   \param const string followNodeName : the name of the node you search
   *   \return the first xmlNode or NULL if not found
   */
  static xmlNodePtr  findFollowNode(const xmlNodePtr  node, const std::string& followNodeName);

  /** \fn static xmlNodePtr  findFollowNode(xmlNodePtr  startNode)
   *   \brief Find the first node with name followNodeName since the startNode node
   *   \param const xmlNodePtr  node : the node you want to begin search
   *   \param const string followNodeName : the name of the node you search
   *   \return the first xmlNode or NULL if not found
   */
  static xmlNodePtr  findFollowNode(const xmlNodePtr  node);

  /** \fn int getNodeChildrenNumber(const xmlNodePtr node)
   *   \brief get the number of children of a parent node given in parameters
   *   \param const xmlNodePtr : the node you want to know how many children it has
   *   \return int : 0 if the node has no child, -1 if the node doesn't exists, a positive number otherwise
   */
  static int getNodeChildrenNumber(const xmlNodePtr node);

private :

  /** \fn template<class T> void string2Vector(const string& strIn, vector<T>& vOut)
   *   \brief get a string and convert its content (separate by spaces) into a vector of T
   *   \param the string to convert
   *   \param a vector<T>, in-out parameter
   */
  template<class T> static void string2Vector(const std::string& strIn, std::vector<T>& vOut)
  {
    std::string strOut = strIn;

    // Remove all '\n' and '\t' in the string
    strOut.erase(std::remove(strOut.begin(), strOut.end(), '\n'), strOut.end());
    strOut.erase(std::remove(strOut.begin(), strOut.end(), '\t'), strOut.end());
    // Remove spaces at the beginning of the string
    //strOut = strOut.substr( strOut.find_first_not_of(' ') );

    std::istringstream iss(strOut);
    T value;

    std::string::size_type pos1, pos2;
    pos1 = strOut.find_first_not_of(' ', 0);
    while (pos1 != std::string::npos)
    {
      iss >> value;
      vOut.push_back(value);
      pos2 = strOut.find_first_of(' ', pos1);
      pos1 = strOut.find_first_not_of(' ', pos2);
    }
  };

  /** \fn static SiconosVector getSiconosRowMatrixValue(const xmlNodePtr  matrixRowNode, int rowSize)
   *   \brief Return a SiconosVector computed from a row of a matrix
   *   \param const xmlNodePtr  matrixRowNode : the matrix row node you want to get in SiconosVector type
   *   \param int rowSize : the size of the row
   *   \return A SiconosVector
   */
  static SimpleVector getSiconosRowMatrixValue(const xmlNodePtr  matrixRowNode, const int& rowSize);

  /** \fn static void setSiconosRowMatrixValue(const xmlNodePtr  matrixRowNode,  SiconosVector v, int colSize)
   *   \brief Set the row describes by matrixRowNode  of a matrix of col size colSize to v vector
   *   \param const xmlNodePtr  matrixRowNode : the matrix row node you want to set
   *   \param SiconosVector v : the new value of the row
   *   \param int colSize : the col size of the matrix who contains the row
   */
  static void setSiconosRowMatrixValue(const xmlNodePtr ,  const SiconosVector &, const unsigned int&);
};

#endif
