/* Siconos-Kernel version 1.1.3, Copyright INRIA 2005-2006.
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
 *   \version 1.1.3.
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

  /** \fn static SimpleVector getSiconosVectorValue(xmlNode * vectorNode)
   *   \brief Return a SimpleVector read from a vector-type node
   *   \param xmlNode*: the vector node to be read
   *   \return A simpleVector
   */
  static SimpleVector getSiconosVectorValue(const xmlNode *);

  /** \fn static vector<int> getIntVector(xmlNode * vectorNode)
   *   \brief Return a vector of int, read from vector-type node
   *   \param xmlNode*: the vector node to be read
   *   \return a vector<int>
   */
  static std::vector<int> getIntVector(const xmlNode *);

  /** \fn static SimpleMatrix getSiconosMatrixValue(xmlNode * siconosMatrixNode)
   *   \brief Return a SimpleMatrix computed from a siconos matrix node
   *   \param siconosMatrixNode : the matrix node you want to get in SimpleMatrix type
   *   \return A Sicono   else cout<<"!!! /!\ the number of dynamical system of the interaction is uneven, a couple of dynamical system in interaction is not complete."<<endl;
  sMatrix
   */
  static SimpleMatrix getSiconosMatrixValue(const xmlNode * siconosMatrixNode);


  //    /** \fn static vector<SiconosVector*> getVectorMemoryValue(const xmlNode * memoryNode)
  //    *   \brief Return a vector of SiconosVector computed from a memory node
  //    *   \param memoryNode : the memory node you want to get in a vector of SiconosVector type
  //    *   \return A  vector of SiconosVector
  //    */
  //    static vector<SiconosVector*> getVectorMemoryValue(const xmlNode * memoryNode);
  //
  //    /** \fn static SiconosMemory getMemoryValue(const xmlNode * memoryNode);
  //    *   \brief Return a SiconosMemory of SiconosVector computed from a memory node
  //    *   \param memoryNode : the memory node you want to get in a vector of SiconosVector type
  //    *   \return SiconosMemory
  //    */
  //    static SiconosMemory getMemoryValue(const xmlNode * memoryNode);



  //-----------------

  /** \fn static bool hasAttributeValue(const xmlNode * node, const std::string attributeName)
   *   \brief Returns the boolean value which shows if an attribut is defined or not
   *   \param node : the node who contents the attribute you want to if it exists
   *   \param attributeName : the attribute of the node you want to know if it exists
   *   \return true if the attribute exists, otherwise false
   */
  static bool hasAttributeValue(const xmlNode * node, const std::string& attributeName);


  /** \fn static string getStringAttributeValue(xmlNode * node, string attributeName)
   *   \brief Return the int value of the attribute attributeName of the node node
   *   \param node : the node who contents the attribute you want
   *   \param attributeName : the attribute of the node you want to have the int value
   *   \return The string value of the attribute attributeName contents in the node node
   */
  static std::string getStringAttributeValue(const xmlNode * node, const std::string& attributeName);

  /** \fn static int getIntegeAttributerValue(xmlNode * node, string attributeName)
   *   \brief Return the int value of the attribute attributeName of the node node
   *   \param node : the node who contents the attribute you want
   *   \param attributeName : the attribute of the node you want to have the int value
   *   \return The int value of the attribute attributeName contents in the node node
   */
  static int getIntegerAttributeValue(const xmlNode * node, const std::string& attributeName);

  /** \fn static double getDoubleAttributeValue(xmlNode * node, string attributeName)
   *   \brief Return the double value of the attribute attributeName of the node node
   *   \param node : the node who contents the attribute you want
   *   \param attributeName : the attribute of the node you want to have the double value
   *   \return The double value of the attribute attributeName contents in the node node
   */
  static double getDoubleAttributeValue(const xmlNode * node, const std::string& attributeName);

  /** \fn static bool getBooleanAttributeValue(xmlNode * node, string attributeName)
   *   \brief Return the boolean value of the attribute attributeName of the node node
   *   \param node : the node who contents the attribute you want
   *   \param attributeName : the attribute of the node you want to have the boolean value
   *   \return The boolean value of the attribute attributeName contents in the node node
   */
  static bool getBooleanAttributeValue(const xmlNode * node, const std::string& attributeName);


  /** \fn static string getStringContentValue(xmlNode * node)
   *   \brief Return the string content of the node node
   *   \param node : the node you want the string content
   *   \return The string value of the content of the node node
   */
  static std::string getStringContentValue(const xmlNode * node);

  /** \fn static int getIntegerContentValue(xmlNode * node)
   *   \brief Return the int content of the node node
   *   \param node : the node you want the int content
   *   \return The int value of the content of the node node
   */
  static int getIntegerContentValue(const xmlNode * node);

  /** \fn static double getDoubleContentValue(xmlNode * node)
   *   \brief Return the double content of the node node
   *   \param node : the node you want the double content
   *   \return The double value of the content of the node node
   */
  static double getDoubleContentValue(const xmlNode * node);

  /** \fn static bool getBooleanContentValue(const xmlNode * node)
   *   \brief Return the boolean content of the node node
   *   \param node : the node you want the double content
   *   \return bool value of the content of the node node
   */
  static bool getBooleanContentValue(const xmlNode * node);


  /** \fn static vector<int> getVectorIntContentValue(const xmlNode * siconosVectorNode)
   *   \brief Return a vector<int>
   *   \param vectorNode : the vector node you want to get in vector<int> type
   *   \return A vector<int>
   */
  static std::vector<int> getVectorIntContentValue(const xmlNode * vectorNode);


  //-----------------


  /** \fn static void setSiconosVectorNodeValue(xmlNode * siconosVectorNode, const SiconosVector& v)
   *   \brief Change values of a siconosVectorNode from a SiconosVector
   * if the vector is greater than VectorMaxSize, it will be saved in an external file
   * else, if the user already used an external file, this file will still be used
   *       else the vector is stored in the xml input/output file
   *   \param siconosVectorNode : the vector node you want to set
   *   \param v : the vector you want to copy the value in the siconosVectorNode
   *   \exception XMLException
   */
  static void setSiconosVectorNodeValue(const xmlNode * siconosVectorNode, const SiconosVector &v);

  /** \fn static void setSiconosMatrixValue(xmlNode * siconosMatrixNode, const SiconosMatrix& m)
   *   \brief Change values of a siconosMatrixNode from a SiconosMatrix
   *   \param siconosMatrixNode : the matrix node you want to set
   *   \param m : the matrix you want to copy the value in the siconosVectorNode
   *   \exception XMLException
   */
  static void setSiconosMatrixNodeValue(const xmlNode * siconosMatrixNode, const SiconosMatrix& m);

  //    /** \fn static void setVectorMemoryValue(xmlNode * memoryNode, vector<SiconosVector*> memory)
  //    *   \brief Change values of a memoryNode from a vector<SiconosVector>
  //    *   \param memoryNode : the memory node you want to set
  //    *   \param memory : the memory you want to copy the value in the memoryNode
  //    *   \exception XMLException
  //    */
  //    static void setVectorMemoryValue(const xmlNode * memoryNode, const vector<SiconosVector*> memory);
  //
  //    /** \fn static void setMemoryValue(const xmlNode * memoryNode, const SiconosMemory & memory)
  //    *   \brief Change values of a memoryNode from a SiconosMemory
  //    *   \param memoryNode : the memory node you want to set
  //    *   \param memory : the memory you want to copy the value in the memoryNode
  //    *   \exception XMLException
  //    */
  //    static void setMemoryValue( xmlNode * memoryNode, const SiconosMemory & memory);

  /** \fn static void setStringAttributeValue(xmlNode * node, char * attributeName, string value)
   *   \brief Set a string value at a node attribute
   *   \param node : the concern node
   *   \param attributeName : the concern attribute you want to set a string value
   *   \param value : the string value to set
   */
  static void setStringAttributeValue(const xmlNode * node, const std::string attributeName, const std::string value);

  /** \fn static void setIntegerAttributeValue(xmlNode * node, char * attributeName, int value)
   *   \brief Set a integer value at a node attribute
   *   \param node : the concern node
   *   \param attributeName : the concern attribute you want to set a integer value
   *   \param value : the integer value to set
   */
  static void setIntegerAttributeValue(const xmlNode * node, const std::string attributeName, const int value);

  /** \fn static void setDoublAttributeeValue(xmlNode * node, char * attributeName, double value)
   *   \brief Set a double value at a node attribute
   *   \param node : the concern node
   *   \param attributeName : the concern attribute you want to set a double value
   *   \param value : the double value to set
   */
  static void setDoubleAttributeValue(const xmlNode * node, const std::string attributeName, const double value);

  /** \fn static void setBooleanAttributeeValue(xmlNode * node, char * attributeName, bool value)
   *   \brief Set a boolean value at a node attribute
   *   \param node : the concern node
   *   \param attributeName : the concern attribute you want to set a boolean value
   *   \param value : the boolean value to set
   */
  static void setBooleanAttributeValue(const xmlNode * node, const std::string attributeName, const bool value);

  /** \fn static void setIntegerContentValue(xmlNode * node, int value)
   *   \brief Set an integer content at a node
   *   \param node : the node you want to set the integer content
   *   \param value : the integer value to set
   */
  static void setIntegerContentValue(const xmlNode * node, const int value);

  /** \fn static void setDoubleContentValue(xmlNode * node, double value)
   *   \brief Set a double content at a node
   *   \param node : the node you want to set the double content
   *   \param value : the double value to set
   */
  static void setDoubleContentValue(const xmlNode * node, const double value);

  /** \fn static void setStringContentValue(xmlNode * node, string value)
   *   \brief Set a string content at a node
   *   \param node : the node you want to set the string content
   *   \param value : the string value to set
   */
  static void setStringContentValue(const xmlNode * node, const std::string value);

  /** \fn void setVectorIntContentValue(const xmlNode * vectorNode, const vector<int> v)
   *   \brief set a vector<int> at a node
   *   \param vectorNode : the vector node you want to set in vector<int> type
   *   \param vector<int> : the vector of int to set
   */
  static void setVectorIntContentValue(const xmlNode * vectorNode, const std::vector<int> v);

  //-----------------

  /** \fn static xmlNode* createMatrixNode(xmlNode*, const string&, const SiconosMatrix& )
   *  \brief creates a new node in the DOM tree to save a SiconosMatrix
   *  \param xmlNode : the root node of the XML object calling this function
   *  \param string : the name of the balise of this SiconosMatrix
   *  \param SiconosMatrix : the SiconosMatrix to save in the XML file
   *  \return xmlNode* : the node created
   */
  static xmlNode* createMatrixNode(xmlNode*, const std::string&, const SiconosMatrix&);

  /** \fn static xmlNode* createVectorNode(xmlNode*, const string&, SiconosVector& )
   *  \brief creates a new node in the DOM tree to save a SiconosVector
   *  \param xmlNode : the root node of the XML object calling this function
   *  \param string : the name of the balise of this SiconosVector
   *  \param SiconosVector& : the SiconosVector to save in the XML file
   *  \return xmlNode* : the node created
   */
  static xmlNode* createVectorNode(xmlNode*, const std::string&, const  SiconosVector&);

  /** \fn static xmlNode* createVectorIntNode(xmlNode*, const string, vector<int> )
   *  \brief creates a new node in the DOM tree to save a SiconosVector
   *  \param xmlNode : the root node of the XML object calling this function
   *  \param string : the name of the balise of this SiconosVector
   *  \param vector<int> : the vector of int to save in the XML file
   *  \return xmlNode* : the node created
   */
  static xmlNode* createVectorIntNode(xmlNode*, const std::string, std::vector<int> v);

  //    /** \fn static xmlNode* createVectorNode(xmlNode*, const string, vector<SiconosVector*> )
  //     *  \brief creates a new node in the DOM tree to save a SiconosVector
  //     *  \param xmlNode : the root node of the XML object calling this function
  //     *  \param string : the name of the balise of this vector
  //     *  \param vector<SiconosVector*> : the vector to save in the XML file
  //     *  \return xmlNode* : the node created
  //     */
  //    static xmlNode* createVectorMemoryNode(xmlNode*, const string, vector<SiconosVector*> );

  //    /** \fn static xmlNode* createSiconosMemoryNode(xmlNode*, const string, SiconosMemory* )
  //     *  \brief creates a new node in the DOM tree to save a SiconosMemory
  //     *  \param xmlNode : the root node of the XML object calling this function
  //     *  \param string : the name of the balise of this SiconosMemory
  //     *  \param SiconosMemory* : the SiconosMemory to save in the XML file
  //     *  \return xmlNode* : the node created
  //     */
  //    static xmlNode* createSiconosMemoryNode(xmlNode*, const string, SiconosMemory* );

  /** \fn static xmlNode* createDoubleNode(xmlNode*, const string, const double)
   *  \brief creates a new node in the DOM tree to save a double value
   *  \param xmlNode : the root node of the XML object calling this function
   *  \param string : the name of the balise of this double value
   *  \param double : the double value to save in the XML file
   *  \return xmlNode* : the node created
   */
  static xmlNode* createDoubleNode(xmlNode*, const std::string, const double);

  /** \fn static xmlNode* createIntegerNode(xmlNode*, const string, const int)
   *  \brief creates a new node in the DOM tree to save an integer value
   *  \param xmlNode : the root node of the XML object calling this function
   *  \param string : the name of the balise of this integer value
   *  \param integer : the integer value to save in the XML file
   *  \return xmlNode* : the node created
   */
  static xmlNode* createIntegerNode(xmlNode*, const std::string, const int);

  /** \fn static xmlNode* createBooleanNode(xmlNode*, const string, const bool)
   *  \brief creates a new node in the DOM tree to save a boolean value
   *  \param xmlNode : the root node of the XML object calling this function
   *  \param string : the name of the balise of this boolean value
   *  \param bool : the boolean value to save in the XML file
   *  \return xmlNode* : the node created
   */
  static xmlNode* createBooleanNode(xmlNode*, const std::string, const bool);

  /** \fn static xmlNode* createStringNode(xmlNode*, const string, const string)
   *  \brief creates a new node in the DOM tree to save a string
   *  \param xmlNode : the root node of the XML object calling this function
   *  \param string : the name of the balise of this string
   *  \param string : the string to save in the XML file
   *  \return xmlNode* : the node created
   */
  static xmlNode* createStringNode(xmlNode*, const std::string, const std::string);

  /** \fn static xmlNode* createSingleNode(xmlNode*, const string)
   *  \brief creates a new node in the DOM tree to save a string
   *  \param xmlNode : the root node of the XML object calling this function
   *  \param string : the name of the balise of this string
   *  \return xmlNode* : the node created
   */
  static xmlNode* createSingleNode(xmlNode*, const std::string);

  /** \fn static xmlNode* createStringAttribute(xmlNode*, const string, const string)
   *  \brief creates a new attribute in the DOM tree to save a string
   *  \param xmlNode : the root node of the XML object calling this function
   *  \param string : the name of the balise of this string
   *  \param string : the string to save in the XML file
   *  \return xmlNode* : the node created
   */
  static void createStringAttribute(xmlNode*, const std::string, const std::string);

  /** \fn static void createBooleanAttribute(xmlNode*, const string, const bool)
   *  \brief adds an attribute to a node in the DOM tree to save a boolean value
   *  \param xmlNode : the root node of the XML object calling this function
   *  \param string : the name of the attribute of this boolean value
   *  \param bool : the boolean value to save in the XML file
   *  \return xmlNode* : the node created
   */
  static void createBooleanAttribute(xmlNode*, const std::string, const bool);

  //-----------------

  /** \fn static xmlNode * findNodeChild(xmlNode * node, string childNodeName)
   *   \brief Find the child node childNodeName of the node 'node'
   *   \param node : the node you want to search a element
   *   \param childNodeName : the name of the node you search
   *   \return the xmlNode or NULL if not found
   */
  static xmlNode * findNodeChild(const xmlNode * node, const std::string& childNodeName);

  /** \fn static xmlNode * findNodeChild(xmlNode * node)
   *   \brief Find the child node childNodeName of the node 'node'
   *   \param node : the node you want to search a element
   *   \return the xmlNode or NULL if not found
   */
  static xmlNode * findNodeChild(const xmlNode * node);

  /** \fn static xmlNode * findFollowNode(xmlNode * startNode, string followNodeName)
   *   \brief Find the first node with name followNodeName since the startNode node
   *   \param const xmlNode * node : the node you want to begin search
   *   \param const string followNodeName : the name of the node you search
   *   \return the first xmlNode or NULL if not found
   */
  static xmlNode * findFollowNode(const xmlNode * node, const std::string& followNodeName);

  /** \fn static xmlNode * findFollowNode(xmlNode * startNode)
   *   \brief Find the first node with name followNodeName since the startNode node
   *   \param const xmlNode * node : the node you want to begin search
   *   \param const string followNodeName : the name of the node you search
   *   \return the first xmlNode or NULL if not found
   */
  static xmlNode * findFollowNode(const xmlNode * node);

  /** \fn int getNodeChildrenNumber(const xmlNode *node)
   *   \brief get the number of children of a parent node given in parameters
   *   \param const xmlNode* : the node you want to know how many children it has
   *   \return int : 0 if the node has no child, -1 if the node doesn't exists, a positive number otherwise
   */
  static int getNodeChildrenNumber(const xmlNode *node);

private :
  /** \fn int findNextfigureIndex(const string s, int start)
   *   \brief string function, which returns the position in the string of the next char different of a space or a tab.
   *   \param string s : the string where we search a character
   *   \param int start : the index in the string where we start the research.
   *   \return int : the position of the next figure, -1 if not found.
   */
  static int findNextfigureIndex(const std::string& s, const int& start);

  /** \fn int findNextBlankIndex(const string s, int start)
   *   \brief string function, which returns the position in the string of the next space char.
   *   \param string s : the string where we search a space character
   *   \param int start : the index in the string where we start the research.
   *   \return int : the position of the next space, -1 if not found
   */
  static int findNextBlankIndex(const std::string& s, const int& start);

  /** \fn vector<double> string2Vector(const string vectorContent, int size)
   *   \brief Return a vector object which values are content in a char array
   *   \param string vectorContent : the node you want to search a element
   *   \param int size : the number of elements contain in vectorContent (the size of the vector to build)
   *   \return vector : the vector computes
   */
  static std::vector<double> string2Vector(const std::string& vectorContent, const int& size);

  /** \fn vector<int> string2Vector(const string vectorContent)
   *   \brief Return a vector object which values are content in a char array
   *   \param string vectorContent : the node you want to search a element
   *   \param int size : the number of elements contain in vectorContent (the size of the vector to build)
   *   \return vector : the vector computes
   */
  static std::vector<int> string2Vector(const std::string& vectorContent);

  /** \fn static SiconosVector getSiconosRowMatrixValue(const xmlNode * matrixRowNode, int rowSize)
   *   \brief Return a SiconosVector computed from a row of a matrix
   *   \param const xmlNode * matrixRowNode : the matrix row node you want to get in SiconosVector type
   *   \param int rowSize : the size of the row
   *   \return A SiconosVector
   */
  static SimpleVector getSiconosRowMatrixValue(const xmlNode * matrixRowNode, const int& rowSize);

  /** \fn static void setSiconosRowMatrixValue(const xmlNode * matrixRowNode,  SiconosVector v, int colSize)
   *   \brief Set the row describes by matrixRowNode  of a matrix of col size colSize to v vector
   *   \param const xmlNode * matrixRowNode : the matrix row node you want to set
   *   \param SiconosVector v : the new value of the row
   *   \param int colSize : the col size of the matrix who contains the row
   */
  static void setSiconosRowMatrixValue(const xmlNode *,  const SiconosVector &, const unsigned int&);
};

#endif
