/* Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
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
#include "SiconosDOMTreeTools.hpp"
#include "ioVector.hpp"
#include "ioMatrix.hpp"
#include "SimpleVector.hpp"
#include "SimpleMatrix.hpp"
#include "RuntimeException.hpp"

using namespace std;

SimpleVector SiconosDOMTreeTools::getSiconosVectorValue(const xmlNodePtr vectorNode)
{
  if (!vectorNode)
    XMLException::selfThrow("SiconosDOMTreeTools - getSiconosVectorValue, node == NULL ");

  // 2 cases:
  //   - read vector from a file
  //   - read in the xml file
  if (xmlHasProp((xmlNodePtr)vectorNode, (xmlChar *)SDTT_VECTORFILE.c_str())) //vector is defined in a extern ascii file
  {
    SimpleVector v(getStringAttributeValue(vectorNode, SDTT_VECTORFILE), true);
    return v;
  }
  else
  {
    //Size
    unsigned int size = getAttributeValue<unsigned int>(vectorNode, SDTT_VECTORSIZE);

    //Content
    xmlChar* tmp = xmlNodeGetContent((xmlNodePtr)vectorNode);
    string vectorContent = (char *)tmp;
    vector<double> tmpV;
    string2Vector(vectorContent, tmpV);
    SimpleVector v(tmpV);
    xmlFree(tmp);

    if (v.size() != size)
    {
      string s("size given in attribute and real size of the loaded vector are different in tag ");
      XMLException::selfThrow("SiconosDOMTreeTools - getSiconosVectorValue : " + s + (char*)vectorNode->name);
    }
    return v;
  }
}

SimpleMatrix SiconosDOMTreeTools::getSiconosMatrixValue(const xmlNodePtr siconosMatrixNode)
{
  if (!siconosMatrixNode)
    XMLException::selfThrow("SiconosDOMTreeTools - getSiconosMatrixValue, node == NULL");

  if (xmlHasProp((xmlNodePtr)siconosMatrixNode, (xmlChar *)SDTT_MATRIXFILE.c_str())) //matrix is defined in a extern ascii file
  {
    SimpleMatrix matrix(getStringAttributeValue(siconosMatrixNode, SDTT_MATRIXFILE), true);
    return matrix;
  }
  else
  {
    //The matrix is precised in the XML DOM Tree
    //number of lines
    unsigned int matrixRowSize = getAttributeValue<unsigned int>(siconosMatrixNode, SDTT_MATRIXROWSIZE);
    //number of columns
    unsigned int matrixColSize = getAttributeValue<unsigned int>(siconosMatrixNode, SDTT_MATRIXCOLSIZE);

    xmlNodePtr node = SiconosDOMTreeTools::findNodeChild(siconosMatrixNode, SDTT_ROW);
    unsigned int i = 0;
    SimpleMatrix matrix(matrixRowSize, matrixColSize);
    SP::SimpleVector v(new SimpleVector(matrixColSize));
    while ((node) && (i < matrixRowSize))
    {
      if (getSiconosRowMatrixValue(node, matrixColSize).size() != matrixColSize)
      {
        string s("A row in the matrix has not the right size in tag ");
        XMLException::selfThrow("SiconosDOMTreeTools - getSiconosMatrixValue : " + s + (char*)node->name);
      }
      *v = getSiconosRowMatrixValue(node, matrixColSize);
      matrix.setRow(i, *v);
      node = SiconosDOMTreeTools::findFollowNode(node, SDTT_ROW);
      i++;
    }
    return matrix;
  }
}


const bool SiconosDOMTreeTools::hasAttributeValue(const xmlNodePtr node, const string& attributeName)
{
  if (xmlHasProp((xmlNodePtr)node, (xmlChar *)attributeName.c_str()))
    return true;
  else
    return false;
}

string SiconosDOMTreeTools::getStringAttributeValue(const xmlNodePtr node, const string& attributeName)
{
  if (!xmlHasProp((xmlNodePtr)node, (xmlChar *)attributeName.c_str()))
    XMLException::selfThrow("SiconosDOMTreeTools - getStringAttributeValue : the attribute " + attributeName + " doesn't exist in tag " + (char*)node->name);
  xmlChar* tmp = xmlGetProp((xmlNodePtr)node, (xmlChar *)(attributeName.c_str()));
  string vOut = (char*)tmp;
  xmlFree(tmp);
  return vOut;
}

string SiconosDOMTreeTools::getStringContentValue(const xmlNodePtr node)
{
  xmlChar* tmp = xmlNodeGetContent((xmlNodePtr)node);
  string vOut = (char*)tmp;
  xmlFree(tmp);
  return vOut;
}

void SiconosDOMTreeTools::setSiconosVectorNodeValue(const xmlNodePtr siconosVectorNode, const SiconosVector& v)
{
  /*
   * if vector size > vectorMaxSize then put the vector in a file linked to the XML file
   */
  if (xmlHasProp((xmlNodePtr)siconosVectorNode, (xmlChar *)SDTT_VECTORFILE.c_str())) //vector is defined in a extern ascii file
  {
    string file = getStringAttributeValue(siconosVectorNode, SDTT_VECTORFILE);
    ioVector io(file, FILE_STORAGE);
    io.write(v);
  }
  else
  {
    //The vector is defined in the XML DOM Tree
    //Size
    unsigned int size = getAttributeValue<unsigned int>(siconosVectorNode, SDTT_VECTORSIZE);

    string vectorName = (char*)siconosVectorNode->name;
    if (size != v.size())
      XMLException::selfThrow("SiconosDOMTreeTools - setSiconosVectorNodeValue : the size of the " + vectorName +
                              " vector you want to save is different of the size defined the Kernel\ncheck the size of your DynamicalSystem");

    SimpleVector sv(v);
    xmlNodeSetContent((xmlNodePtr)siconosVectorNode, (xmlChar *)(sv.toString().c_str()));
  }
}

void SiconosDOMTreeTools::setSiconosMatrixNodeValue(const xmlNodePtr siconosMatrixNode, const SiconosMatrix& matrix)
{
  /*
   * if matrix size > xxx then put the matrix in a file linked to the XML file
   */
  if (xmlHasProp((xmlNodePtr)siconosMatrixNode, (xmlChar *)SDTT_MATRIXFILE.c_str())) //matrix is defined in a extern ascii file
  {
    string file = getStringAttributeValue(siconosMatrixNode, SDTT_MATRIXFILE);
    ioMatrix io(file, FILE_STORAGE);
    io.write(matrix);
  }
  else
  {
    //The matrix is precised in the XML DOM Tree
    //lineSize
    unsigned int matrixColSize = getAttributeValue<unsigned int>(siconosMatrixNode, SDTT_MATRIXCOLSIZE);
    //rowSize
    unsigned int matrixRowSize = getAttributeValue<unsigned int>(siconosMatrixNode, SDTT_MATRIXROWSIZE);

    string matrixName = (char*)siconosMatrixNode->name;

    if (matrixColSize != matrix.size(1))
      XMLException::selfThrow("SiconosDOMTreeTools - setSiconosMatrixValue : the " + matrixName + " matrix col size you want to save is different of the col size defined for it in xml");

    if (matrixRowSize != matrix.size(0))
      XMLException::selfThrow("SiconosDOMTreeTools - setSiconosMatrixValue : the " + matrixName + " matrix row size you want to save is different of the row size defined for it in xml");

    xmlNodePtr node = SiconosDOMTreeTools::findNodeChild(siconosMatrixNode, SDTT_ROW);

    unsigned int i = 0;
    SP::SimpleVector matRow(new SimpleVector(matrix.size(1)));
    while ((node) && (i < matrixRowSize))
    {
      matrix.getRow(i, *matRow);
      setSiconosRowMatrixValue(node, *matRow, matrixColSize);
      node = SiconosDOMTreeTools::findFollowNode(node, SDTT_ROW);
      i++;
    }
  }
}

void SiconosDOMTreeTools::setStringAttributeValue(const xmlNodePtr  node, const string attributeName, const string value)
{
  if (xmlHasProp((xmlNodePtr)node, (xmlChar *)attributeName.c_str()))
    xmlSetProp((xmlNodePtr) node, (xmlChar *)attributeName.c_str(), (xmlChar *)(value.c_str()));
  else
    XMLException::selfThrow("SiconosDOMTreeTools - setStringAttributeValue : the attribute " + attributeName + " doesn't exist in tag " + (char*)node->name);
}


void SiconosDOMTreeTools::setIntegerAttributeValue(const xmlNodePtr  node, const string attributeName, const int value)
{
  if (xmlHasProp((xmlNodePtr)node, (xmlChar *)attributeName.c_str()))
  {
    string stringValue;
    stringstream sstr;
    sstr << value;
    sstr >> stringValue;
    xmlSetProp((xmlNodePtr) node, (xmlChar *)attributeName.c_str(), (xmlChar *) stringValue.c_str());
  }
  else
    XMLException::selfThrow("SiconosDOMTreeTools - setIntegerAttributeValue : the attribute " + attributeName + "doesn't exist in tag " + (char*)node->name);
}

void SiconosDOMTreeTools::setDoubleAttributeValue(const xmlNodePtr  node, const string attributeName, const double value)
{
  if (xmlHasProp((xmlNodePtr)node, (xmlChar *)attributeName.c_str()))
  {
    string stringValue;
    stringstream sstr;
    sstr << value;
    sstr >> stringValue;
    xmlSetProp((xmlNodePtr) node, (xmlChar *)attributeName.c_str(), (xmlChar *) stringValue.c_str());
  }
  else
    XMLException::selfThrow("SiconosDOMTreeTools - setDoubleAttributeValue : the attribute " + attributeName + " doesn't exist in tag " + (char*)node->name);
}


void SiconosDOMTreeTools::setBooleanAttributeValue(const xmlNodePtr  node, const string attributeName, const bool value)
{
  if (xmlHasProp((xmlNodePtr)node, (xmlChar *)attributeName.c_str()))
  {
    string stringValue = "false";
    if (value) stringValue = "true";
    xmlSetProp((xmlNodePtr) node, (xmlChar *)attributeName.c_str(), (xmlChar *) stringValue.c_str());
  }
  else
    XMLException::selfThrow("SiconosDOMTreeTools - setBooleanAttributeValue : the attribute " + attributeName + " doesn't exist in tag " + (char*)node->name);
}


void SiconosDOMTreeTools::setStringContentValue(const xmlNodePtr  node, const string value)
{
  xmlNodeSetContent((xmlNodePtr) node, (xmlChar *)value.c_str());
}


void SiconosDOMTreeTools::setIntegerContentValue(const xmlNodePtr  node, const int value)
{
  string stringValue;
  stringstream sstr;

  sstr << value;
  sstr >> stringValue;
  xmlNodeSetContent((xmlNodePtr) node, (xmlChar *)stringValue.c_str());
}


void SiconosDOMTreeTools::setDoubleContentValue(const xmlNodePtr  node, const double value)
{
  string stringValue;
  stringstream sstr;

  sstr << value;
  sstr >> stringValue;
  xmlNodeSetContent((xmlNodePtr) node, (xmlChar *)stringValue.c_str());
}


void SiconosDOMTreeTools::setVectorIntContentValue(const xmlNodePtr  vectorNode, const vector<int> v)
{
  char element[100];
  string vectorContent = "";
  unsigned int i = 0;

  while (i < v.size())
  {
    strcpy(element, "");
    sprintf(element, "%i", v[i]);
    if (i > 0)
    {
      vectorContent += " ";
      vectorContent += element;
    }
    else vectorContent = element;
    i++;
  }
  xmlNodeSetContent((xmlNodePtr) vectorNode, (xmlChar *)vectorContent.c_str());
}

// -----------------

xmlNodePtr SiconosDOMTreeTools::createMatrixNode(xmlNodePtr rootNode, const string& name, const SiconosMatrix& matrix)
{
  /*
   * \todo if the SiconosMatrix is too big, the SiconosMatrix must be saved in an extern file and only the name of this file must be written in the XML file
   */
  xmlNodePtr node;
  node = xmlNewNode(NULL, BAD_CAST name.c_str());
  //  node = xmlNewChild(rootNode, NULL, BAD_CAST name.c_str(), NULL);

  string col, row;
  stringstream sstr, sstr2;

  sstr << matrix.size(1);
  sstr >> col;
  sstr2 << matrix.size(0);
  sstr2 >> row;

  xmlNodePtr  rowNode;

  xmlNewProp(node, (xmlChar*)(SDTT_MATRIXCOLSIZE.c_str()), (xmlChar*)col.c_str());
  xmlNewProp(node, (xmlChar*)(SDTT_MATRIXROWSIZE.c_str()), (xmlChar*)row.c_str());
  SP::SimpleVector matRow(new SimpleVector(matrix.size(1)));
  for (unsigned int i = 0; i < matrix.size(0); i++)
  {
    rowNode = new xmlNode();
    rowNode = xmlNewNode(NULL, BAD_CAST SDTT_ROW.c_str());
    matrix.getRow(i, *matRow);
    setSiconosRowMatrixValue(rowNode, *matRow, matrix.size(1));
    xmlAddChildList(node, rowNode);
    delete rowNode;
  }
  xmlAddChildList(rootNode, node);

  return node;
}

xmlNodePtr  SiconosDOMTreeTools::createVectorNode(xmlNodePtr  rootNode, const string& name, const  SiconosVector& v)
{
  xmlNodePtr node;

  if (v.size() < VECTOR_MAX_SIZE)
  {
    node = xmlNewNode(NULL, BAD_CAST name.c_str());

    string size;
    stringstream sstr;

    sstr << v.size();
    sstr >> size;

    xmlNewProp(node, (xmlChar*)(SDTT_VECTORSIZE.c_str()), (xmlChar*)size.c_str());

    setSiconosVectorNodeValue(node, v);

    xmlAddChildList(rootNode, node);
  }
  else
  {
    node = xmlNewNode(NULL, BAD_CAST name.c_str());
    string file = name + ".dat"; // \todo : add a time stamp to make this file unique
    xmlNewProp(node, (xmlChar *)SDTT_VECTORFILE.c_str(), (xmlChar*) file.c_str());

    setSiconosVectorNodeValue(node, v);

    xmlAddChildList(rootNode, node);
  }

  return node;
}

xmlNodePtr  SiconosDOMTreeTools::createVectorIntNode(xmlNodePtr  rootNode, const string name, vector<int> v)
{
  xmlNodePtr node;
  node = xmlNewNode(NULL, BAD_CAST name.c_str());

  setVectorIntContentValue(node, v);

  xmlAddChildList(rootNode, node);

  return node;
}

xmlNodePtr  SiconosDOMTreeTools::createDoubleNode(xmlNodePtr  rootNode, const string name, const double d)
{
  string stringValue;
  stringstream sstr;

  sstr << d;
  sstr >> stringValue;

  xmlNodePtr node;
  node = xmlNewNode(NULL, BAD_CAST name.c_str());

  xmlNodeSetContent((xmlNodePtr) node, (xmlChar *)stringValue.c_str());
  xmlAddChildList(rootNode, node);
  return node;
}

xmlNodePtr  SiconosDOMTreeTools::createIntegerNode(xmlNodePtr  rootNode, const string name, const int i)
{
  string stringValue;
  stringstream sstr;

  sstr << i;
  sstr >> stringValue;

  xmlNodePtr node;
  node = xmlNewNode(NULL, BAD_CAST name.c_str());

  xmlNodeSetContent((xmlNodePtr) node, (xmlChar *)stringValue.c_str());
  xmlAddChildList(rootNode, node);
  return node;
}

xmlNodePtr  SiconosDOMTreeTools::createBooleanNode(xmlNodePtr  rootNode, const string name, const bool b)
{
  string stringValue;

  xmlNodePtr node;
  node = xmlNewNode(NULL, BAD_CAST name.c_str());

  if (b) stringValue = "true";
  else stringValue = "false";
  xmlNodeSetContent((xmlNodePtr) node, (xmlChar *)stringValue.c_str());

  xmlAddChildList(rootNode, node);
  return node;
}

xmlNodePtr  SiconosDOMTreeTools::createStringNode(xmlNodePtr  rootNode, const string name, const string s)
{
  xmlNodePtr node;
  node = xmlNewNode(NULL, BAD_CAST name.c_str());

  //  string str = "\"" + s + "\"";
  string str = s;
  xmlNodeSetContent((xmlNodePtr) node, (xmlChar *)str.c_str());

  xmlAddChildList(rootNode, node);
  return node;
}

xmlNodePtr  SiconosDOMTreeTools::createSingleNode(xmlNodePtr  rootNode, const string name)
{
  xmlNodePtr node;
  node = xmlNewNode(NULL, BAD_CAST name.c_str());

  xmlAddChildList(rootNode, node);
  return node;
}

void SiconosDOMTreeTools::createStringAttribute(xmlNodePtr  node, const string name, const string s)
{
  xmlNewProp(node, (xmlChar*)(name.c_str()), (xmlChar*)s.c_str());
}

void SiconosDOMTreeTools::createBooleanAttribute(xmlNodePtr  node, const string s, const bool b)
{
  string stringValue;
  if (b) stringValue = "true";
  else    stringValue = "false";
  xmlNewProp(node, (xmlChar*)(s.c_str()), (xmlChar*)stringValue.c_str());
}

// -----------------

xmlNodePtr  SiconosDOMTreeTools::findNodeChild(const xmlNodePtr  node, const string& childNodeName)
{
  xmlNodePtr childNode = NULL;

  for (childNode = node->children; childNode; childNode = childNode->next)
  {
    if (childNode->type == XML_ELEMENT_NODE)
    {
      if (!xmlStrcmp(childNode->name, (xmlChar *)childNodeName.c_str()))
        return childNode;
    }
  }

  //Not found
  return NULL;
}

xmlNodePtr  SiconosDOMTreeTools::findNodeChild(const xmlNodePtr  node)
{
  xmlNodePtr childNode = NULL;

  for (childNode = node->children; childNode; childNode = childNode->next)
  {
    if (childNode->type == XML_ELEMENT_NODE)
    {
      return childNode;
    }
  }
  //Not found
  return NULL;
}


xmlNodePtr  SiconosDOMTreeTools::findFollowNode(const xmlNodePtr  node, const string& followNodeName)
{
  xmlNodePtr  n = (xmlNodePtr)node->next;
  while (n)
  {
    if (n->type == XML_ELEMENT_NODE)
    {
      if (!xmlStrcmp(n->name, (xmlChar *)followNodeName.c_str()))
        return n;
    }
    n = n->next;
  }

  //Not found
  return NULL;
}

xmlNodePtr  SiconosDOMTreeTools::findFollowNode(const xmlNodePtr  node)
{
  xmlNodePtr  n = (xmlNodePtr)node->next;
  while (n)
  {
    if (n->type == XML_ELEMENT_NODE)
    {
      return n;
    }
    n = n->next;
  }
  //Not found
  return NULL;
}

SimpleVector SiconosDOMTreeTools::getSiconosRowMatrixValue(const xmlNodePtr  siconosMatrixRowNode, const int& colSize)
{
  if (xmlHasProp((xmlNodePtr) siconosMatrixRowNode, (xmlChar *)SDTT_VECTORFILE.c_str())) //row is defined in a extern ascii file
  {
    SimpleVector v(getStringAttributeValue(siconosMatrixRowNode, SDTT_VECTORFILE), true);
    return v;
  }
  else
  {
    //The row is precised in the XML DOM Tree
    //Content
    xmlChar * tmp =  xmlNodeGetContent((xmlNodePtr) siconosMatrixRowNode);
    string vectorContent = (char *)tmp;
    vector<double> tmpV;
    string2Vector(vectorContent, tmpV);
    SimpleVector v(tmpV);
    xmlFree(tmp);
    return v;
  }
}

void SiconosDOMTreeTools::setSiconosRowMatrixValue(const xmlNodePtr  siconosMatrixRowNode, const SiconosVector &v, const unsigned int& colSize)
{
  if (colSize != v.size())
    XMLException::selfThrow("SiconosDOMTreeTools - setSiconosRowMatrixValue : a row size you want to save is different of the matrix size defined in xml");

  char element[100];
  unsigned int i = 0;
  string vectorContent = "";
  while (i < v.size())
  {
    strcpy(element, "");
    sprintf(element, N_DOUBLE_PRECISION, v(i));
    if (i > 0)
    {
      vectorContent += " ";
      vectorContent += element;
    }
    else vectorContent = element;
    i++;
  }
  xmlNodeSetContent((xmlNodePtr)siconosMatrixRowNode, (xmlChar *)(vectorContent.c_str()));
}

int SiconosDOMTreeTools::getNodeChildrenNumber(const xmlNodePtr node)
{
  int res = 0;
  xmlNodePtr n;

  if (!(node) && (node->type != XML_ELEMENT_NODE))
    res = -1;
  else
  {
    n = SiconosDOMTreeTools::findNodeChild((const xmlNodePtr) node);
    //if( n != NULL ) res++;

    while (n)
    {
      n = findFollowNode(n);
      res++;
    }
  }
  return res;
}
