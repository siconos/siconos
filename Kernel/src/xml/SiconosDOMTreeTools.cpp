#include "SiconosDOMTreeTools.h"
#include "check.h"


SimpleVector SiconosDOMTreeTools::getSiconosVectorValue(const xmlNode * siconosVectorNode)
{
  if (siconosVectorNode != NULL)
  {
    if (xmlHasProp((xmlNodePtr)siconosVectorNode, (xmlChar *)SDTT_VECTORFILE.c_str())) //vector is defined in a extern ascii file
    {
      SimpleVector v(getStringAttributeValue(siconosVectorNode, SDTT_VECTORFILE), true);

      return v;
    }

    //The vector is precised in the XML DOM Tree

    //Size
    int size = getIntegerAttributeValue(siconosVectorNode, SDTT_VECTORSIZE);

    //Content
    string vectorContent = (char *)xmlNodeGetContent((xmlNode *)siconosVectorNode);
    SimpleVector v(string2Vector(vectorContent, size));

    if (v.size() == size)
      return v;
    else
    {
      string s("the size given in attribute and the size of the vector effectively loaded are not the same in tag ");
      XMLException::selfThrow("SiconosDOMTreeTools - getSiconosVectorValue : " + s + (char*)siconosVectorNode->name);
    }
  }
  else
  {
    cout << "getSiconosVectorValue - siconosVectorNode == NULL, node not found in the DOM tree, perhaps this attribute is only optional" << endl;
    SimpleVector v;
    return  v;
  }
}


SiconosMatrix SiconosDOMTreeTools::getSiconosMatrixValue(const xmlNode * siconosMatrixNode)
{
  if (siconosMatrixNode != NULL)
  {
    if (xmlHasProp((xmlNodePtr)siconosMatrixNode, (xmlChar *)SDTT_MATRIXFILE.c_str())) //matrix is defined in a extern ascii file
    {
      SiconosMatrix m(getStringAttributeValue(siconosMatrixNode, SDTT_MATRIXFILE), true);
      return m;
    }

    //The matrix is precised in the XML DOM Tree
    //lineSize
    int matrixColSize = getIntegerAttributeValue(siconosMatrixNode, SDTT_MATRIXCOLSIZE);
    //rowSize
    int matrixRowSize = getIntegerAttributeValue(siconosMatrixNode, SDTT_MATRIXROWSIZE);

    xmlNode *node = SiconosDOMTreeTools::findNodeChild(siconosMatrixNode, SDTT_ROW);
    int i = 0;
    SiconosMatrix matrix(matrixRowSize, matrixColSize);
    SimpleVector v;

    while ((node != NULL) && (i < matrixRowSize))
    {
      v = getSiconosRowMatrixValue(node, matrixColSize);
      if (v.size() == matrixColSize)
      {
        matrix.addRow(i, v);
        node = SiconosDOMTreeTools::findFollowNode(node, SDTT_ROW);
        i++;
      }
      else
      {
        string s("A row in the matrix has not the right size in tag ");
        XMLException::selfThrow("SiconosDOMTreeTools - getSiconosMatrixValue : " + s + (char*)node->name);
      }
    }

    return matrix;
  }
  else
  {
    cout << "getSiconosMatrixValue - siconosMatrixNode == NULL, node not found in the DOM tree, perhaps this attribute is only optional" << endl;
    SiconosMatrix m;
    return  m;
  }
}



bool SiconosDOMTreeTools::hasAttributeValue(const xmlNode * node, const string attributeName)
{
  if (xmlHasProp((xmlNodePtr)node, (xmlChar *)attributeName.c_str()))
    return true;
  else
    return false;
}

string SiconosDOMTreeTools::getStringAttributeValue(const xmlNode * node, const string attributeName)
{
  if (xmlHasProp((xmlNodePtr)node, (xmlChar *)attributeName.c_str()))
    return string((char  *)xmlGetProp((xmlNode*)node, (xmlChar *)(attributeName.c_str())));
  else
    XMLException::selfThrow("SiconosDOMTreeTools - getStringAttributeValue : the attribute " + attributeName + " doesn't exist in tag " + (char*)node->name);

}


int SiconosDOMTreeTools::getIntegerAttributeValue(const xmlNode * node, const string attributeName)
{
  if (xmlHasProp((xmlNodePtr)node, (xmlChar *)attributeName.c_str()))
    return atoi((char  *)xmlGetProp((xmlNode *)node, (xmlChar *)attributeName.c_str()));
  else
    XMLException::selfThrow("SiconosDOMTreeTools - getIntegerAttributeValue : the attribute " + attributeName + " doesn't exist in tag " + (char*)node->name);

}


double SiconosDOMTreeTools::getDoubleAttributeValue(const xmlNode * node, const string  attributeName)
{
  double propDoubleValue;
  char * propCharValue;

  if (xmlHasProp((xmlNodePtr)node, (xmlChar *)attributeName.c_str()))
  {
    propCharValue = (char  *)xmlGetProp((xmlNode *)node, (xmlChar *)attributeName.c_str());
    stringstream sstr;
    sstr <<  propCharValue;
    sstr >> propDoubleValue;
    return propDoubleValue;
  }
  else
    XMLException::selfThrow("SiconosDOMTreeTools - getDoubleAttributeValue : the attribute " + attributeName + " doesn't exist in tag " + (char*)node->name);

}


bool SiconosDOMTreeTools::getBooleanAttributeValue(const xmlNode * node, const string attributeName)
{
  if (xmlHasProp((xmlNodePtr)node, (xmlChar *)attributeName.c_str()))
  {
    bool propBooleanValue = false;
    string val = (char  *)xmlGetProp((xmlNode *)node, (xmlChar *)attributeName.c_str());
    if (val == "true") propBooleanValue = true;
    return propBooleanValue;
  }
  else
    XMLException::selfThrow("SiconosDOMTreeTools - getBooleanAttributeValue : the attribute " + attributeName + " doesn't exist in tag " + (char*)node->name);
}



string SiconosDOMTreeTools::getStringContentValue(const xmlNode * node)
{
  return string((char *)xmlNodeGetContent((xmlNode *)node));
}


int SiconosDOMTreeTools::getIntegerContentValue(const xmlNode * node)
{
  return atoi((char *)xmlNodeGetContent((xmlNode *)node));
}


double SiconosDOMTreeTools::getDoubleContentValue(const xmlNode * node)
{
  return atof((char *)xmlNodeGetContent((xmlNode *)node));
}

bool SiconosDOMTreeTools::getBooleanContentValue(const xmlNode * node)
{
  if (strcmp((char *)xmlNodeGetContent((xmlNode *)node), "true") == 0) return true;
  else return false;
}

vector<int> SiconosDOMTreeTools::getVectorIntContentValue(const xmlNode * vectorNode)
{
  return string2Vector((char *)xmlNodeGetContent((xmlNode *)vectorNode));
}


void SiconosDOMTreeTools::setSiconosVectorValue(const xmlNode * siconosVectorNode, const SiconosVector& v)
{
  /*
   * if vector size > vectorMaxSize then put the vector in a file linked to the XML file
   */
  if (xmlHasProp((xmlNodePtr)siconosVectorNode, (xmlChar *)SDTT_VECTORFILE.c_str())) //vector is defined in a extern ascii file
  {
    string file = getStringAttributeValue(siconosVectorNode, SDTT_VECTORFILE);
    v.write(file, /*ASCII*/ /*N_ASCII*/ FILE_STORAGE); //For the moment only ASCII file are managed
  }
  else
  {
    //The vector is defined in the XML DOM Tree
    //Size
    int size = getIntegerAttributeValue(siconosVectorNode, SDTT_VECTORSIZE);

    string vectorName = (char*)siconosVectorNode->name;
    if (size != v.size())
      XMLException::selfThrow("SiconosDOMTreeTools - setSiconosVectorValue : the size of the " + vectorName +
                              " vector you want to save is different of the size defined the Kernel\ncheck the size of your DynamicalSystem");

    SimpleVector sv(v);
    xmlNodeSetContent((xmlNodePtr)siconosVectorNode, (xmlChar *)(sv.toString().c_str()));
  }
}

void SiconosDOMTreeTools::setSiconosVectorValue(const xmlNode * siconosVectorNode, SiconosVector *v)
{
  /*
   * if vector size > vectorMaxSize then put the vector in a file linked to the XML file
   */
  if (xmlHasProp((xmlNodePtr)siconosVectorNode, (xmlChar *)SDTT_VECTORFILE.c_str())) //vector is defined in a extern ascii file
  {
    string file = getStringAttributeValue(siconosVectorNode, SDTT_VECTORFILE);
    v->write(file, /*ASCII*/ /*N_ASCII*/ FILE_STORAGE); //For the moment only ASCII file are managed
  }
  else
  {
    //The vector is defined in the XML DOM Tree
    //Size
    int size = getIntegerAttributeValue(siconosVectorNode, SDTT_VECTORSIZE);

    string vectorName = (char*)siconosVectorNode->name;

    if (size != v->size())
      XMLException::selfThrow("SiconosDOMTreeTools - setSiconosVectorValue : the size of the " + vectorName +
                              " vector you want to save is different of the size defined the Kernel\ncheck the size of your DynamicalSystem");

    xmlNodeSetContent((xmlNodePtr)siconosVectorNode, (xmlChar *)(v->toString().c_str()));
  }
}


void SiconosDOMTreeTools::setSiconosMatrixValue(const xmlNode * siconosMatrixNode, SiconosMatrix matrix)
{
  /*
   * if matrix size > xxx then put the matrix in a file linked to the XML file
   */
  if (xmlHasProp((xmlNodePtr)siconosMatrixNode, (xmlChar *)SDTT_MATRIXFILE.c_str())) //matrix is defined in a extern ascii file
  {
    string file = getStringAttributeValue(siconosMatrixNode, SDTT_MATRIXFILE);
    matrix.write(file, /*ASCII*/ /*N_ASCII*/ FILE_STORAGE); //For the moment only ASCII file are managed
  }
  else
  {
    //The matrix is precised in the XML DOM Tree
    //lineSize
    int matrixColSize = getIntegerAttributeValue(siconosMatrixNode, SDTT_MATRIXCOLSIZE);
    //rowSize
    int matrixRowSize = getIntegerAttributeValue(siconosMatrixNode, SDTT_MATRIXROWSIZE);

    string matrixName = (char*)siconosMatrixNode->name;

    if (matrixColSize != matrix.size(1))
      XMLException::selfThrow("SiconosDOMTreeTools - setSiconosMatrixValue : the " + matrixName + " matrix col size you want to save is different of the col size defined for it in xml");

    if (matrixRowSize != matrix.size(0))
      XMLException::selfThrow("SiconosDOMTreeTools - setSiconosMatrixValue : the " + matrixName + " matrix row size you want to save is different of the row size defined for it in xml");

    xmlNode *node = SiconosDOMTreeTools::findNodeChild(siconosMatrixNode, SDTT_ROW);

    int i = 0;
    while ((node != NULL) && (i < matrixRowSize))
    {
      setSiconosRowMatrixValue(node, matrix.getRow(i), matrixColSize);
      node = SiconosDOMTreeTools::findFollowNode(node, SDTT_ROW);
      i++;
    }
  }
}

void SiconosDOMTreeTools::setSiconosMatrixValue(const xmlNode * siconosMatrixNode, SiconosMatrix *matrix)
{
  /*
   * if matrix size > xxx then put the matrix in a file linked to the XML file
   */
  if (xmlHasProp((xmlNodePtr)siconosMatrixNode, (xmlChar *)SDTT_MATRIXFILE.c_str())) //matrix is defined in a extern ascii file
  {
    string file = getStringAttributeValue(siconosMatrixNode, SDTT_MATRIXFILE);
    matrix->write(file, /*ASCII*/ /*N_ASCII*/ FILE_STORAGE); //For the moment only ASCII file are managed
  }
  else
  {
    //The matrix is precised in the XML DOM Tree
    //lineSize
    int matrixColSize = getIntegerAttributeValue(siconosMatrixNode, SDTT_MATRIXCOLSIZE);
    //rowSize
    int matrixRowSize = getIntegerAttributeValue(siconosMatrixNode, SDTT_MATRIXROWSIZE);

    string matrixName = (char*)siconosMatrixNode->name;

    if (matrixColSize != matrix->size(1))
    {
      matrixColSize = matrix->size(1);
      cout << "SiconosDOMTreeTools - setSiconosMatrixValue : the " + matrixName + " matrix col size you want to save is different of the col size already defined." << endl;
      setIntegerAttributeValue(siconosMatrixNode, SDTT_MATRIXCOLSIZE, matrix->size(1));
    }

    if (matrixRowSize != matrix->size(0))
    {
      matrixRowSize = matrix->size(0);
      cout << "SiconosDOMTreeTools - setSiconosMatrixValue : the " + matrixName + " matrix row size you want to save is different of the row size already defined." << endl;
      setIntegerAttributeValue(siconosMatrixNode, SDTT_MATRIXROWSIZE, matrix->size(0));
    }

    xmlNode *node = SiconosDOMTreeTools::findNodeChild(siconosMatrixNode, SDTT_ROW);

    int i = 0;
    while ((node != NULL) && (i < matrixRowSize))
    {
      setSiconosRowMatrixValue(node, matrix->getRow(i), matrixColSize);
      node = SiconosDOMTreeTools::findFollowNode(node, SDTT_ROW);
      i++;
    }
  }
}


void SiconosDOMTreeTools::setStringAttributeValue(const xmlNode * node, const string attributeName, const string value)
{
  if (xmlHasProp((xmlNodePtr)node, (xmlChar *)attributeName.c_str()))
    xmlSetProp((xmlNode *) node, (xmlChar *)attributeName.c_str(), (xmlChar *)(value.c_str()));
  else
    XMLException::selfThrow("SiconosDOMTreeTools - setStringAttributeValue : the attribute " + attributeName + " doesn't exist in tag " + (char*)node->name);
}


void SiconosDOMTreeTools::setIntegerAttributeValue(const xmlNode * node, const string attributeName, const int value)
{
  if (xmlHasProp((xmlNodePtr)node, (xmlChar *)attributeName.c_str()))
  {
    string stringValue;
    stringstream sstr;
    sstr << value;
    sstr >> stringValue;
    xmlSetProp((xmlNode *) node, (xmlChar *)attributeName.c_str(), (xmlChar *) stringValue.c_str());
  }
  else
    XMLException::selfThrow("SiconosDOMTreeTools - setIntegerAttributeValue : the attribute " + attributeName + "doesn't exist in tag " + (char*)node->name);
}

void SiconosDOMTreeTools::setDoubleAttributeValue(const xmlNode * node, const string attributeName, const double value)
{
  if (xmlHasProp((xmlNodePtr)node, (xmlChar *)attributeName.c_str()))
  {
    string stringValue;
    stringstream sstr;
    sstr << value;
    sstr >> stringValue;
    xmlSetProp((xmlNode *) node, (xmlChar *)attributeName.c_str(), (xmlChar *) stringValue.c_str());
  }
  else
    XMLException::selfThrow("SiconosDOMTreeTools - setDoubleAttributeValue : the attribute " + attributeName + " doesn't exist in tag " + (char*)node->name);
}


void SiconosDOMTreeTools::setBooleanAttributeValue(const xmlNode * node, const string attributeName, const bool value)
{
  if (xmlHasProp((xmlNodePtr)node, (xmlChar *)attributeName.c_str()))
  {
    string stringValue = "false";
    if (value) stringValue = "true";
    xmlSetProp((xmlNode *) node, (xmlChar *)attributeName.c_str(), (xmlChar *) stringValue.c_str());
  }
  else
    XMLException::selfThrow("SiconosDOMTreeTools - setBooleanAttributeValue : the attribute " + attributeName + " doesn't exist in tag " + (char*)node->name);
}


void SiconosDOMTreeTools::setStringContentValue(const xmlNode * node, const string value)
{
  xmlNodeSetContent((xmlNode *) node, (xmlChar *)value.c_str());
}


void SiconosDOMTreeTools::setIntegerContentValue(const xmlNode * node, const int value)
{
  string stringValue;
  stringstream sstr;

  sstr << value;
  sstr >> stringValue;
  xmlNodeSetContent((xmlNode *) node, (xmlChar *)stringValue.c_str());
}


void SiconosDOMTreeTools::setDoubleContentValue(const xmlNode * node, const double value)
{
  string stringValue;
  stringstream sstr;

  sstr << value;
  sstr >> stringValue;
  xmlNodeSetContent((xmlNode *) node, (xmlChar *)stringValue.c_str());
}


void SiconosDOMTreeTools::setVectorIntContentValue(const xmlNode * vectorNode, const vector<int> v)
{
  char element[100];
  string vectorContent = "";
  int i = 0;

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
  xmlNodeSetContent((xmlNode *) vectorNode, (xmlChar *)vectorContent.c_str());
}

// -----------------

xmlNode* SiconosDOMTreeTools::createMatrixNode(xmlNode* rootNode, const string name, SiconosMatrix* matrix)
{
  /*
   * \todo if the SiconosMatrix is too big, the SiconosMatrix must be saved in an extern file and only the name of this file must be written in the XML file
   */
  xmlNode *node;
  node = xmlNewNode(NULL, BAD_CAST name.c_str());
  //  node = xmlNewChild(rootNode, NULL, BAD_CAST name.c_str(), NULL);

  string col, row;
  stringstream sstr, sstr2;

  sstr << matrix->size(1);
  sstr >> col;
  sstr2 << matrix->size(0);
  sstr2 >> row;

  xmlNode* rowNode;

  xmlNewProp(node, (xmlChar*)(SDTT_MATRIXCOLSIZE.c_str()), (xmlChar*)col.c_str());
  xmlNewProp(node, (xmlChar*)(SDTT_MATRIXROWSIZE.c_str()), (xmlChar*)row.c_str());
  for (int i = 0; i < matrix->size(0); i++)
  {
    rowNode = new xmlNode();
    rowNode = xmlNewNode(NULL, BAD_CAST SDTT_ROW.c_str());
    setSiconosRowMatrixValue(rowNode, matrix->getRow(i), matrix->size(1));
    xmlAddChildList(node, rowNode);
  }
  xmlAddChildList(rootNode, node);

  return node;
}

xmlNode* SiconosDOMTreeTools::createVectorNode(xmlNode* rootNode, const string name, SiconosVector* v)
{
  xmlNode *node;

  if (v->size() < VECTOR_MAX_SIZE)
  {
    node = xmlNewNode(NULL, BAD_CAST name.c_str());

    string size;
    stringstream sstr;

    sstr << v->size();
    sstr >> size;

    xmlNewProp(node, (xmlChar*)(SDTT_VECTORSIZE.c_str()), (xmlChar*)size.c_str());

    setSiconosVectorValue(node, v);


    xmlAddChildList(rootNode, node);
  }
  else
  {
    node = xmlNewNode(NULL, BAD_CAST name.c_str());
    string file = name + ".dat"; // \todo : add a time stamp to make this file unique
    xmlNewProp(node, (xmlChar *)SDTT_VECTORFILE.c_str(), (xmlChar*) file.c_str());

    setSiconosVectorValue(node, v);

    xmlAddChildList(rootNode, node);
  }

  return node;
}

xmlNode* SiconosDOMTreeTools::createVectorIntNode(xmlNode* rootNode, const string name, vector<int> v)
{
  xmlNode *node;
  node = xmlNewNode(NULL, BAD_CAST name.c_str());

  setVectorIntContentValue(node, v);

  xmlAddChildList(rootNode, node);

  return node;
}

xmlNode* SiconosDOMTreeTools::createDoubleNode(xmlNode* rootNode, const string name, const double d)
{
  string stringValue;
  stringstream sstr;

  sstr << d;
  sstr >> stringValue;

  xmlNode *node;
  node = xmlNewNode(NULL, BAD_CAST name.c_str());

  xmlNodeSetContent((xmlNode *) node, (xmlChar *)stringValue.c_str());
  xmlAddChildList(rootNode, node);
  return node;
}

xmlNode* SiconosDOMTreeTools::createIntegerNode(xmlNode* rootNode, const string name, const int i)
{
  string stringValue;
  stringstream sstr;

  sstr << i;
  sstr >> stringValue;

  xmlNode *node;
  node = xmlNewNode(NULL, BAD_CAST name.c_str());

  xmlNodeSetContent((xmlNode *) node, (xmlChar *)stringValue.c_str());
  xmlAddChildList(rootNode, node);
  return node;
}

xmlNode* SiconosDOMTreeTools::createBooleanNode(xmlNode* rootNode, const string name, const bool b)
{
  string stringValue;

  xmlNode *node;
  node = xmlNewNode(NULL, BAD_CAST name.c_str());

  if (b) stringValue = "true";
  else stringValue = "false";
  xmlNodeSetContent((xmlNode *) node, (xmlChar *)stringValue.c_str());

  xmlAddChildList(rootNode, node);
  return node;
}

xmlNode* SiconosDOMTreeTools::createStringNode(xmlNode* rootNode, const string name, const string s)
{
  xmlNode *node;
  node = xmlNewNode(NULL, BAD_CAST name.c_str());

  //  string str = "\"" + s + "\"";
  string str = s;
  xmlNodeSetContent((xmlNode *) node, (xmlChar *)str.c_str());

  xmlAddChildList(rootNode, node);
  return node;
}

xmlNode* SiconosDOMTreeTools::createSingleNode(xmlNode* rootNode, const string name)
{
  xmlNode *node;
  node = xmlNewNode(NULL, BAD_CAST name.c_str());

  xmlAddChildList(rootNode, node);
  return node;
}

xmlNode* SiconosDOMTreeTools::createStringAttribute(xmlNode* node, const string name, const string s)
{
  xmlNewProp(node, (xmlChar*)(name.c_str()), (xmlChar*)s.c_str());
}

void SiconosDOMTreeTools::createBooleanAttribute(xmlNode* node, const string s, const bool b)
{
  string stringValue;
  if (b) stringValue = "true";
  else    stringValue = "false";

  xmlNewProp(node, (xmlChar*)(s.c_str()), (xmlChar*)stringValue.c_str());
}

// -----------------

xmlNode * SiconosDOMTreeTools::findNodeChild(const xmlNode * node, const string childNodeName)
{
  xmlNode *childNode = NULL;

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

xmlNode * SiconosDOMTreeTools::findNodeChild(const xmlNode * node)
{
  xmlNode *childNode = NULL;

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


xmlNode * SiconosDOMTreeTools::findFollowNode(const xmlNode * node, const string followNodeName)
{
  xmlNode * n = (xmlNode *)node->next;
  while (n != NULL)
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

xmlNode * SiconosDOMTreeTools::findFollowNode(const xmlNode * node)
{
  xmlNode * n = (xmlNode *)node->next;
  while (n != NULL)
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


int SiconosDOMTreeTools::findNextfigureIndex(const string s, const int start)
{
  int res = -1;
  int cpt = start;
  int len = s.length();

  // Exception : out of range
  if (start >= len)
    XMLException::selfThrow("SiconosDOMTreeTools - findNextBlankIndex : the index given in parameter is greather than the size of the string");

  while ((cpt < len) && (res == -1))
  {
    if (s[cpt] == ' ')
      cpt++;
    else res = cpt;
  }
  return res;
}


int SiconosDOMTreeTools::findNextBlankIndex(const string s, const int start)
{
  int res = -1;
  int cpt = start;
  int len = s.length();

  // Exception : out of range
  if (start >= len)
    XMLException::selfThrow("SiconosDOMTreeTools - findNextBlankIndex : the index given in parameter is greather than the size of the string");

  while ((cpt < len) && (res == -1))
  {
    if (s[cpt] != ' ')
      cpt++;
    else res = cpt;
  }
  return res;
}


vector<double> SiconosDOMTreeTools::string2Vector(string vectorContent, const int size)
{
  vector <double> vect(size);
  int start = 0, nb = 0, end = 0, nestStart = 0;;
  string stmp;

  // suppresses tabs and line breaks.
  for (int i = 0; i < vectorContent.length(); i++)
  {
    if (vectorContent[i] == '\n' || vectorContent[i] == '\t')
    {
      vectorContent[i] = ' ';
    }
  }

  while ((nb < size) && (end != -1))
  {
    start = findNextfigureIndex(vectorContent, end);
    end = findNextBlankIndex(vectorContent, start);
    stmp = vectorContent.substr(start, end - start);
    vect[nb] = atof(stmp.c_str());
    nb++;
  }
  return vect;
}


vector<int> SiconosDOMTreeTools::string2Vector(string vectorContent)
{
  vector <int> vect;
  int start = 0, nb = 0, end = 0, nestStart = 0, size;
  string stmp;

  // suppresses tabs and line breaks.
  for (int i = 0; i < vectorContent.length(); i++)
  {
    if (vectorContent[i] == '\n' || vectorContent[i] == '\t')
    {
      vectorContent[i] = ' ';
    }
  }

  while (end != -1)
  {
    start = findNextfigureIndex(vectorContent, end);
    end = findNextBlankIndex(vectorContent, start);
    stmp = vectorContent.substr(start, end - start);
    vect.push_back(atoi(stmp.c_str()));
    nb++;
  }
  return vect;
}


SimpleVector SiconosDOMTreeTools::getSiconosRowMatrixValue(const xmlNode * siconosMatrixRowNode, const int colSize)
{
  if (xmlHasProp((xmlNodePtr) siconosMatrixRowNode, (xmlChar *)SDTT_VECTORFILE.c_str())) //row is defined in a extern ascii file
  {
    SimpleVector v(getStringAttributeValue(siconosMatrixRowNode, SDTT_VECTORFILE), true);
    return v;
  }

  //The row is precised in the XML DOM Tree
  //Content
  string vectorContent = (char *)xmlNodeGetContent((xmlNode *) siconosMatrixRowNode);

  SimpleVector v(string2Vector(vectorContent.c_str(), colSize));
  return v;
}


void SiconosDOMTreeTools::setSiconosRowMatrixValue(const xmlNode * siconosMatrixRowNode, const SiconosVector &v, int colSize)
{
  if (colSize != v.size())
    XMLException::selfThrow("SiconosDOMTreeTools - setSiconosRowMatrixValue : a row size you want to save is different of the matrix size defined in xml");

  char element[100];
  int i = 0, end = 0;
  string vectorContent = "";
  while (i < v.size())
  {
    strcpy(element, "");
    sprintf(element, /*DOUBLE_PRECISION*/ N_DOUBLE_PRECISION, v(i));
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

int SiconosDOMTreeTools::getNodeChildrenNumber(const xmlNode *node)
{
  int res = 0;
  xmlNode *n;

  if ((node == NULL) && (node->type != XML_ELEMENT_NODE))
  {
    res = -1;
  }
  else
  {
    n = SiconosDOMTreeTools::findNodeChild((const xmlNode*) node);
    //if( n != NULL ) res++;

    while (n != NULL)
    {
      n = findFollowNode(n);
      res++;
    }
  }
  return res;
}
