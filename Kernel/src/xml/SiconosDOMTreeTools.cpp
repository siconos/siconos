
#include "SiconosDOMTreeTools.h"

#include "check.h"

//#include "KernelDefaultConfig.h"


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
    //    OUT("SiconosDOMTreeTools::setSiconosMatrixValue - to Matrix file\n");
    string file = getStringAttributeValue(siconosMatrixNode, SDTT_MATRIXFILE);
    matrix.write(file, /*ASCII*/ /*N_ASCII*/ FILE_STORAGE); //For the moment only ASCII file are managed
  }
  else
  {
    //    OUT("SiconosDOMTreeTools::setSiconosMatrixValue - to XML file\n");
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

    //cout<<"#####################"<<endl;
    //matrix->display();
    //cout<<"   ---"<<endl;
    //cout<<matrixColSize<<" - "<<matrixRowSize<<endl;
    //cout<<"/#####################"<<endl;

    if (matrixColSize != matrix->size(1))
    {
      //XMLException::selfThrow("SiconosDOMTreeTools - setSiconosMatrixValue : the " + matrixName + " matrix col size you want to save is different of the col size defined for it in xml");
      cout << "SiconosDOMTreeTools - setSiconosMatrixValue : the " + matrixName + " matrix col size you want to save is different of the col size already defined." << endl;
      setIntegerAttributeValue(siconosMatrixNode, SDTT_MATRIXCOLSIZE, matrix->size(1));
    }

    if (matrixRowSize != matrix->size(0))
    {
      //XMLException::selfThrow("SiconosDOMTreeTools - setSiconosMatrixValue : the " + matrixName + " matrix row size you want to save is different of the row size defined for it in xml");
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


//void SiconosDOMTreeTools::setVectorMemoryValue(const xmlNode * memoryNode, const vector<SiconosVector*> memory)
//{
//    xmlNode *oldNode=SiconosDOMTreeTools::findNodeChild(memoryNode, SDTT_MEMORY);
//  xmlNode *node;
//  string stringValue;
//  stringstream sstr;
//  node=oldNode;
//
//  int i=0;
//
//  while ((node!=NULL)&&(i<memory.size()))
//  {
//    setSiconosVectorValue(node, *(memory[i]));
//    oldNode=node;
//    node=SiconosDOMTreeTools::findFollowNode(node, SDTT_MEMORY);
//    i++;
//  }
//
//  while (i<memory.size()) //not enought nodes in the DOM tree to save memory
//  {
//    node = xmlNewNode(NULL, BAD_CAST SDTT_MEMORY.c_str());
//    //sstr << *memory[i];
//    //sstr >> stringValue;
//    stringValue = (*(memory[i])).toString();
//    xmlNodeSetContent(node, (const xmlChar *)stringValue.c_str());
//    oldNode->next=node;
//    i++;
//  }
//}
//
//void SiconosDOMTreeTools::setMemoryValue( xmlNode * memoryNode, const SiconosMemory & memory)
//{
//  int i=0;
//  xmlNode *oldNode;
//  xmlNode *node, *parent;
//  string stringValue;
//
//  if( memoryNode != NULL )
//  {
//    oldNode = SiconosDOMTreeTools::findNodeChild(memoryNode, SDTT_MEMORY);
//    node=oldNode;
//
//    while ( (node!=NULL)&&(i<memory.getNbVectorsInMemory()) )
//    {
//      setSiconosVectorValue( node, memory.getSiconosVector(i) );
//      oldNode = node;
//      node = SiconosDOMTreeTools::findFollowNode(node, SDTT_MEMORY);
//      i++;
//    }
//
//    while ( (i<memory.getMemorySize()) && (i<memory.getNbVectorsInMemory()))  //not enought nodes in the DOM tree to save memory
//    {
//      node = xmlNewNode(NULL, BAD_CAST SDTT_MEMORY.c_str());
//      stringValue = memory.getSiconosVector(i)->toString();
//      xmlNodeSetContent(node, (const xmlChar *)stringValue.c_str());
//      xmlAddChild( memoryNode, node );
////      oldNode->next = node;
////      oldNode = node;
//      i++;
//    }
//    setIntegerAttributeValue( memoryNode, SDTT_MEMORYSIZE, memory.getMemorySize() );
//  }
//  else
//  {
//    XMLException::selfThrow("SiconosDOMTreeTools - setMemoryValue : memoryNode == NULL");
//  }
//
//
//}


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
    //    rowNode = xmlNewChild(node, NULL, BAD_CAST SDTT_ROW.c_str(), NULL);
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

    //    v.write(file, /*ASCII*/ /*N_ASCII*/ DefaultFileStorage); //For the moment only ASCII file are managed

    setSiconosVectorValue(node, v);

    xmlAddChildList(rootNode, node);
  }

  return node;
}

xmlNode* SiconosDOMTreeTools::createVectorIntNode(xmlNode* rootNode, const string name, vector<int> v)
{
  xmlNode *node;
  node = xmlNewNode(NULL, BAD_CAST name.c_str());

  //  string size;
  //  stringstream sstr;
  //
  //  sstr << v.size();
  //  sstr >> size;

  setVectorIntContentValue(node, v);

  xmlAddChildList(rootNode, node);

  return node;
}

//xmlNode* SiconosDOMTreeTools::createVectorMemoryNode(xmlNode* rootNode, const string name, vector<SiconosVector*> vect)
//{
//  /*
//   * \todo if the SiconosVectors are too big, the SiconosVectors must be saved in an external files and only the name of these files must be written in the XML file
//   */
//  xmlNode *node;
//
//    /* \todo */
//
//  xmlAddChildList(rootNode, node);
//  return node;
//}

//xmlNode* SiconosDOMTreeTools::createSiconosMemoryNode(xmlNode* rootNode, const string name, SiconosMemory* smem)
//{
//  /*
//   * \todo if the SiconosVectors are too big, the SiconosVectors must be saved in an external files and only the name of these files must be written in the XML file
//   */
//  xmlNode *node;
//
//  /* \todo */
//
//  xmlAddChildList(rootNode, node);
//  return node;
//}

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

    //    n = node->next;
    //    while( n != NULL )
    //    {
    //      if (n->type == XML_ELEMENT_NODE) res++;
    //      n = n->next;
    //    }
  }
  return res;
}

//$Log: SiconosDOMTreeTools.cpp,v $
//Revision 1.58  2005/03/23 15:03:56  jbarbier
//- adaptation to the LMGC90 tags in non smooth dynamical system and strategy
//
//Revision 1.57  2005/03/22 15:55:05  jbarbier
//- class NewtonImpactFriction non smooth law added to the kernel
//
//- xml schema modified for this new class
//- xml schema modified to accept a "joker" for further use of a LMGC90 mechanical plugin
//
//- new test added for the loading/saving of a NewtonImpactFrictionNSL
//
//Revision 1.56  2005/02/15 15:15:33  charlety
//
//_ modified some very slow functions to increase performance
//
//Revision 1.55  2005/01/20 09:05:34  jbarbier
//- configuration file available and usable
//
//- save of vectors and matrices into external files (version 0.1)
//
//Revision 1.54  2005/01/10 17:06:37  jbarbier
//- attribute "size" is now unused in the code
//
//- xml schema v1.2 is in progress
//
//Revision 1.53  2004/12/08 12:49:39  jbarbier
//- changes in the XML Schema, respect of the recommandations of the W3C
//version 1.1
//
//- changes in all balises DS, Relation, NSLaw, OneStepIntegrator, OneStepNSProblem
//in the XML files into specific names like LagrangianNLDS, LinearSystemDS, ...
//for the DS
//
//Revision 1.52  2004/12/06 10:10:34  jbarbier
//- integration of Numerics and use of Numerics on the bouncing ball sample
//
//- Numerics is now needed to run the bouncing ball sample!
//
//Revision 1.51  2004/09/21 11:49:10  jbarbier
//- correction in the XML save for a manual construction of the platform :
//    DS_Concerned of the Interaction
//    DS_Concerned of the Integrator
//
//- test updated for these changes
//
//Revision 1.50  2004/09/14 13:24:54  charlety
//
//_ changes in the interface of SiconosVector
//
//Revision 1.49  2004/09/10 11:26:29  charlety
//
//_ Integration of the new version of the SiconosVector in the platform. the class simpleVector is used mostly to replace old SiconosVector. When a vector can be composite or simple, like the state of a dynamical system, a pointer on SiconosVector is used, and the vector is initialized simple or composite when the system is initialized.
//
//_ All the tests which worked with the previous version of the vector are OK with the new version.
//
//_ Example SICONOS and bouncingBall are OK
//
//_ some comments have still to be adapted to NewSiconosVector .
//
//_ This version of NewSiconosVector could be called 0.9. some details have to be fixed, it will be done before the end of September.
//
//Revision 1.48  2004/08/23 14:30:03  jbarbier
//- All the dynamical systems can be created in a comand program and added to a
//NSDS. The save is OK, but the creation of the boundary conditions is not yet
//finished.
//
//Revision 1.47  2004/08/20 07:34:23  jbarbier
//- creation of Model, NSDS in comand program succeed in creating SiconosModelXML,
//NSDSXML
//
//Revision 1.46  2004/08/11 14:43:45  jbarbier
//- beginning of the mechanism of creation without XML input file of the objects of the platform with the
//creatObjtect methods
//
//- function saveWToXML for Moreau integrator, and same specific functions to save
//M,q and Q,p for LCP and QP
//
//- function to check coherency of the Model
//
//Revision 1.45  2004/08/10 14:51:49  jbarbier
//- functions initialize() of the Lsodar and Adams OneStepIntegrator completed to
//call the function initialize() of the base class
//
//Revision 1.44  2004/07/30 14:37:15  jbarbier
//- saving methods for DynamicalSystemXML and LagrangianNLDSXML
//
//Revision 1.43  2004/07/29 14:25:45  jbarbier
