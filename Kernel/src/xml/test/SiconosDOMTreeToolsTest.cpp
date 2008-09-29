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
#include "SiconosDOMTreeToolsTest.h"
#include <float.h>

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// on place cette classe de test dans le registry
CPPUNIT_TEST_SUITE_REGISTRATION(SiconosDOMTreeToolsTest);




void SiconosDOMTreeToolsTest::setUp()
{
  doc = NULL;
  root = NULL;
  child = NULL;
}

void SiconosDOMTreeToolsTest::tearDown()
{
  //cout<<"tearDown" << endl;
  //  xmlFreeDoc(doc);
  xmlCleanupParser();
  /*  doc = NULL;
    root = NULL;
    child = NULL;*/
}


//______________________________________________________________________________


void SiconosDOMTreeToolsTest::testGetStringContentValue()
{
  doc = xmlParseFile("string.xml");
  root = xmlDocGetRootElement(doc);
  child = root->children;
  child = child->next;

  string s = SiconosDOMTreeTools::getStringContentValue(child);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetStringContentValue : s", s == "cette chaine a un string", true);
  cout << "SiconosDOMTreeToolsTest >>> testGetStringContentValue ............... OK\n ";
}


void SiconosDOMTreeToolsTest::testGetIntegerContentValue()
{
  doc = xmlParseFile("int.xml");
  root = xmlDocGetRootElement(doc);
  child = root->children;
  child = child->next;

  int i = SiconosDOMTreeTools::getIntegerContentValue(child);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetIntegerContentValue : i", i == 254, true);
  cout << "SiconosDOMTreeToolsTest >>> testGetIntegerContentValue .............. OK\n ";
}

void SiconosDOMTreeToolsTest::testGetDoubleContentValue()
{
  doc = xmlParseFile("double.xml");
  root = xmlDocGetRootElement(doc);
  child = root->children;
  child = child->next;

  double d = SiconosDOMTreeTools::getDoubleContentValue(child);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetDoubleContentValue : d", d == 0.054, true);
  cout << "SiconosDOMTreeToolsTest >>> testGetDoubleContentValue ............... OK\n ";
}


void SiconosDOMTreeToolsTest::testGetStringAttributeValue()
{
  doc = xmlParseFile("string_attribute.xml");
  root = xmlDocGetRootElement(doc);
  child = root->children;
  child = child->next;

  string s = SiconosDOMTreeTools::getStringAttributeValue(child, "toto");
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetStringAttributeValue : s", s == "cette chaine a un string", true);
  cout << "SiconosDOMTreeToolsTest >>> testGetStringAttributeValue ............. OK\n ";
}

void SiconosDOMTreeToolsTest::testGetDoubleAttributeValue()
{
  doc = xmlParseFile("double_attribute.xml");
  root = xmlDocGetRootElement(doc);
  child = root->children;
  child = child->next;

  double d = SiconosDOMTreeTools::getDoubleAttributeValue(child, "toto");
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetDoubleAttributeValue : d", d == 0.054, true);
  cout << "SiconosDOMTreeToolsTest >>> testGetDoubleAttributeValue ............. OK\n ";
}


void SiconosDOMTreeToolsTest::testGetIntegerAttributeValue()
{
  doc = xmlParseFile("int_attribute.xml");
  root = xmlDocGetRootElement(doc);
  child = root->children;
  child = child->next;

  int i = SiconosDOMTreeTools::getIntegerAttributeValue(child, "toto");
  CPPUNIT_ASSERT_EQUAL_MESSAGE("getIntegerAttributeValue : d", i == 254, true);
  cout << "SiconosDOMTreeToolsTest >>> getIntegerAttributeValue ................ OK\n ";
}


void SiconosDOMTreeToolsTest::testGetBooleanAttributeValue()
{
  doc = xmlParseFile("bool_attribute.xml");
  root = xmlDocGetRootElement(doc);
  child = root->children;
  child = child->next;

  bool b = SiconosDOMTreeTools::getBooleanAttributeValue(child, "toto");
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetBooleanAttributeValue : b", b == false, true);
  cout << "SiconosDOMTreeToolsTest >>> testGetBooleanAttributeValue ............ OK\n ";
}



void SiconosDOMTreeToolsTest::testGetSiconosVectorValue()
{
  doc = xmlParseFile("vector.xml");
  root = xmlDocGetRootElement(doc);
  child = root->children;
  child = child->next;

  /*SP::SiconosVector/SimpleVector v;
  v = SiconosDOMTreeTools::getSiconosVectorValue(child);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetSiconosVectorValue : v", v(3) == 4, true);

  cout<<"SiconosDOMTreeToolsTest >>> testGetSiconosVectorValue ............... OK\n ";
  }


  void SiconosDOMTreeToolsTest::testGetSiconosMatrixValue()
  {
  doc = xmlParseFile("matrix.xml");
  root = xmlDocGetRootElement(doc);
  child = root->children;
  child = child->next;

  SiconosMatrix m;
  m = SiconosDOMTreeTools::getSiconosMatrixValue(child);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetSiconosMatrixValue : m", m(0, 2) == 3, true);
  cout<<"SiconosDOMTreeToolsTest >>> testGetSiconosMatrixValue ............... OK\n ";
  }

  void SiconosDOMTreeToolsTest::testGetSiconosMatrixValue2()
  {
  doc = xmlParseFile("matrix2.xml");
  root = xmlDocGetRootElement(doc);
  child = root->children;
  child = child->next;

  SiconosMatrix m;
  m = SiconosDOMTreeTools::getSiconosMatrixValue(child);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetSiconosMatrixValue2 : m", m(0, 2) == 3, true);
  cout<<"SiconosDOMTreeToolsTest >>> testGetSiconosMatrixValue2 .............. OK\n ";
  }

  //void SiconosDOMTreeToolsTest::testGetVectorMemoryValue()
  //{
  // doc = xmlParseFile("testMemory.xml");
  // root = xmlDocGetRootElement(doc);
  // child = root->children;
  // child = child->next;
  //
  // vector<SP::SiconosVector> vect;
  // vect = SiconosDOMTreeTools::getVectorMemoryValue(child);
  //
  // SiconosVector v("testVector.dat", true);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetVectorMemoryValue size", vect.size() == 2, true);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetVectorMemoryValue size", *(vect[0]) == v, true);
  //
  // cout<<"SiconosDOMTreeToolsTest >>> testGetVectorMemoryValue ................ OK\n ";
  //}


  void SiconosDOMTreeToolsTest::testSetStringContentValue()
  {
  doc = xmlParseFile("string.xml");
  root = xmlDocGetRootElement(doc);
  child = root->children;
  child = child->next;

  string s = SiconosDOMTreeTools::getStringContentValue(child);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetStringContentValue : s", s == "cette chaine a un string", true);

  string newstring = "week-end be heroes de ta race !";
  SiconosDOMTreeTools::setStringContentValue(child, newstring);
  string s1 = SiconosDOMTreeTools::getStringContentValue(child);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetStringContentValue : newstring", s1 == newstring, true);


  cout<<"SiconosDOMTreeToolsTest >>> testSetStringContentValue ............... OK\n ";
  }


  void SiconosDOMTreeToolsTest::testSetIntegerContentValue()
  {
  doc = xmlParseFile("int.xml");
  root = xmlDocGetRootElement(doc);
  child = root->children;
  child = child->next;

  int i = SiconosDOMTreeTools::getIntegerContentValue(child);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetIntegerContentValue : i", i == 254, true);

  int newint = 18071980;
  SiconosDOMTreeTools::setIntegerContentValue(child, newint);

  i = SiconosDOMTreeTools::getIntegerContentValue(child);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetIntegerContentValue : i", i == newint, true);


  cout<<"SiconosDOMTreeToolsTest >>> testGetIntegerContentValue .............. OK\n ";
  }

  void SiconosDOMTreeToolsTest::testSetDoubleContentValue()
  {
  doc = xmlParseFile("double.xml");
  root = xmlDocGetRootElement(doc);
  child = root->children;
  child = child->next;

  double d = SiconosDOMTreeTools::getDoubleContentValue(child);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetDoubleContentValue : d", d == 0.054, true);

  double newdouble = 1E-3;
  SiconosDOMTreeTools::setDoubleContentValue(child, newdouble);

  d = SiconosDOMTreeTools::getDoubleContentValue(child);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetDoubleContentValue : d", d == newdouble, true);


  cout<<"SiconosDOMTreeToolsTest >>> testSetDoubleContentValue ............... OK\n ";
  }


  void SiconosDOMTreeToolsTest::testSetStringAttributeValue()
  {
  doc = xmlParseFile("string_attribute.xml");
  root = xmlDocGetRootElement(doc);
  child = root->children;
  child = child->next;

  string s = SiconosDOMTreeTools::getStringAttributeValue(child, "toto");

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetStringAttributeValue : s", s == "cette chaine a un string", true);

  string newstring = "week-end be heroes de ta race !";
  SiconosDOMTreeTools::setStringAttributeValue(child, "toto", newstring);

  s = SiconosDOMTreeTools::getStringAttributeValue(child, "toto");
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetStringAttributeValue : s", s == newstring, true);

  cout<<"SiconosDOMTreeToolsTest >>> testSetStringAttributeValue ............. OK\n ";
  }


  void SiconosDOMTreeToolsTest::testSetDoubleAttributeValue()
  {
  doc = xmlParseFile("double_attribute.xml");
  root = xmlDocGetRootElement(doc);
  child = root->children;
  child = child->next;

  double d = SiconosDOMTreeTools::getDoubleAttributeValue(child, "toto");
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetDoubleAttributeValue : d", d == 0.054, true);

  double newdouble = 1E-3;
  SiconosDOMTreeTools::setDoubleAttributeValue(child,"toto", newdouble);

  d = SiconosDOMTreeTools::getDoubleAttributeValue(child, "toto");

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetDoubleContentValue : d", d == newdouble, true);

  cout<<"SiconosDOMTreeToolsTest >>> testSetDoubleAttributeValue ............. OK\n ";
  }


  void SiconosDOMTreeToolsTest::testSetIntegerAttributeValue()
  {
  doc = xmlParseFile("int_attribute.xml");
  root = xmlDocGetRootElement(doc);
  child = root->children;
  child = child->next;

  int i = SiconosDOMTreeTools::getIntegerAttributeValue(child, "toto");
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetIntegerAttributeValue : d", i == 254, true);

  int newint = 18071980;
  SiconosDOMTreeTools::setIntegerAttributeValue(child, "toto", newint);

  i = SiconosDOMTreeTools::getIntegerAttributeValue(child, "toto");
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetIntegerContentValue : i", i == newint, true);


  cout<<"SiconosDOMTreeToolsTest >>> testSetIntegerAttributeValue ............ OK\n ";
  }


  void SiconosDOMTreeToolsTest::testSetSiconosVectorValue()
  {
  doc = xmlParseFile("./vector.xml");
  root = xmlDocGetRootElement(doc);
  child = root->children;
  child = child->next;

  /*SiconosVector*/
  SimpleVector v1;
  v1 = SiconosDOMTreeTools::getSiconosVectorValue(child);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetSiconosVectorValue : v1", v1(3) == 4, true);

  v1(3) = 3.1415;
  SiconosDOMTreeTools::setSiconosVectorNodeValue(child, v1 /*, 5*/);
  /*SiconosVector*/
  SimpleVector vv;
  vv = SiconosDOMTreeTools::getSiconosVectorValue(child);


  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetSiconosVectorValue : vv", v1 == vv, true);

  v1(3) = 4;
  SiconosDOMTreeTools::setSiconosVectorNodeValue(child, v1 /*, 5*/);

  cout << "SiconosDOMTreeToolsTest >>> testSetSiconosVectorValue ............... OK\n ";
}


void SiconosDOMTreeToolsTest::testSetSiconosMatrixValue()
{
  doc = xmlParseFile("matrix.xml");
  root = xmlDocGetRootElement(doc);
  child = root->children;
  child = child->next;

  SiconosMatrix m;
  m = SiconosDOMTreeTools::getSiconosMatrixValue(child);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetSiconosMatrixValue : m", m(0, 2) == 3, true);

  double newVal = -0.2548;
  m(0, 2) = newVal;
  SiconosDOMTreeTools::setSiconosMatrixNodeValue(child, m);

  SiconosMatrix m2;
  m2 = SiconosDOMTreeTools::getSiconosMatrixValue(child);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetSiconosMatrixValue : m2", m2(0, 2) == newVal, true);
  cout << "SiconosDOMTreeToolsTest >>> testSetSiconosMatrixValue ............... OK\n ";
}


void SiconosDOMTreeToolsTest::testSetBooleanAttributeValue()
{
  doc = xmlParseFile("bool_attribute.xml");
  root = xmlDocGetRootElement(doc);
  child = root->children;
  child = child->next;
  bool b = SiconosDOMTreeTools::getBooleanAttributeValue(child, "toto");
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetBooleanAttributeValue : b", b == false, true);

  bool newb = true;
  SiconosDOMTreeTools::setBooleanAttributeValue(child, "toto", newb);

  b = SiconosDOMTreeTools::getBooleanAttributeValue(child, "toto");
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetBooleanAttributeValue : b", b == newb, true);

  cout << "SiconosDOMTreeToolsTest >>> testSetBooleanAttributeValue ............ OK\n ";
}


//void SiconosDOMTreeToolsTest::testSetVectorMemoryValue()
//{
//  // modifier le premier
//
//  doc = xmlParseFile("testMemory.xml");
//  root = xmlDocGetRootElement(doc);
//  child = root->children;
//  child = child->next;
//
//  int i;
//
//  vector<double> vd(10);
//  for (i = 0; i < 10; i++) vd[i] = i+1;
//  SiconosVector *v = new SiconosVector(vd);
//
//  vector<SiconosVector*> vect;
//  vect = SiconosDOMTreeTools::getVectorMemoryValue(child);
////
////  //cout<<"PRECISION : "<<DBL_MAX_EXP<<" "<<DBL_MANT_DIG<<" "<<FLT_MANT_DIG<<endl;
////  cout<<"testSetVectorMemoryValue vect.size = "<<vect.size()<<endl;
////  delete vect[0];
////  vect[0] = v;
////  SiconosDOMTreeTools::setVectorMemoryValue(child, vect);
////  const vector<SiconosVector*> vect2 = SiconosDOMTreeTools::getVectorMemoryValue(child);
////  cout <<"v "<<(*v);
//
//  cout<<"SiconosDOMTreeToolsTest >>> testSetVectorMemoryValue ................ OK\n ";
//  v->write("testVector.dat");
//  delete v;
//}

void SiconosDOMTreeToolsTest::testCreateMatrixNode()
{
  doc = xmlParseFile("matrix.xml");
  root = xmlDocGetRootElement(doc);
  child = root->children;
  child = child->next;

  SiconosMatrix mat("matrix.dat", true);
  xmlNode* node = SiconosDOMTreeTools::createMatrixNode(root, "TestMatrix", mat);

  SiconosMatrix matRes;
  matRes = SiconosDOMTreeTools::getSiconosMatrixValue(node);

  CPPUNIT_ASSERT_EQUAL(mat.size(0) == matRes.size(0), true);
  CPPUNIT_ASSERT_EQUAL(mat.size(1) == matRes.size(1), true);
  CPPUNIT_ASSERT_EQUAL(mat == matRes, true);

  cout << "SiconosDOMTreeToolsTest >>> testCreateMatrixNode ............... OK\n ";
}

void SiconosDOMTreeToolsTest::testCreateVectorNode()
{
  doc = xmlParseFile("vector.xml");
  root = xmlDocGetRootElement(doc);
  child = root->children;
  child = child->next;

  /*SiconosVector*/
  SimpleVector vect("vector.dat", true);

  xmlNode* node = SiconosDOMTreeTools::createVectorNode(root, "TestVector", vect);

  /*SiconosVector*/
  SimpleVector vectRes;
  vectRes = SiconosDOMTreeTools::getSiconosVectorValue(node);
  CPPUNIT_ASSERT_EQUAL(vect.size(), vectRes.size());
  CPPUNIT_ASSERT_EQUAL(vect == vectRes, true);

  cout << "SiconosDOMTreeToolsTest >>> testCreateVectorNode ............... OK\n ";
}

void SiconosDOMTreeToolsTest::testCreateDoubleNode()
{
  doc = xmlParseFile("vector.xml");
  root = xmlDocGetRootElement(doc);
  child = root->children;
  child = child->next;

  double d = 3.1415;

  xmlNode* node = SiconosDOMTreeTools::createDoubleNode(root, "TestDouble", d);

  double res = SiconosDOMTreeTools::getDoubleContentValue(node);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testCreateDoubleNode : d", res == d, true);
  cout << "SiconosDOMTreeToolsTest >>> testCreateDoubleNode ............... OK\n ";

  //  xmlSaveFile("totti.xml", doc);
}

void SiconosDOMTreeToolsTest::testCreateIntegerNode()
{
  doc = xmlParseFile("vector.xml");
  root = xmlDocGetRootElement(doc);
  child = root->children;
  child = child->next;

  int i = 3;

  xmlNode* node = SiconosDOMTreeTools::createIntegerNode(root, "TestInteger", i);

  int res = SiconosDOMTreeTools::getIntegerContentValue(node);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testCreateIntegerNode : d", res == i, true);
  cout << "SiconosDOMTreeToolsTest >>> testCreateIntegerNode ............... OK\n ";
}

void SiconosDOMTreeToolsTest::testCreateStringNode()
{
  doc = xmlParseFile("vector.xml");
  root = xmlDocGetRootElement(doc);
  child = root->children;
  child = child->next;

  //string s = "\"totti\"";
  string s = "totti";

  xmlNode* node = SiconosDOMTreeTools::createStringNode(root, "TestString", s);

  string res = SiconosDOMTreeTools::getStringContentValue(node);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testCreateStringNode : d", res == s, true);
  cout << "SiconosDOMTreeToolsTest >>> testCreateStringNode ............... OK\n ";
}

void SiconosDOMTreeToolsTest::testCreateBooleanNode()
{
  doc = xmlParseFile("vector.xml");
  root = xmlDocGetRootElement(doc);
  child = root->children;
  child = child->next;

  bool b = 3.1415;

  xmlNode* node = SiconosDOMTreeTools::createBooleanNode(root, "TestBoolean", b);

  bool res = SiconosDOMTreeTools::getBooleanContentValue(node);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testCreateBooleanNode : d", res == b, true);
  cout << "SiconosDOMTreeToolsTest >>> testCreateBooleanNode ............... OK\n ";
}

//-------------------------------

void SiconosDOMTreeToolsTest::testGetNodeChildrenNumber()
{
  doc = xmlParseFile("testMemory.xml");
  root = xmlDocGetRootElement(doc);
  child = root->children;
  child = child->next;

  int nb = SiconosDOMTreeTools::getNodeChildrenNumber(root);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetNodeChildrenNumber : nb", nb == 2, true);
  cout << "SiconosDOMTreeToolsTest >>> testGetNodeChildrenNumber ............... OK\n ";
}

void SiconosDOMTreeToolsTest::testGetStringAttributeValueException()
{
  doc = xmlParseFile("string_attribute.xml");
  root = xmlDocGetRootElement(doc);
  child = root->children;
  child = child->next;

  string s = SiconosDOMTreeTools::getStringAttributeValue(child, "badAttribute"); // this attribute does not exist
}


void SiconosDOMTreeToolsTest::testGetIntegerAttributeValueException()
{
  doc = xmlParseFile("int_attribute.xml");
  root = xmlDocGetRootElement(doc);
  child = root->children;
  child = child->next;

  int i ;
  i = SiconosDOMTreeTools::getIntegerAttributeValue(child, "badAttribute");
}

void SiconosDOMTreeToolsTest::testGetDoubleAttributeValueException()
{
  doc = xmlParseFile("double_attribute.xml");
  root = xmlDocGetRootElement(doc);
  child = root->children;
  child = child->next;

  double d;
  d = SiconosDOMTreeTools::getDoubleAttributeValue(child, "badAttribute");
}


void SiconosDOMTreeToolsTest::testGetBooleanAttributeValueException()
{
  doc = xmlParseFile("bool_attribute.xml");
  root = xmlDocGetRootElement(doc);
  child = root->children;
  child = child->next;

  bool b;
  b = SiconosDOMTreeTools::getBooleanAttributeValue(child, "badAttribute");
}


void SiconosDOMTreeToolsTest::testSetStringAttributeValueException()
{
  doc = xmlParseFile("string_attribute.xml");
  root = xmlDocGetRootElement(doc);
  child = root->children;
  child = child->next;

  string newstring = "week-end be heroes de ta race !";
  SiconosDOMTreeTools::setStringAttributeValue(child, "BadAttribute", newstring);

}


void SiconosDOMTreeToolsTest::testSetDoubleAttributeValueException()
{
  doc = xmlParseFile("double_attribute.xml");
  root = xmlDocGetRootElement(doc);
  child = root->children;
  child = child->next;

  double newdouble = 1E-3;
  SiconosDOMTreeTools::setDoubleAttributeValue(child, "BadAttribute", newdouble);

}


void SiconosDOMTreeToolsTest::testSetIntegerAttributeValueException()
{
  doc = xmlParseFile("int_attribute.xml");
  root = xmlDocGetRootElement(doc);
  child = root->children;
  child = child->next;

  int newint = 18071980;
  SiconosDOMTreeTools::setIntegerAttributeValue(child, "BadAttribute", newint);
}


void SiconosDOMTreeToolsTest::testSetBooleanAttributeValueException()
{
  doc = xmlParseFile("bool_attribute.xml");
  root = xmlDocGetRootElement(doc);
  child = root->children;
  child = child->next;

  bool newb = true;
  SiconosDOMTreeTools::setBooleanAttributeValue(child, "BadAttribute", newb);
}


void SiconosDOMTreeToolsTest::testGetVectorValueBadSizeException()
{
  // in vector_bad_size, the attribute "size" = 3, and the vector contains 5 values.
  doc = xmlParseFile("vector_bad_size.xml");
  root = xmlDocGetRootElement(doc);
  child = root->children;
  child = child->next;

  /*SiconosVector*/
  SimpleVector v;
  v = SiconosDOMTreeTools::getSiconosVectorValue(child);
}

//void SiconosDOMTreeToolsTest::testString2Vector()
//{
//  string s = " ";
//  vector<double> V;
//
//  V = SiconosDOMTreeTools::string2Vector(s, 0);
//
//  cout<<"SiconosDOMTreeToolsTest >>> testString2Vector ............... OK\n ";
//}

