#include "SiconosObject.h"
#include <iostream>

ObjectSiconos::ObjectSiconos(int aatt, string ttype):
  // Construtors of the Base Class
  BaseObjectSiconos(2 * aatt),
  // Constructors of the Built-in types
  type(ttype), att(aatt),
  // Constructors of the object members (composition)
  compositionObject(5 * aatt)
{
  cout << "ObjectSiconos -- Constructor with minimal data" << endl;
  // ObjectSiconos(); // Dangerous Fix the order ??
  // Initialization of pointers on external objects
  objectsiconosxml = 0;
  externalObject = new ExternalObject(aatt);
}

CompositionObject::CompositionObject(int ii): i(ii)
{
};


ExternalObject::ExternalObject(int ii): i(ii)
{
}

BaseObjectSiconos::BaseObjectSiconos(int ii): i(ii), type("BASE_TYPE")
{
  cout << "BaseObjectSiconos -- default Constructor" << endl;
}
