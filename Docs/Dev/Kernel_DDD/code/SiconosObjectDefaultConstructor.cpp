#include "SiconosObject.h"
#include <iostream>


ObjectSiconos::ObjectSiconos():
  // Construtors of the Base Class
  BaseObjectSiconos(),
  // Constructors of the Built-in types
  type("OBJECTSICONOS_TYPE"), att(0),
  // Constructors of the object members (composition)
  compositionObject()
{
  cout << "ObjectSiconos -- default Constructor" << endl;

  // Initialization of pointers on external objects
  objectsiconosxml = 0;
  externalObject = 0;

  // Do we need to call the new operator ?
  // externalObject = new ExternalObject();
}

BaseObjectSiconos::BaseObjectSiconos(): i(0), type("BASE_TYPE")
{
  cout << "BaseObjectSiconos -- default Constructor" << endl;
}

CompositionObject::CompositionObject(): i(0)
{
};

ExternalObject::ExternalObject(): i(0)
{
}

