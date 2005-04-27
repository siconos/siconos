#include "SiconosObject.h"
#include <iostream>

ObjectSiconos::ObjectSiconos(ObjectSiconosXML *objectsiconosxml):
  // Construtors of the Base Class
  BaseObjectSiconos(),
  // Constructors of the Built-in types
  type("OBJETSICONOS_TYPE"), att(),
  // Constructors of the object members (composition)
  compositionObject()
{
  cout << "ObjectSiconos -- Constructor with minimal data" << endl;
  // ObjectSiconos(); // Dangerous Fix the order ??
  // Initialization of pointers on external objects

  objectsiconosxml = oxml;
  // Loading of the attribute from the ObjectXML
  fillObjectWithObjectXML();
  // Create from factory if needed (new) and  link downwards the External objects
  linkObjectXML();

}
