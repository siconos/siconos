
















#include "SiconosObject.h"
#include <iostream>

ObjectSiconos::ObjectSiconos(int aatt, string ttype): BaseObjectSiconos(), type(ttype), att(aatt), compositionObject()
{
  cout << "ObjectSiconos -- Constructor with minimal data" << endl;
  ObjectSiconos(); // Fix the order ??
  externalObject = new ExternalObject(aatt);
}

CompositionObject::CompositionObject(int ii)
{
  i = ii;
};


ExternalObject::ExternalObject(int ii)
{
  i = ii;
}

