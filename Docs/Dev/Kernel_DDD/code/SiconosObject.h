#include <string>
#include <vector>


using namespace std;

class BaseObjectSiconos
{
public:
  BaseObjectSiconos();
  void display();
private :
  int i;
  string type;
};


class CompositionObject
{
public:
  CompositionObject();
  CompositionObject(int ii);
private :
  int i;
};

class ExternalObject
{
public:
  ExternalObject();
  ExternalObject(int ii);
private :
  int i;
};

class ObjectSiconosXML
{
public :
  ObjectSiconosXML();

private :
  // Built-in types attributes (int, char, float, ..)
  string type;
  int att;
};


class ObjectSiconos : BaseObjectSiconos
{
public :
  ObjectSiconos();
  ObjectSiconos(int aatt, string ttype);
  void display();

private :
  // Built-in types attributes (int, char, float, ..)
  string type;
  int att;


  // Objects members (Composition)
  CompositionObject compositionObject;


  // Pointer on external objects
  ObjectSiconosXML *objectsiconosxml;
  ExternalObject *externalObject;

};

