#include <string>
#include <vector>


using namespace std;

class BaseObjectSiconos
{
public:
  BaseObjectSiconos();
  BaseObjectSiconos(int ii);
  void display();
  void fillObjectWithObjectXML();
  void linkObjectXML();
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
  ObjectSiconos(ObjectSiconosXML *oxml);
  void display();
  void fillObjectWithObjectXML();
  void linkObjectXML();

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

