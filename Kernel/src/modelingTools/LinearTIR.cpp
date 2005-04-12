
#include "LinearTIR.h"
#include "check.h"

LinearTIR::LinearTIR(): Relation()
{
  this->relationType = LINEARTIRELATION;//"LinearTIR";
}

LinearTIR::LinearTIR(RelationXML* relxml): Relation(relxml)
{
  this->relationType = LINEARTIRELATION;//"LinearTIR";
}

LinearTIR::~LinearTIR()
{}


SiconosMatrix* LinearTIR::getCPtr(void)
{
  return &this->C;
}

SiconosMatrix* LinearTIR::getDPtr(void)
{
  return &this->D;
}

SiconosMatrix* LinearTIR::getEPtr(void)
{
  return &this->E;
}

SiconosVector* LinearTIR::getAPtr(void)
{
  return &this->a;
}


void LinearTIR::computeOutput()
{
  IN("LinearTIR::computeOutput\n");

  vector<DynamicalSystem*> vDS = this->interaction->getDynamicalSystems();
  if (vDS.size() == 2)
  {
    /*
     *\WARNING to do with CompositeVector
     */

    //    SiconosVector x(*(vDS[0]->getXPtr()), false);
    //    x.add(*( vDS[1]->getXPtr()));
  }
  else if (vDS.size() == 1)
  {
    /*
     *\WARNING to do with CompositeVector
     */
    //    SiconosVector x(*(vDS[0]->getXPtr()), false);
  }
  else RuntimeException::selfThrow("The interaction doesn't contain the right number of Dynamical Systems");

  /*
   * NOT TERMINATED... SEE LAGRANGIANLINEARRELATION
   */

  RuntimeException::selfThrow("LinearTIR::computeOutput not yet implemented");


  OUT("LinearTIR::computeOutput\n");
}



void LinearTIR::fillRelationWithRelationXML()
{
  Relation::fillRelationWithRelationXML();
  OUT("LinearTIR::fillRelationWithRelationXML\n");
  if (this->relationxml != NULL)
  {
    this->C = (static_cast<LinearTIRXML*>(this->relationxml))->getC();
    this->D = (static_cast<LinearTIRXML*>(this->relationxml))->getD();
    this->E = (static_cast<LinearTIRXML*>(this->relationxml))->getE();
    this->a = (static_cast<LinearTIRXML*>(this->relationxml))->getA();
  }
}

void LinearTIR::display() const
{
  cout << "---------------------------------------------------" << endl;
  cout << "____ data of the LinearTIR " << endl;
  cout << "| C " << endl;
  this->C.display();
  cout << "| D " << endl;
  this->D.display();
  cout << "| E " << endl;
  this->E.display();
  cout << "| a " << endl;
  this->a.display();
  cout << "____________________________" << endl;
  cout << "---------------------------------------------------" << endl;
}

void LinearTIR::saveRelationToXML()
{
  Relation::saveRelationToXML();
  OUT("LinearTIR::saveRelationToXML\n");
  if (this->relationxml != NULL)
  {
    //    this->display();

    (static_cast<LinearTIRXML*>(this->relationxml))->setC(&(this->C));
    (static_cast<LinearTIRXML*>(this->relationxml))->setD(&(this->D));
    (static_cast<LinearTIRXML*>(this->relationxml))->setE(&(this->E));
    (static_cast<LinearTIRXML*>(this->relationxml))->setA(&(this->a));
  }
}

void LinearTIR::createRelation(LinearTIRXML * relationXML,
                               SiconosMatrix* C, SiconosMatrix* D,
                               SiconosMatrix* E, SiconosVector* a)//, Interaction * interaction)
{
  if (relationXML != NULL)
  {
    this->init();
    this->relationxml = relationXML;
    this->relationType = LINEARTIRELATION;//"LinearTIR";
    this->fillRelationWithRelationXML();
  }
  else
  {
    this->C = *C;
    this->D = *D;
    this->E = *E;
    this->a = *a;
  }
}


LinearTIR* LinearTIR::convert(Relation *r)
{
  cout << "LinearTIR::convert (Relation *r)" << endl;
  LinearTIR* ltir = dynamic_cast<LinearTIR*>(r);
  return ltir;
}

