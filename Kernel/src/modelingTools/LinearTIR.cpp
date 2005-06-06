#include "LinearTIR.h"
using namespace std;


LinearTIR::LinearTIR(): Relation(), C(NULL), D(NULL), E(NULL), a(NULL),
  isCAllocatedIn(false), isDAllocatedIn(false), isEAllocatedIn(false),
  isAAllocatedIn(false)
{
  relationType = LINEARTIRELATION;
}

LinearTIR::LinearTIR(RelationXML* relxml): Relation(relxml), C(NULL), D(NULL), E(NULL), a(NULL),
  isCAllocatedIn(true), isDAllocatedIn(true), isEAllocatedIn(true),
  isAAllocatedIn(true)
{
  relationType = LINEARTIRELATION;
  if (relationxml != NULL)
  {
    int row = ((static_cast<LinearTIRXML*>(relationxml))->getC()).size(0);
    int col = ((static_cast<LinearTIRXML*>(relationxml))->getC()).size(1);
    C = new SiconosMatrix(row, col);
    *C = (static_cast<LinearTIRXML*>(relationxml))->getC();
    row = ((static_cast<LinearTIRXML*>(relationxml))->getD()).size(0);
    col = ((static_cast<LinearTIRXML*>(relationxml))->getD()).size(1);
    D = new SiconosMatrix(row, col);
    *D = (static_cast<LinearTIRXML*>(relationxml))->getD();
    row = ((static_cast<LinearTIRXML*>(relationxml))->getE()).size(0);
    col = ((static_cast<LinearTIRXML*>(relationxml))->getE()).size(1);
    E = new SiconosMatrix(row, col);
    *E = (static_cast<LinearTIRXML*>(relationxml))->getE();
    a = new SimpleVector(((static_cast<LinearTIRXML*>(relationxml))->getA()).size());
    *a = (static_cast<LinearTIRXML*>(relationxml))->getA();
  }
  else RuntimeException::selfThrow("LinearTIR::xml constructor, xml file=NULL");
}

LinearTIR::LinearTIR(SiconosMatrix* newC, SiconosMatrix* newD,
                     SiconosMatrix* newE, SiconosVector* newA):
  Relation(), C(NULL), D(NULL), E(NULL), a(NULL), isCAllocatedIn(true),
  isDAllocatedIn(true), isEAllocatedIn(true), isAAllocatedIn(true)
{
  relationType = LINEARTIRELATION;
  C = new SiconosMatrix(newC->size(0), newC->size(1));
  D = new SiconosMatrix(newD->size(0), newD->size(1));
  E = new SiconosMatrix(newE->size(0), newE->size(1));
  a = new SimpleVector(newA->size());
  *C = *newC;
  *D = *newD;
  *E = *newE;
  *a = *newA;
}

LinearTIR::~LinearTIR()
{
  if (isCAllocatedIn) delete C;
  C = NULL;
  if (isDAllocatedIn) delete D;
  D = NULL;
  if (isEAllocatedIn) delete E;
  E = NULL;
  if (isAAllocatedIn) delete a;
  a = NULL;
}

void LinearTIR::computeOutput()
{
  IN("LinearTIR::computeOutput\n");

  vector<DynamicalSystem*> vDS = interaction->getDynamicalSystems();
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

void LinearTIR::display() const
{
  cout << "---------------------------------------------------" << endl;
  cout << "____ data of the LinearTIR " << endl;
  cout << "| C " << endl;
  if (C != NULL) C->display();
  else cout << "->NULL" << endl;
  cout << "| D " << endl;
  if (D != NULL) D->display();
  else cout << "->NULL" << endl;
  cout << "| E " << endl;
  if (E != NULL) E->display();
  else cout << "->NULL" << endl;
  cout << "| a " << endl;
  if (a != NULL) a->display();
  else cout << "->NULL" << endl;
  cout << "____________________________" << endl;
  cout << "---------------------------------------------------" << endl;
}

void LinearTIR::saveRelationToXML()
{
  OUT("LinearTIR::saveRelationToXML\n");
  if (relationxml != NULL)
  {
    (static_cast<LinearTIRXML*>(relationxml))->setC(*C);
    (static_cast<LinearTIRXML*>(relationxml))->setD(*D);
    (static_cast<LinearTIRXML*>(relationxml))->setE(*E);
    (static_cast<LinearTIRXML*>(relationxml))->setA(*a);
  }
}

LinearTIR* LinearTIR::convert(Relation *r)
{
  cout << "LinearTIR::convert (Relation *r)" << endl;
  LinearTIR* ltir = dynamic_cast<LinearTIR*>(r);
  return ltir;
}

