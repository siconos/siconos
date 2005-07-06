#include "LinearTIR.h"
using namespace std;


// Default constructor
LinearTIR::LinearTIR():
  Relation(), C(NULL), D(NULL), F(NULL), e(NULL), B(NULL), a(NULL)
{
  relationType = LINEARTIRELATION;
}

// xml constructor
LinearTIR::LinearTIR(RelationXML* relxml):
  Relation(relxml), C(NULL), D(NULL), F(NULL), e(NULL), B(NULL), a(NULL)
{
  relationType = LINEARTIRELATION;
  if (relationxml != NULL)
  {
    LinearTIRXML * lTIRxml = static_cast<LinearTIRXML *>(relationxml);
    isAllocatedIn.resize(6, false);

    // get size of vector y
    unsigned int size = interaction->getNInteraction();
    if (lTIRxml->hasC())
    {
      C = new SiconosMatrix(size, size);
      isAllocatedIn[0] = true;
      *C = lTIRxml->getC();
    }
    else
      RuntimeException::selfThrow("LinearTIR:: xml constructor: input matrix C is missing in xml file ");

    if (lTIRxml->hasD())
    {
      D = new SiconosMatrix(size, size);
      isAllocatedIn[1] = true;
      *D = lTIRxml->getD();
    }
    if (lTIRxml->hasF())
    {
      unsigned int uSize = lTIRxml->getF().size(1);
      F = new SiconosMatrix(size, uSize);
      isAllocatedIn[2] = true;
      *F = lTIRxml->getF();
    }
    if (lTIRxml->hasE())
    {
      e = new SimpleVector(size);
      isAllocatedIn[3] = true;
      *e = lTIRxml->getE();
    }
    if (lTIRxml->hasB())
    {
      B = new SiconosMatrix(size, size);
      isAllocatedIn[4] = true;
      *B = lTIRxml->getB();
    }
    else
      RuntimeException::selfThrow("LinearTIR:: xml constructor: input matrix B is missing in xml file ");

    if (lTIRxml->hasA())
    {
      a = new SimpleVector(size);
      isAllocatedIn[5] = true;
      *a = lTIRxml->getA();
    }
  }
  else RuntimeException::selfThrow("LinearTIR::xml constructor, xml file=NULL");
}

// Minimum data (C, B) constructor
LinearTIR::LinearTIR(const SiconosMatrix& newC, const SiconosMatrix& newB):
  Relation(), C(NULL), D(NULL), F(NULL), e(NULL), B(NULL), a(NULL)
{
  relationType = LINEARTIRELATION;
  // get size of vector y
  unsigned int size = interaction->getNInteraction();
  if (newC.size(0) != size || newC.size(1) != size || newB.size(0) != size || newB.size(1) != size)
    RuntimeException::selfThrow("LinearTIR:: constructor from data, inconsistent size with y vector for input vector or matrix");

  C = new SiconosMatrix(newC.size(0), newC.size(1));
  B = new SiconosMatrix(newB.size(0), newB.size(1));
  *C = newC;
  *B = newB;
  isAllocatedIn.resize(6, false);
  isAllocatedIn[0] = true;
  isAllocatedIn[4] = true;
}

// Constructor from a complete set of data
LinearTIR::LinearTIR(const SiconosMatrix& newC, const SiconosMatrix& newD,
                     const SiconosMatrix& newF, const SimpleVector& newE,
                     const SiconosMatrix& newB, const SimpleVector& newA):
  Relation(), C(NULL), D(NULL), F(NULL), e(NULL), B(NULL), a(NULL)
{
  relationType = LINEARTIRELATION;

  // get size of vector y
  unsigned int size = interaction->getNInteraction();

  if (newC.size(0) != size || newC.size(1) != size || newD.size(0) != size || newD.size(1) != size || newF.size(0) != size
      || newE.size() != size || newB.size(0) != size || newB.size(1) != size || newA.size() != size)
    RuntimeException::selfThrow("LinearTIR:: constructor from data, inconsistent size with y vector for input vector or matrix");

  C = new SiconosMatrix(size, size);
  D = new SiconosMatrix(size, size);
  // \todo check newF.size(1) = size of u and that the ds type is well adapted to this kind of relation
  F = new SiconosMatrix(size, newF.size(1));
  e = new SimpleVector(size);
  B = new SiconosMatrix(size, size);
  a = new SimpleVector(size);
  *C = newC;
  *D = newD;
  *F = newF;
  *e = newE;
  *B = newB;
  *a = newA;
  isAllocatedIn.resize(6, true);
}

LinearTIR::~LinearTIR()
{
  if (isAllocatedIn[0])
  {
    delete C;
    C = NULL;
  }
  if (isAllocatedIn[1])
  {
    delete D;
    D = NULL;
  }
  if (isAllocatedIn[2])
  {
    delete F;
    F = NULL;
  }
  if (isAllocatedIn[3])
  {
    delete e;
    e = NULL;
  }
  if (isAllocatedIn[4])
  {
    delete B;
    B = NULL;
  }
  if (isAllocatedIn[5])
  {
    delete a;
    a = NULL;
  }
}

// setters

void LinearTIR::setC(const SiconosMatrix& newValue)
{
  unsigned int n = interaction->getNInteraction();
  if (newValue.size(0) != n || newValue.size(1) != n)
    RuntimeException::selfThrow("LinearTIR - setC: inconsistent dimensions with problem size for input matrix C");

  if (C == NULL)
  {
    C = new SiconosMatrix(n, n);
    isAllocatedIn[0] = true;
  }
  *C = newValue;
}

void LinearTIR::setCPtr(SiconosMatrix *newPtr)
{
  if (isAllocatedIn[0]) delete C;
  C = newPtr;
  isAllocatedIn[0] = false;
}

void LinearTIR::setD(const SiconosMatrix& newValue)
{
  unsigned int n = interaction->getNInteraction();
  if (newValue.size(0) != n || newValue.size(1) != n)
    RuntimeException::selfThrow("LinearTIR - setD: inconsistent dimensions with problem size for input matrix D");

  if (D == NULL)
  {
    D = new SiconosMatrix(n, n);
    isAllocatedIn[1] = true;
  }
  *D = newValue;
}

void LinearTIR::setDPtr(SiconosMatrix *newPtr)
{
  if (isAllocatedIn[1])  delete D;
  D = newPtr;
  isAllocatedIn[1] = false;
}

void LinearTIR::setF(const SiconosMatrix& newValue)
{
  unsigned int n = interaction->getNInteraction();
  if (newValue.size(0) != n)
    RuntimeException::selfThrow("LinearTIR - setF: inconsistent dimensions with problem size for input matrix F");

  unsigned int uSize = newValue.size(1);
  if (F == NULL)
  {
    F = new SiconosMatrix(n, uSize);
    isAllocatedIn[2] = true;
  }
  else if (uSize != F->size(1))
    RuntimeException::selfThrow("LinearTIR - setF: inconsistent dimensions with problem size for input matrix F");

  *F = newValue;
}

void LinearTIR::setFPtr(SiconosMatrix *newPtr)
{
  if (isAllocatedIn[2]) delete F;
  F = newPtr;
  isAllocatedIn[2] = false;
}

void LinearTIR::setE(const SimpleVector& newValue)
{
  unsigned int n = interaction->getNInteraction();
  if (newValue.size() != n)
    RuntimeException::selfThrow("LinearTIR - setE: inconsistent dimensions with problem size for input vector e");

  if (e == NULL)
  {
    e = new SimpleVector(n);
    isAllocatedIn[3] = true;
  }
  *e = newValue;
}

void LinearTIR::setEPtr(SimpleVector *newPtr)
{
  if (isAllocatedIn[3]) delete e;
  e = newPtr;
  isAllocatedIn[3] = false;
}

void LinearTIR::setB(const SiconosMatrix& newValue)
{
  unsigned int n = interaction->getNInteraction();
  if (newValue.size(0) != n || newValue.size(1) != n)
    RuntimeException::selfThrow("LinearTIR - setB: inconsistent dimensions with problem size for input matrix B");

  if (B == NULL)
  {
    B = new SiconosMatrix(n, n);
    isAllocatedIn[4] = true;
  }
  *B = newValue;
}

void LinearTIR::setBPtr(SiconosMatrix *newPtr)
{
  if (isAllocatedIn[4]) delete B;
  B = newPtr;
  isAllocatedIn[4] = false;
}

void LinearTIR::setA(const SimpleVector& newValue)
{
  unsigned int n = interaction->getNInteraction();
  if (newValue.size() != n)
    RuntimeException::selfThrow("LinearTIR - setA: inconsistent dimensions with problem size for input vector a");

  if (a == NULL)
  {
    a = new SimpleVector(n);
    isAllocatedIn[5] = true;
  }
  *a = newValue;
}

void LinearTIR::setAPtr(SimpleVector *newPtr)
{
  if (isAllocatedIn[5]) delete a;
  a = newPtr;
  isAllocatedIn[5] = false;
}

void LinearTIR::computeOutput()
{
  IN("LinearTIR::computeOutput\n");
  /*
    vector<DynamicalSystem*> vDS = interaction->getDynamicalSystems();
    if (vDS.size() == 2)
    {
    \todo
    }
    else if (vDS.size() == 1)
    {
    \todo
    }
    else RuntimeException::selfThrow("The interaction doesn't contain the right number of Dynamical Systems");
  */
  RuntimeException::selfThrow("LinearTIR::computeOutput not yet implemented");
  OUT("LinearTIR::computeOutput\n");
}

void LinearTIR::computeInput()
{
  RuntimeException::selfThrow("LinearTIR::computeInput not yet implemented");
}

void LinearTIR::display() const
{
  cout << " ===== Linear Time Invariant relation display ===== " << endl;
  cout << "| C " << endl;
  if (C != NULL) C->display();
  else cout << "->NULL" << endl;
  cout << "| D " << endl;
  if (D != NULL) D->display();
  else cout << "->NULL" << endl;
  cout << "| F " << endl;
  if (F != NULL) F->display();
  else cout << "->NULL" << endl;
  cout << "| e " << endl;
  if (e != NULL) e->display();
  else cout << "->NULL" << endl;
  cout << "| B " << endl;
  if (B != NULL) B->display();
  else cout << "->NULL" << endl;
  cout << "| a " << endl;
  if (a != NULL) a->display();
  else cout << "->NULL" << endl;
  cout << " ================================================== " << endl;
}

void LinearTIR::saveRelationToXML()
{
  OUT("LinearTIR::saveRelationToXML\n");
  if (relationxml != NULL)
  {
    LinearTIRXML * lTIRxml = (static_cast<LinearTIRXML*>(relationxml));
    lTIRxml->setC(*C);
    lTIRxml->setD(*D);
    lTIRxml->setF(*F);
    lTIRxml->setE(*e);
    lTIRxml->setB(*B);
    lTIRxml->setA(*a);
  }
}

LinearTIR* LinearTIR::convert(Relation *r)
{
  cout << "LinearTIR::convert (Relation *r)" << endl;
  LinearTIR* ltir = dynamic_cast<LinearTIR*>(r);
  return ltir;
}

