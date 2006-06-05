/* Siconos-Kernel version 1.2.0, Copyright INRIA 2005-2006.
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
#include "LinearTIR.h"
using namespace std;


// xml constructor
LinearTIR::LinearTIR(RelationXML* relxml, Interaction * inter):
  Relation(relxml, inter), C(NULL), D(NULL), F(NULL), e(NULL), B(NULL)
{
  relationType = LINEARTIRELATION;
  if (relationxml != NULL)
  {
    LinearTIRXML * lTIRxml = static_cast<LinearTIRXML *>(relationxml);
    isAllocatedIn.resize(5, false);
    unsigned int sizeY, sizeX, size; // size of vector y and of vector x

    if (inter != NULL)
    {
      // get size of vector y from linked interaction
      size = interaction->getInteractionSize();
    }
    // === Output ===

    // get matrices values.
    // The only required variables are C and B, other are optional.
    // If ouput plug-in is provided ...
    if (!isOutputPlugged) // ie connected to default plug-in
    {
      if (lTIRxml->hasC())
      {
        sizeY = lTIRxml->getC().size(0);
        sizeX = lTIRxml->getC().size(1);
      }
      else
        RuntimeException::selfThrow("LinearTIR:: xml constructor: input matrix C is missing in xml file ");

      if (inter != NULL && size != sizeY)
        RuntimeException::selfThrow("LinearTIR:: xml constructor, inconsistent size between C and y vector");

      C = new SimpleMatrix(lTIRxml->getC());
      isAllocatedIn[0] = true;

      if (lTIRxml->hasD())
      {
        if (lTIRxml->getD().size(0) != sizeY || lTIRxml->getD().size(0) != sizeY)
          RuntimeException::selfThrow("LinearTIR:: xml constructor, inconsistent size between D and C");
        D = new SimpleMatrix(lTIRxml->getD());
        isAllocatedIn[1] = true;
      }
      if (lTIRxml->hasF())
      {
        if (lTIRxml->getF().size(0) != sizeY)
          RuntimeException::selfThrow("LinearTIR:: xml constructor, inconsistent size between F and C");
        F = new SimpleMatrix(lTIRxml->getF());
        isAllocatedIn[2] = true;
      }
      if (lTIRxml->hasE())

      {
        if (lTIRxml->getE().size() != sizeY)
          RuntimeException::selfThrow("LinearTIR:: xml constructor, inconsistent size between e and C");
        e = new SimpleVector(lTIRxml->getE());
        isAllocatedIn[3] = true;
      }
    }
    else if (lTIRxml->hasC() || lTIRxml->hasD() || lTIRxml->hasF() || lTIRxml->hasE())
      cout << "Warning: LinearTIR xml constructor, you give plug-in function and matrices values for output definition -> conflict. Plug-in will be used." << endl;

    // === Input (lambda/R) ===

    if (!isInputPlugged) // ie input connected to default plug-in
    {
      if (lTIRxml->hasB())
      {
        unsigned int sizeB = lTIRxml->getB().size(1);
        if (inter != NULL && size != sizeB)
          RuntimeException::selfThrow("LinearTIR:: xml constructor, inconsistent size between B and y vector");

        if (lTIRxml->getB().size(0) != sizeX || (C != NULL && sizeB != C->size(0)))
          RuntimeException::selfThrow("LinearTIR:: xml constructor, inconsistent size between B and ds vector");
        B = new SimpleMatrix(lTIRxml->getB());
        isAllocatedIn[4] = true;
      }
      else
        RuntimeException::selfThrow("LinearTIR:: xml constructor: input matrix B is missing in xml file ");
    }
    else if (lTIRxml->hasB())
      cout << "Warning: LinearTIR xml constructor, you give plug-in function and matrices values for input definition -> conflict. Plug-in will be used." << endl;
  }
  else RuntimeException::selfThrow("LinearTIR::xml constructor, xml file=NULL");
}

// Minimum data (C, B) constructor
LinearTIR::LinearTIR(const SiconosMatrix& newC, const SiconosMatrix& newB, Interaction* inter):
  Relation(inter), C(NULL), D(NULL), F(NULL), e(NULL), B(NULL)
{
  isOutputPlugged = false;
  isInputPlugged  = false;
  relationType = LINEARTIRELATION;
  unsigned int sizeY, sizeX;
  sizeY = newC.size(0);
  sizeX = newC.size(1);
  if (newB.size(0) != sizeX || newB.size(1) != sizeY)
    RuntimeException::selfThrow("LinearTIR:: constructor from data, inconsistent size between C and B");

  if (inter != NULL)
  {
    // get size of vector y
    unsigned int size = interaction->getInteractionSize();
    if (sizeY != size)
      RuntimeException::selfThrow("LinearTIR:: constructor from data, inconsistent size with y vector for input vector or matrix");
  }

  C = new SimpleMatrix(sizeY, sizeX);
  B = new SimpleMatrix(sizeX, sizeY);
  *C = newC;
  *B = newB;
  isAllocatedIn.resize(5, false);
  isAllocatedIn[0] = true;
  isAllocatedIn[4] = true;
}

// Constructor from a complete set of data
LinearTIR::LinearTIR(const SiconosMatrix& newC, const SiconosMatrix& newD,
                     const SiconosMatrix& newF, const SimpleVector& newE,
                     const SiconosMatrix& newB, Interaction* inter):
  Relation(inter), C(NULL), D(NULL), F(NULL), e(NULL), B(NULL)
{
  relationType = LINEARTIRELATION;
  isOutputPlugged = false;
  isInputPlugged  = false;
  unsigned int sizeY, sizeX; // size of vector y and of vector x

  sizeY = newC.size(0);
  sizeX = newC.size(1);

  if (inter != NULL)
  {
    unsigned int size = interaction->getInteractionSize();
    if (size != sizeY)
      RuntimeException::selfThrow("LinearTIR:: constructor from data, inconsistent size between C and y vector");
  }

  C = new SimpleMatrix(sizeY, sizeX);
  *C = newC;

  if (newD.size(0) != sizeY || newD.size(1) != sizeY)
    RuntimeException::selfThrow("LinearTIR:: constructor from data, inconsistent size between C and D");
  D = new SimpleMatrix(sizeY, sizeY);
  *D = newD;

  // \todo check newF.size(1) = size of u and that the ds type is well adapted to this kind of relation
  if (newF.size(0) != sizeY)
    RuntimeException::selfThrow("LinearTIR:: constructor from data, inconsistent size between F and C");
  F = new SimpleMatrix(sizeY, newF.size(1));
  *F = newF;

  if (newE.size(0) != sizeY)
    RuntimeException::selfThrow("LinearTIR:: constructor from data, inconsistent size between e and C");
  e = new SimpleVector(sizeY);
  *e = newE;

  if (newB.size(0) != sizeX || newB.size(1) != sizeY)
    RuntimeException::selfThrow("LinearTIR:: constructor from data, inconsistent size between C and B");
  B = new SimpleMatrix(sizeX, sizeY);
  *B = newB;

  isAllocatedIn.resize(5, true);
}

// Copy constructor (inter is optional)
LinearTIR::LinearTIR(const Relation & newLTIR, Interaction* inter):
  Relation(newLTIR, inter), C(NULL), D(NULL), F(NULL), e(NULL), B(NULL)
{
  if (relationType != LINEARTIRELATION)
    RuntimeException::selfThrow("LinearTIR:: copy constructor, inconsistent relation types for copy");

  // Warning: the interaction link is not copyed!!!

  const LinearTIR * ltir = static_cast<const LinearTIR*>(&newLTIR);
  isAllocatedIn.resize(5, false);

  // Since this is a copy, we suppose that various sizes of members of newLTIR are consistent alltogether
  // -> no more tests on that subject.

  // === Output ===
  if (!isOutputPlugged)
  {
    C = new SimpleMatrix(ltir->getC());
    isAllocatedIn[0] = true; // C
    if (ltir->getDPtr() != NULL)
    {
      D = new SimpleMatrix(ltir->getD());
      isAllocatedIn[1] = true;
    }
    if (ltir->getFPtr() != NULL)
    {
      F = new SimpleMatrix(ltir->getF());
      isAllocatedIn[2] = true;
    }
    if (ltir->getEPtr() != NULL)
    {
      e = new SimpleVector(ltir->getE());
      isAllocatedIn[3] = true;
    }
  }
  else
    cout << "Warning: LinearTIR copy constructor, original relations uses plug-in function for output definition." << endl;

  // === Input ===
  if (!isInputPlugged)
  {
    B = new SimpleMatrix(ltir->getB());
    isAllocatedIn[4] = true; // B
  }
  else
    cout << "Warning: LinearTIR copy constructor, original relations uses plug-in function for input definition." << endl;

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
}

// setters

void LinearTIR::setC(const SiconosMatrix& newValue)
{
  isOutputPlugged = false;
  unsigned int sizeY;
  if (interaction != NULL)
  {
    sizeY = interaction->getInteractionSize();
    if (newValue.size(0) != sizeY)
      RuntimeException::selfThrow("LinearTIR - setC: inconsistent dimensions with problem size for input matrix C");
  }

  if (C == NULL)
  {
    C = new SimpleMatrix(newValue);
    isAllocatedIn[0] = true;
  }
  else
  {
    if (newValue.size(1) == C->size(1))
      *C = newValue;
    else
      RuntimeException::selfThrow("LinearTIR - setC: inconsistent dimensions with problem size for input matrix C");
  }
}

void LinearTIR::setCPtr(SiconosMatrix *newPtr)
{
  isOutputPlugged = false;
  if (isAllocatedIn[0]) delete C;
  if (interaction != NULL)
  {
    unsigned int sizeY = interaction->getInteractionSize();
    if (newPtr->size(0) != sizeY)
      RuntimeException::selfThrow("LinearTIR - setCPtr: inconsistent dimensions with problem size for input matrix C");
  }

  C = newPtr;
  isAllocatedIn[0] = false;
}

void LinearTIR::setD(const SiconosMatrix& newValue)
{
  if (isOutputPlugged)
    cout << " LinearTIR:setD warning: before this set, output was a plug-in, do not forget to set other required operators (C ...)" << endl;
  isOutputPlugged = false;
  unsigned int sizeY = newValue.size(0);
  if (sizeY != newValue.size(1))
    RuntimeException::selfThrow("LinearTIR - setD:  D is not square!!");

  if (interaction != NULL)
  {
    unsigned int size = interaction->getInteractionSize();
    if (size != sizeY)
      RuntimeException::selfThrow("LinearTIR - setD: inconsistent dimensions with problem size for input matrix D");
  }

  if (D == NULL)
  {
    D = new SimpleMatrix(newValue);
    isAllocatedIn[1] = true;
  }
  else
  {
    if (sizeY == D->size(0))
      *D = newValue;
    else
      RuntimeException::selfThrow("LinearTIR - setD: inconsistent dimensions with problem size for input matrix D");
  }
}

void LinearTIR::setDPtr(SiconosMatrix *newPtr)
{
  if (isOutputPlugged)
    cout << " LinearTIR:setDPtr warning: before this set, output was a plug-in, do not forget to set other required operators (C ...)" << endl;
  isOutputPlugged = false;
  if (isAllocatedIn[1])  delete D;
  if (interaction != NULL)
  {
    unsigned int sizeY = interaction->getInteractionSize();
    if (newPtr->size(0) != sizeY || newPtr->size(1) != sizeY)
      RuntimeException::selfThrow("LinearTIR - setDPtr: inconsistent dimensions with problem size for input matrix D");
  }
  D = newPtr;
  isAllocatedIn[1] = false;
}

void LinearTIR::setF(const SiconosMatrix& newValue)
{
  if (isOutputPlugged)
    cout << " LinearTIR:setF warning: before this set, output was a plug-in, do not forget to set other required operators (C ...)" << endl;
  isOutputPlugged = false;
  unsigned int sizeY = newValue.size(0);
  if (interaction != NULL)
  {
    unsigned int size = interaction->getInteractionSize();
    if (size != sizeY)
      RuntimeException::selfThrow("LinearTIR - setF: inconsistent dimensions with problem size for input matrix F");
  }

  if (F == NULL)
  {
    F = new SimpleMatrix(newValue);
    isAllocatedIn[2] = true;
  }
  else
  {
    if (newValue.size(1) == F->size(1))
      *F = newValue;
    else
      RuntimeException::selfThrow("LinearTIR - setF: inconsistent dimensions with problem size for input matrix F");
  }
}

void LinearTIR::setFPtr(SiconosMatrix *newPtr)
{
  if (isOutputPlugged)
    cout << " LinearTIR:setFPtr warning: before this set, output was a plug-in, do not forget to set other required operators (C ...)" << endl;
  isOutputPlugged = false;
  if (isAllocatedIn[2]) delete F;
  if (interaction != NULL)
  {
    unsigned int sizeY = interaction->getInteractionSize();
    if (newPtr->size(0) != sizeY)
      RuntimeException::selfThrow("LinearTIR - setFPtr: inconsistent dimensions with problem size for input matrix F");
  }
  F = newPtr;
  isAllocatedIn[2] = false;
}

void LinearTIR::setE(const SimpleVector& newValue)
{
  if (isOutputPlugged)
    cout << " LinearTIR:setE warning: before this set, output was a plug-in, do not forget to set other required operators (C ...)" << endl;
  isOutputPlugged = false;
  unsigned int sizeY = newValue.size();
  if (interaction != NULL)
  {
    unsigned int size = interaction->getInteractionSize();
    if (size != sizeY)
      RuntimeException::selfThrow("LinearTIR - setE: inconsistent dimensions with problem size for input vector e");
  }

  if (e == NULL)
  {
    e = new SimpleVector(newValue);
    isAllocatedIn[3] = true;
  }
  else
  {
    if (sizeY == e->size())
      *e = newValue;
    else
      RuntimeException::selfThrow("LinearTIR - setE: inconsistent dimensions with problem size for input vector e");
  }
}

void LinearTIR::setEPtr(SimpleVector *newPtr)
{
  if (isOutputPlugged)
    cout << " LinearTIR:setEPtr warning: before this set, output was a plug-in, do not forget to set other required operators (C ...)" << endl;
  isOutputPlugged = false;
  if (isAllocatedIn[3]) delete e;
  if (interaction != NULL)
  {
    unsigned int sizeY = interaction->getInteractionSize();
    if (newPtr->size() != sizeY)
      RuntimeException::selfThrow("LinearTIR - setEPtr: inconsistent dimensions with problem size for input matrix E");
  }
  e = newPtr;
  isAllocatedIn[3] = false;
}

void LinearTIR::setB(const SiconosMatrix& newValue)
{
  isInputPlugged = false;
  unsigned int sizeY = newValue.size(1);
  unsigned int sizeX = newValue.size(0);

  if (interaction != NULL)
  {
    unsigned int size = interaction->getInteractionSize();
    if (size != sizeY)
      RuntimeException::selfThrow("LinearTIR - setB: inconsistent dimensions with problem size for input matrix B");
  }

  if (B == NULL)
  {
    B = new SimpleMatrix(newValue);
    isAllocatedIn[4] = true;
  }
  else
  {
    if (sizeX == B->size(0))
      *B = newValue;
    else
      RuntimeException::selfThrow("LinearTIR - setB: inconsistent dimensions with problem size for input matrix B");
  }
}

void LinearTIR::setBPtr(SiconosMatrix *newPtr)
{
  isInputPlugged = false;
  if (isAllocatedIn[4]) delete B;
  if (interaction != NULL)
  {
    unsigned int sizeY = interaction->getInteractionSize();
    if (newPtr->size(1) != sizeY)
      RuntimeException::selfThrow("LinearTIR - setBPtr: inconsistent dimensions with problem size for input matrix B");
  }
  B = newPtr;
  isAllocatedIn[4] = false;
}

void LinearTIR::getCBlockDSPtr(DynamicalSystem * ds, SiconosMatrix& CBlock) const
{
  unsigned int k = 0;

  DSSet vDS = interaction ->getDynamicalSystems();
  DSIterator itDS = vDS.begin();

  // look for ds
  while (*itDS != ds && itDS != vDS.end())
  {
    k += (*itDS)->getN();
    itDS++;
  }

  // check dimension
  if ((*itDS)->getN() != CBlock.size(1))
    RuntimeException::selfThrow("LinearTIR - getCBlockDSPtr: inconsistent sizes between CBlock and DS");

  // get block
  unsigned int l = k + (*itDS)->getN() - 1;
  vector<unsigned int> index_list(4);
  index_list[0] = 0;
  index_list[1] = C->size(0) - 1;
  index_list[2] = k;
  index_list[3] = l;
  C->getBlock(index_list, CBlock);
}

void LinearTIR::getCBlockDSPtr(const int& DSNumber, SiconosMatrix& CBlock) const
{
  unsigned int k = 0;

  DSSet vDS = interaction ->getDynamicalSystems();
  DSIterator itDS = vDS.begin();

  // look for DS number DSNumber ...
  while ((*itDS)->getNumber() != DSNumber && itDS != vDS.end())
  {
    k += (*itDS)->getN();
    itDS++;
  }

  // check dimension
  if ((*itDS)->getN() != CBlock.size(1))
    RuntimeException::selfThrow("LinearTIR - getCBlockDSPtr: inconsistent sizes between CBlock and DS");

  // get block
  unsigned int l = k + (*itDS)->getN() - 1;
  vector<unsigned int> index_list(4);
  index_list[0] = 0;
  index_list[1] = C->size(0) - 1;
  index_list[2] = k;
  index_list[3] = l;
  C->getBlock(index_list, CBlock);
}

void LinearTIR::getBBlockDSPtr(DynamicalSystem* ds, SiconosMatrix& BBlock) const
{
  unsigned int k = 0;

  DSSet vDS = interaction ->getDynamicalSystems();
  DSIterator itDS = vDS.begin();

  while ((*itDS) != ds && itDS != vDS.end())
  {
    k += (*itDS)->getN();
    itDS++;
  }
  // check dimension
  if ((*itDS)->getN() != BBlock.size(0))
    RuntimeException::selfThrow("LinearTIR - getBBlockDSPtr: inconsistent sizes between BBlock and DS");

  // get block
  unsigned int l = k + (*itDS)->getN() - 1;
  vector<unsigned int> index_list(4);
  index_list[0] = k;
  index_list[1] = l;
  index_list[2] = 0;
  index_list[3] = B->size(1) - 1;

  B->getBlock(index_list, BBlock);
}

void LinearTIR::getBBlockDSPtr(const int& DSNumber, SiconosMatrix& BBlock) const
{
  unsigned int k = 0;

  DSSet vDS = interaction ->getDynamicalSystems();
  DSIterator itDS = vDS.begin();

  while ((*itDS)->getNumber() != DSNumber && itDS != vDS.end())
  {
    k += (*itDS)->getN();
    itDS++;
  }
  // check dimension
  if ((*itDS)->getN() != BBlock.size(0))
    RuntimeException::selfThrow("LinearTIR - getBBlockDSPtr: inconsistent sizes between BBlock and DS");

  // get block
  unsigned int l = k + (*itDS)->getN() - 1;
  vector<unsigned int> index_list(4);
  index_list[0] = k;
  index_list[1] = l;
  index_list[2] = 0;
  index_list[3] = B->size(1) - 1;

  B->getBlock(index_list, BBlock);
}

void LinearTIR::computeOutput(const double& time)
{
  if (!isOutputPlugged)
  {
    DSSet vDS = interaction->getDynamicalSystems();
    BlockVector *xTmp = new BlockVector();
    BlockVector *uTmp = new BlockVector();
    DSIterator it;
    for (it = vDS.begin(); it != vDS.end(); it++)
    {
      // Put x and u of each DS into a block
      // Warning: use copy constructors, no link between pointers
      if (((*it)->getType() != LDS) && ((*it)->getType() != LITIDS))
        RuntimeException::selfThrow("LinearTIR - computeOutput: not yet implemented for DS type " + (*it)->getType());

      xTmp->add((*it)->getX());
      if ((*it)->getUPtr() != NULL)
        uTmp->add(*((*it)->getUPtr())) ;
    }

    SimpleVector *y = interaction->getYPtr(0);
    SimpleVector *lambda = interaction->getLambdaPtr(0);

    // compute y
    *y = *C * *xTmp;

    if (D != NULL)
      *y += *D * *lambda;

    if (F != NULL)
      *y += *F * *uTmp;

    if (e != NULL)
      *y += *e;

    // \todo update y, yDot ... depending on the relative degree.

    // free memory
    delete xTmp;
    delete uTmp;
  }
  else
    Relation::computeOutput(time);
}
void LinearTIR::computeFreeOutput(const double& time)
{
  if (!isOutputPlugged)
  {
    DSSet vDS = interaction->getDynamicalSystems();
    BlockVector *xTmp = new BlockVector();
    BlockVector *uTmp = new BlockVector();
    DSIterator it;

    for (it = vDS.begin(); it != vDS.end(); it++)
    {
      // Put xFree and u of each DS into a block
      // Warning: use copy constructors, no link between pointers
      if (((*it)->getType() != LDS) && ((*it)->getType() != LITIDS))
        RuntimeException::selfThrow("LinearTIR - computeFreeOutput: not yet implemented for DS type " + (*it)->getType());

      xTmp->add((*it)->getXFree());
      if ((*it)->getUPtr() != NULL)
        uTmp->add(*((*it)->getUPtr())) ;
    }

    SimpleVector *yFree = interaction->getYPtr(0);
    // warning : yFree is saved in y !!

    // compute yFree
    *yFree = *C * *xTmp;

    if (F != NULL)
      *yFree += *F * *uTmp ;

    if (e != NULL)
      *yFree += *e;

    // \todo update y, yDot ... depending on the relative degree.

    // free memory
    delete xTmp;
    delete uTmp;
  }
  else
    Relation::computeFreeOutput(time);
}

void LinearTIR::computeInput(const double& time)
{
  if (!isOutputPlugged)
  {
    DSSet vDS = interaction->getDynamicalSystems();
    DSIterator it;
    BlockVector *r = new BlockVector();
    for (it = vDS.begin(); it != vDS.end(); it++)
    {
      // Put r of each DS into a block
      // Warning: use addPtr -> link between pointers
      bool isComp = (*it)->getRPtr()->isBlock();
      if (isComp)
      {
        BlockVector * tmp = static_cast<BlockVector*>((*it)->getRPtr());
        r->addPtr(tmp->getVectorPtr(0));
        r->addPtr(tmp->getVectorPtr(1));
      }
      else
        r->addPtr(static_cast<SimpleVector*>((*it)->getRPtr()));
    }

    SimpleVector *lambda = interaction->getLambdaPtr(0);

    *r += *B * *lambda;
    delete r;
  }
  else
    Relation::computeInput(time);
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
  cout << " ================================================== " << endl;
}

void LinearTIR::saveRelationToXML() const
{
  if (relationxml != NULL)
  {
    LinearTIRXML * lTIRxml = (static_cast<LinearTIRXML*>(relationxml));
    lTIRxml->setC(*C);
    lTIRxml->setD(*D);
    lTIRxml->setF(*F);
    lTIRxml->setE(*e);
    lTIRxml->setB(*B);
  }
}

LinearTIR* LinearTIR::convert(Relation *r)
{
  LinearTIR* ltir = dynamic_cast<LinearTIR*>(r);
  return ltir;
}

// Default (private) constructor
LinearTIR::LinearTIR():
  Relation(), C(NULL), D(NULL), F(NULL), e(NULL), B(NULL)
{
  relationType = LINEARTIRELATION;
  isAllocatedIn.resize(5, false);
}

