//$Id: SiconosMemory.cpp,v 1.12 2005/02/11 17:36:03 charlety Exp $
#include "SiconosMemory.h"
#include "check.h"

SiconosMemory::SiconosMemory()
{
  IN("SiconosMemory()\n");

  this->memorySize = 0;
  this->nbVectorsInMemory = 0;
  this->vectorMemory.resize(memorySize);
  this->memoryXML = NULL;

  OUT("SiconosMemory()\n");
}

SiconosMemory::SiconosMemory(int memorySize)
{
  IN("SiconosMemory(int memorySize)\n");

  this->memorySize = memorySize;
  this->nbVectorsInMemory = 0;
  this->vectorMemory.clear();
  this->vectorMemory.resize(this->memorySize);
  for (int i = 0; i < this->memorySize; i++)
    (this->vectorMemory)[i] = new SimpleVector();
  this->memoryXML = NULL;

  OUT("SiconosMemory(int memorySize)\n");
}

SiconosMemory::SiconosMemory(int memorySize, SiconosMemoryXML *memoryXML)
{
  IN("SiconosMemory(int memorySize, SiconosMemoryXML *memoryXML)\n");

  this->memorySize = memorySize;
  this->memoryXML = memoryXML;

  this->nbVectorsInMemory = 0;
  this->vectorMemory.clear();
  this->vectorMemory.resize(this->memorySize);
  for (int i = 0; i < this->memorySize; i++)
    (this->vectorMemory)[i] = new SimpleVector();

  if (memoryXML != NULL) this->setVectorMemory(memorySize, memoryXML->getSiconosMemoryVector());

  OUT("SiconosMemory(int memorySize, SiconosMemoryXML *memoryXML)\n");
}

SiconosMemory::SiconosMemory(SiconosMemoryXML *memoryXML)
{
  IN("SiconosMemory(SiconosMemoryXML *memoryXML)\n");
  this->memorySize = 0;
  this->nbVectorsInMemory = 0;
  this->vectorMemory.resize(this->memorySize);
  this->memoryXML = NULL;

  this->setVectorMemory(memoryXML->getSiconosMemorySize(), memoryXML->getSiconosMemoryVector());
  this->memoryXML = memoryXML;

  OUT("SiconosMemory(SiconosMemoryXML *memoryXML)\n");
}


SiconosMemory::SiconosMemory(vector<SiconosVector*> V)
{
  IN("SiconosMemory(vector<SiconosVector*> V)\n");

  this->memorySize = V.size();
  this->nbVectorsInMemory = this->memorySize;
  this->vectorMemory.clear();
  this->vectorMemory.resize(this->memorySize);
  for (int i = 0; i < this->memorySize; i++)
  {
    (this->vectorMemory)[i] = new SimpleVector();
    *(this->vectorMemory[i]) = *(V[i]);
  }

  OUT("SiconosMemory(vector<SiconosVector*> V)\n");
}

SiconosMemory::SiconosMemory(int memorySize, vector<SiconosVector*> V)
{
  IN("SiconosMemory(int memorySize, vector<SiconosVector*> V)\n");

  if (memorySize < V.size())
  {
    // exception : out of range
    SiconosMemoryException::selfThrow("SiconosMemory(int memorySize, vector<SiconosVector*> V) : V.size > memorySize");
  }
  else
  {
    this->memorySize = memorySize;
    this->nbVectorsInMemory = V.size();
    this->vectorMemory.clear();
    this->vectorMemory.resize(this->memorySize);
    for (int i = 0; i < this->memorySize; i++)
    {
      (this->vectorMemory)[i] = new SimpleVector();
      if (i < this->nbVectorsInMemory)
        *(this->vectorMemory[i]) = *(V[i]);
    }
  }
  this->memoryXML = NULL;

  OUT("SiconosMemory(int memorySize, vector<SiconosVector*> V)\n");
}

SiconosMemory::SiconosMemory(const SiconosMemory&  source)
{
  IN("SiconosMemory(const SiconosMemory&  source) \n");
  int i;

  this->memorySize = source.memorySize;
  this->nbVectorsInMemory = source.nbVectorsInMemory;

  this->vectorMemory.clear();
  this->vectorMemory.resize(this->memorySize);
  for (i = 0; i < this->memorySize; i++)
    if (i < this->nbVectorsInMemory)
      (this->vectorMemory)[i] = new SimpleVector(*(source.vectorMemory[i]));
    else
      (this->vectorMemory)[i] = new SimpleVector();
  this->memoryXML = NULL;

  OUT("SiconosMemory(const SiconosMemory&  source) \n");
}


SiconosMemory::~SiconosMemory()
{
  IN("~SiconosMemory()\n");

  for (int i = 0; i < this->memorySize; i++)
    delete(this->vectorMemory)[i];
  this->vectorMemory.clear();
  this->memorySize = 0;
  this->nbVectorsInMemory = 0;

  OUT("~SiconosMemory()\n");
}


/* ********************************************** */

void SiconosMemory::setVectorMemory(vector<SiconosVector*> V)
{
  IN("SiconosMemory::setVectorMemory(vector<SiconosVector*> V)\n");

  for (int i = 0; i < this->memorySize; i++)
    delete(this->vectorMemory)[i];

  this->memorySize = V.size();
  this->nbVectorsInMemory = this->memorySize;
  this->vectorMemory.clear();
  this->vectorMemory.resize(this->memorySize);
  for (int i = 0; i < this->memorySize; i++)
  {
    (this->vectorMemory)[i] = new SimpleVector();
    *(this->vectorMemory[i]) = *(V[i]);
  }

  OUT("SiconosMemory::setVectorMemory(vector<SiconosVector*> V)\n");
}

void SiconosMemory::setVectorMemory(int mSize, vector<SiconosVector*> V)
{
  IN("SiconosMemory::setVectorMemory(int mSize, vector<SiconosVector*> V)\n");

  if (mSize < V.size())
  {
    // exception : out of range
    //SiconosMemoryException::selfThrow("setVectorMemory(int mSize, vector<SiconosVector*> V) : V.size > mSize");
    cout << "setVectorMemory(int mSize, vector<SiconosVector*> V) : V.size > mSize.\nThe saved SiconosVector which position in the vector of SiconosVector is greater than mSize, will be lost." << endl;

    /*
     * the node in the DOM tree which are unused will be deleted
     */
    if (this->memoryXML != NULL)
      this->memoryXML->deleteUnusedMemoryNodes(mSize);
    else SiconosMemoryException::selfThrow("setVectorMemory(int mSize, vector<SiconosVector*> V) : V.size > mSize");
  }
  else
  {
    int i;

    if (this->memorySize != mSize)
    {
      for (i = 0; i < this->memorySize; i++)
        delete(this->vectorMemory)[i];
      this->memorySize = mSize;

      this->vectorMemory.clear();
      this->vectorMemory.resize(this->memorySize);
      for (i = 0; i < this->memorySize; i++)
        (this->vectorMemory)[i] = new SimpleVector();
    }

    this->nbVectorsInMemory = V.size();
    for (i = 0; i < this->memorySize; i++)
    {
      if (i < this->nbVectorsInMemory)
        *(this->vectorMemory[i]) = *(V[i]);
    }
  }

  OUT("SiconosMemory::setVectorMemory(int mSize, vector<SiconosVector*> V)\n");
}


SiconosVector* SiconosMemory::getSiconosVector(int index) const
{
  IN("SiconosMemory::getSiconosVector(int index) const\n");

  if ((index < 0) || (index >= this->memorySize))
  {
    // exception : out of range
    SiconosMemoryException::selfThrow("getSiconosVector(int index) : index > memorySize");
  }
  else if (index >= this->nbVectorsInMemory)
  {
    // Warning : not significant value
    cout << "WARNING --- SiconosMemory::getSiconosVector(int index) : index > nbVectorsInMemory" << endl;
    OUT("SiconosMemory::getSiconosVector(int index) const\n");
    return NULL;
  }
  else
  {
    OUT("SiconosMemory::getSiconosVector(int index) const\n");
    return (this->vectorMemory)[index];
  }
}


void SiconosMemory::swap(SiconosVector* v)
{
  IN("SiconosMemory::swap(SiconosVector* v)\n");
  int i;
  SiconosVector* tmp;

  // we add always in head of vectorMemory
  // direct swapping
  if (this->nbVectorsInMemory < this->memorySize)
  {
    this->nbVectorsInMemory ++;
  }
  else if (this->memorySize <= 0)
  {
    // exception : out of range
    SiconosMemoryException::selfThrow("swap(SiconosVector* v) : memorySize <= 0");
  }

  tmp = (this->vectorMemory)[this->nbVectorsInMemory - 1];

  for (i = this->nbVectorsInMemory - 1; i >= 0; i--)
  {
    if (i == 0)
    {
      // where we put v (we copy it)
      (this->vectorMemory)[i] = tmp;
      *(this->vectorMemory)[i] = *v;
    }
    else
    {
      // we move the old vectors in memory
      (this->vectorMemory)[i] = (this->vectorMemory)[i - 1];
    }
  }

  OUT("SiconosMemory::swap(SiconosVector* v)\n");
}



void SiconosMemory::display() const
{
  cout << "| memorySize : " << this->memorySize << endl;
  cout << "| nbVectorsInMemory : " << this->nbVectorsInMemory << endl;
  cout << "| vectorMemory size : " << this->vectorMemory.size() << endl;
  for (int i = 0; i < nbVectorsInMemory; i++)
    //for (int i = 0; i < this->memorySize; i++)
  {
    cout << "Memory " << i << " >>> " << (this->vectorMemory)[i] << " | ";
    this->vectorMemory[i]->display();
  }
}


SiconosMemory& SiconosMemory::operator = (const SiconosMemory& source)
{
  IN("SiconosMemory::operator = (const SiconosMemory& source)\n");

  int i;

  for (i = 0; i < this->memorySize; i++)
    delete(this->vectorMemory)[i];

  this->memorySize = source.memorySize;

  this->vectorMemory.clear();
  this->vectorMemory.resize(this->memorySize);
  for (i = 0; i < this->memorySize; i++)
    (this->vectorMemory)[i] = new SimpleVector();

  this->nbVectorsInMemory = source.nbVectorsInMemory;
  for (i = 0; i < this->memorySize; i++)
  {
    if (i < this->nbVectorsInMemory)
      *(this->vectorMemory[i]) = *(source.vectorMemory[i]);
  }
  this->memoryXML = source.memoryXML;

  OUT("SiconosMemory::operator = (const SiconosMemory& source)\n");
  return *this;
}


void SiconosMemory::fillMemoryWithMemoryXML()
{
  IN("SiconosMemory::fillMemoryWithMemoryXML()\n");
  int size = this->memoryXML->getSiconosMemorySize();
  vector<SiconosVector*> v = this->memoryXML->getSiconosMemoryVector();
  this->setVectorMemory(size, v);
  OUT("SiconosMemory::fillMemoryWithMemoryXML()\n");
}


//$Log: SiconosMemory.cpp,v $
//Revision 1.12  2005/02/11 17:36:03  charlety
//
//_ little "inspection of code"
//_ basic getters and setters passed inline
//_ getters functions passed const
//
//Revision 1.11  2004/09/10 11:26:22  charlety
//
//_ Integration of the new version of the SiconosVector in the platform. the class simpleVector is used mostly to replace old SiconosVector. When a vector can be composite or simple, like the state of a dynamical system, a pointer on SiconosVector is used, and the vector is initialized simple or composite when the system is initialized.
//
//_ All the tests which worked with the previous version of the vector are OK with the new version.
//
//_ Example SICONOS and bouncingBall are OK
//
//_ some comments have still to be adapted to NewSiconosVector .
//
//_ This version of NewSiconosVector could be called 0.9. some details have to be fixed, it will be done before the end of September.
//
//Revision 1.10  2004/08/05 14:58:54  jbarbier
//- correction of the SiconosMemory::setVectorMemory(int mSize,
//vector<SiconosVector*> V) method when there's no memoryXML object
//
//Revision 1.9  2004/08/04 11:03:22  jbarbier
//- about the SiconosMemory : when a SiconosMemory has a maxSize greater than the
//number of steps in memory required by an integrator, the oldest SiconosVector
//are deleted
//
//- the way to initialize the SiconosMemory by the integrator has been updated to
//match with these changes
//
//Revision 1.8  2004/07/28 14:13:48  charlety
//
//_ add of function to get SiconosMemory objects from XML (xMemory, xDotMemory in DynamicalSystem, etc.)
//
//Revision 1.7  2004/07/28 08:23:23  charlety
//
//_ Found a bug in a constructor of SiconosMemory
//_ try to get a SiconosMemory from XML successfull
//
//Revision 1.6  2004/07/16 06:23:31  charlety
//
//_ modification of the operator = of SiconosVector : the left operand can now be a composite vector
//
//_ added functions to load / save a SiconosMemory in an XML file
//
//Revision 1.5  2004/07/13 12:44:53  jbarbier
//- integration of the SiconosMemory to the plateform
//- files SiconosMemoryXML added
//
//Revision 1.4  2004/07/12 11:22:54  charlety
//
//_ the state vector x of the dynamical system is now plugged to q and the velocity when this system is a Lagrangian one.
//
//_ A function to force a SiconosVector to be composite has been written. This is a temporary solution. We should use the operator = and the constructor by copy instead. This problem will be fixed later in the summer.
//
//Revision 1.3  2004/07/09 11:14:53  charlety
//
//_ Added a constructor by copy and an operator = in class SiconosMemory
//_ the getters on memory in DynamicalSystems return now some pointers
//
//Revision 1.2  2004/07/08 09:13:43  charlety
//
//_ creation of a Siconos exception dedicated to the class SiconosMemory
//_ new tests for SiconosMemory
//
//Revision 1.1  2004/07/07 13:53:13  charlety
//
//_ First version of the memory object
//_ some little bugs corrected otherwhere
//
