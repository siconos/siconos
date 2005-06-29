#include "NewSiconosVector.h"
using namespace std;

SiconosVector::SiconosVector()
{}

SiconosVector::~SiconosVector()
{}

SiconosVector::SiconosVector(const SiconosVector &)
{
  SiconosVectorException::selfThrow("SiconosVector::copy constructor, not allowed for abstract class");
}

