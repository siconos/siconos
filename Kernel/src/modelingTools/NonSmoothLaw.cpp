#include "NonSmoothLaw.h"


NonSmoothLaw::NonSmoothLaw()
{
  this->nslawxml = NULL;
}

NonSmoothLaw::NonSmoothLaw(NonSmoothLawXML* nslawxml)
{
  this->nslawxml = nslawxml;
}

NonSmoothLaw::~NonSmoothLaw()
{}


void NonSmoothLaw::fillNonSmoothLawWithNonSmoothLawXML()
{
  IN("NonSmoothLaw::fillNonSmoothLawWithNonSmoothLawXML\n");
  if (this->nslawxml == NULL)
    RuntimeException::selfThrow("NonSmoothLaw::fillNonSmoothLawWithNonSmoothLawXML - The NonSmoothLawXML object doesn't exists");
  OUT("NonSmoothLaw::fillNonSmoothLawWithNonSmoothLawXML\n");
}

void NonSmoothLaw::saveNonSmoothLawToXML()
{
  IN("NonSmoothLaw::saveNonSmoothLawToXML\n");
  if (this->nslawxml == NULL)
    RuntimeException::selfThrow("NonSmoothLaw::saveNonSmoothLawToXML - The NonSmoothLawXML object doesn't exists");
  OUT("NonSmoothLaw::saveNonSmoothLawToXML\n");
}

void NonSmoothLaw::display() const
{
  IN("NonSmoothLaw::display\n");
  cout << "---------------------------------------------------" << endl;
  cout << "____ data of the NonSmoothLaw " << endl;
  cout << " NonSmoothLaw Type :" << this->nsLawType << endl;
  cout << "____________________________" << endl;
  cout << "---------------------------------------------------" << endl;
  OUT("NonSmoothLaw::display\n");
}

