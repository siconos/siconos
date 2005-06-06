#include "NonSmoothLaw.h"
using namespace std;

// Constructors
// warning -> this is an abstract class, so constructors are usefull only for
// calls in derived classes constructors
NonSmoothLaw::NonSmoothLaw(): nsLawType("none"), nslawxml(NULL)
{}

NonSmoothLaw::NonSmoothLaw(NonSmoothLawXML* newNsLawXml): nsLawType("none"), nslawxml(newNsLawXml)
{}

NonSmoothLaw::~NonSmoothLaw()
{}

void NonSmoothLaw::saveNonSmoothLawToXML()
{
  IN("NonSmoothLaw::saveNonSmoothLawToXML\n");
  if (nslawxml == NULL)
    RuntimeException::selfThrow("NonSmoothLaw::saveNonSmoothLawToXML - The NonSmoothLawXML object doesn't exists");
  OUT("NonSmoothLaw::saveNonSmoothLawToXML\n");
}

void NonSmoothLaw::display() const
{
  IN("NonSmoothLaw::display\n");
  cout << "---------------------------------------------------" << endl;
  cout << "____ data of the NonSmoothLaw " << endl;
  cout << " NonSmoothLaw Type :" << nsLawType << endl;
  cout << "____________________________" << endl;
  cout << "---------------------------------------------------" << endl;
  OUT("NonSmoothLaw::display\n");
}

