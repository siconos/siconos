//$Id: NonSmoothLaw.cpp,v 1.4 2005/02/11 17:36:02 charlety Exp $
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

//$Log: NonSmoothLaw.cpp,v $
//Revision 1.4  2005/02/11 17:36:02  charlety
//
//_ little "inspection of code"
//_ basic getters and setters passed inline
//_ getters functions passed const
//
//Revision 1.3  2004/09/14 13:49:54  jbarbier
//- files added in sample/ to run run the main_siconos test program
//
//- all the platform can now be saved in an XML file when it is created manually
//
//Revision 1.2  2004/08/11 14:43:45  jbarbier
//- beginning of the mechanism of creation without XML input file of the objects of the platform with the
//creatObjtect methods
//
//- function saveWToXML for Moreau integrator, and same specific functions to save
//M,q and Q,p for LCP and QP
//
//- function to check coherency of the Model
//
//Revision 1.1  2004/07/06 14:54:48  acary
//Renaming NSLaw into NonSmoothLaw
//Renaming RelayNSLaw into RelayNSL
//Renaming CCNSLaw into ComplementarityConditionNSL
//Renaming NewtonImpactLaw into NewtonImpactLawNSL
//
//Revision 1.9  2004/06/30 13:35:55  acary
//Formalization and Computation of the LCP with the NewtonImpactLawNSL
//
//Revision 1.8  2004/06/30 09:12:20  acary
//Added tag CVS and Macro IN and OUT
//