//$Id: ComplementarityConditionNSL.cpp,v 1.9 2005/03/08 12:41:36 jbarbier Exp $
#include "ComplementarityConditionNSL.h"

#include "check.h"


ComplementarityConditionNSL::ComplementarityConditionNSL(): NonSmoothLaw()
{
  this->nsLawType = COMPLEMENTARITYCONDITIONNSLAW;
}

ComplementarityConditionNSL::ComplementarityConditionNSL(NonSmoothLawXML* nslawxml): NonSmoothLaw(nslawxml)
{
  this->nsLawType = COMPLEMENTARITYCONDITIONNSLAW;
}

ComplementarityConditionNSL::~ComplementarityConditionNSL()
{}

bool ComplementarityConditionNSL::isVerified(void) const
{
  bool res = false;

  // to do

  return res;
}

void ComplementarityConditionNSL::fillNonSmoothLawWithNonSmoothLawXML()
{
  IN("ComplementarityConditionNSL::fillNonSmoothLawWithNonSmoothLawXML\n");
  NonSmoothLaw::fillNonSmoothLawWithNonSmoothLawXML();
  if (this->nslawxml != NULL)
  {}
  else RuntimeException::selfThrow("ComplementarityConditionNSL::fillNonSmoothLawWithNonSmoothLawXML - ComplementarityConditionNSLXML object not exists");
  OUT("ComplementarityConditionNSL::fillNonSmoothLawWithNonSmoothLawXML\n");
}

void ComplementarityConditionNSL::saveNonSmoothLawToXML()
{
  IN("ComplementarityConditionNSL::saveNonSmoothLawToXML\n");
  NonSmoothLaw::saveNonSmoothLawToXML();
  if (this->nslawxml != NULL)
  {}
  else RuntimeException::selfThrow("ComplementarityConditionNSL::saveNonSmoothLawToXML - ComplementarityConditionNSLXML object not exists");
  OUT("ComplementarityConditionNSL::saveNonSmoothLawToXML\n");
}

void ComplementarityConditionNSL::createNonSmoothLaw(ComplementarityConditionNSLXML * nslawXML)//, Interaction * interaction)
{
  if (nslawXML != NULL)
  {
    this->nslawxml = nslawXML;
    this->nsLawType = COMPLEMENTARITYCONDITIONNSLAW;
    this->fillNonSmoothLawWithNonSmoothLawXML();
  }
  else
  {}
}

ComplementarityConditionNSL* ComplementarityConditionNSL::convert(NonSmoothLaw* nsl)
{
  cout << "ComplementarityConditionNSL::convert (NonSmoothLaw* nsl)" << endl;
  ComplementarityConditionNSL* ccnsl = dynamic_cast<ComplementarityConditionNSL*>(nsl);
  return ccnsl;
}


//$Log: ComplementarityConditionNSL.cpp,v $
//Revision 1.9  2005/03/08 12:41:36  jbarbier
//- constant variables files modified :
//Some constants added in SiconosConst
//
//all global tag of the modeling tools are in XMLTagsName, other tags are specific to an XML class
//
//Revision 1.8  2005/03/07 13:17:18  jbarbier
//- new test : Ball2D, with a ball moving in a 2D system
//
//- another constant variables moved/refactored in XMLTagsName
//- making uniform the name of the constant variables
//
//Revision 1.7  2005/02/11 17:35:54  charlety
//
//_ little "inspection of code"
//_ basic getters and setters passed inline
//_ getters functions passed const
//
//Revision 1.6  2005/01/31 16:26:19  charlety
//
//_ Added a method named "convert" to classes which inherits from another. This is necessary for Python interface, in order to be able to use down-casting mechanism.
//
//Revision 1.5  2004/09/03 14:41:40  jbarbier
//- new functions to create the boundary condition of the dynamical systems
//- new functions to add an interaction to a NSDS
//- new functions to create the relation and the non-smooth law of an interaciton
//
//Revision 1.4  2004/08/12 11:55:14  jbarbier
//- new methods createModel, createNSDS, createStrategy, ...
//they now allow to make the link with upper objects of the platform
//it will be used for the creation of the platform without XML input file
//
//- the createModel method is finished but the attributes of the other objects
//of the platform are missing for the conctruction
//
//Revision 1.3  2004/08/11 14:43:45  jbarbier
//- beginning of the mechanism of creation without XML input file of the objects of the platform with the
//creatObjtect methods
//
//- function saveWToXML for Moreau integrator, and same specific functions to save
//M,q and Q,p for LCP and QP
//
//- function to check coherency of the Model
//
//Revision 1.2  2004/07/27 09:32:43  jbarbier
//- functions createNSLaw for complementarityConditionNSL, newtonImpactLawNSL and RelayNSL
//
//Revision 1.1  2004/07/06 14:54:48  acary
//Renaming NSLaw into NonSmoothLaw
//Renaming RelayNSLaw into RelayNSL
//Renaming CCNSLaw into ComplementarityConditionNSL
//Renaming NewtonImpactLaw into NewtonImpactLawNSL
//
//Revision 1.10  2004/06/30 09:12:17  acary
//Added tag CVS and Macro IN and OUT
//
