//$Id: RelayNSL.cpp,v 1.12 2005/03/08 12:41:36 jbarbier Exp $
#include "RelayNSL.h"

#include "check.h"

RelayNSL::RelayNSL(): NonSmoothLaw()
{
  c = 0.0;
  d = 0.0;
  this->nsLawType = RELAYNSLAW;
}

RelayNSL::RelayNSL(NonSmoothLawXML* nslawxml): NonSmoothLaw(nslawxml)
{
  c = 0.0;
  d = 0.0;
  this->nsLawType = RELAYNSLAW;
}

RelayNSL::RelayNSL(double c, double d)
{
  this->c = c;
  this->d = d;
  this->nsLawType = RELAYNSLAW;
}

RelayNSL::~RelayNSL()
{}

bool RelayNSL::isVerified(void) const
{
  bool res = false;

  // to do

  return res;
}

void RelayNSL::fillNonSmoothLawWithNonSmoothLawXML()
{
  IN("RelayNSL::fillNonSmoothLawWithNonSmoothLawXML\n");
  NonSmoothLaw::fillNonSmoothLawWithNonSmoothLawXML();
  if (this->nslawxml != NULL)
  {
    this->c = (static_cast<RelayNSLXML*>(this->nslawxml))->getC();
    this->d = (static_cast<RelayNSLXML*>(this->nslawxml))->getD();

    //    this->display();
  }
  else RuntimeException::selfThrow("RelayNSL::fillNonSmoothLawWithNonSmoothLawXML - object NonSmoothLawXML does not exist");
  OUT("RelayNSL::fillNonSmoothLawWithNonSmoothLawXML\n");
}

void RelayNSL::display() const
{
  cout << "------------------------------------" << endl;
  cout << "____ data of the RelayNSL" << endl;
  cout << "| c : " << this->c << endl;
  cout << "| d : " << this->d << endl;
  cout << "____________________________" << endl;
  cout << "------------------------------------" << endl;
}

void RelayNSL::saveNonSmoothLawToXML()
{
  IN("RelayNSL::saveNonSmoothLawToXML\n");
  static_cast<RelayNSLXML*>(this->nslawxml)->setC(this->c);
  static_cast<RelayNSLXML*>(this->nslawxml)->setD(this->d);
  //NonSmoothLaw::fillNonSmoothLawWithNonSmoothLawXML();
  OUT("RelayNSL::saveNonSmoothLawToXML\n");
}

void RelayNSL::createNonSmoothLaw(RelayNSLXML * nslawXML, double c, double d)//, Interaction * interaction)
{
  if (nslawXML != NULL)
  {
    this->nslawxml = nslawXML;
    this->nsLawType = RELAYNSLAW;
    this->fillNonSmoothLawWithNonSmoothLawXML();
  }
  else
  {
    this->c = c;
    this->d = d;
  }
}


RelayNSL* RelayNSL::convert(NonSmoothLaw* nsl)
{
  cout << "RelayNSL::convert (NonSmoothLaw* nsl)" << endl;
  RelayNSL* rnsl = dynamic_cast<RelayNSL*>(nsl);
  return rnsl;
}

//$Log: RelayNSL.cpp,v $
//Revision 1.12  2005/03/08 12:41:36  jbarbier
//- constant variables files modified :
//Some constants added in SiconosConst
//
//all global tag of the modeling tools are in XMLTagsName, other tags are specific to an XML class
//
//Revision 1.11  2005/03/07 13:17:20  jbarbier
//- new test : Ball2D, with a ball moving in a 2D system
//
//- another constant variables moved/refactored in XMLTagsName
//- making uniform the name of the constant variables
//
//Revision 1.10  2005/02/11 17:36:02  charlety
//
//_ little "inspection of code"
//_ basic getters and setters passed inline
//_ getters functions passed const
//
//Revision 1.9  2005/01/31 16:26:24  charlety
//
//_ Added a method named "convert" to classes which inherits from another. This is necessary for Python interface, in order to be able to use down-casting mechanism.
//
//Revision 1.8  2004/09/21 11:49:09  jbarbier
//- correction in the XML save for a manual construction of the platform :
//    DS_Concerned of the Interaction
//    DS_Concerned of the Integrator
//
//- test updated for these changes
//
//Revision 1.7  2004/09/16 11:35:24  jbarbier
//- save of the TimeDiscretisation in a XML file in manual creation of the
//platform which was forgotten is now available.
//
//- the save of the platform's data can be done when the platform is created with
//an XML input file and completed with dynmical systems, interactions, one-step
//non smooth problem and one-step integrator.
//
//Revision 1.6  2004/09/03 14:41:47  jbarbier
//- new functions to create the boundary condition of the dynamical systems
//- new functions to add an interaction to a NSDS
//- new functions to create the relation and the non-smooth law of an interaciton
//
//Revision 1.5  2004/08/17 15:12:43  jbarbier
//- methods createDynamicalSystem, createBoundaryCondition, createInteraction,
//createRelation and createNSLaw completed with the required attributes
//
//Revision 1.4  2004/08/12 11:55:18  jbarbier
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
//Revision 1.11  2004/06/30 09:12:20  acary
//Added tag CVS and Macro IN and OUT
//
