//$Id: LinearEC.cpp,v 1.6 2005/03/15 14:44:03 jbarbier Exp $
#include "LinearEC.h"

LinearEC::LinearEC(): EqualityConstraint()
{
  this->type = LINEAREC;
}

LinearEC::LinearEC(EqualityConstraintXML *ecxml): EqualityConstraint(ecxml)
{
  this->type = LINEAREC;
}

LinearEC::~LinearEC()
{}

void LinearEC::createEqualityConstraint(EqualityConstraintXML *ecXML ,
                                        int number,  SiconosMatrix *G,
                                        vector<DSInputOutput*> *dsioVector)
{
  if (ecXML != NULL)
  {
    this->ecXML = ecXML;
    this->type = NLINEAREC;
    this->fillEqualityConstraintWithEqualityConstraintXML();
  }
  else
  {
    this->ecXML = NULL;
    this->type = NLINEAREC;
    this->number = number;
    this->G = *G;
    this->dsioVector = *dsioVector;
  }
}

//$Log: LinearEC.cpp,v $
//Revision 1.6  2005/03/15 14:44:03  jbarbier
//- pySiconos.i edited to remove local paths
//
//- checkCoherency checks whether the DSInputOutputs and EqualityConstraints have unique numbers
//
//Revision 1.5  2005/03/15 09:57:47  jbarbier
//- EqualityConstraint save OK
//
//Revision 1.4  2005/03/14 16:05:27  jbarbier
//- manual creation of DSInputOutput saving OK
//
//- in progress for EqualityConstraint
//
//Revision 1.3  2005/03/11 15:06:20  jbarbier
//- save to XML methods of EqualityConstraint and DSInputOutput added
//
//- XML loading process modified : Model loads NSDS, then NSDS loads the DynamicalSystems, EqualityConstraints, Interactions; Modle loads Strategy, then Strategy loads TimeDiscretisation, then the Integrators, then the OneStepNSProblem
//
//Revision 1.2  2005/03/09 15:30:32  jbarbier
//- add of LagrangianEC class
//
//- in progress : implementation of the EqualityConstraint and DSInputOutput - create methods
//
//Revision 1.1  2005/01/17 10:56:25  jbarbier
//- classes EqualityConstraint and DSInputOutput added with inherited classes
//
//- classes EqualityConstraintXML and DSInputOutputXML added with inherited classes
//