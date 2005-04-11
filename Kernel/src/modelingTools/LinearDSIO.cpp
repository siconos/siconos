#include "LinearDSIO.h"
#include "check.h"

LinearDSIO::LinearDSIO(): DSInputOutput()
{
  this->dsioType = LINEARDSIO;
}

LinearDSIO::LinearDSIO(DSInputOutputXML* dsioxml): DSInputOutput(dsioxml)
{
  this->dsioType = LINEARDSIO;
}

LinearDSIO::~LinearDSIO()
{}


void LinearDSIO::fillDSInputOutputWithDSInputOutputXML()
{
  if (this->dsioxml != NULL)
  {
    this->number = this->dsioxml->getNumber();
    this->A = static_cast<LinearDSIOXML*>(this->dsioxml)->getA();
    this->B = static_cast<LinearDSIOXML*>(this->dsioxml)->getB();
  }
  else RuntimeException::selfThrow("DSInputOutput::fillDSInputOutputWithDSInputOutputXML - object DSInputOutputXML does not exist");
}

void LinearDSIO::saveDSInputOutputToXML()
{
  if (this->dsioxml != NULL)
  {
    /*
     * these attributes are only required for LagrangianNonLinear DSInputOutput !
     */
    static_cast<LinearDSIOXML*>(this->dsioxml)->setA(&(this->A));
    static_cast<LinearDSIOXML*>(this->dsioxml)->setB(&(this->B));
  }
  else RuntimeException::selfThrow("DSInputOutput::saveDSInputOutputToXML - object DSInputOutputXML does not exist");
}

void LinearDSIO::createDSInputOutput(DSInputOutputXML * dsioXML, int number,
                                     SiconosMatrix *A, SiconosMatrix *B)
{
  if (dsioXML != NULL)
  {
    //    this->init();
    this->dsioxml = dsioXML;
    this->dsioType = LINEARDSIO;
    this->fillDSInputOutputWithDSInputOutputXML();
  }
  else
  {
    //this->dsioxml = dsioXML;
    this->dsioType = LINEARDSIO;
    this->number = number;
    this->A = *A;
    this->B = *B;
  }
}

