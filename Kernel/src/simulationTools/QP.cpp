
#include "QP.h"
#include "QPXML.h"
#include "check.h"

QP::QP(): OneStepNSProblem()
{
  this->nspbType = QP_OSNSP;
}

QP::QP(OneStepNSProblemXML* osnspbxml): OneStepNSProblem(osnspbxml)
{
  this->nspbType = QP_OSNSP;
}


QP::~QP()
{}


SiconosMatrix* QP::getQPtr(void)
{
  return &this->Q;
}

SimpleVector* QP::getPPtr(void)
{
  return &this->p;
}



void QP::formalise(double time)
{
  //OUT("QP::formaliseOneStepNSProblem\n");
}


void QP::compute(void)
{
  //OUT("QP::computeOneStepNSProblem\n");
}


void QP::fillNSProblemWithNSProblemXML()
{
  OUT("QP::fillNSProblemWithNSProblemXML");
  OneStepNSProblem::fillNSProblemWithNSProblemXML();
  if (this->onestepnspbxml != NULL)
  {
    this->Q = (static_cast<QPXML*>(this->onestepnspbxml))->getQ();
    this->p = (static_cast<QPXML*>(this->onestepnspbxml))->getP();

    //    this->display();
  }
  else RuntimeException::selfThrow("QP::fillNSProblemWithNSProblemXML - the OneStepNSProblemXML object does not exist");
}

void QP::display() const
{
  cout << "------------------------------------------------------" << endl;
  cout << "____ data of the DynamicalSystem read from a XML file" << endl;
  cout << "| Q " << endl;
  this->Q.display();
  cout << "| p " << endl ;
  this->p.display();
  cout << "____________________________" << endl;
  cout << "------------------------------------------------------" << endl;
}

void QP::saveNSProblemToXML()
{
  OUT("QP::saveNSProblemToXML");
  OneStepNSProblem::saveNSProblemToXML();
  if (this->onestepnspbxml != NULL)
  {
    /*
     * the Q et p of the LCP can be saved by calling the saveQToXML() and savePToXML()
     */
    //(static_cast<QPXML*>(this->onestepnspbxml))->setQ( &(this->Q) );
    //(static_cast<QPXML*>(this->onestepnspbxml))->setP( &(this->p) );
  }
  else RuntimeException::selfThrow("QP::saveNSProblemToXML - the OneStepNSProblemXML object does not exist");
}

void QP::savePToXML()
{
  IN("QP::savePToXML\n");
  if (this->onestepnspbxml != NULL)
  {
    (static_cast<QPXML*>(this->onestepnspbxml))->setP(&(this->p));
  }
  else RuntimeException::selfThrow("QP::savePToXML - OneStepNSProblemXML object not exists");
  OUT("QP::savePToXML\n");
}

void QP::saveQToXML()
{
  IN("QP::saveQToXML\n");
  if (this->onestepnspbxml != NULL)
  {
    (static_cast<QPXML*>(this->onestepnspbxml))->setQ(&(this->Q));
  }
  else RuntimeException::selfThrow("QP::saveQToXML - OneStepNSProblemXML object not exists");
  OUT("QP::saveQToXML\n");
}

void QP::createOneStepNSProblem(OneStepNSProblemXML * osnspbXML, Strategy * strategy)
{
  if (osnspbXML != NULL)
  {
    this->onestepnspbxml = osnspbXML;
    this->nspbType = QP_OSNSP;

    this->fillNSProblemWithNSProblemXML();
  }
  else
  {
    this->strategy = strategy;
    this->nspbType = QP_OSNSP;
    this->fillInteractionVector();
  }
  this->init();
}


QP* QP::convert(OneStepNSProblem* osnsp)
{
  cout << "QP::convert (DynamicalSystem* osnsp)" << endl;
  QP* qp = dynamic_cast<QP*>(osnsp);
  return qp;
}

