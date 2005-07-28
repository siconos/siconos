
#include "QP.h"
using namespace std;


QP::QP(OneStepNSProblemXML* osnspbxml, Strategy* newStrat):
  OneStepNSProblem(osnspbxml, newStrat), Q(NULL), p(NULL),
  isQAllocatedIn(true), isPAllocatedIn(true)
{
  nspbType = QP_OSNSP;
  if (onestepnspbxml != NULL)
  {
    QPXML * xmlqp = (static_cast<QPXML*>(onestepnspbxml));
    int size = (xmlqp->getP()).size();
    n = size;
    Q = new SiconosMatrix(size, size);
    p = new SimpleVector(size);
    if (xmlqp->hasQ())
      *Q = (static_cast<QPXML*>(onestepnspbxml))->getQ();
    if (xmlqp->hasP())
      *p = (static_cast<QPXML*>(onestepnspbxml))->getP();
  }
  else RuntimeException::selfThrow("QP::xml constructor, xml file=NULL");
}

QP::~QP()
{
  if (isQAllocatedIn)
  {
    delete Q;
    Q = NULL;
  }
  if (isPAllocatedIn)
  {
    delete p;
    p = NULL;
  }
}

void QP::compute(const double& time)
{
  RuntimeException::selfThrow("QP::compute not yet implemented");

}

void QP::display() const
{
  cout << "------------------------------------------------------" << endl;
  cout << "____ data of the DynamicalSystem read from a XML file" << endl;
  cout << "| Q " << endl;
  if (Q != NULL) Q->display();
  else cout << "->NULL" << endl;
  cout << "| p " << endl ;
  if (p != NULL) p->display();
  else cout << "->NULL" << endl;
  cout << "____________________________" << endl;
  cout << "------------------------------------------------------" << endl;
}

void QP::saveNSProblemToXML()
{
  OUT("QP::saveNSProblemToXML");
  OneStepNSProblem::saveNSProblemToXML();
  if (onestepnspbxml != NULL)
  {}
  else RuntimeException::selfThrow("QP::saveNSProblemToXML - the OneStepNSProblemXML object does not exist");
}

void QP::savePToXML()
{
  IN("QP::savePToXML\n");
  if (onestepnspbxml != NULL)
  {
    (static_cast<QPXML*>(onestepnspbxml))->setP(*p);
  }
  else RuntimeException::selfThrow("QP::savePToXML - OneStepNSProblemXML object not exists");
  OUT("QP::savePToXML\n");
}

void QP::saveQToXML()
{
  IN("QP::saveQToXML\n");
  if (onestepnspbxml != NULL)
  {
    (static_cast<QPXML*>(onestepnspbxml))->setQ(*Q);
  }
  else RuntimeException::selfThrow("QP::saveQToXML - OneStepNSProblemXML object not exists");
  OUT("QP::saveQToXML\n");
}

QP* QP::convert(OneStepNSProblem* osnsp)
{
  cout << "QP::convert (DynamicalSystem* osnsp)" << endl;
  QP* qp = dynamic_cast<QP*>(osnsp);
  return qp;
}

// Default private constructor
QP::QP(): OneStepNSProblem(), Q(NULL), p(NULL),
  isQAllocatedIn(false), isPAllocatedIn(false)
{
  nspbType = QP_OSNSP;
}
