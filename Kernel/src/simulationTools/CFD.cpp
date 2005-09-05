#include "CFD.h"
using namespace std;

CFD::CFD(): OneStepNSProblem(), nCfd(0), w(NULL), z(NULL), M(NULL), q(NULL),
  isWAllocatedIn(false), isZAllocatedIn(false), isMAllocatedIn(false), isQAllocatedIn(false)
{
  nspbType = CFD_OSNSP;
}

CFD::CFD(OneStepNSProblemXML* osnspbxml, Strategy * newStrat):
  OneStepNSProblem(osnspbxml, newStrat), nCfd(0), w(NULL), z(NULL), M(NULL), q(NULL),
  isWAllocatedIn(true), isZAllocatedIn(true), isMAllocatedIn(true), isQAllocatedIn(true)
{
  nspbType = CFD_OSNSP;
  if (osnspbxml != NULL)
  {
    // no getter-xml for nCfd ...
    nCfd = ((static_cast<CFDXML*>(onestepnspbxml))->getQ()).size();
    n = nCfd;
    w = new SimpleVector(nCfd);
    z = new SimpleVector(nCfd);
    M = new SiconosMatrix(nCfd, nCfd);
    q = new SimpleVector(nCfd);
    *M = (static_cast<CFDXML*>(onestepnspbxml))->getM();
    *q = (static_cast<CFDXML*>(onestepnspbxml))->getQ();
  }
  else RuntimeException::selfThrow("CFD: xml constructor, xml file=NULL");
}

CFD::~CFD()
{
  if (isWAllocatedIn)
  {
    delete w;
    w = NULL;
  }
  if (isZAllocatedIn)
  {
    delete z;
    z = NULL;
  }
  if (isMAllocatedIn)
  {
    delete M;
    M = NULL;
  }
  if (isQAllocatedIn)
  {
    delete q;
    q = NULL;
  }
}

void CFD::preCFD(const double& time)
{
  IN("CFD::formalize(void)\n");
  OUT("CFD::formalize(void)\n");

  // formalisations specific to CFD problem
  // ...
}


void CFD::compute(const double& time)
{
  IN("CFD::compute(void)\n");
  OUT("CFD::compute(void)\n");
}


void CFD::computeM(void)
{
  OUT("CFD::computeM(void)\n");
}

void CFD::computeQ(const double& time)
{
  IN("CFD::computeQ(void)\n");
  OUT("CFD::computeQ(void)\n");
}

void CFD::display() const
{
  cout << "------------------------------------------------------" << endl;
  cout << "____ data of the CFD " << endl;
  cout << "| nCfd : " << nCfd << endl;
  cout << "| CFD Matrix M  : " << endl;
  if (M != NULL) M->display();
  else cout << "-> NULL" << endl;
  cout << "| CFD vector q : " << endl;
  if (q != NULL) q->display();
  else cout << "-> NULL" << endl;
  cout << "____________________________" << endl;
  cout << "------------------------------------------------------" << endl;
}

void CFD::saveNSProblemToXML()
{
  IN("CFD::saveNSProblemToXML\n");
  OneStepNSProblem::saveNSProblemToXML();
  if (onestepnspbxml != NULL)
  {
    (static_cast<CFDXML*>(onestepnspbxml))->setM(*M);
    (static_cast<CFDXML*>(onestepnspbxml))->setQ(*q);
  }
  else RuntimeException::selfThrow("CFD::saveNSProblemToXML - OneStepNSProblemXML object not exists");
  OUT("CFD::saveNSProblemToXML\n");
}

void CFD::saveMToXML()
{
  IN("CFD::saveMToXML\n");
  if (onestepnspbxml != NULL)
  {
    (static_cast<CFDXML*>(onestepnspbxml))->setM(*M);
  }
  else RuntimeException::selfThrow("CFD::saveMToXML - OneStepNSProblemXML object not exists");
  OUT("CFD::saveMToXML\n");
}

void CFD::saveQToXML()
{
  IN("CFD::saveQToXML\n");
  if (onestepnspbxml != NULL)
  {
    (static_cast<CFDXML*>(onestepnspbxml))->setQ(*q);
  }
  else RuntimeException::selfThrow("CFD::saveQToXML - OneStepNSProblemXML object not exists");
  OUT("CFD::saveQToXML\n");
}

CFD* CFD::convert(OneStepNSProblem* osnsp)
{
  cout << "CFD::convert (OneStepNSProblem* osnsp)" << endl;
  CFD* lcp = dynamic_cast<CFD*>(osnsp);
  return lcp;
}
