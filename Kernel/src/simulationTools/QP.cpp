/* Siconos version 1.0, Copyright INRIA 2005.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
*/

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
