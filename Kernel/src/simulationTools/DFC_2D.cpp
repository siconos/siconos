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
#include "DFC_2D.h"
using namespace std;

DFC_2D::DFC_2D(): OneStepNSProblem(), nDfc_2D(0), w(NULL), z(NULL), M(NULL), q(NULL),
  isWAllocatedIn(false), isZAllocatedIn(false), isMAllocatedIn(false), isQAllocatedIn(false)
{
  nspbType = DFC_2D_OSNSP;
}

DFC_2D::DFC_2D(OneStepNSProblemXML* osnspbxml, Strategy * newStrat):
  OneStepNSProblem(osnspbxml, newStrat), nDfc_2D(0), w(NULL), z(NULL), M(NULL), q(NULL),
  isWAllocatedIn(true), isZAllocatedIn(true), isMAllocatedIn(true), isQAllocatedIn(true)
{
  nspbType = DFC_2D_OSNSP;
  if (osnspbxml != NULL)
  {
    // no getter-xml for nDfc_2D ...
    nDfc_2D = ((static_cast<DFC_2DXML*>(onestepnspbxml))->getQ()).size();
    n = nDfc_2D;
    w = new SimpleVector(nDfc_2D);
    z = new SimpleVector(nDfc_2D);
    M = new SiconosMatrix(nDfc_2D, nDfc_2D);
    q = new SimpleVector(nDfc_2D);
    *M = (static_cast<DFC_2DXML*>(onestepnspbxml))->getM();
    *q = (static_cast<DFC_2DXML*>(onestepnspbxml))->getQ();
  }
  else RuntimeException::selfThrow("DFC_2D: xml constructor, xml file=NULL");
}

DFC_2D::~DFC_2D()
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

void DFC_2D::preDFC_2D(const double& time)
{
  IN("DFC_2D::formalize(void)\n");
  OUT("DFC_2D::formalize(void)\n");

  // formalisations specific to DFC_2D problem
  // ...
}


void DFC_2D::compute(const double& time)
{
  IN("DFC_2D::compute(void)\n");
  OUT("DFC_2D::compute(void)\n");
}


void DFC_2D::computeM(void)
{
  OUT("DFC_2D::computeM(void)\n");
}

void DFC_2D::computeQ(const double& time)
{
  IN("DFC_2D::computeQ(void)\n");
  OUT("DFC_2D::computeQ(void)\n");
}

void DFC_2D::display() const
{
  cout << "------------------------------------------------------" << endl;
  cout << "____ data of the DFC_2D " << endl;
  cout << "| nDfc_2D : " << nDfc_2D << endl;
  cout << "| DFC_2D Matrix M  : " << endl;
  if (M != NULL) M->display();
  else cout << "-> NULL" << endl;
  cout << "| DFC_2D vector q : " << endl;
  if (q != NULL) q->display();
  else cout << "-> NULL" << endl;
  cout << "____________________________" << endl;
  cout << "------------------------------------------------------" << endl;
}

void DFC_2D::saveNSProblemToXML()
{
  IN("DFC_2D::saveNSProblemToXML\n");
  OneStepNSProblem::saveNSProblemToXML();
  if (onestepnspbxml != NULL)
  {
    (static_cast<DFC_2DXML*>(onestepnspbxml))->setM(*M);
    (static_cast<DFC_2DXML*>(onestepnspbxml))->setQ(*q);
  }
  else RuntimeException::selfThrow("DFC_2D::saveNSProblemToXML - OneStepNSProblemXML object not exists");
  OUT("DFC_2D::saveNSProblemToXML\n");
}

void DFC_2D::saveMToXML()
{
  IN("DFC_2D::saveMToXML\n");
  if (onestepnspbxml != NULL)
  {
    (static_cast<DFC_2DXML*>(onestepnspbxml))->setM(*M);
  }
  else RuntimeException::selfThrow("DFC_2D::saveMToXML - OneStepNSProblemXML object not exists");
  OUT("DFC_2D::saveMToXML\n");
}

void DFC_2D::saveQToXML()
{
  IN("DFC_2D::saveQToXML\n");
  if (onestepnspbxml != NULL)
  {
    (static_cast<DFC_2DXML*>(onestepnspbxml))->setQ(*q);
  }
  else RuntimeException::selfThrow("DFC_2D::saveQToXML - OneStepNSProblemXML object not exists");
  OUT("DFC_2D::saveQToXML\n");
}

DFC_2D* DFC_2D::convert(OneStepNSProblem* osnsp)
{
  cout << "DFC_2D::convert (OneStepNSProblem* osnsp)" << endl;
  DFC_2D* lcp = dynamic_cast<DFC_2D*>(osnsp);
  return lcp;
}
