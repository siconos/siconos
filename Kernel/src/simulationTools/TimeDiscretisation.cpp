/* Siconos-Kernel version 1.1.2, Copyright INRIA 2005-2006.
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
#include "TimeDiscretisation.h"
using namespace std;

// --- CONSTRUCTORS ---

// IO Constructors -> XML
TimeDiscretisation::TimeDiscretisation(TimeDiscretisationXML * tdXML, Strategy* str):
  h(0.0), nSteps(0), tk(NULL), hMin(0.0), hMax(0.0), constant(0), timeDiscretisationXML(tdXML),
  k(0), strategy(str), isTkAllocatedIn(false)
{
  if (strategy == NULL) RuntimeException::selfThrow("TimeDiscretisation::xml constructor - Strategy=NULL!");
  if (timeDiscretisationXML != NULL)
  {
    // --- Check what are the given data ---
    bool hasNSteps = timeDiscretisationXML->hasN();
    bool hasH = timeDiscretisationXML->hasH();
    bool hasTk = timeDiscretisationXML->hasTk();

    // --- Read the data ---
    if (hasH) h = timeDiscretisationXML->getH();
    if (hasNSteps) nSteps = timeDiscretisationXML->getN();
    if (timeDiscretisationXML->hasHMin()) hMin = timeDiscretisationXML->getHMin();
    if (timeDiscretisationXML->hasHMax()) hMax = timeDiscretisationXML->getHMax();
    constant = timeDiscretisationXML->isConstant();

    // --- Check if given data are coherent and switch to the right case ---
    int tdCase = checkTimeDiscretisationCase(hasTk, hasH, hasNSteps, hasT());
    // --- Compute the corresponding time discretisation scheme ---
    if (tdCase == 1 || tdCase == 2)
    {
      tk = new SimpleVector(timeDiscretisationXML->getTk().size());
      isTkAllocatedIn = true;
      *tk = timeDiscretisationXML->getTk();
    }
    computeTimeDiscretisation(tdCase);
  }
  else RuntimeException::selfThrow("TimeDiscretisation: xml constructor - TimeDiscretisationXML = NULL");
}

// --- Straightforward constructors ---

// Provide tk and strategy
TimeDiscretisation::TimeDiscretisation(SimpleVector *newTk, Strategy* str):
  h(0), nSteps(0), tk(NULL), hMin(0), hMax(0), constant(1),
  timeDiscretisationXML(NULL), k(0), strategy(str), isTkAllocatedIn(true)
{
  if (strategy == NULL) RuntimeException::selfThrow("TimeDiscretisation:: data constructor - Strategy=NULL!");
  // Allocate memory for tk and fill it
  tk = new SimpleVector(newTk->size());
  *tk = *newTk;
  int tdCase = checkTimeDiscretisationCase(1, 0, 0, hasT());
  computeTimeDiscretisation(tdCase);
}

// Provide h and nSteps, calculate tk
TimeDiscretisation::TimeDiscretisation(const double& newH, const int& newNSteps, Strategy* str):
  h(newH), nSteps(newNSteps), tk(NULL), hMin(newH), hMax(newH), constant(1),
  timeDiscretisationXML(NULL), k(0), strategy(str), isTkAllocatedIn(false)
{
  if (strategy == NULL) RuntimeException::selfThrow("TimeDiscretisation::data constructor - Strategy=NULL!");
  int tdCase = checkTimeDiscretisationCase(0, 1, 1, hasT());
  computeTimeDiscretisation(tdCase);
}

// Provide nSteps, calculate h and tk
TimeDiscretisation::TimeDiscretisation(const int& newNSteps, Strategy* str):
  h(0), nSteps(newNSteps), tk(NULL), hMin(0), hMax(0), constant(1),
  timeDiscretisationXML(NULL), k(0), strategy(str), isTkAllocatedIn(false)
{
  if (strategy == NULL) RuntimeException::selfThrow("TimeDiscretisation::data constructor - Strategy=NULL!");
  int tdCase = checkTimeDiscretisationCase(0, 0, 1, hasT());
  computeTimeDiscretisation(tdCase);
}

// Provide h, calculate nSteps and tk
TimeDiscretisation::TimeDiscretisation(const double& newH, Strategy* str):
  h(newH), nSteps(0), tk(NULL), hMin(newH), hMax(newH), constant(1),
  timeDiscretisationXML(NULL), k(0), strategy(str), isTkAllocatedIn(false)
{
  if (strategy == NULL) RuntimeException::selfThrow("TimeDiscretisation::data constructor - Strategy=NULL!");
  int tdCase = checkTimeDiscretisationCase(0, 1, 0, hasT());
  computeTimeDiscretisation(tdCase);
}

// --- Destructor ---
TimeDiscretisation::~TimeDiscretisation()
{
  if (isTkAllocatedIn)
  {
    delete tk;
    tk = NULL;
  }
}

// --- Constructors related functions ---
const int TimeDiscretisation::checkTimeDiscretisationCase(const bool& hasTk, const bool& hasH, const bool& hasNSteps, const bool& hasT) const
{
  int tdCase;
  if (hasTk && ! hasH && ! hasNSteps && ! hasT) tdCase = 1;
  else if (hasTk && ! hasH && ! hasNSteps &&  hasT) tdCase = 2;
  else if (! hasTk && hasH &&  hasNSteps && ! hasT) tdCase = 3;
  else if (! hasTk && hasH &&  ! hasNSteps &&  hasT) tdCase = 4;
  else if (! hasTk && ! hasH &&  hasNSteps &&  hasT) tdCase = 5;
  else if (hasTk &&  hasH &&  hasNSteps &&  hasT) tdCase = 6;
  else tdCase = 0;
  return(tdCase);
}

void TimeDiscretisation::computeTimeDiscretisation(const int& tdCase)
{
  // load time min value from the model
  double t0 = getT0();
  double T;
  double tol = 1e-15;
  switch (tdCase)
  {
  case 1: // Given: tk - Compute h, nSteps and T
    //assert((*tk)(0) == t0);
    if ((*tk)(0) != t0) RuntimeException::selfThrow("TimeDiscretisation::constructor - t0 and tk(0) are different");
    nSteps = tk->size() - 1;
    setT((*tk)(nSteps));
    h = (*tk)(1) - (*tk)(0);
    break;
  case 2: // Given: tk, T - Compute h, nSteps
    T = getT();
    if ((*tk)(0) != t0) RuntimeException::selfThrow("TimeDiscretisation::constructor - t0 and tk(0) are different");
    if (tk->getValue(tk->size() - 1) != T) RuntimeException::selfThrow("TimeDiscretisation::constructor - T and tk(N) are different");
    nSteps = tk->size() - 1;
    h = (*tk)(1) - (*tk)(0);
    break;
  case 3: // Given: nSteps, h - Compute T, tk
    setT(t0 + nSteps * h);
    tk = new SimpleVector(nSteps + 1);
    isTkAllocatedIn = true;
    (*tk)(0) = t0;
    (*tk)(nSteps) = getT();
    for (int i = 1; i < nSteps; i++)(*tk)(i) = t0 + i * h;
    break;
  case 4: // Given: h, T - Compute nSteps, tk
    T  = getT();
    nSteps = (int)floor((T - t0) / h);
    if ((t0 + (h * nSteps) - T) >= tol) cout << " Warning, (nSteps x step size) seems to differ from final T?" << endl;
    tk = new SimpleVector(nSteps + 1);
    isTkAllocatedIn = true;
    (*tk)(0) = t0;
    (*tk)(nSteps) = getT();
    for (int i = 1; i < nSteps; i++)(*tk)(i) = t0 + i * h;
    break;
  case 5: // Given: nSteps, T - Compute h, tk
    T  = getT();
    h = (T - t0) / ((double)nSteps);
    tk = new SimpleVector(nSteps + 1);
    isTkAllocatedIn = true;
    (*tk)(0) = t0;
    (*tk)(nSteps) = getT();
    for (int i = 1; i < nSteps; i++)(*tk)(i) = t0 + i * h;
    break;
  case 6: // all the data are known. Just check they are coherent.
    if ((*tk)(0) != t0) RuntimeException::selfThrow("TimeDiscretisation::constructor - incompatible data");
    if (tk->getValue(tk->size() - 1) != T) RuntimeException::selfThrow("TimeDiscretisation::constructor - incompatible data");
    if ((t0 + h * nSteps) != T)  RuntimeException::selfThrow("TimeDiscretisation::constructor - incompatible data");
    if (h != (*tk)(1) - (*tk)(0)) RuntimeException::selfThrow("TimeDiscretisation::constructor - incompatible data");
    if (nSteps != (int)(tk->size()) - 1)  RuntimeException::selfThrow("TimeDiscretisation::constructor - incompatible data");
    break;
  default:
    RuntimeException::selfThrow("TimeDiscretisation::constructor - wrong scheme, wrong number of data");
  }
}

// --- Functions to manage t0 and T (data of the model)
const double TimeDiscretisation::getT0() const
{
  return strategy->getModelPtr()->getT0();
}

const bool TimeDiscretisation::hasT() const
{
  // When the model is filled, T is set to negative value if not given by user.
  if (getT() < 0) return(0);
  else return (1);
}

const double TimeDiscretisation::getT() const
{
  return strategy->getModelPtr()->getFinalT();
}

inline void TimeDiscretisation::setT(const double& newValue)
{
  strategy->getModelPtr()->setFinalT(newValue);
}

// --- Other functions ---

void TimeDiscretisation::display() const
{
  cout << "-----------------------------------------------------" << endl;
  cout << "____ data of the TimeDiscretisation " << endl;
  cout << "| h : " << h << endl;
  cout << "| nSteps : " << nSteps << endl;
  cout << "| tk " << endl;
  if (tk != NULL) tk->display();
  else cout << "-> NULL" << endl;
  cout << "| hMin : " << hMin << endl;
  cout << "| hMax : " << hMax << endl;
  cout << "| constant : " << constant << endl;
  cout << "-----------------------------------------------------" << endl << endl;
}

// --- XML functions ---

void TimeDiscretisation::saveTimeDiscretisationToXML()
{
  IN("TimeDiscretisation::saveTimeDiscretisationToXML\n");
  if (timeDiscretisationXML != NULL)
  {
    timeDiscretisationXML->setH(h);
    timeDiscretisationXML->setN(nSteps);
    timeDiscretisationXML->setTkNode(*tk);
    timeDiscretisationXML->setHMin(hMin);
    timeDiscretisationXML->setHMax(hMax);
  }
  else RuntimeException::selfThrow("TimeDiscretisation::saveTimeDiscretisationToXML - TimeDiscretisationXML object not exists");
  OUT("TimeDiscretisation::saveTimeDiscretisationToXML\n");
}

// --- Default (private) constructor ---
TimeDiscretisation::TimeDiscretisation(): h(0.0), nSteps(0), tk(NULL), hMin(0.0), hMax(0.0), constant(0),
  timeDiscretisationXML(NULL), k(0), strategy(NULL), isTkAllocatedIn(false)
{}
