/* Siconos-Kernel version 1.1.4, Copyright INRIA 2005-2006.
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

// PRIVATE METHODS

// --- Default constructor ---
TimeDiscretisation::TimeDiscretisation(): h(0.0), nSteps(0), tk(NULL), hMin(0.0), hMax(0.0), constant(0),
  timeDiscretisationXML(NULL), k(0), strategy(NULL), isTkAllocatedIn(false), tdCase(0)
{}

void TimeDiscretisation::checkCase(const bool& hasTk, const bool& hasH, const bool& hasNSteps, const bool& hasT)
{
  if (hasTk && hasT) tdCase = 1;     // Case of constructor I
  else if (hasTk && !hasT) tdCase = 2; // Case of constructor I
  else if (hasH &&  hasNSteps &&  hasT) tdCase = 3; // Constructor II
  else if (hasH &&  hasNSteps && !hasT) tdCase = 4; // Constructor II
  else if (!hasH &&  hasNSteps &&  hasT) tdCase = 5; // Constructor III
  else if (!hasH &&  hasNSteps && !hasT) tdCase = 6; // Constructor III
  else if (hasH && !hasNSteps &&  hasT) tdCase = 7; // Constructor IV
  else if (hasH && !hasNSteps && !hasT) tdCase = 8; // Constructor IV
  else tdCase = 0;
}

void TimeDiscretisation::compute(const int& tdCase)
{
  // load time min value from the model
  double t0 = getT0();
  double T;
  double tol = 1e-12;
  switch (tdCase)
  {
  case 1: // Given: tk, T - Compute h, nSteps
    // check tk(0) = t0
    if (fabs((*tk)(0) - t0) > tol) RuntimeException::selfThrow("TimeDiscretisation::compute - t0 and tk(0) are different");
    nSteps = tk->size() - 1;
    // get final T from Model and check tk(nSteps) = T
    T = getT();
    if (fabs((*tk)(nSteps) - T) > tol) RuntimeException::selfThrow("TimeDiscretisation::compute - T and tk(N) are different");
    // compute h
    h = (*tk)(1) - (*tk)(0);
    break;
  case 2: // Given: tk - Compute h, nSteps and T
    //assert((*tk)(0) == t0);
    if (fabs((*tk)(0) - t0) > tol) RuntimeException::selfThrow("TimeDiscretisation::compute - t0 and tk(0) are different");
    nSteps = tk->size() - 1;
    // set T in the Model
    setT((*tk)(nSteps));
    // compute h
    h = (*tk)(1) - (*tk)(0);
    break;
  case 3:// Given: h, nSteps, T => too many inputs
    RuntimeException::selfThrow("TimeDiscretisation::compute - redundant input values. Only two are needed among T, h and nSteps.");
  case 4: // Given: h, nSteps - Compute T, tk
    tk = new SimpleVector(nSteps + 1);
    isTkAllocatedIn = true;
    // set T in the Model and in tk
    setT(t0 + nSteps * h);
    // compute tk
    (*tk)(0) = t0;
    for (unsigned int i = 1; i < nSteps; ++i)(*tk)(i) = t0 + i * h;
    break;
  case 5: // Given: nSteps, T - Compute h, tk
    // get T from Model and compute h
    T  = getT();
    h = (T - t0) / nSteps;
    // compute tk
    tk = new SimpleVector(nSteps + 1);
    isTkAllocatedIn = true;
    (*tk)(0) = t0;
    (*tk)(nSteps) = getT();
    for (unsigned int i = 1; i < nSteps; ++i)(*tk)(i) = t0 + i * h;
    break;
  case 6: // Given: nSteps => T is missing
    RuntimeException::selfThrow("TimeDiscretisation::compute - T is required in Model, since only nSteps is provided.");
  case 7: // Given: h, T - Compute nSteps, tk
    T  = getT();
    nSteps = (unsigned int)ceil((T - t0) / h);
    // Compute tk
    tk = new SimpleVector(nSteps + 1);
    isTkAllocatedIn = true;
    (*tk)(0) = t0;
    (*tk)(nSteps) = getT();
    for (unsigned int i = 1; i < nSteps; i++)(*tk)(i) = t0 + i * h;
    break;
  case 8: // Given: h => T is missing
    RuntimeException::selfThrow("TimeDiscretisation::constructor - T is required in Model, since only h is provided.");
  default:
    RuntimeException::selfThrow("TimeDiscretisation::compute - wrong scheme, wrong number of input data");
  }

  // compute hMin/hMax. hMin may differ from h when h and T are fixed and (T-t0)/nSteps is not an integer.
  hMin = (*tk)(nSteps) - (*tk)(nSteps - 1);
  hMax = h;
}

// PUBLIC METHODS

// --- CONSTRUCTORS ---

// IO Constructors -> XML
TimeDiscretisation::TimeDiscretisation(TimeDiscretisationXML * tdXML, Strategy* str):
  h(0.0), nSteps(0), tk(NULL), hMin(0.0), hMax(0.0), constant(0), timeDiscretisationXML(tdXML),
  k(0), strategy(str), isTkAllocatedIn(false), tdCase(0)
{
  if (strategy == NULL) RuntimeException::selfThrow("TimeDiscretisation::xml constructor - Strategy=NULL!");
  strategy->setTimeDiscretisationPtr(this);

  if (timeDiscretisationXML != NULL)
  {
    // --- Check what are the given data ---
    bool hasNSteps = timeDiscretisationXML->hasN();
    bool hasH = timeDiscretisationXML->hasH();
    bool hasTk = timeDiscretisationXML->hasTk();

    // Eliminate cases with too many inputs
    if ((hasTk && hasH) || (hasTk && hasNSteps))
      RuntimeException::selfThrow("TimeDiscretisation: xml constructor - Too many input data, some of them are useless.");

    // --- Read the data ---
    if (hasH) h = timeDiscretisationXML->getH();
    if (hasNSteps) nSteps = timeDiscretisationXML->getN();
    if (timeDiscretisationXML->hasHMin()) hMin = timeDiscretisationXML->getHMin();
    if (timeDiscretisationXML->hasHMax()) hMax = timeDiscretisationXML->getHMax();
    constant = timeDiscretisationXML->isConstant();

    // --- Check if given data are coherent and switch to the right case ---
    checkCase(hasTk, hasH, hasNSteps, hasT());
    // --- Compute the corresponding time discretisation scheme ---
    if (tdCase == 1 || tdCase == 2)
    {
      tk = new SimpleVector(timeDiscretisationXML->getTk().size());
      isTkAllocatedIn = true;
      *tk = timeDiscretisationXML->getTk();
    }
    compute(tdCase);
  }
  else RuntimeException::selfThrow("TimeDiscretisation: xml constructor - TimeDiscretisationXML = NULL");
}

// --- Straightforward constructors ---

// Provide tk and strategy (I)
TimeDiscretisation::TimeDiscretisation(SimpleVector *newTk, Strategy* str):
  h(0), nSteps(0), tk(NULL), hMin(0), hMax(0), constant(1),
  timeDiscretisationXML(NULL), k(0), strategy(str), isTkAllocatedIn(true), tdCase(0)
{
  if (strategy == NULL) RuntimeException::selfThrow("TimeDiscretisation:: data constructor - Strategy=NULL!");
  strategy->setTimeDiscretisationPtr(this);
  // Allocate memory for tk and fill it
  tk = new SimpleVector(newTk->size());
  *tk = *newTk;
  checkCase(1, 0, 0, hasT());
  compute(tdCase);
}

// Provide h and nSteps, calculate tk (II)
TimeDiscretisation::TimeDiscretisation(const double& newH, const unsigned int& newNSteps, Strategy* str):
  h(newH), nSteps(newNSteps), tk(NULL), hMin(newH), hMax(newH), constant(1),
  timeDiscretisationXML(NULL), k(0), strategy(str), isTkAllocatedIn(false), tdCase(0)
{
  if (strategy == NULL) RuntimeException::selfThrow("TimeDiscretisation::data constructor - Strategy=NULL!");
  strategy->setTimeDiscretisationPtr(this);
  checkCase(0, 1, 1, hasT());
  compute(tdCase);
}

// Provide nSteps, calculate h and tk (III)
TimeDiscretisation::TimeDiscretisation(const unsigned int& newNSteps, Strategy* str):
  h(0), nSteps(newNSteps), tk(NULL), hMin(0), hMax(0), constant(1),
  timeDiscretisationXML(NULL), k(0), strategy(str), isTkAllocatedIn(false), tdCase(0)
{
  if (strategy == NULL) RuntimeException::selfThrow("TimeDiscretisation::data constructor - Strategy=NULL!");
  strategy->setTimeDiscretisationPtr(this);
  checkCase(0, 0, 1, hasT());
  compute(tdCase);
}

// Provide h, calculate nSteps and tk (IV)
TimeDiscretisation::TimeDiscretisation(const double& newH, Strategy* str):
  h(newH), nSteps(0), tk(NULL), hMin(newH), hMax(newH), constant(1),
  timeDiscretisationXML(NULL), k(0), strategy(str), isTkAllocatedIn(false), tdCase(0)
{
  if (strategy == NULL) RuntimeException::selfThrow("TimeDiscretisation::data constructor - Strategy=NULL!");
  strategy->setTimeDiscretisationPtr(this);
  checkCase(0, 1, 0, hasT());
  compute(tdCase);
}

// --- Destructor ---
TimeDiscretisation::~TimeDiscretisation()
{
  if (isTkAllocatedIn) delete tk;
  tk = NULL;
}

void TimeDiscretisation::setH(const double& newH)
{
  h = newH;
  // depending on initial way of construction, the whole time discretisation
  // is re-built or not.
  if (tdCase == 1 || tdCase == 2 || tdCase == 5 || tdCase == 6)
    checkCase(0, 1, 0, hasT());
  // else nothing, we keep the same tdCase value. This corresponds to the cases where h was an input.

  compute(tdCase);
}

void TimeDiscretisation::setNSteps(const unsigned int& newNSteps)
{
  nSteps = newNSteps;
  // depending on initial way of construction, the whole time discretisation
  // is re-built or not.
  if (tdCase == 1 || tdCase == 2 || tdCase == 7 || tdCase == 8)
    checkCase(0, 0, 1, hasT());
  // else nothing, we keep the same tdCase value. This corresponds to the cases where nSteps was an input.

  compute(tdCase);
}

void TimeDiscretisation::setTk(const SimpleVector& newValue)
{
  cout << " /!\\ TimeDiscretisation::setTk - Warning: you set a new tk vector, this will change nSteps and h values /!\\ " << endl;

  if (tk->size() != newValue.size())
  {
    if (isTkAllocatedIn) delete tk;
    tk = new SimpleVector(newValue);
    isTkAllocatedIn = true;
  }
  else
    *tk = newValue;

  // Update other time discretisation values
  checkCase(1, 0, 0, hasT());
  compute(tdCase);
}

void TimeDiscretisation::setTkPtr(SimpleVector *newPtr)
{
  cout << " /!\\ TimeDiscretisation::setTk - Warning: you set a new tk vector, this will change nSteps and h values /!\\ " << endl;

  if (isTkAllocatedIn) delete tk;
  tk = newPtr;
  isTkAllocatedIn = false;
  // Update other time discretisation values
  checkCase(1, 0, 0, hasT());
  compute(tdCase);
}

// --- Functions to manage t0 and T (data of the model)
const double TimeDiscretisation::getT0() const
{
  return strategy->getModelPtr()->getT0();
}

void TimeDiscretisation::setT0(const double& newValue)
{
  if (strategy->getModelPtr() == NULL)
    RuntimeException::selfThrow("TimeDiscretisation::setT0, no Model linked to this time discretisation");
  (*tk)(0) = newValue;
  strategy->getModelPtr()->t0 = newValue;
  if (hasT() && newValue >= getT())
    RuntimeException::selfThrow("TimeDiscretisation::setT0, input value for t0 greater than T (final time).");

  // recompute time discretisation values with new t0.
  compute(tdCase);
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
  (*tk)(nSteps) = newValue;
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
  if (timeDiscretisationXML != NULL)
  {
    timeDiscretisationXML->setH(h);
    timeDiscretisationXML->setN(nSteps);
    timeDiscretisationXML->setTkNode(*tk);
    timeDiscretisationXML->setHMin(hMin);
    timeDiscretisationXML->setHMax(hMax);
  }
  else RuntimeException::selfThrow("TimeDiscretisation::saveTimeDiscretisationToXML - TimeDiscretisationXML object not exists");
}

