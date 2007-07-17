/* Siconos-Kernel version 2.1.1, Copyright INRIA 2005-2006.
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
#include "TimeDiscretisationXML.h"
#include "RuntimeException.h"
#include "Model.h"
using namespace std;

// PRIVATE METHODS

// --- Default constructor ---
TimeDiscretisation::TimeDiscretisation(): h(0.0), nSteps(0), tk(NULL), hMin(0.0), hMax(0.0), constant(0),
  timeDiscretisationXML(NULL), model(NULL), isTkAllocatedIn(false), tdCase(0), isUpToDate(false)
{}

void TimeDiscretisation::setCase(bool hasTk, bool hasH, bool hasNSteps)
{
  if (hasTk) tdCase = 1; // Case of constructor I
  else if (hasH &&  hasNSteps) tdCase = 2; // Constructor II
  else if (!hasH &&  hasNSteps) tdCase = 3; // Constructor III
  else if (hasH && !hasNSteps) tdCase = 4; // Constructor IV
  else tdCase = 0;
  isUpToDate = false;
}

// --- CONSTRUCTORS ---

// IO Constructors -> XML
TimeDiscretisation::TimeDiscretisation(TimeDiscretisationXML * tdXML, Model* m):
  h(0.0), nSteps(0), tk(NULL), hMin(0.0), hMax(0.0), constant(0), timeDiscretisationXML(tdXML),
  model(m), isTkAllocatedIn(false), tdCase(0), isUpToDate(false)
{
  if (timeDiscretisationXML == NULL)
    RuntimeException::selfThrow("TimeDiscretisation: xml constructor - TimeDiscretisationXML = NULL");

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
  if (hasTk) setTk(timeDiscretisationXML->getTk());

  constant = timeDiscretisationXML->isConstant();

  // --- set tdCase value.
  setCase(hasTk, hasH, hasNSteps);
}

// --- Straightforward constructors ---

// Provide tk and model (I)
TimeDiscretisation::TimeDiscretisation(SiconosVector *newTk, Model* m):
  h(0), nSteps(0), tk(NULL), hMin(0), hMax(0), constant(1),
  timeDiscretisationXML(NULL), model(m), isTkAllocatedIn(true), tdCase(0), isUpToDate(false)
{
  // Allocate memory for tk and fill it
  tk = new SimpleVector(*newTk);

  setCase(1, 0, 0);
}

// Provide h and nSteps, calculate tk (II)
TimeDiscretisation::TimeDiscretisation(double newH, unsigned int newNSteps, Model* m):
  h(newH), nSteps(newNSteps), tk(NULL), hMin(newH), hMax(newH), constant(1),
  timeDiscretisationXML(NULL), model(m), isTkAllocatedIn(false), tdCase(0), isUpToDate(false)
{
  setCase(0, 1, 1);
}

// Provide nSteps, calculate h and tk (III)
TimeDiscretisation::TimeDiscretisation(unsigned int newNSteps, Model* m):
  h(0), nSteps(newNSteps), tk(NULL), hMin(0), hMax(0), constant(1),
  timeDiscretisationXML(NULL), model(m), isTkAllocatedIn(false), tdCase(0), isUpToDate(false)
{
  setCase(0, 0, 1);
}

// Provide h, calculate nSteps and tk (IV)
TimeDiscretisation::TimeDiscretisation(double newH, Model* m):
  h(newH), nSteps(0), tk(NULL), hMin(newH), hMax(newH), constant(1),
  timeDiscretisationXML(NULL), model(m), isTkAllocatedIn(false), tdCase(0), isUpToDate(false)
{
  setCase(0, 1, 0);
}

// --- Destructor ---
TimeDiscretisation::~TimeDiscretisation()
{
  if (isTkAllocatedIn) delete tk;
  tk = NULL;
}

void TimeDiscretisation::setTk(const SiconosVector& newValue)
{
  if (tk == NULL)
  {
    tk = new SimpleVector(newValue);
    isTkAllocatedIn = true;
  }
  else
  {
    if (tk->size() != newValue.size())
    {
      if (isTkAllocatedIn) delete tk;
      tk = new SimpleVector(newValue);
      isTkAllocatedIn = true;
    }
    else
      *tk = newValue;
  }

  // Update tdCase
  setCase(1, 0, 0);
}

void TimeDiscretisation::setTkPtr(SiconosVector *newPtr)
{
  if (isTkAllocatedIn) delete tk;
  tk = newPtr;
  isTkAllocatedIn = false;
  // Update tdCase
  setCase(1, 0, 0);
}

const bool TimeDiscretisation::hasT() const
{
  // When the model is filled, T is set to negative value if not given by user.
  if (model->getFinalT() < 0) return false;
  else return true;
}

void TimeDiscretisation::initialize()
{

  if (model == NULL)
    RuntimeException::selfThrow("TimeDiscretisation::initialize - Not linked to a valid Model (NULL pointer).");

  if (!isUpToDate) // Initialization only if something as changed in the object (or at the first call ...)
  {
    // load time min value from the model
    double t0 = model->getT0();
    double T;
    double tol = 1e-12;
    switch (tdCase)
    {
    case 1: // Given: tk, possibly T - Compute h, nSteps
      // check tk(0) = t0
      if (fabs((*tk)(0) - t0) > tol) RuntimeException::selfThrow("TimeDiscretisation::initialize - t0 and tk(0) are different");
      nSteps = tk->size() - 1;
      if (hasT()) // T from Model, check tk(nSteps) == T.
      {
        T = model->getFinalT();
        if (fabs((*tk)(nSteps) - T) > tol) RuntimeException::selfThrow("TimeDiscretisation::initialize - T and tk(N) are different");
      }
      else model->setFinalT((*tk)(nSteps)); // Set T in the Model with tk(nSteps).
      // compute h
      h = (*tk)(1) - (*tk)(0);
      break;
    case 2: // Given: h, nSteps - Compute T, tk
      // Check if T is given in the Model, exception if yes.
      if (hasT())
        RuntimeException::selfThrow("TimeDiscretisation::initialize - redundant input values. Only two are needed among T, h and nSteps.");
      // Build tk
      if (isTkAllocatedIn) delete tk;
      tk = new SimpleVector(nSteps + 1);
      isTkAllocatedIn = true;
      for (unsigned int i = 0; i < nSteps + 1; ++i)(*tk)(i) = t0 + i * h;
      // set T in the Model
      model->setFinalT((*tk)(nSteps));
      break;
    case 3: // Given: nSteps, T - Compute h, tk
      // Check if T is given in the Model, exception if no.
      if (!hasT())
        RuntimeException::selfThrow("TimeDiscretisation::initialize - T is required in Model, since only nSteps is provided.");
      // get T from Model and compute h
      T  = model->getFinalT();
      h = (T - t0) / nSteps;
      // compute tk
      if (isTkAllocatedIn) delete tk;
      tk = new SimpleVector(nSteps + 1);
      isTkAllocatedIn = true;
      for (unsigned int i = 0; i < nSteps; ++i)(*tk)(i) = t0 + i * h;
      (*tk)(nSteps) = T;
      break;
    case 4: // Given: h, T - Compute nSteps, tk
      // Check if T is given in the Model, exception if no.
      if (!hasT())
        RuntimeException::selfThrow("TimeDiscretisation::initialize - T is required in Model, since only h is provided.");
      T  = model->getFinalT();
      nSteps = (unsigned int)ceil((T - t0) / h);
      // Compute tk
      if (isTkAllocatedIn) delete tk;
      tk = new SimpleVector(nSteps + 1);
      isTkAllocatedIn = true;
      for (unsigned int i = 0; i < nSteps; i++)(*tk)(i) = t0 + i * h;
      (*tk)(nSteps) = T;
      break;
    default:
      RuntimeException::selfThrow("TimeDiscretisation::initialize - wrong scheme, wrong number of input data");
    }

    // compute hMin/hMax. hMin may differ from h when h and T are fixed and (T-t0)/nSteps is not an integer.
    hMin = (*tk)(nSteps) - (*tk)(nSteps - 1);
    hMax = h;
    isUpToDate = true;
  }
}
// --- Other functions ---

void TimeDiscretisation::display() const
{
  cout << "====> Time Disretisation :" << endl;
  cout << " time step size: " << h << ", number of time steps : " << nSteps << endl;
  cout << "====" << endl;
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

