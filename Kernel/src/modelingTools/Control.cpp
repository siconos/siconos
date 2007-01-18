/* Siconos-Kernel version 2.0.1, Copyright INRIA 2005-2006.
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
#include "Control.h"

using namespace std;

void Control::initParameter(const std::string id)
{
  if (parametersList[id] == NULL)
  {
    parametersList[id] = new SimpleVector(1);
    string alloc = "parameter_for_" + id;
    isAllocatedIn[alloc] = true;
  }
}

// ===== CONSTRUCTORS =====

// Default constructor (protected)
Control::Control(const string type): type(type), uSize(0), u(NULL), T(NULL), nsds(NULL), ds(NULL),
  computeUPtr(NULL), computeTPtr(NULL)
{}

// copy constructor
Control::Control(const Control& newDS): type(type), uSize(0), u(NULL), T(NULL), nsds(NULL), ds(NULL),
  computeUPtr(NULL), computeTPtr(NULL)
{}

// --- Destructor ---
Control::~Control()
{
  if (isAllocatedIn["u"]) delete u;
  u = NULL;
  if (isAllocatedIn["T"]) delete T;
  T = NULL;

  map<string, SiconosVector*>::iterator it;
  for (it = parametersList.begin(); it != parametersList.end(); ++it)
  {
    string alloc = "parameter_for_" + it->first;
    if (isAllocatedIn[alloc]) delete it->second;
  }
  parametersList.clear();
}

void  Control::setUSize(const unsigned int newUSize)
{
  if (isAllocatedIn["u"]) delete u;
  uSize = newUSize;
  u = new SimpleVector(uSize);
  isAllocatedIn["u"] = true;
  u->zero();
}

// Three steps to set u:
//  - Check if uSize has been given (default value=0 in all constructors)
//  - Allocate memory for u, if necessary
//  - Set value for u
void Control::setU(const SiconosVector& newValue)
{
  if (uSize == 0 || newValue.size() != uSize)
    RuntimeException::selfThrow("Control::setU - inconsistent sizes between u input and uSize - Maybe you forget to set uSize?");
  if (u != NULL)
    *u = newValue;

  else
  {
    if (newValue.isBlock())
      u = new BlockVector(newValue);
    else
      u = new SimpleVector(newValue);
    isAllocatedIn["u"] = true;
  }
  isPlugin["u"] = false;
}

void Control::setUPtr(SiconosVector* newPtr)
{
  if (uSize == 0 || newPtr->size() != uSize)
    RuntimeException::selfThrow("Control::setUPtr - inconsistent sizes between u input and uSize - Maybe you forget to set uSize?");
  // check dimensions ...

  if (isAllocatedIn["u"]) delete u;
  u = newPtr;
  isAllocatedIn["u"] = false;
  isPlugin["u"] = false;
}

void Control::setT(const SiconosMatrix& newValue)
{
  // check dimensions ...
  if (uSize == 0 || newValue.size(1) != uSize || newValue.size(0) != ds->getDim())
    RuntimeException::selfThrow("Control::setT - inconsistent sizes between T input, uSize and/or n - Maybe you forget to set n or uSize?");

  if (T != NULL)
    *T = newValue;
  else
  {
    T = new SimpleMatrix(newValue);
    isAllocatedIn["T"] = true;
  }
  isPlugin["T"] = false;
}

void Control::setTPtr(SiconosMatrix *newPtr)
{
  // check dimensions ...
  if (uSize == 0 || newPtr->size(1) != uSize || newPtr->size(0) != ds->getDim())
    RuntimeException::selfThrow("Control::setTPtr - inconsistent sizes between T input, uSize and/or n - Maybe you forget to set n or uSize?");

  if (isAllocatedIn["T"]) delete T;
  T = newPtr;
  isAllocatedIn["T"] = false;
  isPlugin["T"] = false;
}


void Control::setComputeUFunction(const string pluginPath, const string functionName)
{
  // since u is not allocated by default, memory must be reserved for it
  if (uSize == 0)
    RuntimeException::selfThrow("Control::setComputeUFunction - uSize is equal to 0 - Maybe you forget to set it?");

  if (u == NULL)
  {
    u = new SimpleVector(uSize);
    isAllocatedIn["u"] = true;
    u->zero();
  }

  initParameter("u");

  computeUPtr = NULL;
  cShared.setFunction(&computeUPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  pluginNames["u"] = plugin + ":" + functionName;
  isPlugin["u"] = true;
}

void Control::setComputeTFunction(const string pluginPath, const string functionName)
{
  // since T is not allocated by default, memory must be reserved for it
  if (uSize == 0)
    RuntimeException::selfThrow("Control::setComputeUFunction - uSize is equal to 0 - Maybe you forget to set it?");

  if (T == NULL)
  {
    T = new SimpleMatrix(ds->getDim(), uSize);
    isAllocatedIn["T"] = true;
    T->zero();
  }

  initParameter("T");

  computeTPtr = NULL;
  cShared.setFunction(&computeTPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  pluginNames["T"] = plugin + ":" + functionName;
  isPlugin["T"] = true;
}

void Control::setParameters(const std::map<string, SiconosVector*>& newMap)
{
  // copy!!

  map<string, SiconosVector*>::const_iterator it;
  for (it = newMap.begin(); it != newMap.end(); ++it)
  {
    parametersList[it->first] = new SimpleVector(*(it->second));
    string alloc = "parameter_for_" + it->first;
    isAllocatedIn[alloc] = true;
  }
}

void Control::setParameter(const SiconosVector& newValue, const string id)
{
  parametersList[id] = new SimpleVector(newValue);
  string alloc = "parameter_for_" + id;
  isAllocatedIn[alloc] = true;
}

void Control::setParameterPtr(SiconosVector *newPtr, const string id)
{
  parametersList[id] = newPtr;
  string alloc = "parameter_for_" + id;
  isAllocatedIn[alloc] = false;
}

void Control::computeU(const double time)
{
  if (isPlugin["u"])
  {
    if (computeUPtr == NULL)
      RuntimeException::selfThrow("Control::computeU() is not linked to a plugin function");
    if (u == NULL)
      RuntimeException::selfThrow("Control::computeU(), warning: u = NULL");

    SiconosVector* param = parametersList["u"];
    SiconosVector *x = ds->getXPtr();
    unsigned int n = ds->getDim();
    computeUPtr(uSize, n, time, &((*x)(0)), &(*u)(0), &(*param)(0));
  }
  // else nothing!
}

void Control::computeU(const double time, SiconosVector* xx)
{
  if (isPlugin["u"])
  {
    if (computeUPtr == NULL)
      RuntimeException::selfThrow("Control::computeU() is not linked to a plugin function");
    if (u == NULL)
      RuntimeException::selfThrow("Control::computeU(), warning: u = NULL");
    SiconosVector* param = parametersList["u"];
    unsigned int n = ds->getDim();
    computeUPtr(uSize, n, time, &(*xx)(0), &(*u)(0), &(*param)(0));
  }
  // else nothing!
}

void Control::computeT()
{
  if (isPlugin["T"])
  {
    if (computeTPtr == NULL)
      RuntimeException::selfThrow("Control::computeT() is not linked to a plugin function");
    if (T == NULL)
      RuntimeException::selfThrow("Control::computeT(), warning: T = NULL");
    SiconosVector* param = parametersList["T"];
    SiconosVector *x = ds->getXPtr();
    unsigned int n = ds->getDim();
    computeTPtr(uSize, n, &((*x)(0)), &(*T)(0, 0), &(*param)(0));
  }
  // else nothing!
}

void Control::display() const
{
  cout << " ===== Control display =====" << endl;
  cout << " ===========================" << endl;
}

