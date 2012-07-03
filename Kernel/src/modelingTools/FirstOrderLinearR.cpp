/* Siconos-Kernel, Copyright INRIA 2005-2011.
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
* Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
*/
#include "FirstOrderLinearR.hpp"
#include "LinearRXML.hpp"
#include "Interaction.hpp"
#include "Plugin.hpp"

using namespace std;
using namespace RELATION;



FirstOrderLinearR::FirstOrderLinearR():
  FirstOrderR(LinearR)
{
  ;
}

// xml constructor
FirstOrderLinearR::FirstOrderLinearR(SP::RelationXML relxml):
  FirstOrderR(relxml, LinearR)
{

}

// Constructor with C and B plug-in names
FirstOrderLinearR::FirstOrderLinearR(const string& CName, const string& BName):
  FirstOrderR(LinearR)
{
  // Warning: we cannot allocate memory for C/D matrix since no interaction
  // is connected to the relation. This will be done during initialize.
  // We only set the name of the plugin-function and connect it to the user-defined function.
  _pluginJachx->setComputeFunction(CName);
  _pluginJacLg->setComputeFunction(BName);
  //  Plugin::setFunction(_plunginJachx,CName);
  //Plugin::setFunction(_pluginJacLg,BName);
  //  Jach[0].reset(new PluggedMatrix(CName));
  //  Jacg[0].reset(new PluggedMatrix(BName));
}

// Constructor from a complete set of data (plugin)
FirstOrderLinearR::FirstOrderLinearR(const string& CName, const string& DName, const string& FName, const string& EName, const string& BName): FirstOrderR(LinearR)
{
  _pluginJachx->setComputeFunction(CName);
  _pluginJachlambda->setComputeFunction(DName);
  _pluginJacLg->setComputeFunction(BName);
  _pluginf->setComputeFunction(FName);
  _plugine->setComputeFunction(EName);

  //   Plugin::setFunction(_plunginJachx,CName);
  //   Plugin::setFunction(_pluginJachlambda,DName);
  //   Plugin::setFunction(_pluginJacLg,BName);

  //   Plugin::setFunction(_pluginf,FName);
  //   Plugin::setFunction(_plugine,EName);



}

// Minimum data (C, B as pointers) constructor
FirstOrderLinearR::FirstOrderLinearR(SP::SiconosMatrix newC, SP::SiconosMatrix newB):
  FirstOrderR(LinearR)
{
  _jachx = newC;
  _jacglambda = newB;
}

// // Constructor from a complete set of data
FirstOrderLinearR::FirstOrderLinearR(SP::SiconosMatrix  newC, SP::SiconosMatrix newD, SP::SiconosMatrix newF, SP::SiconosVector newE, SP::SiconosMatrix newB):
  FirstOrderR(LinearR)
{
  _jachx = newC;
  _jacglambda = newB;
  _jachlambda = newD;
  _F = newF;
  _e = newE;
}

// Minimum data (C, B as matrices) constructor
/*

FirstOrderLinearR::FirstOrderLinearR(const SiconosMatrix& newC, const SiconosMatrix& newB):
FirstOrderR(LinearR)
{
assert(false&&"FirstOrderLinearR::FirstOrderLinearR copy matrix ???\n");
C = createSPtrSiconosMatrix(newC);
// C.reset(new SiconosMatrix(newC));
B.reset(new SiconosMatrix(newB));
}
// Constructor from a complete set of data (matrices)
FirstOrderLinearR::FirstOrderLinearR(const SiconosMatrix& newC, const SiconosMatrix& newD, const SiconosMatrix& newF, const SiconosVector& newE, const SiconosMatrix& newB):
FirstOrderR(LinearR)
{
C.reset(new SiconosMatrix(newC));
B.reset(newB);
D.reset(newD);
F.reset(newF);
e.reset(newE);

}*/

void FirstOrderLinearR::initialize(Interaction& inter)
{
  // Note: do not call FirstOrderR::initialize to avoid jacobianH and jacobianG allocation.

  // Update data member (links to DS variables)
  if (!_jachx)
    RuntimeException::selfThrow("FirstOrderLinearR::initialize() C is null and is a required input.");
  if (!_jacglambda)
    RuntimeException::selfThrow("FirstOrderLinearR::initialize() B is null and is a required input.");

  // Check if various operators sizes are consistent.
  // Reference: interaction.
  unsigned int sizeY = inter.getSizeOfY();
  unsigned int sizeX = inter.getSizeOfDS();
  unsigned int sizeZ = inter.getSizez();

  // The initialization of each matrix/vector depends on the way the Relation was built ie if the matrix/vector
  // was read from xml or not
  if (_jachx->size(0) == 0)
    _jachx->resize(sizeX, sizeY);
  else
    assert((_jachx->size(0) == sizeY && _jachx->size(1) == sizeX) && "FirstOrderLinearR::initialize , inconsistent size between C and Interaction.");


  if (_jacglambda->size(0) == 0)
    _jacglambda->resize(sizeY, sizeX);
  else
    assert((_jacglambda->size(1) == sizeY && _jacglambda->size(0) == sizeX) && "FirstOrderLinearR::initialize , inconsistent size between B and interaction.");

  // C and B are the minimum inputs. The others may remain null.

  if (_jachlambda)
  {
    if (_jachlambda->size(0) == 0)
      _jachlambda->resize(sizeY, sizeY);
    else
      assert((_jachlambda->size(0) == sizeY || _jachlambda->size(1) == sizeY) && "FirstOrderLinearR::initialize , inconsistent size between C and D.");
  }

  if (!_F && _pluginf->fPtr)
  {
    _F.reset(new SimpleMatrix(sizeY, sizeZ));
  }

  if (_F)
  {
    if (_F->size(0) == 0)
      _F->resize(sizeY, sizeZ);
    else
      assert(((_F->size(0) == sizeY) && (_F->size(1) == sizeZ)) && "FirstOrderLinearR::initialize , inconsistent size between C and F.");
  }
  if (!_e && _plugine->fPtr)
  {
    _e.reset(new SiconosVector(sizeY));
  }

  if (_e)
  {
    if (_e->size() == 0)
      _e->resize(sizeY);
    else
      assert(_e->size() == sizeY && "FirstOrderLinearR::initialize , inconsistent size between C and e.");
  }
}

// // setters




void FirstOrderLinearR::computeC(const double time, Interaction& inter)
{
  if (_pluginJachx->fPtr)
  {
    SiconosVector workZ(*inter.data(z));
    ((FOMatPtr1)(_pluginJachx->fPtr))(time, _jachx->size(0), _jachx->size(1), &(*_jachx)(0, 0), workZ.size(), &(workZ)(0));
    // Copy data that might have been changed in the plug-in call.
    *inter.data(z) = workZ;
  }
}

void FirstOrderLinearR::computeD(const double time, Interaction& inter)
{
  if (_pluginJachlambda->fPtr)
  {
    SiconosVector workZ = *inter.data(z);
    ((FOMatPtr1)(_pluginJachlambda->fPtr))(time, _jachlambda->size(0), _jachlambda->size(1), &(*_jachlambda)(0, 0), workZ.size(), &(workZ)(0));
    // Copy data that might have been changed in the plug-in call.
    *inter.data(z) = workZ;
  }
}

void FirstOrderLinearR::computeF(const double time, Interaction& inter)
{
  if (_pluginf->fPtr)
  {
    SiconosVector workZ = *inter.data(z);
    ((FOMatPtr1)(_pluginf->fPtr))(time, _F->size(0), _F->size(1), &(*_F)(0, 0), workZ.size(), &(workZ)(0));
    // Copy data that might have been changed in the plug-in call.
    *inter.data(z) = workZ;
  }
}

void FirstOrderLinearR::computeE(const double time, Interaction& inter)
{

  if (_plugine->fPtr)
  {
    SiconosVector workZ = *inter.data(z);
    ((FOVecPtr) _plugine->fPtr)(time, _e->size(), &(*_e)(0), workZ.size(), &(workZ)(0));
    // Copy data that might have been changed in the plug-in call.
    *inter.data(z) = workZ;
  }
}

void FirstOrderLinearR::computeb(const double time, Interaction& inter)
{
  if (_pluginJacLg->fPtr)
  {
    SiconosVector workZ = *inter.data(z);
    ((FOMatPtr1) _pluginJacLg->fPtr)(time, _jacglambda->size(0), _jacglambda->size(1), &(*_jacglambda)(0, 0), workZ.size(), &(workZ)(0));
    // Copy data that might have been changed in the plug-in call.
    *inter.data(z) = workZ;
  }
}

void FirstOrderLinearR::computeh(const double time, Interaction& inter)
{
  computeOutput(time, inter, 0);
}

void FirstOrderLinearR::computeg(const double time, Interaction& inter)
{
  computeInput(time, inter, 0);
}
void FirstOrderLinearR::computeOutput(const double time, Interaction& inter, unsigned int level)
{
  computeC(time, inter);
  computeD(time, inter);
  computeF(time, inter);
  computeE(time, inter);

  // Note that the second argument remains unamed since it is not used: for first order systems, we always compute
  // y[0]

  // We get y and lambda of the interaction (pointers)
  SiconosVector& y = *inter.y(0);
  SiconosVector& lambda = *inter.lambda(0);

  // compute y
  if (_jachx)
    prod(*_jachx, *inter.data(x), y);
  else
    y.zero();

  if (_jachlambda)
    prod(*_jachlambda, lambda, y, false);

  if (_e)
    y += *_e;

  if (_F)
    prod(*_F, *inter.data(z), y, false);
}

void FirstOrderLinearR::computeInput(double time, Interaction& inter, unsigned int level)
{
  computeb(time, inter);

  // We get lambda of the interaction (pointers)
  SiconosVector& lambda = *inter.lambda(level);
  prod(*_jacglambda, lambda, *inter.data(r), false);
}

void FirstOrderLinearR::display() const
{
  cout << " ===== Linear Time Invariant relation display ===== " << endl;
  cout << "| C " << endl;
  if (_jachx) _jachx->display();
  else cout << "->NULL" << endl;
  cout << "| D " << endl;
  if (_jachlambda) _jachlambda->display();
  else cout << "->NULL" << endl;
  cout << "| F " << endl;
  if (_F) _F->display();
  else cout << "->NULL" << endl;
  cout << "| e " << endl;
  if (_e) _e->display();
  else cout << "->NULL" << endl;
  cout << "| B " << endl;
  if (_jacglambda) _jacglambda->display();
  else cout << "->NULL" << endl;
  cout << " ================================================== " << endl;
}

void FirstOrderLinearR::setComputeEFunction(const std::string& pluginPath, const std::string& functionName)
{
  FirstOrderR::setComputeEFunction(pluginPath, functionName);
}
void FirstOrderLinearR::setComputeEFunction(FOVecPtr ptrFunct)
{
  _plugine->setComputeFunction((void*)ptrFunct);
}


void FirstOrderLinearR::saveRelationToXML() const
{
  if (!_relationxml)
    RuntimeException::selfThrow("FirstOrderLinearR::saveRelationToXML, no yet implemented.");

  SP::LinearRXML folrXML = (boost::static_pointer_cast<LinearRXML>(_relationxml));
  folrXML->setC(*_jachx);
  folrXML->setD(*_jachlambda);
  folrXML->setF(*_F);
  folrXML->setE(*_e);
  folrXML->setB(*_jacglambda);
}

FirstOrderLinearR* FirstOrderLinearR::convert(Relation *r)
{
  return dynamic_cast<FirstOrderLinearR*>(r);
}

