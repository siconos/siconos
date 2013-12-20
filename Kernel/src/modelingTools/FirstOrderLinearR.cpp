/* Siconos-Kernel, Copyright INRIA 2005-2012.
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

#include "BlockVector.hpp"

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
FirstOrderLinearR::FirstOrderLinearR(const std::string& Cname, const std::string& Bname):
  FirstOrderR(LinearR)
{
  // Warning: we cannot allocate memory for C/D matrix since no interaction
  // is connected to the relation. This will be done during initialize.
  // We only set the name of the plugin-function and connect it to the user-defined function.
  _pluginJachx->setComputeFunction(Cname);
  _pluginJacLg->setComputeFunction(Bname);
}

// Constructor from a complete set of data (plugin)
FirstOrderLinearR::FirstOrderLinearR(const std::string& Cname, const std::string& Dname, const std::string& Fname, const std::string& Ename, const std::string& Bname): FirstOrderR(LinearR)
{
  _pluginJachx->setComputeFunction(Cname);
  _pluginJachlambda->setComputeFunction(Dname);
  _pluginJacLg->setComputeFunction(Bname);
  _pluginf->setComputeFunction(Fname);
  _plugine->setComputeFunction(Ename);
}

// Minimum data (C, B as pointers) constructor
FirstOrderLinearR::FirstOrderLinearR(SP::SiconosMatrix C, SP::SiconosMatrix B):
  FirstOrderR(LinearR)
{
  _jachx = C;
  _jacglambda = B;
}

// // Constructor from a complete set of data
FirstOrderLinearR::FirstOrderLinearR(SP::SiconosMatrix C, SP::SiconosMatrix D, SP::SiconosMatrix F, SP::SiconosVector E, SP::SiconosMatrix B):
  FirstOrderR(LinearR)
{
  _jachx = C;
  _jacglambda = B;
  _jachlambda = D;
  _F = F;
  _e = E;
}

void FirstOrderLinearR::initialize(Interaction& inter)
{
  // Note: do not call FirstOrderR::initialize to avoid jacobianH and jacobianG allocation.

  // get interesting size
  unsigned int sizeY = inter.getSizeOfY();
  unsigned int sizeX = inter.getSizeOfDS();
  unsigned int sizeZ = inter.getSizez();
  // Update data member (links to DS variables)
  if (!_jachx)
    _jachx.reset(new SimpleMatrix(sizeY, sizeX));
  if (!_jacglambda)
    _jacglambda.reset(new SimpleMatrix(sizeX, sizeY));

  // Check if various operators sizes are consistent.
  // Reference: interaction.

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




void FirstOrderLinearR::computeC(double time, Interaction& inter)
{
  if (_pluginJachx->fPtr)
  {
    SiconosVector workZ(*inter.data(z));
    ((FOMatPtr1)(_pluginJachx->fPtr))(time, _jachx->size(0), _jachx->size(1), &(*_jachx)(0, 0), workZ.size(), &(workZ)(0));
    // Copy data that might have been changed in the plug-in call.
    *inter.data(z) = workZ;
  }
}

void FirstOrderLinearR::computeD(double time, Interaction& inter)
{
  if (_pluginJachlambda->fPtr)
  {
    SiconosVector workZ = *inter.data(z);
    ((FOMatPtr1)(_pluginJachlambda->fPtr))(time, _jachlambda->size(0), _jachlambda->size(1), &(*_jachlambda)(0, 0), workZ.size(), &(workZ)(0));
    // Copy data that might have been changed in the plug-in call.
    *inter.data(z) = workZ;
  }
}

void FirstOrderLinearR::computeF(double time, Interaction& inter)
{
  if (_pluginf->fPtr)
  {
    SiconosVector workZ = *inter.data(z);
    ((FOMatPtr1)(_pluginf->fPtr))(time, _F->size(0), _F->size(1), &(*_F)(0, 0), workZ.size(), &(workZ)(0));
    // Copy data that might have been changed in the plug-in call.
    *inter.data(z) = workZ;
  }
}

void FirstOrderLinearR::computeE(double time, Interaction& inter)
{

  if (_plugine->fPtr)
  {
    SiconosVector workZ = *inter.data(z);
    ((FOVecPtr) _plugine->fPtr)(time, _e->size(), &(*_e)(0), workZ.size(), &(workZ)(0));
    // Copy data that might have been changed in the plug-in call.
    *inter.data(z) = workZ;
  }
}

void FirstOrderLinearR::computeb(double time, Interaction& inter)
{
  if (_pluginJacLg->fPtr)
  {
    SiconosVector workZ = *inter.data(z);
    ((FOMatPtr1) _pluginJacLg->fPtr)(time, _jacglambda->size(0), _jacglambda->size(1), &(*_jacglambda)(0, 0), workZ.size(), &(workZ)(0));
    // Copy data that might have been changed in the plug-in call.
    *inter.data(z) = workZ;
  }
}

void FirstOrderLinearR::computeh(double time, Interaction& inter)
{
  computeOutput(time, inter, 0);
}

void FirstOrderLinearR::computeg(double time, Interaction& inter)
{
  computeInput(time, inter, 0);
}
void FirstOrderLinearR::computeOutput(double time, Interaction& inter, unsigned int level)
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
  std::cout << " ===== Linear Time Invariant relation display ===== " <<std::endl;
  std::cout << "| C " <<std::endl;
  if (_jachx) _jachx->display();
  else std::cout << "->NULL" <<std::endl;
  std::cout << "| D " <<std::endl;
  if (_jachlambda) _jachlambda->display();
  else std::cout << "->NULL" <<std::endl;
  std::cout << "| F " <<std::endl;
  if (_F) _F->display();
  else std::cout << "->NULL" <<std::endl;
  std::cout << "| e " <<std::endl;
  if (_e) _e->display();
  else std::cout << "->NULL" <<std::endl;
  std::cout << "| B " <<std::endl;
  if (_jacglambda) _jacglambda->display();
  else std::cout << "->NULL" <<std::endl;
  std::cout << " ================================================== " <<std::endl;
}

void FirstOrderLinearR::saveRelationToXML() const
{
  if (!_relationxml)
    RuntimeException::selfThrow("FirstOrderLinearR::saveRelationToXML, no yet implemented.");

  SP::LinearRXML folrXML = (std11::static_pointer_cast<LinearRXML>(_relationxml));
  folrXML->setC(*_jachx);
  folrXML->setD(*_jachlambda);
  folrXML->setF(*_F);
  folrXML->setE(*_e);
  folrXML->setB(*_jacglambda);
}

