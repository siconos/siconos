/* Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
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
  _plunginJacxh->setComputeFunction(CName);
  _pluginJacLg->setComputeFunction(BName);
  //  Plugin::setFunction(_plunginJacxh,CName);
  //Plugin::setFunction(_pluginJacLg,BName);
  //  JacH[0].reset(new PluggedMatrix(CName));
  //  JacG[0].reset(new PluggedMatrix(BName));
}

// Constructor from a complete set of data (plugin)
FirstOrderLinearR::FirstOrderLinearR(const string& CName, const string& DName, const string& FName, const string& EName, const string& BName): FirstOrderR(LinearR)
{
  _plunginJacxh->setComputeFunction(CName);
  _pluginJacLh->setComputeFunction(DName);
  _pluginJacLg->setComputeFunction(BName);
  _pluginf->setComputeFunction(FName);
  _plugine->setComputeFunction(EName);

  //   Plugin::setFunction(_plunginJacxh,CName);
  //   Plugin::setFunction(_pluginJacLh,DName);
  //   Plugin::setFunction(_pluginJacLg,BName);

  //   Plugin::setFunction(_pluginf,FName);
  //   Plugin::setFunction(_plugine,EName);



}

// Minimum data (C, B as pointers) constructor
FirstOrderLinearR::FirstOrderLinearR(SP::SiconosMatrix newC, SP::SiconosMatrix newB):
  FirstOrderR(LinearR)
{
  JacXH = newC;
  JacLG = newB;
}

// // Constructor from a complete set of data
FirstOrderLinearR::FirstOrderLinearR(SP::SiconosMatrix  newC, SP::SiconosMatrix newD, SP::SiconosMatrix newF, SP::SiconosVector newE, SP::SiconosMatrix newB):
  FirstOrderR(LinearR)
{
  JacXH = newC;
  JacLG = newB;
  JacLH = newD;
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

void FirstOrderLinearR::initialize(SP::Interaction inter)
{
  assert(inter && "FirstOrderLinearR::initialize failed. No Interaction linked to the present relation.");
  _interaction = inter;

  // Note: do not call FirstOrderR::initialize to avoid jacobianH and jacobianG allocation.

  // Update data member (links to DS variables)
  initDSLinks();
  if (!JacXH)
    RuntimeException::selfThrow("FirstOrderLinearR::initialize() C is null and is a required input.");
  if (!JacLG)
    RuntimeException::selfThrow("FirstOrderLinearR::initialize() B is null and is a required input.");

  // Check if various operators sizes are consistent.
  // Reference: interaction.
  unsigned int sizeY = interaction()->getSizeOfY();
  unsigned int sizeX = interaction()->getSizeOfDS();
  unsigned int sizeZ = interaction()->getSizeZ();

  // The initialization of each matrix/vector depends on the way the Relation was built ie if the matrix/vector
  // was read from xml or not
  if (JacXH->size(0) == 0)
    JacXH->resize(sizeX, sizeY);
  else
    assert((JacXH->size(0) == sizeY && JacXH->size(1) == sizeX) && "FirstOrderLinearR::initialize , inconsistent size between C and Interaction.");


  if (JacLG->size(0) == 0)
    JacLG->resize(sizeY, sizeX);
  else
    assert((JacLG->size(1) == sizeY && JacLG->size(0) == sizeX) && "FirstOrderLinearR::initialize , inconsistent size between B and interaction.");

  // C and B are the minimum inputs. The others may remain null.

  if (JacLH)
  {
    if (JacLH->size(0) == 0)
      JacLH->resize(sizeY, sizeY);
    else
      assert((JacLH->size(0) == sizeY || JacLH->size(1) == sizeY) && "FirstOrderLinearR::initialize , inconsistent size between C and D.");
  }

  if (_F)
  {
    if (_F->size(0) == 0)
      _F->resize(sizeY, sizeZ);
    else
      assert(((_F->size(0) != sizeY) && (_F->size(1) != sizeZ)) && "FirstOrderLinearR::initialize , inconsistent size between C and F.");
  }
  if (!_e && _plugine->fPtr)
  {
    unsigned int sizeY = interaction()->getSizeOfY();
    _e.reset(new SimpleVector(sizeY));
  }

  if (_e)
  {
    if (_e->size() == 0)
      _e->resize(sizeY);
    else
      assert(_e->size() == sizeY && "FirstOrderLinearR::initialize , inconsistent size between C and e.");
  }

  _workZ.reset(new SimpleVector(sizeZ));
}

// // setters




void FirstOrderLinearR::computeC(double time)
{
  if (JacXH)
  {
    if (_plunginJacxh->fPtr)
    {
      unsigned int sizeY = interaction()->getSizeOfY();
      unsigned int sizeX = interaction()->getSizeOfDS();
      unsigned int sizeZ = interaction()->getSizeZ();
      *_workZ = *data[z];
      ((FOMatPtr1)(_plunginJacxh->fPtr))(time, sizeY, sizeX, &(*JacXH)(0, 0), sizeZ, &(*_workZ)(0));
      // Copy data that might have been changed in the plug-in call.
      *data[z] = *_workZ;
    }
  }
}

void FirstOrderLinearR::computeD(double time)
{
  if (JacLH)
  {
    if (_pluginJacLh->fPtr)
    {
      unsigned int sizeY = interaction()->getSizeOfY();
      unsigned int sizeZ = interaction()->getSizeZ();
      *_workZ = *data[z];
      ((FOMatPtr1)(_pluginJacLh->fPtr))(time, sizeY, sizeY, &(*JacLH)(0, 0), sizeZ, &(*_workZ)(0));
      // Copy data that might have been changed in the plug-in call.
      *data[z] = *_workZ;
    }
    // else nothing
  }
}

void FirstOrderLinearR::computeF(double time)
{
  if (_F && _pluginf->fPtr)
  {
    unsigned int sizeY = interaction()->getSizeOfY();
    unsigned int sizeZ = interaction()->getSizeZ();
    *_workZ = *data[z];
    ((FOMatPtr1)(_pluginf->fPtr))(time, sizeY, sizeZ, &(*_F)(0, 0), sizeZ, &(*_workZ)(0));
    // Copy data that might have been changed in the plug-in call.
    *data[z] = *_workZ;
  }
}

void FirstOrderLinearR::computeE(double time)
{

  if (_e && _plugine->fPtr)
  {
    unsigned int sizeY = interaction()->getSizeOfY();
    unsigned int sizeZ = interaction()->getSizeZ();
    *_workZ = *data[z];
    ((FOVecPtr) _plugine->fPtr)(time, sizeY, &(*_e)(0), sizeZ, &(*_workZ)(0));
    // Copy data that might have been changed in the plug-in call.
    *data[z] = *_workZ;
  }
}

void FirstOrderLinearR::computeb(double time)
{
  if (JacLG && _pluginJacLg->fPtr)
  {
    unsigned int sizeY = interaction()->getSizeOfY();
    unsigned int sizeX = interaction()->getSizeOfDS();
    unsigned int sizeZ = interaction()->getSizeZ();
    *_workZ = *data[z];
    ((FOMatPtr1) _pluginJacLg->fPtr)(time, sizeX, sizeY, &(*JacLG)(0, 0), sizeZ, &(*_workZ)(0));
    // Copy data that might have been changed in the plug-in call.
    *data[z] = *_workZ;
  }
}

void FirstOrderLinearR::computeh(double time)
{
  computeOutput(time, 0);
}

void FirstOrderLinearR::computeG(double time)
{
  computeInput(time, 0);
}
void FirstOrderLinearR::computeOutput(double time, unsigned int)
{
  computeC(time);
  computeD(time);
  computeF(time);
  computeE(time);

  // Note that the second argument remains unamed since it is not used: for first order systems, we always compute
  // y[0]

  // We get y and lambda of the interaction (pointers)
  SP::SiconosVector y = interaction()->y(0);
  SP::SiconosVector lambda = interaction()->lambda(0);

  // compute y
  if (JacXH)
    prod(*JacXH, *data[x], *y);
  else
    y->zero();

  if (JacLH)
    prod(*JacLH, *lambda, *y, false);

  if (_e)
    *y += *_e;

  if (_F)
    prod(*_F, *data[z], *y, false);
}

void FirstOrderLinearR::computeInput(double time, unsigned int level)
{
  computeb(time);

  // We get lambda of the interaction (pointers)
  SP::SiconosVector lambda = interaction()->lambda(level);
  prod(*JacLG, *lambda, *data[r], false);
}

void FirstOrderLinearR::display() const
{
  cout << " ===== Linear Time Invariant relation display ===== " << endl;
  cout << "| C " << endl;
  if (JacXH) JacXH->display();
  else cout << "->NULL" << endl;
  cout << "| D " << endl;
  if (JacLH) JacLH->display();
  else cout << "->NULL" << endl;
  cout << "| F " << endl;
  if (_F) _F->display();
  else cout << "->NULL" << endl;
  cout << "| e " << endl;
  if (_e) _e->display();
  else cout << "->NULL" << endl;
  cout << "| B " << endl;
  if (JacLG) JacLG->display();
  else cout << "->NULL" << endl;
  cout << " ================================================== " << endl;
}

// void FirstOrderLinearR::setComputeEFunction(const std::string& pluginPath, const std::string& functionName){
//   FirstOrderR::setComputeEFunction(pluginPath,functionName);

// }
void FirstOrderLinearR::setComputeEFunction(FOVecPtr ptrFunct)
{
  _plugine->setComputeFunction((void*)ptrFunct);
}


void FirstOrderLinearR::saveRelationToXML() const
{
  if (!relationxml)
    RuntimeException::selfThrow("FirstOrderLinearR::saveRelationToXML, no yet implemented.");

  SP::LinearRXML folrXML = (boost::static_pointer_cast<LinearRXML>(relationxml));
  folrXML->setC(*JacXH);
  folrXML->setD(*JacLH);
  folrXML->setF(*_F);
  folrXML->setE(*_e);
  folrXML->setB(*JacLG);
}

FirstOrderLinearR* FirstOrderLinearR::convert(Relation *r)
{
  return dynamic_cast<FirstOrderLinearR*>(r);
}

