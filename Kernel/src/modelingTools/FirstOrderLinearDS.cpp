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
#include "FirstOrderLinearDS.hpp"
#include "FirstOrderLinearDSXML.hpp"
#include "Plugin.hpp"

using namespace std;
typedef void (*computeAfct)(double, unsigned int, unsigned int, double*, unsigned int, double*);

// --- Constructors ---

// From xml file
FirstOrderLinearDS::FirstOrderLinearDS(SP::DynamicalSystemXML dsXML)
  : FirstOrderNonLinearDS(dsXML)
{

  // pointer to xml
  SP::FirstOrderLinearDSXML foldsxml = (boost::static_pointer_cast <FirstOrderLinearDSXML>(dsXML));
  _pluginb.reset(new PluggedObject());
  _pluginA.reset(new PluggedObject());

  // Check if f is given as a plug-in in xml input file.
  if (foldsxml->hasF() || foldsxml->hasJacobianfx())
    RuntimeException::selfThrow("FirstOrderLinearDS - xml constructor, you give a f or its jacobian as a plug-in for a FirstOrderLinearDS -> set rather A and b plug-in.");

  string plugin;
  // A
  if (foldsxml->hasA())
  {
    if (foldsxml->isAPlugin())
    {
      plugin = foldsxml->getAPlugin();
      _A.reset(new SimpleMatrix(_n, _n));
      setComputeAFunction(SSL::getPluginName(plugin), SSL::getPluginFunctionName(plugin));
    }
    else
      _A.reset(new SimpleMatrix(foldsxml->getA()));
  }

  // b
  if (foldsxml->hasB())
  {
    if (foldsxml->isBPlugin())
    {
      _b.reset(new SimpleVector(_n));
      plugin = foldsxml->getBPlugin();
      setComputebFunction(SSL::getPluginName(plugin), SSL::getPluginFunctionName(plugin));
    }
    else
      _b.reset(new SimpleVector(foldsxml->getBVector()));
  }

  checkDynamicalSystem();
}

// From a minimum set of data, A and b connected to a plug-in
FirstOrderLinearDS::FirstOrderLinearDS(SP::SiconosVector newX0, const string& APlugin, const string& bPlugin):
  FirstOrderNonLinearDS(newX0)
{

  _pluginb.reset(new PluggedObject());
  _pluginA.reset(new PluggedObject());
  _pluginA->setComputeFunction(APlugin);
  _pluginb->setComputeFunction(bPlugin);

  _f.reset(new SimpleVector(getDim()));
  _A.reset(new SimpleMatrix(getDim(), getDim()));
  _b.reset(new SimpleVector(getDim()));

  checkDynamicalSystem();
}
FirstOrderLinearDS::FirstOrderLinearDS(const SiconosVector& newX0, const std::string& APlugin, const std::string& bPlugin):
  FirstOrderNonLinearDS(createSPtrSiconosVector((SiconosVector&) newX0))
{
  //  SP::SiconosVector sp_aux = createSPtrSiconosVector ((SiconosVector&) newX0);
  //  FirstOrderNonLinearDS::FirstOrderNonLinearDS();

  _pluginb.reset(new PluggedObject());
  _pluginA.reset(new PluggedObject());
  _pluginA->setComputeFunction(APlugin);
  _pluginb->setComputeFunction(bPlugin);

  _f.reset(new SimpleVector(getDim()));
  _A.reset(new SimpleMatrix(getDim(), getDim()));
  _b.reset(new SimpleVector(getDim()));


  checkDynamicalSystem();

}

// From a minimum set of data, A from a given matrix

FirstOrderLinearDS::FirstOrderLinearDS(SP::SiconosVector newX0, SP::SiconosMatrix newA):
  FirstOrderNonLinearDS(newX0)
{
  _f.reset(new SimpleVector(getDim()));
  _pluginb.reset(new PluggedObject());
  _pluginA.reset(new PluggedObject());
  if ((newA->size(0) != _n) || (newA->size(1) != _n))
    RuntimeException::selfThrow("FirstOrderLinearDS - constructor(number,x0,A): inconsistent dimensions with problem size for input matrix A");

  _A = newA;
  checkDynamicalSystem();
}

// From a minimum set of data, A from a given matrix
FirstOrderLinearDS::FirstOrderLinearDS(SP::SiconosVector newX0, SP::SiconosMatrix newA, SP::SiconosVector newB):
  FirstOrderNonLinearDS(newX0)
{
  _pluginb.reset(new PluggedObject());
  _pluginA.reset(new PluggedObject());
  if ((newA->size(0) != _n) || (newA->size(1) != _n))
    RuntimeException::selfThrow("FirstOrderLinearDS - constructor(x0,A,b): inconsistent dimensions with problem size for input matrix A");

  if (newB->size() != _n)
    RuntimeException::selfThrow("FirstOrderLinearDS - constructor(x0,A,b): inconsistent dimensions with problem size for input vector b ");

  _A = newA;
  _b = newB;
  _f.reset(new SimpleVector(getDim()));

  checkDynamicalSystem();
}

// Copy constructor
FirstOrderLinearDS::FirstOrderLinearDS(const FirstOrderLinearDS & FOLDS): FirstOrderNonLinearDS(FOLDS)
{
  _A.reset(new SimpleMatrix(*(FOLDS.A())));

  if (_b)
    _b.reset(new SimpleVector(*(FOLDS.b())));

  if (Type::value(FOLDS) == Type::FirstOrderLinearDS)
  {
    _pluginA.reset(new PluggedObject(*(FOLDS.getPluginA())));
    _pluginb.reset(new PluggedObject(*(FOLDS.getPluginB())));
  }
}

bool FirstOrderLinearDS::checkDynamicalSystem() // useless ...?
{
  bool output = DynamicalSystem::checkDynamicalSystem();
  if (!output) cout << "FirstOrderLinearDS Warning: your dynamical system seems to be uncomplete (check = false)" << endl;
  return output;
}

void FirstOrderLinearDS::initRhs(double time)
{
  computeRhs(time); // If necessary, this will also compute A and b.
  if (! _jacxRhs)  // if not allocated with a set or anything else
  {
    if (_A && ! _M)  // if M is not defined, then A = jacobianRhsx, no memory allocation for that one.
      _jacxRhs = _A;
    else if (_A && _M)
      _jacxRhs.reset(new SimpleMatrix(_n, _n));
    // else no allocation, jacobian is equal to 0.
  }
  computeJacobianRhsx(time);
}

void FirstOrderLinearDS::updatePlugins(double time)
{
  computeA(time);
  if (_b)
    computeb(time);
}
/*
void FirstOrderLinearDS::setA(const Plugged_Matrix_FTime& newValue)
{
  assert(newValue.size(0)==n&&"FirstOrderLinearDS - setA: inconsistent dimensions with problem size for input matrix A.");
  assert(newValue.size(1)==n&&"FirstOrderLinearDS - setA: inconsistent dimensions with problem size for input matrix A.");

  if( ! A )
    A.reset(new Plugged_Matrix_FTime(newValue));
  else
    *A = newValue;
}

void FirstOrderLinearDS::setB(const Plugged_Vector_FTime& newValue)
{
  assert(newValue.size()==n&&"FirstOrderLinearDS - setB: inconsistent dimensions with problem size for input vector b");

  if( ! b )
    b.reset(new Plugged_Vector_FTime(newValue));
  else
    *b = newValue;
}
*/
void FirstOrderLinearDS::setComputeAFunction(const string& pluginPath, const string& functionName)
{
  _pluginA->setComputeFunction(pluginPath, functionName);
  //   Plugin::setFunction(&_APtr, pluginPath, functionName);
  //   SSL::buildPluginName(pluginNameAPtr,pluginPath,functionName);
}

void FirstOrderLinearDS::setComputeAFunction(LDSPtrFunction fct)
{
  _pluginA->setComputeFunction((void*)fct);
  //  _APtr=fct;
}
void FirstOrderLinearDS::setComputebFunction(const string& pluginPath, const string& functionName)
{
  //  Plugin::setFunction(&_bPtr, pluginPath, functionName);
  _pluginb->setComputeFunction(pluginPath, functionName);
  if (!_b)
    _b.reset(new SimpleVector(getDim()));
  //  SSL::buildPluginName(pluginNamebPtr,pluginPath,functionName);
}

void FirstOrderLinearDS::setComputebFunction(LDSPtrFunction fct)
{
  _pluginb->setComputeFunction((void*)fct);
  //  _bPtr = fct;
}


void FirstOrderLinearDS::computeA(const double time)
{
  if (_A && _pluginA->fPtr)
  {
    ((computeAfct)_pluginA->fPtr)(time, _n, _n, &(*_A)(0, 0), _z->size(), &(*_z)(0));
  }
}

void FirstOrderLinearDS::computeb(const double time)
{
  if (_b && _pluginb->fPtr)
    ((LDSPtrFunction)_pluginb->fPtr)(time, _n, &(*_b)(0), _z->size(), &(*_z)(0));
}
/*This function is called only by Lsodar and eventDriven*/
void FirstOrderLinearDS::computeRhs(const double time, const bool)
{
  // second argument is useless at the time - Used in derived classes
  // compute A=jacobianfx

  *_x[1] = * _r;

  if (_A)
  {
    computeA(time);
    prod(*_A, *_x[0], *_x[1], false);
  }

  // compute and add b if required
  if (_b)
  {
    computeb(time);
    *_x[1] += *_b;
  }

  if (_M)
  {
    // allocate invM at the first call of the present function
    if (! _invM)
      _invM.reset(new SimpleMatrix(*_M));

    _invM->PLUForwardBackwardInPlace(*_x[1]);
  }
}

void FirstOrderLinearDS::computeJacobianRhsx(const double time, const bool)
{
  computeA(time);

  if (_M && _A)
  {
    *_jacxRhs = *_A;
    // copy M into invM for LU-factorisation, at the first call of this function.
    if (! _invM)
      _invM.reset(new SimpleMatrix(*_M));
    // solve MjacobianRhsx = A
    _invM->PLUForwardBackwardInPlace(*_jacxRhs);
  }
  // else jacobianRhsx = A, pointers equality.

}

void FirstOrderLinearDS::display() const
{
  cout << "=== Linear system display, " << _number << endl;
  cout << "=============================" << endl;
}

void FirstOrderLinearDS::saveSpecificDataToXML()
{
  if (!_dsxml)
    RuntimeException::selfThrow("FirstOrderLinearDS::saveDSToXML - The DynamicalSystemXML object doesn't exists");
  boost::static_pointer_cast<FirstOrderLinearDSXML>(_dsxml)->setA(*_A);

  // b
  if (_b)
  {
    if (!(boost::static_pointer_cast <FirstOrderLinearDSXML>(_dsxml))->isBPlugin())
      boost::static_pointer_cast<FirstOrderLinearDSXML>(_dsxml)->setB(*_b);
  }

  else RuntimeException::selfThrow("FirstOrderLinearDS::saveDSToXML - The DynamicalSystemXML object doesn't exists");
}

FirstOrderLinearDS* FirstOrderLinearDS::convert(DynamicalSystem* ds)
{
  FirstOrderLinearDS* lsds = dynamic_cast<FirstOrderLinearDS*>(ds);
  return lsds;
}
void FirstOrderLinearDS::computef(double time)
{
  updatePlugins(time);
  prod(*_A, *_x[0], *_f);
  if (_b)
  {
    computeb(time);
    *_f += *_b;
  }
}

void FirstOrderLinearDS::computef(double time, SP::SiconosVector x2)
{
  //  RuntimeException::selfThrow("FirstOrderLinearDS::computeF - Must not be used");
  updatePlugins(time);
  prod(*_A, *x2, *_f);
  if (_b)
  {
    *_f += *_b;
  }
}

void FirstOrderLinearDS::setA(const SiconosMatrix& newA)
{
  if (_A)
    *_A = newA;
  else
    _A.reset(new SimpleMatrix(newA));
}

void FirstOrderLinearDS::zeroPlugin()
{
  if (_pluginM)
    _pluginM.reset(new PluggedObject());
  if (_pluginA)
    _pluginA.reset(new PluggedObject());
  if (_pluginb)
    _pluginb.reset(new PluggedObject());
}
