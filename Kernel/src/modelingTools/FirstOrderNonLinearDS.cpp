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
#include "FirstOrderNonLinearDS.h"
#include "FirstOrderNonLinearDSXML.h"
#include "BlockVector.h"
#include "Plugin.hpp"
#include "PluginTypes.hpp"

using namespace std;

// ===== CONSTRUCTORS =====

// From a minimum set of data
FirstOrderNonLinearDS::FirstOrderNonLinearDS(SP::SiconosVector newX0): DynamicalSystem(DS::FONLDS, newX0->size())
{
  // == Initial conditions ==
  _x0 = newX0;

  // == Current state ==
  // x is composed of two blocks of size n, x[0] = \f$ x \f$ and x[1]=\f$ \dot x \f$.
  // x[0] initialized with x0.

  _x[0].reset(new SimpleVector(*_x0));
  _x[1].reset(new SimpleVector(_n));

  //mG
  mG_alpha.reset(new SimpleVector(_n));
  mResidur.reset(new SimpleVector(_n));
  mXp.reset(new SimpleVector(getDim()));
  mXq.reset(new SimpleVector(getDim()));
  _workFree.reset(new SimpleVector(getDim()));
  mfold.reset(new SimpleVector(getDim()));

  // == r ==

  _r.reset(new SimpleVector(_n));

  checkDynamicalSystem();
}
FirstOrderNonLinearDS::FirstOrderNonLinearDS(const SiconosVector& newX0): DynamicalSystem(DS::FONLDS, newX0.size())
{
  // == Initial conditions ==
  _x0 = createSPtrSiconosVector((SiconosVector&)newX0);

  // == Current state ==
  // x is composed of two blocks of size n, x[0] = \f$ x \f$ and x[1]=\f$ \dot x \f$.
  // x[0] initialized with x0.

  _x[0].reset(new SimpleVector(*_x0));
  _x[1].reset(new SimpleVector(_n));

  //mG
  mG_alpha.reset(new SimpleVector(_n));
  mResidur.reset(new SimpleVector(_n));
  mXp.reset(new SimpleVector(getDim()));
  mXq.reset(new SimpleVector(getDim()));
  _workFree.reset(new SimpleVector(getDim()));
  mfold.reset(new SimpleVector(getDim()));

  // == r ==

  _r.reset(new SimpleVector(_n));

  checkDynamicalSystem();

}

// From XML file
FirstOrderNonLinearDS::FirstOrderNonLinearDS(SP::DynamicalSystemXML dsXML):
  DynamicalSystem(dsXML)
{
  /*  // -- FONLDS xml object --
  SP::FirstOrderNonLinearDSXML fonlds = boost::static_pointer_cast <FirstOrderNonLinearDSXML>(dsxml);

  // === Initial conditions ===
  // Warning: n is set thanks to x0 size
  if ( ! fonlds->hasX0())
    RuntimeException::selfThrow("FirstOrderNonLinearDS:: xml constructor, x0 is a required input");

  x0.reset(new SimpleVector(fonlds->getX0()));

  n = x0->size();

  // === Current state (optional input) ===
  // x is composed of two blocks of size n, (*x)[0] = \f$ x \f$ and (*x)[1]=\f$ \dot x \f$.

  if ( fonlds->hasX())
    x[0].reset(new SimpleVector(fonlds->getX()));
  else // (*x)[0] initialize with x0.
    x[0].reset(new SimpleVector(*x0));
  // build and initialize right-hand side
  x[1].reset(new SimpleVector(n));
  // r

  r.reset(new SimpleVector(n));

  string plugin;

  // f and jacobianXF are required for DynamicalSystem but not for derived class.
  // Then we can not set exception if they are not given.
  if(fonlds->hasM())
    {
      if ( fonlds->isMPlugin())
  {
    plugin = fonlds->getMPlugin();
    setComputeMFunction(SSL::getPluginName( plugin ), SSL::getPluginFunctionName( plugin ));
  }
      else // This means that M is constant
  {
          M.reset(new PMJF(fonlds->getMMatrix()));
    if(M->size(0)!=n || M->size(1)!=n)
      RuntimeException::selfThrow("FirstOrderNonLinearDS:: xml constructor, M size differs from n!");
  }
    }

  if(fonlds->hasF())
    {
      if(fonlds->isFPlugin())
  {
    plugin = fonlds->getFPlugin();
    setComputeFFunction(SSL::getPluginName( plugin ), SSL::getPluginFunctionName( plugin ));
  }
      else
  {
    if(fonlds->getFVector().size()!=n)
      RuntimeException::selfThrow("FirstOrderNonLinearDS:: xml constructor, f size differs from n!");

          mf.reset(new PVF(fonlds->getFVector()));
  }
    }

  if(fonlds->hasJacobianXF())
    {
      if ( fonlds->isJacobianXFPlugin())
  {
    plugin = fonlds->getJacobianXFPlugin();
    setComputeJacobianXFFunction(SSL::getPluginName( plugin ), SSL::getPluginFunctionName( plugin ));
  }
      else // This means that jacobianXF is constant
  {
          jacobianXF.reset(new PMJF(fonlds->getJacobianXFMatrix()));
    if(jacobianXF->size(0)!=n || jacobianXF->size(1)!=n)
      RuntimeException::selfThrow("FirstOrderNonLinearDS:: xml constructor, jacobianXF size differs from n!");
  }
    }

  // Memory
  if ( fonlds->hasXMemory())
    xMemory.reset(new SiconosMemory( fonlds->getXMemoryXML() ));

    checkDynamicalSystem();*/
}

// From a minimum set of data
FirstOrderNonLinearDS::FirstOrderNonLinearDS(const SiconosVector& newX0, const string& fPlugin, const string& jacobianXFPlugin):
  DynamicalSystem(DS::FONLDS, newX0.size())
{
  // == Initial conditions ==
  _x0.reset(new SimpleVector(newX0));

  // == Current state ==
  // x is composed of two blocks of size n, x[0] = \f$ x \f$ and x[1]=\f$ \dot x \f$.
  // x[0] initialized with x0.

  _x[0].reset(new SimpleVector(*_x0));
  _x[1].reset(new SimpleVector(_n));
  //mG
  mG_alpha.reset(new SimpleVector(_n));
  mResidur.reset(new SimpleVector(_n));
  mXp.reset(new SimpleVector(getDim()));
  mXq.reset(new SimpleVector(getDim()));
  _workFree.reset(new SimpleVector(getDim()));
  _r.reset(new SimpleVector(getDim()));
  mfold.reset(new SimpleVector(getDim()));

  // == r ==

  _r.reset(new SimpleVector(_n));

  // == f and its jacobian ==
  // Allocation and link with the plug-in

  Plugin::setFunction(&_computeFPtr, SSL::getPluginName(fPlugin), SSL::getPluginFunctionName(fPlugin));
  //  _pluginJacXF->setComputeFunction
  Plugin::setFunction(&computeJacobianXFPtr, SSL::getPluginName(jacobianXFPlugin), SSL::getPluginFunctionName(jacobianXFPlugin));
  _pluginNameComputeFPtr = fPlugin;
  pluginNameComputeJacobianXFPtr = jacobianXFPlugin;

  checkDynamicalSystem();
}
void FirstOrderNonLinearDS::preparStep()
{
  mXp->zero();
  _r->zero();
};
bool FirstOrderNonLinearDS::checkDynamicalSystem()
{
  DynamicalSystem::checkDynamicalSystem();
  bool output = DynamicalSystem::checkDynamicalSystem();
  if (!output) cout << "FirstOrderNonLinearDS Warning: your dynamical system seems to be uncomplete (check = false)" << endl;
  return output;
}


/*void FirstOrderNonLinearDS::setM(const PMJF& newValue)
{
  assert(newValue.size(0)==n&&"FirstOrderNonLinearDS - setM: inconsistent dimensions with problem size for input matrix M.");
  assert(newValue.size(1)==n&&"FirstOrderNonLinearDS - setM: inconsistent dimensions with problem size for input matrix M.");

  if( ! M )
    M.reset(new PMJF(newValue));
  else
    *M = newValue;
    }*/

void FirstOrderNonLinearDS::setInvM(const SiconosMatrix& newValue)
{
  if (newValue.size(0) != _n || newValue.size(1) != _n)
    RuntimeException::selfThrow("FirstOrderNonLinearDS::setInvM: inconsistent dimensions with problem size for input matrix.");

  if (! _invM)
    _invM.reset(new SimpleMatrix(_n, _n));
  *_invM = newValue;
}

void FirstOrderNonLinearDS::setInvMPtr(SP::SiconosMatrix newPtr)
{
  _invM = newPtr;
}
/*
void FirstOrderNonLinearDS::setF(const PVF& newValue)
{
  assert(newValue.size()==n&&"FirstOrderNonLinearDS - setF: inconsistent dimensions with problem size for input vector f");

  if( ! mf )
    mf.reset(new PVF(newValue));
  else
    *mf = newValue;
    }*/
/*
void FirstOrderNonLinearDS::setJacobianXF(const PMJF& newValue)
{
 assert(newValue.size(0)==n&&"FirstOrderNonLinearDS - setJacobianXF: inconsistent dimensions with problem size for input matrix M.");
 assert(newValue.size(1)==n&&"FirstOrderNonLinearDS - setJacobianXF: inconsistent dimensions with problem size for input matrix M.");

 if( ! jacobianXF )
   jacobianXF.reset(new PMJF(newValue));
 else
   *jacobianXF = newValue;
}
*/
void FirstOrderNonLinearDS::initRhs(double time)
{
  // compute initial values for f and jacobianXF, initialize right-hand side.
  computeRhs(time); // this will compute, if required, f and M.

  if (! _jacXRhs)  // if not allocated with a set or anything else
  {
    if (_jacobianXF && ! _M)  // if M is not defined, then jacobianXF = jacobianXRhs, no memory allocation for that one.
      _jacXRhs = _jacobianXF;
    else if (_jacobianXF && _M)
      _jacXRhs.reset(new SimpleMatrix(_n, _n));

    // else no allocation, jacobian is equal to 0.
  }
  computeJacobianXRhs(time);
}

void FirstOrderNonLinearDS::updatePlugins(double time)
{
  computeM(time);
  computeF(time);
  computeJacobianXF(time);
}

void FirstOrderNonLinearDS::initialize(const string& simulationType, double time, unsigned int sizeOfMemory)
{
  // reset x to x0 and r to zero.
  _r->zero();
  *(_x[0]) = *_x0;

  // If z has not been set, we initialize it with a null vector of size 1, since z is required in plug-in functions call.
  if (! _z)
    _z.reset(new SimpleVector(1));

  // Initialize memory vectors
  initMemory(sizeOfMemory);

  updatePlugins(time);
  *mfold = *mf;

  if (simulationType == "EventDriven")
  {
    // Rhs and its jacobian
    initRhs(time);
  }
}

// ===== MEMORY MANAGEMENT FUNCTIONS =====

void FirstOrderNonLinearDS::initMemory(unsigned int steps)
{
  DynamicalSystem::initMemory(steps);

  if (steps == 0)
    cout << "Warning : FirstOrderNonLinearDS::initMemory with size equal to zero" << endl;
  else
    _rMemory.reset(new SiconosMemory(steps));
}

void FirstOrderNonLinearDS::swapInMemory()
{
  _xMemory->swap(_x[0]);
  _rMemory->swap(_r);
  *mfold = *mf;
}

// ===== COMPUTE PLUGINS FUNCTIONS =====

void FirstOrderNonLinearDS::setComputeMFunction(const string& pluginPath, const string& functionName)
{
  Plugin::setFunction(&pluginComputeM, pluginPath, functionName);
  SSL::buildPluginName(pluginNamePluginComputeM, pluginPath, functionName);
}

void FirstOrderNonLinearDS::setComputeMFunction(FPtr1 fct)
{
  pluginComputeM = fct;
}
void FirstOrderNonLinearDS::setComputeFFunction(const string& pluginPath, const string& functionName)
{
  Plugin::setFunction(&_computeFPtr, pluginPath, functionName);
  SSL::buildPluginName(_pluginNameComputeFPtr, pluginPath, functionName);
}

void FirstOrderNonLinearDS::setComputeFFunction(FPtr1 fct)
{
  _computeFPtr = fct;
}

void FirstOrderNonLinearDS::setComputeJacobianXFFunction(const string& pluginPath, const string& functionName)
{
  Plugin::setFunction(&computeJacobianXFPtr, pluginPath, functionName);
  SSL::buildPluginName(pluginNameComputeJacobianXFPtr, pluginPath, functionName);
}

void FirstOrderNonLinearDS::setComputeJacobianXFFunction(FPtr1 fct)
{
  computeJacobianXFPtr = fct;
}

void FirstOrderNonLinearDS::computeM(double time)
{
  // second argument is useless at the time - Used in derived classes
  if (pluginComputeM)
  {
    (pluginComputeM)(time, _n, &((*(_x[0]))(0)), &(*_M)(0, 0), _z->size(), &(*_z)(0));
  }
}

void FirstOrderNonLinearDS::computeM(double time, SP::SiconosVector x2)
{
  // second argument is useless at the time - Used in derived classes
  if (pluginComputeM)
  {
    assert(x2->size() == _n && "FirstOrderNonLinearDS::computeM(t,x) x size does not fit with the system size.");
    (pluginComputeM)(time, _n, &((*x2)(0)), &(*_M)(0, 0), _z->size(), &(*_z)(0));
  }
  // else nothing!
}

void FirstOrderNonLinearDS::computeF(double time)
{
  if (_computeFPtr)
    (_computeFPtr)(time, _n, &((*(_x[0]))(0)) , &(*mf)(0), _z->size(), &(*_z)(0));
}

void FirstOrderNonLinearDS::computeF(double time, SP::SiconosVector x2)
{
  if (_computeFPtr)
    (_computeFPtr)(time, _n, &((*x2)(0)) , &(*mf)(0), _z->size(), &(*_z)(0));
  // else nothing!
}

void FirstOrderNonLinearDS::computeJacobianXF(double time, bool)
{
  // second argument is useless at the time - Used in derived classes
  if (computeJacobianXFPtr)
    (computeJacobianXFPtr)(time, _n, &((*(_x[0]))(0)), &(*_jacobianXF)(0, 0), _z->size(), &(*_z)(0));
}

void FirstOrderNonLinearDS::computeJacobianXF(double time, SP::SiconosVector x2)
{
  // second argument is useless at the time - Used in derived classes
  if (computeJacobianXFPtr)
    (computeJacobianXFPtr)(time, _n, &((*x2)(0)), &(*_jacobianXF)(0, 0), _z->size(), &(*_z)(0));
}

void FirstOrderNonLinearDS::computeRhs(double time, bool)
{
  // second argument is useless at the time - Used in derived classes

  // compute rhs = M-1*( f + r ).

  *_x[1] = *_r; // Warning: p update is done in Interactions/Relations

  if (mf)
  {
    computeF(time);
    *(_x[1]) += *mf;
  }

  if (_M)
  {
    // allocate invM at the first call of the present function
    if (! _invM)
      _invM.reset(new SimpleMatrix(*_M));
    _invM->PLUForwardBackwardInPlace(*_x[1]);
  }
}

void FirstOrderNonLinearDS::computeJacobianXRhs(double time, bool)
{
  // second argument is useless at the time - Used in derived classes

  // compute jacobian of rhs according to x, = M-1(jacobianXF + jacobianX(T.u))
  // At the time, second term is set to zero.
  computeJacobianXF(time);
  // solve M*jacobianXRhS = jacobianXF
  if (_M && _jacobianXF)
  {
    *_jacXRhs = *_jacobianXF;
    // copy _M into _invM for LU-factorisation, at the first call of this function.
    if (! _invM)
      _invM.reset(new SimpleMatrix(*_M));

    _invM->PLUForwardBackwardInPlace(*_jacXRhs);
  }
  // else jacobianXRhs = jacobianXF, pointers equality set in initRhs

}

// ===== XML MANAGEMENT FUNCTIONS =====

void FirstOrderNonLinearDS::saveSpecificDataToXML()
{
  // -- FirstOrderNonLinearDS  xml object --
  SP::FirstOrderNonLinearDSXML fonlds = boost::static_pointer_cast <FirstOrderNonLinearDSXML>(_dsxml);
  // --- other data ---
  if (!fonlds)
    RuntimeException::selfThrow("FirstOrderNonLinearDS::saveSpecificDataToXML - The DynamicalSystemXML object doesn't exists");

  if (_computeFPtr)
    fonlds->setFPlugin(_pluginNameComputeFPtr);
  if (computeJacobianXFPtr)
    fonlds->setJacobianXFPlugin(pluginNamePluginComputeM);
}

// ===== MISCELLANEOUS ====

void FirstOrderNonLinearDS::display() const
{
  cout << " =====> First Order Non Linear DS (number: " << _number << ")." << endl;
  cout << "- n (size) : " << _n << endl;
  cout << "- x " << endl;
  if (_x[0]) _x[0]->display();
  else cout << "-> NULL" << endl;
  cout << "- x0 " << endl;
  if (_x0) _x0->display();
  else cout << "-> NULL" << endl;
  cout << "- M: " << endl;
  if (_M) _M->display();
  else cout << "-> NULL" << endl;
  cout << " ============================================" << endl;
}

void FirstOrderNonLinearDS::resetNonSmoothPart()
{
  _r->zero();
}

/*must be remove, replace by the RelativeConvergenceCriteron of the simulation*/
/*double FirstOrderNonLinearDS::dsConvergenceIndicator()
{
    double dsCvgIndic;
  // Velocity is used to calculate the indicator.
  SP::SiconosVector diff(new SimpleVector(x[0]->size()));
  // Compute difference between present and previous Newton steps
  SP::SiconosVector valRef = workV[NewtonSave];
  *diff =  *(x[0]) - *valRef;
  if (valRef->norm2()!=0)
    dsCvgIndic= diff->norm2()/(valRef->norm2());
  else
    dsCvgIndic= diff->norm2();
    return (dsCvgIndic);
    }*/

FirstOrderNonLinearDS* FirstOrderNonLinearDS::convert(DynamicalSystem* ds)
{
  FirstOrderNonLinearDS* fonlds = dynamic_cast<FirstOrderNonLinearDS*>(ds);
  return fonlds;
}

