/* Siconos-Kernel, Copyright INRIA 2005-2010.
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
#include "FirstOrderNonLinearDS.hpp"
#include "FirstOrderNonLinearDSXML.hpp"
#include "BlockVector.hpp"
#include "Plugin.hpp"
#include "PluginTypes.hpp"

using namespace std;

// ===== CONSTRUCTORS =====

// From a minimum set of data
FirstOrderNonLinearDS::FirstOrderNonLinearDS(SP::SiconosVector newX0):
  DynamicalSystem(newX0->size())
{
  zeroPlugin();
  // == Initial conditions ==
  _x0 = newX0;

  // == Current state ==
  // x is composed of two blocks of size n, x[0] = \f$ x \f$ and x[1]=\f$ \dot x \f$.
  // x[0] initialized with x0.

  _x[0].reset(new SimpleVector(*_x0));
  _x[1].reset(new SimpleVector(_n));

  //mG
  _g_alpha.reset(new SimpleVector(_n));
  _residur.reset(new SimpleVector(_n));
  _xp.reset(new SimpleVector(getDim()));
  _xq.reset(new SimpleVector(getDim()));
  _workFree.reset(new SimpleVector(getDim()));
  _fold.reset(new SimpleVector(getDim()));

  // == r ==

  _r.reset(new SimpleVector(_n));

  checkDynamicalSystem();
}
FirstOrderNonLinearDS::FirstOrderNonLinearDS(const SiconosVector& newX0):
  DynamicalSystem(newX0.size())
{
  zeroPlugin();
  // == Initial conditions ==
  _x0 = createSPtrSiconosVector((SiconosVector&)newX0);

  // == Current state ==
  // x is composed of two blocks of size n, x[0] = \f$ x \f$ and x[1]=\f$ \dot x \f$.
  // x[0] initialized with x0.

  _x[0].reset(new SimpleVector(*_x0));
  _x[1].reset(new SimpleVector(_n));

  //mG
  _g_alpha.reset(new SimpleVector(_n));
  _residur.reset(new SimpleVector(_n));
  _xp.reset(new SimpleVector(getDim()));
  _xq.reset(new SimpleVector(getDim()));
  _workFree.reset(new SimpleVector(getDim()));
  _fold.reset(new SimpleVector(getDim()));

  // == r ==

  _r.reset(new SimpleVector(_n));

  checkDynamicalSystem();

}

// From XML file
FirstOrderNonLinearDS::FirstOrderNonLinearDS(SP::DynamicalSystemXML dsXML):
  DynamicalSystem(dsXML)
{
  zeroPlugin();
  // -- Type::FirstOrderNonLinearDS xml object --
  SP::FirstOrderNonLinearDSXML fonlds = boost::static_pointer_cast <FirstOrderNonLinearDSXML>(dsXML);

  // === Initial conditions ===
  // Warning: n is set thanks to x0 size
  if (! fonlds->hasX0())
    RuntimeException::selfThrow("FirstOrderNonLinearDS:: xml constructor, x0 is a required input");

  _x0.reset(new SimpleVector(fonlds->getX0()));

  _n = _x0->size();

  // === Current state (optional input) ===
  // x is composed of two blocks of size n, (*x)[0] = \f$ x \f$ and (*x)[1]=\f$ \dot x \f$.

  if (fonlds->hasx())
    _x[0].reset(new SimpleVector(fonlds->getx()));
  else // (*x)[0] initialize with x0.
    _x[0].reset(new SimpleVector(*_x0));
  // build and initialize right-hand side
  _x[1].reset(new SimpleVector(_n));
  // r

  _r.reset(new SimpleVector(_n));

  string plugin;

  // f and jacobianfx are required for DynamicalSystem but not for derived class.
  // Then we can not set exception if they are not given.
  if (fonlds->hasM())
  {
    if (fonlds->isMPlugin())
    {
      plugin = fonlds->getMPlugin();
      setComputeMFunction(SSL::getPluginName(plugin), SSL::getPluginFunctionName(plugin));
    }
    else // This means that M is constant
    {
      _M.reset(new SimpleMatrix(fonlds->getMMatrix()));
      if (_M->size(0) != _n || _M->size(1) != _n)
        RuntimeException::selfThrow("FirstOrderNonLinearDS:: xml constructor, M size differs from n!");
    }
  }

  if (fonlds->hasF())
  {
    if (fonlds->isFPlugin())
    {
      plugin = fonlds->getFPlugin();
      setComputeFFunction(SSL::getPluginName(plugin), SSL::getPluginFunctionName(plugin));
    }
    else
    {
      if (fonlds->getFVector().size() != _n)
        RuntimeException::selfThrow("FirstOrderNonLinearDS:: xml constructor, f size differs from n!");

      _f.reset(new SimpleVector(fonlds->getFVector()));
    }
  }

  if (fonlds->hasJacobianfx())
  {
    if (fonlds->isJacobianfxPlugin())
    {
      plugin = fonlds->getJacobianfxPlugin();
      setComputeJacobianfxFunction(SSL::getPluginName(plugin), SSL::getPluginFunctionName(plugin));
    }
    else // This means that jacobianfx is constant
    {
      _jacobianfx.reset(new SimpleMatrix(fonlds->getJacobianfxMatrix()));
      if (_jacobianfx->size(0) != _n || _jacobianfx->size(1) != _n)
        RuntimeException::selfThrow("FirstOrderNonLinearDS:: xml constructor, jacobianfx size differs from n!");
    }
  }

  // Memory
  if (fonlds->hasXMemory())
    _xMemory.reset(new SiconosMemory(fonlds->getXMemoryXML()));

  checkDynamicalSystem();
}

// From a minimum set of data
FirstOrderNonLinearDS::FirstOrderNonLinearDS(const SiconosVector& newX0, const string& fPlugin, const string& jacobianfxPlugin):
  DynamicalSystem(newX0.size())
{
  zeroPlugin();
  // == Initial conditions ==
  _x0.reset(new SimpleVector(newX0));

  // == Current state ==
  // x is composed of two blocks of size n, x[0] = \f$ x \f$ and x[1]=\f$ \dot x \f$.
  // x[0] initialized with x0.

  _x[0].reset(new SimpleVector(*_x0));
  _x[1].reset(new SimpleVector(_n));
  //mG
  _g_alpha.reset(new SimpleVector(_n));
  _residur.reset(new SimpleVector(_n));
  _xp.reset(new SimpleVector(getDim()));
  _xq.reset(new SimpleVector(getDim()));
  _workFree.reset(new SimpleVector(getDim()));
  _r.reset(new SimpleVector(getDim()));
  _fold.reset(new SimpleVector(getDim()));

  // == r ==

  _r.reset(new SimpleVector(_n));

  // == f and its jacobian ==
  // Allocation and link with the plug-in
  _pluginf->setComputeFunction(fPlugin);
  _pluginJacxf->setComputeFunction(jacobianfxPlugin);
  //  Plugin::setFunction(&computeJacobianfxPtr, SSL::getPluginName( jacobianfxPlugin ),SSL::getPluginFunctionName( jacobianfxPlugin ));
  //  _pluginNameComputeFPtr = fPlugin;
  //  pluginNameComputeJacobianfxPtr = jacobianfxPlugin;

  checkDynamicalSystem();
}
void FirstOrderNonLinearDS::zeroPlugin()
{
  // DynamicalSystem::zeroPlugin();
  _pluginf.reset(new PluggedObject());
  _pluginJacxf.reset(new PluggedObject());
  _pluginM.reset(new PluggedObject());
}
void FirstOrderNonLinearDS::preparStep()
{
  _xp->zero();
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

  if( ! _f )
    _f.reset(new PVF(newValue));
  else
    *_f = newValue;
    }*/
/*
void FirstOrderNonLinearDS::setJacobianfx(const PMJF& newValue)
{
 assert(newValue.size(0)==n&&"FirstOrderNonLinearDS - setJacobianfx: inconsistent dimensions with problem size for input matrix M.");
 assert(newValue.size(1)==n&&"FirstOrderNonLinearDS - setJacobianfx: inconsistent dimensions with problem size for input matrix M.");

 if( ! jacobianfx )
   jacobianfx.reset(new PMJF(newValue));
 else
   *jacobianfx = newValue;
}
*/
void FirstOrderNonLinearDS::initRhs(double time)
{
  // compute initial values for f and jacobianfx, initialize right-hand side.
  computeRhs(time); // this will compute, if required, f and M.

  if (! _jacxRhs)  // if not allocated with a set or anything else
  {
    if (_jacobianfx && ! _M)  // if M is not defined, then jacobianfx = jacobianRhsx, no memory allocation for that one.
      _jacxRhs = _jacobianfx;
    else if (_jacobianfx && _M)
      _jacxRhs.reset(new SimpleMatrix(_n, _n));

    // else no allocation, jacobian is equal to 0.
  }
  computeJacobianRhsx(time);
}

void FirstOrderNonLinearDS::updatePlugins(double time)
{
  computeM(time);
  computef(time);
  computeJacobianfx(time);
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
  if (_f)
    *_fold = *_f;

  if (simulationType == "EventDriven")
  {
    // Rhs and its jacobian ==> the right is to put in initOSNA of EventDriven
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
  *_fold = *_f;
}

// ===== COMPUTE PLUGINS FUNCTIONS =====

void FirstOrderNonLinearDS::setComputeMFunction(const string& pluginPath, const string& functionName)
{
  _pluginM->setComputeFunction(pluginPath, functionName);
}

void FirstOrderNonLinearDS::setComputeMFunction(FPtr1 fct)
{
  _pluginM->setComputeFunction((void *)fct);
}
void FirstOrderNonLinearDS::setComputeFFunction(const string& pluginPath, const string& functionName)
{
  _pluginf->setComputeFunction(pluginPath, functionName);
}

void FirstOrderNonLinearDS::setComputeFFunction(FPtr1 fct)
{
  _pluginf->setComputeFunction((void *)fct);
}

void FirstOrderNonLinearDS::setComputeJacobianfxFunction(const string& pluginPath, const string& functionName)
{
  _pluginJacxf->setComputeFunction(pluginPath, functionName);
}

void FirstOrderNonLinearDS::setComputeJacobianfxFunction(FPtr1 fct)
{
  _pluginJacxf->setComputeFunction((void *)fct);
}

void FirstOrderNonLinearDS::computeM(double time)
{
  // second argument is useless at the time - Used in derived classes
  if (_pluginM->fPtr)
  {
    ((FNLDSPtrfct)_pluginM->fPtr)(time, _n, &((*(_x[0]))(0)), &(*_M)(0, 0), _z->size(), &(*_z)(0));
  }
}

void FirstOrderNonLinearDS::computeM(double time, SP::SiconosVector x2)
{
  // second argument is useless at the time - Used in derived classes
  if (_pluginM->fPtr)
  {
    assert(x2->size() == _n && "FirstOrderNonLinearDS::computeM(t,x) x size does not fit with the system size.");
    ((FNLDSPtrfct)_pluginM->fPtr)(time, _n, &((*x2)(0)), &(*_M)(0, 0), _z->size(), &(*_z)(0));
  }
  // else nothing!
}

void FirstOrderNonLinearDS::computef(double time)
{
  if (_pluginf->fPtr)
    ((FNLDSPtrfct)_pluginf->fPtr)(time, _n, &((*(_x[0]))(0)) , &(*_f)(0), _z->size(), &(*_z)(0));
}

void FirstOrderNonLinearDS::computef(double time, SP::SiconosVector x2)
{
  if (_pluginf->fPtr)
    ((FNLDSPtrfct)_pluginf->fPtr)(time, _n, &((*x2)(0)) , &(*_f)(0), _z->size(), &(*_z)(0));
  // else nothing!
}

void FirstOrderNonLinearDS::computeJacobianfx(double time, bool)
{
  // second argument is useless at the time - Used in derived classes
  if (_pluginJacxf->fPtr)
    ((FNLDSPtrfct)_pluginJacxf->fPtr)(time, _n, &((*(_x[0]))(0)), &(*_jacobianfx)(0, 0), _z->size(), &(*_z)(0));
}

void FirstOrderNonLinearDS::computeJacobianfx(double time, SP::SiconosVector x2)
{
  // second argument is useless at the time - Used in derived classes
  if (_pluginJacxf->fPtr)
    ((FNLDSPtrfct)_pluginJacxf->fPtr)(time, _n, &((*x2)(0)), &(*_jacobianfx)(0, 0), _z->size(), &(*_z)(0));
}

void FirstOrderNonLinearDS::computeRhs(double time, bool)
{
  // second argument is useless at the time - Used in derived classes

  // compute rhs = M-1*( f + r ).

  *_x[1] = *_r; // Warning: p update is done in Interactions/Relations

  if (_f)
  {
    computef(time);
    *(_x[1]) += *_f;
  }

  if (_M)
  {
    // allocate invM at the first call of the present function
    if (! _invM)
      _invM.reset(new SimpleMatrix(*_M));
    _invM->PLUForwardBackwardInPlace(*_x[1]);
  }
}

void FirstOrderNonLinearDS::computeJacobianRhsx(double time, bool)
{
  // second argument is useless at the time - Used in derived classes

  // compute jacobian of rhs according to x, = M-1(jacobianfx + jacobianX(T.u))
  // At the time, second term is set to zero.
  computeJacobianfx(time);
  // solve M*jacobianXRhS = jacobianfx
  if (_M && _jacobianfx)
  {
    *_jacxRhs = *_jacobianfx;
    // copy _M into _invM for LU-factorisation, at the first call of this function.
    if (! _invM)
      _invM.reset(new SimpleMatrix(*_M));

    _invM->PLUForwardBackwardInPlace(*_jacxRhs);
  }
  // else jacobianRhsx = jacobianfx, pointers equality set in initRhs

}

// ===== XML MANAGEMENT FUNCTIONS =====

void FirstOrderNonLinearDS::saveSpecificDataToXML()
{
  // -- FirstOrderNonLinearDS  xml object --
  SP::FirstOrderNonLinearDSXML fonlds = boost::static_pointer_cast <FirstOrderNonLinearDSXML>(_dsxml);
  // --- other data ---
  if (!fonlds)
    RuntimeException::selfThrow("FirstOrderNonLinearDS::saveSpecificDataToXML - The DynamicalSystemXML object doesn't exists");

  if (_pluginf->fPtr)
    fonlds->setFPlugin(_pluginf->getPluginName());
  if (_pluginJacxf->fPtr)
    fonlds->setJacobianfxPlugin(_pluginJacxf->getPluginName());
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

