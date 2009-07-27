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

using namespace std;

// ===== CONSTRUCTORS =====

// From a minimum set of data
FirstOrderNonLinearDS::FirstOrderNonLinearDS(const SiconosVector& newX0): DynamicalSystem(DS::FONLDS, newX0.size())
{
  // == Initial conditions ==
  x0.reset(new SimpleVector(newX0));

  // == Current state ==
  // x is composed of two blocks of size n, x[0] = \f$ x \f$ and x[1]=\f$ \dot x \f$.
  // x[0] initialized with x0.

  x[0].reset(new SimpleVector(*x0));
  x[1].reset(new SimpleVector(n));

  //mG
  mG_alpha.reset(new SimpleVector(n));
  mResidur.reset(new SimpleVector(n));
  mXp.reset(new SimpleVector(getDim()));
  mXq.reset(new SimpleVector(getDim()));
  mXfree.reset(new SimpleVector(getDim()));
  mfold.reset(new SimpleVector(getDim()));

  // == r ==

  r.reset(new SimpleVector(n));

  checkDynamicalSystem();
}

// From XML file
FirstOrderNonLinearDS::FirstOrderNonLinearDS(SP::DynamicalSystemXML dsXML):
  DynamicalSystem(dsXML)
{
  // -- FONLDS xml object --
  SP::FirstOrderNonLinearDSXML fonlds = boost::static_pointer_cast <FirstOrderNonLinearDSXML>(dsxml);

  // === Initial conditions ===
  // Warning: n is set thanks to x0 size
  if (! fonlds->hasX0())
    RuntimeException::selfThrow("FirstOrderNonLinearDS:: xml constructor, x0 is a required input");

  x0.reset(new SimpleVector(fonlds->getX0()));

  n = x0->size();

  // === Current state (optional input) ===
  // x is composed of two blocks of size n, (*x)[0] = \f$ x \f$ and (*x)[1]=\f$ \dot x \f$.

  if (fonlds->hasX())
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
  if (fonlds->hasM())
  {
    if (fonlds->isMPlugin())
    {
      plugin = fonlds->getMPlugin();
      setComputeMFunction(SSL::getPluginName(plugin), SSL::getPluginFunctionName(plugin));
    }
    else // This means that M is constant
    {
      M.reset(new PMJF(fonlds->getMMatrix()));
      if (M->size(0) != n || M->size(1) != n)
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
      if (fonlds->getFVector().size() != n)
        RuntimeException::selfThrow("FirstOrderNonLinearDS:: xml constructor, f size differs from n!");

      mf.reset(new PVF(fonlds->getFVector()));
    }
  }

  if (fonlds->hasJacobianXF())
  {
    if (fonlds->isJacobianXFPlugin())
    {
      plugin = fonlds->getJacobianXFPlugin();
      setComputeJacobianXFFunction(SSL::getPluginName(plugin), SSL::getPluginFunctionName(plugin));
    }
    else // This means that jacobianXF is constant
    {
      jacobianXF.reset(new PMJF(fonlds->getJacobianXFMatrix()));
      if (jacobianXF->size(0) != n || jacobianXF->size(1) != n)
        RuntimeException::selfThrow("FirstOrderNonLinearDS:: xml constructor, jacobianXF size differs from n!");
    }
  }

  // Memory
  if (fonlds->hasXMemory())
    xMemory.reset(new SiconosMemory(fonlds->getXMemoryXML()));

  checkDynamicalSystem();
}

// From a minimum set of data
FirstOrderNonLinearDS::FirstOrderNonLinearDS(const SiconosVector& newX0, const string& fPlugin, const string& jacobianXFPlugin):
  DynamicalSystem(DS::FONLDS, newX0.size())
{
  // == Initial conditions ==
  x0.reset(new SimpleVector(newX0));

  // == Current state ==
  // x is composed of two blocks of size n, x[0] = \f$ x \f$ and x[1]=\f$ \dot x \f$.
  // x[0] initialized with x0.

  x[0].reset(new SimpleVector(*x0));
  x[1].reset(new SimpleVector(n));
  //mG
  mG_alpha.reset(new SimpleVector(n));
  mResidur.reset(new SimpleVector(n));
  mXp.reset(new SimpleVector(getDim()));
  mXq.reset(new SimpleVector(getDim()));
  mXfree.reset(new SimpleVector(getDim()));
  r.reset(new SimpleVector(getDim()));
  mfold.reset(new SimpleVector(getDim()));

  // == r ==

  r.reset(new SimpleVector(n));

  // == f and its jacobian ==
  // Allocation and link with the plug-in
  setComputeFFunction(SSL::getPluginName(fPlugin), SSL::getPluginFunctionName(fPlugin));
  setComputeJacobianXFFunction(SSL::getPluginName(jacobianXFPlugin), SSL::getPluginFunctionName(jacobianXFPlugin));
  checkDynamicalSystem();
}
void FirstOrderNonLinearDS::preparStep()
{
  mXp->zero();
  r->zero();
};
bool FirstOrderNonLinearDS::checkDynamicalSystem()
{
  DynamicalSystem::checkDynamicalSystem();
  bool output = DynamicalSystem::checkDynamicalSystem();
  if (!output) cout << "FirstOrderNonLinearDS Warning: your dynamical system seems to be uncomplete (check = false)" << endl;
  return output;
}


void FirstOrderNonLinearDS::setM(const PMJF& newValue)
{
  assert(newValue.size(0) == n && "FirstOrderNonLinearDS - setM: inconsistent dimensions with problem size for input matrix M.");
  assert(newValue.size(1) == n && "FirstOrderNonLinearDS - setM: inconsistent dimensions with problem size for input matrix M.");

  if (! M)
    M.reset(new PMJF(newValue));
  else
    *M = newValue;
}

void FirstOrderNonLinearDS::setInvM(const SiconosMatrix& newValue)
{
  if (newValue.size(0) != n || newValue.size(1) != n)
    RuntimeException::selfThrow("FirstOrderNonLinearDS::setInvM: inconsistent dimensions with problem size for input matrix.");

  if (! invM)
    invM.reset(new SimpleMatrix(n, n));
  *invM = newValue;
}

void FirstOrderNonLinearDS::setInvMPtr(SP::SiconosMatrix newPtr)
{
  invM = newPtr;
}

void FirstOrderNonLinearDS::setF(const PVF& newValue)
{
  assert(newValue.size() == n && "FirstOrderNonLinearDS - setF: inconsistent dimensions with problem size for input vector f");

  if (! mf)
    mf.reset(new PVF(newValue));
  else
    *mf = newValue;
}

void FirstOrderNonLinearDS::setJacobianXF(const PMJF& newValue)
{
  assert(newValue.size(0) == n && "FirstOrderNonLinearDS - setJacobianXF: inconsistent dimensions with problem size for input matrix M.");
  assert(newValue.size(1) == n && "FirstOrderNonLinearDS - setJacobianXF: inconsistent dimensions with problem size for input matrix M.");

  if (! jacobianXF)
    jacobianXF.reset(new PMJF(newValue));
  else
    *jacobianXF = newValue;
}

void FirstOrderNonLinearDS::initRhs(double time)
{
  // compute initial values for f and jacobianXF, initialize right-hand side.
  computeRhs(time); // this will compute, if required, f and M.

  if (! jacobianXRhs)  // if not allocated with a set or anything else
  {
    if (jacobianXF && ! M)  // if M is not defined, then jacobianXF = jacobianXRhs, no memory allocation for that one.
      jacobianXRhs = jacobianXF;
    else if (jacobianXF && M)
      jacobianXRhs.reset(new PMJF(n, n));

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
  r->zero();
  *(x[0]) = *x0;

  // If z has not been set, we initialize it with a null vector of size 1, since z is required in plug-in functions call.
  if (! z)
    z.reset(new SimpleVector(1));

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
    rMemory.reset(new SiconosMemory(steps));
}

void FirstOrderNonLinearDS::swapInMemory()
{
  xMemory->swap(x[0]);
  rMemory->swap(r);
  *mfold = *mf;
}

// ===== COMPUTE PLUGINS FUNCTIONS =====

void FirstOrderNonLinearDS::setComputeMFunction(const string& pluginPath, const string& functionName)
{
  if (!M)
    M.reset(new PMJF(n, n));
  M->setComputeFunction(pluginPath, functionName);
}

void FirstOrderNonLinearDS::setComputeMFunction(FPtr1 fct)
{
  if (!M)
    M.reset(new PMJF(n, n));
  M->setComputeFunction(fct);
}

void FirstOrderNonLinearDS::setComputeFFunction(const string& pluginPath, const string& functionName)
{
  if (! mf)
    mf.reset(new PVF(n));
  mf->setComputeFunction(pluginPath, functionName);
}

void FirstOrderNonLinearDS::setComputeFFunction(FPtr1 fct)
{
  if (! mf)
    mf.reset(new PVF(n));
  mf->setComputeFunction(fct);
}

void FirstOrderNonLinearDS::setComputeJacobianXFFunction(const string& pluginPath, const string& functionName)
{
  if (!jacobianXF)
    jacobianXF.reset(new PMJF(n, n));
  jacobianXF->setComputeFunction(pluginPath, functionName);
}

void FirstOrderNonLinearDS::setComputeJacobianXFFunction(FPtr1 fct)
{
  if (!jacobianXF)
    jacobianXF.reset(new PMJF(n, n));
  jacobianXF->setComputeFunction(fct);
}

void FirstOrderNonLinearDS::computeM(double time)
{
  // second argument is useless at the time - Used in derived classes
  if (M->isPlugged())
  {
    if (!M->fPtr)
      RuntimeException::selfThrow("FirstOrderNonLinearDS::computeM() is not linked to a plugin function");
    (M->fPtr)(time, n, &((*(x[0]))(0)), &(*M)(0, 0), z->size(), &(*z)(0));
  }
  // else nothing!
}

void FirstOrderNonLinearDS::computeM(double time, SP::SiconosVector x2)
{
  // second argument is useless at the time - Used in derived classes
  if (M->isPlugged())
  {
    if (!M->fPtr)
      RuntimeException::selfThrow("FirstOrderNonLinearDS::computeM() is not linked to a plugin function");
    if (x2->size() != n)
      RuntimeException::selfThrow("FirstOrderNonLinearDS::computeM(t,x) x size does not fit with the system size.");

    (M->fPtr)(time, n, &((*x2)(0)), &(*M)(0, 0), z->size(), &(*z)(0));
  }
  // else nothing!
}

void FirstOrderNonLinearDS::computeF(double time)
{
  if (mf->isPlugged())
  {
    if (!mf->fPtr)
      RuntimeException::selfThrow("FirstOrderNonLinearDS::computeF() f is not linked to a plugin function");
    (mf->fPtr)(time, n, &((*(x[0]))(0)) , &(*mf)(0), z->size(), &(*z)(0));
  }
  // else nothing!
}

void FirstOrderNonLinearDS::computeF(double time, SP::SiconosVector x2)
{
  if (mf->isPlugged())
  {
    if (!mf->fPtr)
      RuntimeException::selfThrow("FirstOrderNonLinearDS::computeF() f is not linked to a plugin function");
    if (x2->size() != n)
      RuntimeException::selfThrow("FirstOrderNonLinearDS::computeF(t,x) x size does not fit with the system size.");
    (mf->fPtr)(time, n, &((*x2)(0)) , &(*mf)(0), z->size(), &(*z)(0));
  }
  // else nothing!
}

void FirstOrderNonLinearDS::computeJacobianXF(double time, bool)
{
  // second argument is useless at the time - Used in derived classes
  if (jacobianXF->isPlugged())
  {
    if (!jacobianXF->fPtr)
      RuntimeException::selfThrow("FirstOrderNonLinearDS::computeJacobianXF() is not linked to a plugin function");
    (jacobianXF->fPtr)(time, n, &((*(x[0]))(0)), &(*jacobianXF)(0, 0), z->size(), &(*z)(0));
  }
  // else nothing!
}

void FirstOrderNonLinearDS::computeJacobianXF(double time, SP::SiconosVector x2)
{
  // second argument is useless at the time - Used in derived classes
  if (jacobianXF->isPlugged())
  {
    if (!jacobianXF->fPtr)
      RuntimeException::selfThrow("FirstOrderNonLinearDS::computeJacobianXF() is not linked to a plugin function");
    if (x2->size() != n)
      RuntimeException::selfThrow("FirstOrderNonLinearDS::computeJacobianXF(t,x) x size does not fit with the system size.");

    (jacobianXF->fPtr)(time, n, &((*x2)(0)), &(*jacobianXF)(0, 0), z->size(), &(*z)(0));
  }
  // else nothing!
}

void FirstOrderNonLinearDS::computeRhs(double time, bool)
{
  // second argument is useless at the time - Used in derived classes

  // compute rhs = M-1*( f + r ).

  *x[1] = *r; // Warning: p update is done in Interactions/Relations

  if (mf)
  {
    computeF(time);
    *(x[1]) += *mf;
  }

  if (M)
  {
    // allocate invM at the first call of the present function
    if (! invM)
      invM.reset(new SimpleMatrix(*M));
    invM->PLUForwardBackwardInPlace(*x[1]);
  }
}

void FirstOrderNonLinearDS::computeJacobianXRhs(double time, bool)
{
  // second argument is useless at the time - Used in derived classes

  // compute jacobian of rhs according to x, = M-1(jacobianXF + jacobianX(T.u))
  // At the time, second term is set to zero.
  computeJacobianXF(time);
  // solve M*jacobianXRhS = jacobianXF
  if (M && jacobianXF)
  {
    *jacobianXRhs = *jacobianXF;
    // copy M into invM for LU-factorisation, at the first call of this function.
    if (! invM)
      invM.reset(new SimpleMatrix(*M));

    invM->PLUForwardBackwardInPlace(*jacobianXRhs);
  }
  // else jacobianXRhs = jacobianXF, pointers equality set in initRhs

}

// ===== XML MANAGEMENT FUNCTIONS =====

void FirstOrderNonLinearDS::saveSpecificDataToXML()
{
  // -- FirstOrderNonLinearDS  xml object --
  SP::FirstOrderNonLinearDSXML fonlds = boost::static_pointer_cast <FirstOrderNonLinearDSXML>(dsxml);
  // --- other data ---
  if (!fonlds)
    RuntimeException::selfThrow("FirstOrderNonLinearDS::saveSpecificDataToXML - The DynamicalSystemXML object doesn't exists");

  if (mf->isPlugged())
    fonlds->setFPlugin(mf->getPluginName());
  if (jacobianXF->isPlugged())
    fonlds->setJacobianXFPlugin(jacobianXF->getPluginName());
}

// ===== MISCELLANEOUS ====

void FirstOrderNonLinearDS::display() const
{
  cout << " =====> First Order Non Linear DS (number: " << number << ")." << endl;
  cout << "- n (size) : " << n << endl;
  cout << "- x " << endl;
  if (x[0]) x[0]->display();
  else cout << "-> NULL" << endl;
  cout << "- x0 " << endl;
  if (x0) x0->display();
  else cout << "-> NULL" << endl;
  cout << "- M: " << endl;
  if (M) M->display();
  else cout << "-> NULL" << endl;
  cout << " ============================================" << endl;
}

void FirstOrderNonLinearDS::resetNonSmoothPart()
{
  r->zero();
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

