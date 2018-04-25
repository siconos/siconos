/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/

#include "MatrixIntegrator.hpp"
#include "SiconosAlgebra.hpp"
#include "FirstOrderLinearTIDS.hpp"
#include "EventDriven.hpp"
#include "SubPluggedObject.hpp"
#include "LsodarOSI.hpp"
#include "TimeDiscretisation.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "EventsManager.hpp"

//#define DEBUG_WHERE_MESSAGES
#define DEBUG_NOCOLOR
#define DEBUG_STDOUT
#define DEBUG_MESSAGES
#include <debug.h>

MatrixIntegrator::MatrixIntegrator(const DynamicalSystem& ds, const NonSmoothDynamicalSystem& nsds,
                                   const  TimeDiscretisation & td, SP::SiconosMatrix E): _E(E)
{
  commonInit(ds, nsds, td);
  _mat.reset(new SimpleMatrix(*E));
  _mat->zero();
}

MatrixIntegrator::MatrixIntegrator(const DynamicalSystem& ds, const NonSmoothDynamicalSystem& nsds,
                                   const TimeDiscretisation & td, SP::PluggedObject plugin, const unsigned int p):
  _plugin(plugin)
{
  commonInit(ds, nsds, td);
  unsigned int n = ds.n();
  _mat.reset(new SimpleMatrix(n, p, 0));
  _spo.reset(new SubPluggedObject(*_plugin, n, p));
  std11::static_pointer_cast<FirstOrderLinearDS>(_DS)->setPluginB(_spo);
  _isConst = false;
}

MatrixIntegrator::MatrixIntegrator(const DynamicalSystem& ds, const NonSmoothDynamicalSystem& nsds, const  TimeDiscretisation & td)
{
  unsigned int n = ds.n();
  _mat.reset(new SimpleMatrix(n, n, 0));
  commonInit(ds, nsds, td);
}

void MatrixIntegrator::commonInit(const DynamicalSystem& ds, const NonSmoothDynamicalSystem& nsds, const TimeDiscretisation & td)
{
  DEBUG_BEGIN("MatrixIntegrator::commonInit(...)\n");
  _TD.reset(new TimeDiscretisation(td));

  DEBUG_EXPR(ds.display(););


  Type::Siconos dsType = Type::value(ds);
  if (dsType == Type::FirstOrderLinearTIDS)
  {
     _DS.reset(new FirstOrderLinearTIDS(static_cast<const FirstOrderLinearTIDS&>(ds)));
     _isConst = _TD->hConst();
  }
  else if (dsType == Type::FirstOrderLinearDS)
  {
     const FirstOrderLinearDS& cfolds = static_cast<const FirstOrderLinearDS&>(ds);
     _DS.reset(new FirstOrderLinearDS(cfolds));
     if (_E)
     {
       DEBUG_PRINT("Matrix integration of the constant (possibly time variant) part\n");
       if (cfolds.getPluginB()->isPlugged())
       {
         DEBUG_PRINT("B is plugged\n");
         std11::static_pointer_cast<FirstOrderLinearDS>(_DS)->setPluginB(cfolds.getPluginB());
       }
       _isConst = (_TD->hConst()) && (cfolds.hasConstantB()) ? true : false;
     }
     else
     {
       DEBUG_PRINT("Matrix integration of the linear (possibly time variant) part\n");
       if (cfolds.getPluginA()->isPlugged())
       {
         DEBUG_PRINT("A is plugged\n");
         std11::static_pointer_cast<FirstOrderLinearDS>(_DS)->setPluginA(cfolds.getPluginA());
       }
       _isConst = (_TD->hConst()) && (cfolds.hasConstantA()) ? true : false;
     }
     //_isConst = (_TD->hConst()) && !(cfolds.getPluginA()->isPlugged()) && !(cfolds.getPluginB()->isPlugged()) ? true : false;
  }

  _DS->setNumber(9999999);
  DEBUG_EXPR(_DS->display(););
  // integration stuff
  _nsds.reset(new NonSmoothDynamicalSystem());
  _nsds->sett0(nsds.t0());
  _nsds->setT(nsds.finalT());


  _OSI.reset(new LsodarOSI());
  _nsds->insertDynamicalSystem(_DS);
  _sim.reset(new EventDriven(_nsds, _TD, 0));
  _sim->associate(_OSI, _DS);
  _sim->setName("Matrix integrator simulation");
  //change tolerance
  _OSI->setTol(1, 10 * MACHINE_PREC, 5 * MACHINE_PREC);
  DEBUG_END("MatrixIntegrator::commonInit(...)\n");
}

void MatrixIntegrator::integrate()
{
  DEBUG_BEGIN("MatrixIntegrator::integrate()\n");
  DEBUG_EXPR(std::cout << name() << std::endl ; );
  SiconosVector& x0 = *_DS->x0();
  SiconosVector& x = *_DS->x();

  SP::SiconosVector Ecol = static_cast<FirstOrderLinearDS&>(*_DS).b();
 
  
  if (!Ecol && _E)
  {
    DEBUG_PRINT("!Ecol && _E. we reset Ecol\n");
    Ecol.reset(new SiconosVector(_DS->n(), 0));
    static_cast<FirstOrderLinearDS&>(*_DS).setbPtr(Ecol);
  }
  
  unsigned int p = _mat->size(1);
  for (unsigned int i = 0; i < p; i++)
  {
    
    x0.zero();
    if (_E)
      _E->getCol(i, *Ecol);
    else if (_plugin)
      _spo->setIndex(i);
    else
    {
      x0(i) = 1;
    }
    
    DEBUG_EXPR(if (_E) _E->display(););
    DEBUG_EXPR(Ecol->display(););
    DEBUG_EXPR(x0.display(););
    DEBUG_PRINTF("integration interval = [%g,%g]", _sim->getTk(), _sim->getTkp1());
    //Reset LsodarOSI
    //_OSI->setIsInitialized(false);
    _DS->resetToInitialState();
    _sim->setIstate(1);
    _sim->advanceToEvent();
    _mat->setCol(i, x);
  }

  _sim->processEvents();
  //_DS->resetToInitialState();

  DEBUG_EXPR(_mat->display(););
  DEBUG_EXPR(_DS->display(););
  DEBUG_END("MatrixIntegrator::integrate()\n");
}
