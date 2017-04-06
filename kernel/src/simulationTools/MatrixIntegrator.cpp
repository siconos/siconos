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
#include "Model.hpp"
#include "EventDriven.hpp"
#include "SubPluggedObject.hpp"
#include "LsodarOSI.hpp"
#include "TimeDiscretisation.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "EventsManager.hpp"

MatrixIntegrator::MatrixIntegrator(const DynamicalSystem& ds, const Model& m, SP::SiconosMatrix E): _E(E)
{
  commonInit(ds, m);
  _mat.reset(new SimpleMatrix(*E));
  _mat->zero();
}

MatrixIntegrator::MatrixIntegrator(const DynamicalSystem& ds, const Model& m, SP::PluggedObject plugin, const unsigned int p):
  _plugin(plugin)
{
  commonInit(ds, m);
  unsigned int n = ds.n();
  _mat.reset(new SimpleMatrix(n, p, 0));
  _spo.reset(new SubPluggedObject(*_plugin, n, p));
  std11::static_pointer_cast<FirstOrderLinearDS>(_DS)->setPluginB(_spo);
  _isConst = false;
}

MatrixIntegrator::MatrixIntegrator(const DynamicalSystem& ds, const Model& m)
{
  unsigned int n = ds.n();
  _mat.reset(new SimpleMatrix(n, n, 0));
  commonInit(ds, m);
}

void MatrixIntegrator::commonInit(const DynamicalSystem& ds, const Model& m)
{
  _TD.reset(new TimeDiscretisation(m.simulation()->eventsManager()->timeDiscretisation()));
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
     // std11::static_pointer_cast<FirstOrderLinearDS>(_DS)->zeroPlugin();
     if (cfolds.getPluginA()->isPlugged())
     {
       std11::static_pointer_cast<FirstOrderLinearDS>(_DS)->setPluginA(cfolds.getPluginA());
     }
     _isConst = (_TD->hConst()) && !(cfolds.getPluginA()->isPlugged()) ? true : false;
  }

  // integration stuff
  _model.reset(new Model(m.t0(), m.finalT()));
  _OSI.reset(new LsodarOSI());
  _model->nonSmoothDynamicalSystem()->insertDynamicalSystem(_DS);
  _sim.reset(new EventDriven(_TD, 0));
  _sim->prepareIntegratorForDS(_OSI, _DS, _model, m.t0());
  _model->setSimulation(_sim);
  _model->initialize();

  //change tolerance
  _OSI->setTol(1, 10 * MACHINE_PREC, 5 * MACHINE_PREC);

}

void MatrixIntegrator::integrate()
{
  SiconosVector& x = *_DS->x();
  SP::SiconosVector Ecol = static_cast<FirstOrderLinearDS&>(*_DS).b();
  if (!Ecol && _E)
  {
    Ecol.reset(new SiconosVector(_DS->n(), 0));
    static_cast<FirstOrderLinearDS&>(*_DS).setbPtr(Ecol);
  }
  unsigned int p = _mat->size(1);
  for (unsigned int i = 0; i < p; i++)
  {
    x.zero();
    if (_E)
      _E->getCol(i, *Ecol);
    else if (_plugin)
      _spo->setIndex(i);
    else
      x(i) = 1;
    //Reset LsodarOSI
    _sim->setIstate(1);
    _sim->advanceToEvent();
    _mat->setCol(i, x);
  }
  _sim->processEvents();
}
