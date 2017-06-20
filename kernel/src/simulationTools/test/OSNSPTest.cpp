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
#include "OSNSPTest.hpp"
#include "EventsManager.hpp"

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(OSNSPTest);


void OSNSPTest::setUp()
{
  _A.reset(new SimpleMatrix(_n, _n, 0));
  _b.reset(new SiconosVector(_n, 0));
  _x0.reset(new SiconosVector(_n, 0));
}

void OSNSPTest::init()
{
  _DS.reset(new FirstOrderLinearTIDS(_x0, _A, _b));
  _TD.reset(new TimeDiscretisation(_t0, _h));
  _model.reset(new Model(_t0, _T));
  _osi.reset(new EulerMoreauOSI(_theta));
  _model->nonSmoothDynamicalSystem()->insertDynamicalSystem(_DS);
  _sim.reset(new TimeStepping(_TD, 0));
  _sim->setNonSmoothDynamicalSystemPtr(_model->nonSmoothDynamicalSystem());
  _sim->prepareIntegratorForDS(_osi, _DS, _model, _t0);
  _model->setSimulation(_sim);
  _model->initialize();
}

void OSNSPTest::tearDown()
{}

void OSNSPTest::testAVI()
{
  std::cout << "===========================================" <<std::endl;
  std::cout << " ===== OSNSP tests start ... ===== " <<std::endl;
  std::cout << "===========================================" <<std::endl;
  std::cout << "------- implicit Twisting relation  -------" <<std::endl;
  _h = 1e-1;
  _T = 20.0;
  double G = 10.0;
  double beta = .3;
  _A->zero();
  (*_A)(0, 1) = 1.0;
  _x0->zero();
  (*_x0)(0) = 10.0;
  (*_x0)(1) = 10.0;
  SP::SimpleMatrix B(new SimpleMatrix(_n, _n, 0));
  SP::SimpleMatrix C(new SimpleMatrix(_n, _n));
  (*B)(1, 0) = G;
  (*B)(1, 1) = G*beta;
  C->eye();
  SP::FirstOrderLinearTIR rel(new FirstOrderLinearTIR(C, B));
  SP::SimpleMatrix H(new SimpleMatrix(4, 2));
  (*H)(0, 0) = 1.0;
  (*H)(1, 0) = -_h/2.0;
  (*H)(2, 0) = -1.0;
  (*H)(3, 0) = _h/2.0;
  (*H)(1, 1) = 1.0;
  (*H)(3, 1) = -1.0;
  SP::SiconosVector K(new SiconosVector(4));
  (*K)(0) = -1.0;
  (*K)(1) = -1.0;
  (*K)(2) = -1.0;
  (*K)(3) = -1.0;
  SP::NonSmoothLaw nslaw(new NormalConeNSL(_n, H, K));
  _DS.reset(new FirstOrderLinearTIDS(_x0, _A, _b));
  _TD.reset(new TimeDiscretisation(_t0, _h));
  _model.reset(new Model(_t0, _T));
  SP::Interaction inter(new Interaction(nslaw, rel));
  _osi.reset(new EulerMoreauOSI(_theta));
  _model->nonSmoothDynamicalSystem()->insertDynamicalSystem(_DS);
  _model->nonSmoothDynamicalSystem()->link(inter, _DS);
  _sim.reset(new TimeStepping(_TD));
  _sim->setNonSmoothDynamicalSystemPtr(_model->nonSmoothDynamicalSystem());
  _sim->prepareIntegratorForDS(_osi, _DS, _model, _t0);
  SP::AVI osnspb(new AVI());
  _sim->insertNonSmoothProblem(osnspb);
  _model->setSimulation(_sim);
  _model->initialize();
  SimpleMatrix dataPlot((unsigned)ceil((_T - _t0) / _h) + 10, 5);
  SiconosVector& xProc = *_DS->x();
  SiconosVector& lambda = *inter->lambda(0);
  unsigned int k = 0;
  dataPlot(0, 0) = _t0;
  dataPlot(0, 1) = (*_x0)(0);
  dataPlot(0, 2) = (*_x0)(1);
  dataPlot(0, 3) = -1.0;
  dataPlot(0, 4) = -1.0;
  while (_sim->hasNextEvent())
  {
    _sim->computeOneStep();
    k++;
    dataPlot(k, 0) = _sim->nextTime();
    dataPlot(k, 1) = xProc(0);
    dataPlot(k, 2) = xProc(1);
    dataPlot(k, 3) = lambda(0);
    dataPlot(k, 4) = lambda(1);
    _sim->nextStep();
  }
  std::cout <<std::endl <<std::endl;
  dataPlot.resize(k, dataPlot.size(1));
  ioMatrix::write("testAVI.dat", "ascii", dataPlot, "noDim");
  // Reference Matrix
  SimpleMatrix dataPlotRef(dataPlot);
  dataPlotRef.zero();
  ioMatrix::read("testAVI.ref", "ascii", dataPlotRef);
  SP::SiconosVector err(new SiconosVector(dataPlot.size(1)));
  (dataPlot - dataPlotRef).normInfByColumn(err);
  err->display();

  double maxErr = err->getValue(0) > err->getValue(1) ? (err->getValue(0) > err->getValue(2) ? err->getValue(0) : err->getValue(2)) : (err->getValue(1) > err->getValue(2) ? err->getValue(1) : err->getValue(2));

  std::cout << "------- Integration Ok, error = " << maxErr << " -------" <<std::endl;
  if (maxErr > _tol)
  {
    dataPlot.display();
    (dataPlot - dataPlotRef).display();
  }
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAVI : ",  maxErr < _tol, true);
}
