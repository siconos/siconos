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
#include "ZOHTest.hpp"
using namespace std;

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(ZOHTest);


void ZOHTest::setUp()
{
  _t0 = 0;
  _T = 10;
  _h = 0.1;
  _n = 2;
  _tol = 1e-12;
  // Download a Model from Template.xml file
  _A.reset(new SimpleMatrix(_n, _n, 0));
  _b.reset(new SimpleVector(_n, 0));
  _x0.reset(new SimpleVector(_n, 0));
}

void ZOHTest::init()
{
  _DS.reset(new FirstOrderLinearTIDS(_x0, _A, _b));
  _TD.reset(new TimeDiscretisation(_t0, _h));
  _model.reset(new Model(_t0, _T));
  _model->nonSmoothDynamicalSystem()->insertDynamicalSystem(_DS);
  _sim.reset(new TimeStepping(_TD, 0));
  _ZOH.reset(new ZeroOrderHold(_DS));
  _sim->insertIntegrator(_ZOH);
  _model->initialize(_sim);
}

void ZOHTest::tearDown()
{}

void ZOHTest::testMatrixExp0()
{
  cout << "===========================================" << endl;
  cout << " ===== ZOH tests start ... ===== " << endl;
  cout << "===========================================" << endl;
  cout << "------- Compute matrix exponential of the identity matrix -------" << endl;
  _A->eye();
  init();
  _sim->computeOneStep();
  _sim->nextStep();
  SP::SiconosMatrix tmpM(new SimpleMatrix(_n, _n, 0));
  tmpM->eye();
  *tmpM = (*tmpM) * exp(_h);
  const SiconosMatrix& Phi = _ZOH->getPhi(*_DS);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testMatrixExp0 : ", Phi.size(0) == _n, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testMatrixExp0 : ", Phi.size(1) == _n, true);
  double diff = (*tmpM - Phi).normInf();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testMatrixExp0 : ", diff < _tol, true);
  cout << "------- First computation ok, error = " << diff << " -------" << endl;
  cout << endl << endl;
}

void ZOHTest::testMatrixExp1()
{
  cout << "===========================================" << endl;
  cout << " ===== ZOH tests start ... ===== " << endl;
  cout << "===========================================" << endl;
  cout << "------- Compute matrix exponential of a upper triangular matrix -------" << endl;
  _A->zero();
  (*_A)(0, 1) = 1;
  init();
  _sim->computeOneStep();
  _sim->nextStep();
  SP::SiconosMatrix tmpM(new SimpleMatrix(_n, _n, 0));
  tmpM->eye();
  (*tmpM)(0, 1) = _h;
  const SiconosMatrix& Phi = _ZOH->getPhi(*_DS);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testMatrixExp1 : ", Phi.size(0) == _n, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testMatrixExp1 : ", Phi.size(1) == _n, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testMatrixExp1 : ", (*tmpM - Phi).normInf() < _tol, true);
  cout << "------- Second computation ok, error = " << (*tmpM - Phi).normInf() << " -------" << endl;
  cout << endl << endl;
}

void ZOHTest::testMatrixIntegration1()
{
  cout << "===========================================" << endl;
  cout << " ===== ZOH tests start ... ===== " << endl;
  cout << "===========================================" << endl;
  cout << "------- Integrate an oscillator -------" << endl;
  _A->zero();
  (*_A)(0, 1) = 1;
  (*_A)(1, 0) = -1;
  _x0->zero();
  (*_x0)(0) = 1;
  init();
  SimpleMatrix dataPlot(ceil((_T - _t0) / _h) + 10, 3);
  SiconosVector& xProc = *_DS->x();
  unsigned int k = 0;
  dataPlot(0, 0) = _t0;
  dataPlot(0, 1) = (*_x0)(0);
  dataPlot(0, 2) = (*_x0)(1);
  while (_sim->nextTime() < _T)
  {
    _sim->computeOneStep();
    k++;
    dataPlot(k, 0) = _sim->nextTime();
    dataPlot(k, 1) = xProc(0);
    dataPlot(k, 2) = xProc(1);
    _sim->nextStep();
  }
  dataPlot.display();
  cout << endl << endl;
  ioMatrix io("testMatrixIntegration1.dat", "ascii");
  dataPlot.resize(k, 3);
  io.write(dataPlot, "noDim");
  // Reference Matrix
  SimpleMatrix dataPlotRef(dataPlot);
  dataPlotRef.zero();
  //magic line to compute the following:
  // python -c "import numpy as np; t = np.linspace(0, 9.9, 100); np.savetxt('testMatrixIntegration1.ref', np.transpose([t, np.cos(t), -np.sin(t)]))" && sed -i "1i100 3" testMatrixIntegration1.ref
  ioMatrix ref("testMatrixIntegration1.ref", "ascii");
  ref.read(dataPlotRef);
  cout << "------- Integration Ok, error = " << (dataPlot - dataPlotRef).normInf() << " -------" << endl;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testMatrixExp1 : ", (dataPlot - dataPlotRef).normInf() < _tol, true);
}

void ZOHTest::testMatrixIntegration2()
{
  cout << "===========================================" << endl;
  cout << " ===== ZOH tests start ... ===== " << endl;
  cout << "===========================================" << endl;
  cout << "------- Integrate x \\in -sgn(x)  -------" << endl;
  _A->zero();
  _x0->zero();
  (*_x0)(0) = 1;
  (*_x0)(1) = -1;
  SP::SiconosMatrix B(new SimpleMatrix(_n, _n));
  SP::SiconosMatrix C(new SimpleMatrix(_n, _n));
  B->eye();
  C->eye();
  SP::FirstOrderLinearTIR rel(new FirstOrderLinearTIR(C, B));
  SP::SimpleMatrix D(new SimpleMatrix(_n, _n, 0));
  rel->setDPtr(D);
  SP::NonSmoothLaw nslaw(new RelayNSL(_n));
  _DS.reset(new FirstOrderLinearTIDS(_x0, _A, _b));
  _TD.reset(new TimeDiscretisation(_t0, _h));
  _model.reset(new Model(_t0, _T));
  std::string id = "interaction for control";
  SP::Interaction inter(new Interaction(id, _DS, _n, _n, nslaw, rel));
  SP::NonSmoothDynamicalSystem myNSDS(new NonSmoothDynamicalSystem(_DS));
  myNSDS->insertInteraction(inter, true);
  _model->setNonSmoothDynamicalSystemPtr(myNSDS);
  _sim.reset(new TimeStepping(_TD, 1));
  _ZOH.reset(new ZeroOrderHold(_DS));
  _sim->insertIntegrator(_ZOH);
  SP::Relay osnspb(new Relay());
  _sim->insertNonSmoothProblem(osnspb);
  _model->initialize(_sim);
  SimpleMatrix dataPlot(ceil((_T - _t0) / _h) + 10, 5);
  SiconosVector& xProc = *_DS->x();
  SiconosVector& lambda = *inter->lambda(0);
  unsigned int k = 0;
  dataPlot(0, 0) = _t0;
  dataPlot(0, 1) = (*_x0)(0);
  dataPlot(0, 2) = (*_x0)(1);
  dataPlot(0, 3) = 0;
  dataPlot(0, 4) = 0;
  while (_sim->nextTime() < _T)
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
  dataPlot.display();
  cout << endl << endl;
  ioMatrix io("testMatrixIntegration2.dat", "ascii");
  dataPlot.resize(k, 5);
  io.write(dataPlot, "noDim");
  // Reference Matrix
  SimpleMatrix dataPlotRef(dataPlot);
  dataPlotRef.zero();
  ioMatrix ref("testMatrixIntegration2.ref", "ascii");
  ref.read(dataPlotRef);
  cout << "------- Integration Ok, error = " << (dataPlot - dataPlotRef).normInf() << " -------" << endl;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testMatrixExp2 : ", (dataPlot - dataPlotRef).normInf() < _tol, true);
}

void ZOHTest::testMatrixIntegration3()
{
  cout << "===========================================" << endl;
  cout << " ===== ZOH tests start ... ===== " << endl;
  cout << "===========================================" << endl;
  cout << "------- Integrate Orlov's controller  -------" << endl;
  _h = .001;
  _A->zero();
  (*_A)(0, 1) = 1;
  _x0->zero();
  (*_x0)(0) = 1;
  (*_x0)(1) = 1;
  SP::SiconosMatrix B(new SimpleMatrix(_n, _n, 0));
  SP::SiconosMatrix C(new SimpleMatrix(_n, _n));
  (*B)(1, 0) = 2;
  (*B)(1, 1) = 1;
  C->eye();
  SP::FirstOrderLinearTIR rel(new FirstOrderLinearTIR(C, B));
  SP::SimpleMatrix D(new SimpleMatrix(_n, _n, 0));
  rel->setDPtr(D);
  SP::NonSmoothLaw nslaw(new RelayNSL(_n));
  _DS.reset(new FirstOrderLinearTIDS(_x0, _A, _b));
  _TD.reset(new TimeDiscretisation(_t0, _h));
  _model.reset(new Model(_t0, _T));
  std::string id = "interaction for control";
  SP::Interaction inter(new Interaction(id, _DS, _n, _n, nslaw, rel));
  SP::NonSmoothDynamicalSystem myNSDS(new NonSmoothDynamicalSystem(_DS));
  myNSDS->insertInteraction(inter, true);
  _model->setNonSmoothDynamicalSystemPtr(myNSDS);
  _sim.reset(new TimeStepping(_TD, 1));
  _ZOH.reset(new ZeroOrderHold(_DS));
  _sim->insertIntegrator(_ZOH);
  SP::Relay osnspb(new Relay());
  _sim->insertNonSmoothProblem(osnspb);
  _model->initialize(_sim);
  SimpleMatrix dataPlot(ceil((_T - _t0) / _h) + 10, 7);
  SiconosVector& xProc = *_DS->x();
  SiconosVector& lambda = *inter->lambda(0);
  SimpleVector sampledControl(_n);
  unsigned int k = 0;
  dataPlot(0, 0) = _t0;
  dataPlot(0, 1) = (*_x0)(0);
  dataPlot(0, 2) = (*_x0)(1);
  dataPlot(0, 3) = 0;
  dataPlot(0, 4) = 0;
  dataPlot(0, 5) = 0;
  dataPlot(0, 6) = 0;
  while (_sim->nextTime() < _T)
  {
    _sim->computeOneStep();
    prod(*B, lambda, sampledControl, true);
    k++;
    dataPlot(k, 0) = _sim->nextTime();
    dataPlot(k, 1) = xProc(0);
    dataPlot(k, 2) = xProc(1);
    dataPlot(k, 3) = sampledControl(0);
    dataPlot(k, 4) = sampledControl(1);
    dataPlot(k, 5) = lambda(0);
    dataPlot(k, 6) = lambda(1);
    _sim->nextStep();
  }
  dataPlot.display();
  cout << endl << endl;
  ioMatrix io("testMatrixIntegration3.dat", "ascii");
  dataPlot.resize(k, 7);
  io.write(dataPlot, "noDim");
  // Reference Matrix
  SimpleMatrix dataPlotRef(dataPlot);
  dataPlotRef.zero();
  ioMatrix ref("testMatrixIntegration3.ref", "ascii");
  ref.read(dataPlotRef);
  cout << "------- Integration Ok, error = " << (dataPlot - dataPlotRef).normInf() << " -------" << endl;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testMatrixExp3 : ", (dataPlot - dataPlotRef).normInf() < _tol, true);
}

void ZOHTest::testMatrixIntegration4()
{
  cout << "===========================================" << endl;
  cout << " ===== ZOH tests start ... ===== " << endl;
  cout << "===========================================" << endl;
  cout << "------- Integrate Orlov's controller  -------" << endl;
  _h = .001;
  _A->zero();
  (*_A)(0, 1) = 1;
  _x0->zero();
  (*_x0)(0) = 1;
  (*_x0)(1) = 1;
  SP::SiconosMatrix B(new SimpleMatrix(_n, _n, 0));
  SP::SiconosMatrix C(new SimpleMatrix(_n, _n));
  (*B)(1, 0) = 2;
  (*B)(1, 1) = 1;
  C->eye();
  SP::FirstOrderLinearTIR rel(new FirstOrderLinearTIR(C, B));
  SP::SimpleMatrix D(new SimpleMatrix(_n, _n, 0));
  rel->setDPtr(D);
  SP::NonSmoothLaw nslaw(new RelayNSL(_n));
  _DS.reset(new FirstOrderLinearDS(_x0, _A, _b));
  _TD.reset(new TimeDiscretisation(_t0, _h));
  _model.reset(new Model(_t0, _T));
  std::string id = "interaction for control";
  SP::Interaction inter(new Interaction(id, _DS, _n, _n, nslaw, rel));
  SP::NonSmoothDynamicalSystem myNSDS(new NonSmoothDynamicalSystem(_DS));
  myNSDS->insertInteraction(inter, true);
  _model->setNonSmoothDynamicalSystemPtr(myNSDS);
  _sim.reset(new TimeStepping(_TD, 1));
  _ZOH.reset(new ZeroOrderHold(_DS));
  _sim->insertIntegrator(_ZOH);
  SP::Relay osnspb(new Relay());
  _sim->insertNonSmoothProblem(osnspb);
  _model->initialize(_sim);
  SimpleMatrix dataPlot(ceil((_T - _t0) / _h) + 10, 7);
  SiconosVector& xProc = *_DS->x();
  SiconosVector& lambda = *inter->lambda(0);
  SimpleVector sampledControl(_n);
  unsigned int k = 0;
  dataPlot(0, 0) = _t0;
  dataPlot(0, 1) = (*_x0)(0);
  dataPlot(0, 2) = (*_x0)(1);
  dataPlot(0, 3) = 0;
  dataPlot(0, 4) = 0;
  dataPlot(0, 5) = 0;
  dataPlot(0, 6) = 0;
  while (_sim->nextTime() < _T)
  {
    _sim->computeOneStep();
    prod(*B, lambda, sampledControl, true);
    k++;
    dataPlot(k, 0) = _sim->nextTime();
    dataPlot(k, 1) = xProc(0);
    dataPlot(k, 2) = xProc(1);
    dataPlot(k, 3) = sampledControl(0);
    dataPlot(k, 4) = sampledControl(1);
    dataPlot(k, 5) = lambda(0);
    dataPlot(k, 6) = lambda(1);
    _sim->nextStep();
  }
  dataPlot.display();
  cout << endl << endl;
  ioMatrix io("testMatrixIntegration4.dat", "ascii");
  dataPlot.resize(k, 7);
  io.write(dataPlot, "noDim");
  // Reference Matrix
  SimpleMatrix dataPlotRef(dataPlot);
  dataPlotRef.zero();
  ioMatrix ref("testMatrixIntegration4.ref", "ascii");
  ref.read(dataPlotRef);
  cout << "------- Integration Ok, error = " << (dataPlot - dataPlotRef).normInf() << " -------" << endl;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testMatrixExp4 : ", (dataPlot - dataPlotRef).normInf() < _tol, true);
}
