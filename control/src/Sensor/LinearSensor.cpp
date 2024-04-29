/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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
#include "LinearSensor.hpp"

#include "DynamicalSystem.hpp"
#include "SiconosException.hpp"
#include "SiconosMatrixVectorOp.hpp"  // mat-vec prod
#include "SiconosVector.hpp"
#include "SiconosMatrix.hpp"
siconos::control::LinearSensor::LinearSensor(
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds)
    : ControlSensor(SensorType::Linear, ds) {}

siconos::control::LinearSensor::LinearSensor(
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds,
    std::shared_ptr<siconos::algebra::SiconosMatrix> matC,
    std::shared_ptr<siconos::algebra::SiconosMatrix> matD)
    : ControlSensor(SensorType::Linear, ds), _matC(matC), _matD(matD) {}

void siconos::control::LinearSensor::initialize(
    const siconos::modeling::NonSmoothDynamicalSystem& nsds) {
  // Call initialize of base class
  ControlSensor::initialize(nsds);

  // consistency checks
  if (!_matC) {
    THROW_EXCEPTION("LinearSensor::initialize - no C matrix was given");
  }

  auto colC = _matC->size(1);
  auto rowC = _matC->size(0);
  // What happen here if we have more than one DS ?
  // This may be unlikely to happen.
  //  _DS = _model->nonSmoothDynamicalSystem()->dynamicalSystemNumber(0);
  if (colC != _DS->n()) {
    THROW_EXCEPTION(
        " LinearSensor::initialize - The number of column of the C matrix must be equal to "
        "the length of x");
  }
  if (_matD) {
    auto rowD = _matD->size(0);
    if (rowC != rowD) {
      THROW_EXCEPTION("C and D must have the same number of rows");
    }
  }

  // --- Get the values ---
  // -> saved in a matrix data
  // -> event
  _storedY = std::make_shared<siconos::algebra::SiconosVector>(rowC);
  //  (_data[_eSensor])["StoredY"] = storedY;
  // set the dimension of the output
  *_storedY = siconos::algebra::prod(*_matC, *_DSx);
}

void siconos::control::LinearSensor::capture() {
  *_storedY = siconos::algebra::prod(*_matC, *_DSx);
  // untested
  if (_matD)
    //    *_storedY += siconos::algebra::prod(*_matD, *_DS->z());
    //  _dataPlot->setSubRow(_k, 1, _storedY);
    _k++;

  if (_delay > 0) {
    _bufferY.push_back(_storedY);
  }
}
void siconos::control::LinearSensor::setC(const siconos::algebra::SiconosMatrix& C) {
  *_matC = C;
}

void siconos::control::LinearSensor::setD(const siconos::algebra::SiconosMatrix& D) {
  *_matD = D;
}
