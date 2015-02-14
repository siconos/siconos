/* Siconos-Kernel, Copyright INRIA 2005-2015
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

#include <iostream>

#include <FirstOrderLinearTIDS.hpp>
#include <TimeStepping.hpp>
#include <SiconosVector.hpp>
#include <SimpleMatrix.hpp>
#include <NormalConeNSL.hpp>
#include <AVI.hpp>

#include "Twisting.hpp"

#include "ControlSensor.hpp"
#include "ActuatorFactory.hpp"

Twisting::Twisting(SP::ControlSensor sensor, double gain, double beta, double hControl):
  CommonSMC(TWISTING, sensor)
{
  _u.reset(new SiconosVector(2));
  if (beta <= 0.0 || beta >= 1.0)
  {
    std::cout << "Twisting constructor: beta is not in (0, 1)" << std::endl;
  }

  _B.reset(new SimpleMatrix(2, 2));
  (*_B)(1, 0) = gain;
  (*_B)(1, 1) = gain*beta;

  setNSdata(hControl);
}

Twisting::Twisting(SP::ControlSensor sensor, double hControl):
  CommonSMC(TWISTING, sensor)
{
  _u.reset(new SiconosVector(2));
  setNSdata(hControl);
}

void Twisting::setNSdata(double hControl)
{
  SP::SimpleMatrix H(new SimpleMatrix(4, 2));
  (*H)(0, 0) = 1.0;
  (*H)(1, 0) = -hControl/2.0;
  (*H)(2, 0) = -1.0;
  (*H)(3, 0) = hControl/2.0;
  (*H)(1, 1) = 1.0;
  (*H)(3, 1) = -1.0;

  SP::SiconosVector K(new SiconosVector(4));
  (*K)(0) = -1.0;
  (*K)(1) = -1.0;
  (*K)(2) = -1.0;
  (*K)(3) = -1.0;

  _nsLawSMC.reset(new NormalConeNSL(2, H, K));
  _OSNSPB_SMC.reset(new AVI());
  _numericsSolverId = SICONOS_AVI_CAOFERRIS;
}

void Twisting::initialize(const Model& m)
{
  // basic check
  if (!_nsLawSMC || !_OSNSPB_SMC)
  {
    RuntimeException::selfThrow("Twisting::initialize - nslaw or osnsp not set. If you used the constructor with only the ControlSensor as argument, you need to manually call setNSdata");
  }
  CommonSMC::initialize(m);
}


Twisting::~Twisting()
{
}

void Twisting::actuate()
{

  *(_DS_SMC->x()) = _sensor->y();

  _simulationSMC->computeOneStep();
  _simulationSMC->nextStep();

  // discontinous part
  *_us = *_lambda;
  *_u = *_us;
  _indx++;

}

AUTO_REGISTER_ACTUATOR(TWISTING, Twisting)
