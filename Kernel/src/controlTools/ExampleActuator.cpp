/* Siconos-Kernel version 1.3.0, Copyright INRIA 2005-2006.
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
#include "ExampleActuator.h"
#include "DynamicalSystem.h"
#include "ActuatorFactory.h"
#include "ioMatrix.h"
using namespace std;
using namespace ActuatorFactory;

ExampleActuator::ExampleActuator(): Actuator()
{}

ExampleActuator::ExampleActuator(const std::string& name, TimeDiscretisation* t): Actuator(name, t)
{}

ExampleActuator::~ExampleActuator()
{}

void ExampleActuator::initialize()
{
  // Call initialize of base class
  Actuator::initialize();
}

void ExampleActuator::actuate()
{
  cout << "Actuator action ... " << endl;
  DSIterator itDS;
  SiconosVector * myZ = new SimpleVector(3);
  (*myZ)(0) = 12;
  (*myZ)(1) = 132;
  (*myZ)(2) = 212;

  for (itDS = allDS.begin(); itDS != allDS.end(); ++itDS)
    (*itDS)->setZPtr(myZ);

}

ExampleActuator* ExampleActuator::convert(Actuator* s)
{
  ExampleActuator* sp = dynamic_cast<ExampleActuator*>(s);
  return sp;
}

AUTO_REGISTER_ACTUATOR("ExampleActuator", ExampleActuator);


