/* Siconos-Kernel version 1.3.0, Copyright INRIA 2005-2011.
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
#include "DelayedSampledActuator.hpp"
#include "DynamicalSystem.hpp"
#include "ActuatorFactory.hpp"

#include "Model.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "SensorX.hpp"

extern "C" {
#include "ControlLaw.h"
}

using namespace std;
using namespace ActuatorFactory;


DelayedSampledActuator::DelayedSampledActuator(): Actuator()
{}

DelayedSampledActuator::DelayedSampledActuator(int type, SP::TimeDiscretisation t): Actuator(type, t)
{}

DelayedSampledActuator::~DelayedSampledActuator()
{}

void DelayedSampledActuator::initialize()
{
  // Call initialize of base class
  Actuator::initialize();

  //ajout du system dynamique dans les attributs de notre objet
  this->addDynamicalSystemPtr(model()->nonSmoothDynamicalSystem()->dynamicalSystem(0));
}

void DelayedSampledActuator::actuate()
{
  DSIterator itDS;
  SensorsIterator itS;
  int nDof = model()->nonSmoothDynamicalSystem()->dynamicalSystem(0)->getDim();
  int ncont = 0;
  double t = model()->currentTime();
  SiconosVector * myZ(new SiconosVector(nDof));
  SiconosVector * state(new SiconosVector(2 * nDof));

  // ici récupération de la valeur du Sensor
  itS = getSensors()->begin();
  //  EventsContainer capteurEvents = (*itS)->getEvents();
  //  DataSet data;
  SP::Event event  = (*itS)->event();
  //  Event * event = *(capteurEvents.begin());
  DataSet * data(new DataSet((*itS)->getData()));
  *state = *(((*data)[event])["StoredX"]);
  //for(unsigned int i=0; i<nDof; i++)
  //cout << (*state)(i) <<" ";
  //cout << endl;

  // et là calcul des couples correspondant
  controlLaw(&t, &((*state)(0)), &((*state)(0)) + nDof, &nDof, &ncont, &((*myZ)(0)));

  //mises à jour de la variable gardant les couples pour le systeme dynamique: Z
  itDS = dynamicalSystems()->begin();
  (*itDS)->setz(*myZ);

  delete myZ;
  delete state;
}

DelayedSampledActuator* DelayedSampledActuator::convert(Actuator* s)
{
  DelayedSampledActuator* sp = dynamic_cast<DelayedSampledActuator*>(s);
  return sp;
}

AUTO_REGISTER_ACTUATOR(2, DelayedSampledActuator);


