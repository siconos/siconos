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
  int nDof = model()->nonSmoothDynamicalSystem()->dynamicalSystem(0)->dimension();
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


