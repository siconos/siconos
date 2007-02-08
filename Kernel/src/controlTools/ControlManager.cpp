/* Siconos-Kernel version 2.0.1, Copyright INRIA 2005-2006.
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

#include "ControlManager.h"
#include "Model.h"
#include "Sensor.h"
#include "SensorFactory.h"

using namespace std;

ControlManager::ControlManager(): model(NULL)
{}

ControlManager::ControlManager(Model* m): model(m)
{}

ControlManager::~ControlManager()
{
  SensorsIterator it;
  for (it = allSensors.begin(); it != allSensors.end(); ++it)
    if ((*it) != NULL) delete(*it);
  allSensors.clear();
}

void ControlManager::initialize()
{
}

void ControlManager::addSensor(const string& type, TimeDiscretisation* t)
{
  if (t->getModelPtr() != model)
    RuntimeException::selfThrow("ControlManager::addSensor(...,timeDiscretisation) failed. The Model linked to the controlManager is different from the one of the TimeDiscretisation.");

  SensorFactory::Registry& regSensor(SensorFactory::Registry::get()) ;
  allSensors.insert(regSensor.instantiate("SensorPosition", t));
}

void ControlManager::display() const
{
  cout << "ControlManager " ;
  if (model != NULL)
    cout << "linked to model named: " << model->getTitle() << "." << endl;
  else
    cout << "not linked to a model." << endl;
}
