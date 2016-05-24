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
#include "SensorX.hpp"
#include "SensorFactory.hpp"
#include "DynamicalSystem.hpp"
#include "Model.hpp"
#include "TimeDiscretisation.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "LagrangianDS.hpp"

using namespace std;
using namespace SensorFactory;

SensorX::SensorX(): Sensor()
{}

SensorX::SensorX(int name, SP::TimeDiscretisation t): Sensor(name, t)
{}

SensorX::~SensorX()
{
  storedX.reset();
}

void SensorX::initialize()
{
  // Call initialize of base class
  Sensor::initialize();


  //On n'as pas besoin de garder tout les relevés capteurs, on ne garde donc qu'un seul evenement
  //pour y associer notre vecteur de données.

  //Comme on veut récuperer un vecteur a un temps donné, on créer une copie dans un autre vecteur
  storedX.reset(new SiconosVector(model()->nonSmoothDynamicalSystem()->dynamicalSystem(0)->n()));
  (_data[_eSensor])["StoredX"] = storedX;
}

void SensorX::capture()
{
  //capture du vecteur d'état
  storedX = model()->nonSmoothDynamicalSystem()->dynamicalSystem(0)->x();
}

SensorX* SensorX::convert(Sensor* s)
{
  SensorX* sp = dynamic_cast<SensorX*>(s);
  return sp;
}

AUTO_REGISTER_SENSOR(2, SensorX);


