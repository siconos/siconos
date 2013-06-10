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
  storedX.reset(new SiconosVector(model()->nonSmoothDynamicalSystem()->dynamicalSystem(0)->getN()));
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


