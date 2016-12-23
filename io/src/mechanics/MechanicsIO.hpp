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

#ifndef MechanicsIO_hpp
#define MechanicsIO_hpp

#include <MechanicsFwd.hpp>
#include <SiconosPointers.hpp>
#include <SiconosFwd.hpp>

class MechanicsIO
{
protected:

  template<typename T, typename G>
  SP::SimpleMatrix visitAllVerticesForVector(const G& graph) const;

  template<typename T, typename G>
  SP::SiconosVector visitAllVerticesForDouble(const G& graph) const;

public:
  /** default constructor
   */
  MechanicsIO() {};

  /** get all positions: translation (x,y,z) + orientation quaternion
   * (qw, qx, qy, qz)
   * \param model the model
   * \return a SP::SimpleMatrix where the columns are
             id, x, y, z, qw, qx, qy, qz
     id is the DynamicalSystem number + 1
   */
  SP::SimpleMatrix positions(const Model& model) const;

  /** get all velocities: translation (xdot, ydot, zdot) + orientation velocities
      ox, oy, oz
   * \param model the model
      \return a matrix where the columns are id, xdot, ydot, zdot,
      ox, oy, oz
      id is the DynamicalSystem number + 1
  */
  SP::SimpleMatrix velocities(const Model& model) const;

  /** get the coordinates of all contact points, normals, reactions and velocities
   * \param model the model
   * \param index_set the index set number.
      \return a matrix where the columns are mu x y z, nx, ny, nz, rx, ry, rz, vx, vy, vz, ox, oy, oz, id
  */
  SP::SimpleMatrix contactPoints(const Model& model, unsigned int index_set=1) const;

  /** get the domain of each contact point
   * \param model the model
      \return a matrix where the columns are domain, id
  */
  SP::SimpleMatrix domains(const Model& model) const;
};


#endif
