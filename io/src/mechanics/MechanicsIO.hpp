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

#ifndef MechanicsIO_hpp
#define MechanicsIO_hpp

#include <MechanicsFwd.hpp>
#ifdef HAVE_SICONOS_MECHANISMS
#include <MechanismsFwd.hpp>
#endif
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
   * \param nsds current nonsmooth dynamical system
   * \return a SP::SimpleMatrix where the columns are
   *          id, x, y, z, qw, qx, qy, qz
   *   id is the DynamicalSystem number + 1
   */
  SP::SimpleMatrix positions(const NonSmoothDynamicalSystem& nsds) const;

  /** get all velocities: translation (xdot, ydot, zdot) + orientation velocities
   * ox, oy, oz
   * \param nsds current nonsmooth dynamical system
   *   \return a matrix where the columns are id, xdot, ydot, zdot,
   *   ox, oy, oz
   * id is the DynamicalSystem number + 1
   */
  SP::SimpleMatrix velocities(const NonSmoothDynamicalSystem& nsds) const;

  /** get the coordinates of all contact points, normals, reactions and velocities
   * \param nsds current nonsmooth dynamical system
   * \param index_set the index set number.
   \return a matrix where the columns are mu x y z, nx, ny, nz, rx, ry, rz, vx, vy, vz, ox, oy, oz, id
  */
  SP::SimpleMatrix contactPoints(const NonSmoothDynamicalSystem& nsds, unsigned int index_set=1) const;

  /** get the contact information that is the ds linked by the interaction
   * \param nsds current nonsmooth dynamical system
   * \param index_set the index set number.
   \return a matrix where the columns are interaction id, ds1 number, ds2 number, static object number (if possible)
  */

  SP::SimpleMatrix contactInfo(const NonSmoothDynamicalSystem& nsds, unsigned int index_set=1) const;

  /** get the dissipation values  of all contact points
   * \param nsds current nonsmooth dynamical system
   * \param index_set the index set number.
   * \param omega the value of the weigth for the weight in the computaion of the contact work
   *  by default omega =1/2 and the contact work corresponds to the theoretical formula
   *        1/2 (v^+ + v^-)^\top p
   * otherwise it corresponds to v_{k+omega} p 
   * \param tol double for the computation of contact status
   \return a matrix where the columns are id, normal contact work, tangent contact work, friction dissipation, contact status
  */
  SP::SimpleMatrix contactContactWork(const NonSmoothDynamicalSystem& nsds,
				      unsigned int index_set=1,
				      double omega = 0.5) const;

  /** get the domain of each contact point
   * \param nsds current nonsmooth dynamical system
   * \return a matrix where the columns are domain, id
  */
  SP::SimpleMatrix domains(const NonSmoothDynamicalSystem& nsds) const;
};


#endif
