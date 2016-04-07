/* Siconos-IO, Copyright INRIA 2005-2011.
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
   * \return a SP::SimpleMatrix where the columns are x, y, z, qw, qx, qy, qz 
   */
  SP::SimpleMatrix positions(const Model& model) const;

  /** get all velocities: translation (xdot, ydot, zdot) + orientation velocities 
      ox, oy, oz
   * \param model the model
      \return a matrix where the columns are xdot, ydot, zdot,
      ox, oy, oz
  */
  SP::SimpleMatrix velocities(const Model& model) const;

  /** get the coordinates of all contact points, normals, reactions and velocities
   * \param model the model
      \return a matrix where the columns are mu x y z, nx, ny, nz, rx, ry, rz, vx, vy, vz, ox, oy, oz, id
  */
  SP::SimpleMatrix contactPoints(const Model& model) const;
};


#endif
