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

/*! \file Rope.h

  Rope class
*/
#pragma once

#include <vector>

#include "MechanicalProperties.h"
#include "Point.h"

namespace siconos::fem::cable {

class Pylon;

/** Describes a  discrete piece of cable, between two pylons.
 *
 */
class Rope {
 public:
  /** Build a rope span (a piece of cable between two pylons)

      \param a_pile1 first pylon supporting the rope
      \param a_pile2 second pylon
      \param a_tol tolerance used to perform (stop) Newton-Raphson iterations
      \param n_max max. number of iterations used in Newton-Raphson
  */
  Rope(const Pylon &a_pile0, const Pylon &a_pile1, double a_tol, int n_max);

  Rope(Rope &&) = default;

  Rope(const Rope &) = default;

  virtual ~Rope() noexcept = default;

  /** Compute initial profile of the Rope using Catenary equations

      Update positions, internal forces and tension in the rope.

      \param[in] a_meca mechanical properties of the rope span
      \param[in] nb_nodes number of nodes uses to compute profile with catenary equations
      \param[in] a_T0 initial tension applied at the beginning of the rope
      \para[in] m a_R0 support reaction at the beginning of the rope
  */
  void compute(const class MechanicalProperties &a_meca, int nb_nodes, double a_T0,
               const Point &a_R0);

  int computeNbNodes(int nb_elem, double L);

  /** Computes the profile of the rope
   * \param[out] a_q nodes coordinates
   * \param[out] a_R reaction forces at nodes
   * \param[out] a_TS tension at nodes
   * \param[in] q_offset offset (index) in a_q vector
   * \param[in] a_reverse true when 'down' ropeway is concerned
   */
  int computeMesh(std::vector<Point> &a_q, std::vector<Point> &a_R, std::vector<double> &a_TS,
                  int q_offset, bool a_reverse = false);

  double get_T0();
  double get_LastT();
  Point get_LastR();
  double get_L();
  const Point &supportReaction() const;
  const MechanicalProperties &get_meca() const;
  const Pylon &left_pylon() const;

  void to_json(ojson &j);

 private:
  Point ropeway_inc;

  /** Coordinates of rope nodes */
  std::vector<Point> nodes_coords_ = {};

  /** Internal forces (at each node) [x,y,z]-> [H,V,B] */
  std::vector<Point> internal_forces_ = {};

  /** Tension (at each node)*/
  std::vector<double> TS_ = {};

  /** Mechanical properties of the rope (EA, rho ...) */
  MechanicalProperties meca;

  /** First pylon supporting the rope */
  const Pylon &left_pylon_;

  /** Second pylon supporting the rope */
  const Pylon &right_pylon_;

  /** Support reaction at contact [H,V,B]*/
  Point support_reaction_;

  /** True if the rope is the last one in the line */
  bool m_last_{false};

  /** level of tolerance used in get_adm_1C to compute ropeway_inc */
  double tol_ = 1e-7;

  /** max number of iterations used in get_adm_1C to compute ropeway_inc */
  int max_iter_ = 50;

  /** number of nodes used to discretize the rope */
  int number_of_nodes_{0};

  Rope() = delete;
  Rope &operator=(const Rope &) = delete;
  Rope &operator=(Rope &&) = delete;

  /**
      Modifies cable_inc such that the cable equation is satisfied (using tolerance tol and max
     number of iteration)

      \param a_meca object which handles the geometric and material description of the cable
      \param bc a vector of two Pylon objects
      \return [lenght, slope_y, slope_z]
  */
  Point get_adm_1C(const MechanicalProperties &a_meca, const std::vector<Pylon> &bc);

  /**
     \returns an initial guess for the length and the slopes
     \param bc a vector of two Pylons
   */
  Point guess(const std::vector<Pylon> &bc);

  /**
     computes the residual equation and the jacobian of a fixed-fixed cable with imposed
     initial tension

     \param a_meca object which handles the geometric and material description of the cable
     \param bc a vector of two Pylons (end points of the rope)
     \param cable_inc = [L, etaY, etaZ], with L
     the unstretched length of the cable, etaY the initial y-slope and etaZ the initial z-slope
     \param[out] residu the resulting residual
     \param[out] J the resulting jacobian
   */
  void cable_eq(const MechanicalProperties &a_meca, const std::vector<Pylon> &bc,
                const Point &cable_inc, Point &residu, std::vector<std::vector<double>> &J);

  /** Computes initial profile of the cable

     \param[in] a__meca mechanical properties of the cable
     \param[in] cable_inc [L, etaY, etaZ], L is the unstretched length of the cable,
     etaY the initial y-slope and etaZ the initial z-slope
     \param[in] nb_nodes number of nodes
     \param[in, out] a_q nodes coordinates
     \param[in, out] a_R internal forces at each node
     \param[in, out] a_TS tension at each node
     \param[in] q_offset
     \param[in] a_reverse

  */
  void get_profile_1C(const MechanicalProperties &a_meca, const Point &cable_inc, int nb_nodes,
                      std::vector<Point> &a_q, std::vector<Point> &a_R,
                      std::vector<double> &a_TS, int q_offset = 0, bool a_reverse = false);
};
}  // namespace siconos::fem::cable
