/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2021 INRIA.
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

/*! \file FirstOrderType1R.hpp
\brief non linear relations, with y depending on dynamical systems state and r
on lambda.
 */

#ifndef FirstOrderType1R_H
#define FirstOrderType1R_H

#include "FirstOrderR.hpp"

/** Pointer to function for plug-in for operators related to output and its
 * gradients.*/
typedef void (*Type1Ptr)(unsigned int, double *, unsigned int, double *,
                         unsigned int, double *);

/** FirstOrder Non Linear Relation.
 *
 * Derived from FirstOrderR - See this class for more comments.
 *
 *  Relation for First Order Dynamical Systems, with:
 * \rststar
 * .. math::
 *
 *    y &= h(x,z)\\
 *    r &= g(\lambda,z)
 *
 * \endrststar
 *
 * Operators (and their corresponding plug-in):
- h: saved in Interaction as y (plug-in: output[0])
- \f$ \nabla_x h \f$: jacobianH[0] ( output[1] )
- g: saved in DS as r ( input[0])
- \f$ \nabla_\lambda g \f$: jacobianG[0] ( input[1] )
 *
 */
class FirstOrderType1R : public FirstOrderR {
protected:
  // serialization hooks
  ACCEPT_SERIALIZATION(FirstOrderType1R);

public:
  /** default constructor */
  FirstOrderType1R() : FirstOrderR(RELATION::Type1R){};

  /** build from plugin for \f$h(x,z)\f$ and \f$g(\lambda, z)\f$
   *  \param pluginh the plugin to compute h
   *  \param pluging the plugin to compute g
   */
  FirstOrderType1R(const std::string &pluginh, const std::string &pluging);

  /** build from plugin for \f$h(x,z)\f$,\f$g(\lambda, z) \f$ and their
   * gradients \param pluginh the plugin to compute h \param pluging the plugin
   * to compute g \param pluginJachx the plugin to compute \f$\nabla_x h\f$
   *  \param pluginJacglambda the plugin to compute \f$\nabla_{\lambda} g\f$
   */
  FirstOrderType1R(const std::string &pluginh, const std::string &pluging,
                   const std::string &pluginJachx,
                   const std::string &pluginJacglambda);

  /** destructor */
  virtual ~FirstOrderType1R() noexcept = default;

  /** initialize the relation (check sizes, memory allocation ...)
   * \param inter the interaction that owns this relation
   */
  void initialize(Interaction &inter) override;

  /** check sizes of the relation specific operators.
   * \param inter an Interaction using this relation
   */
  inline void checkSize(Interaction &inter) override {};

  /** to compute the output y = h(t,x,...) of the Relation
      \param time current time value
      \param x coordinates of the dynamical systems involved in the relation
      \param z user defined parameters (optional)
      \param y the resulting vector
  */
  void computeh(double time, const BlockVector &x, BlockVector &z,
                SiconosVector &y);

  /** to compute the nonsmooth input r = g(t,x,...) of the Relation
      \param time current time value
      \param lambda interaction \f$\lambda\f$ vector
      \param z user defined parameters (optional)
      \param r the resulting vector
  */
  void computeg(double time, const SiconosVector &lambda, BlockVector &z,
                BlockVector &r);

  /** to compute \f$ C = \nabla_x h \f$
      \param time current time value
      \param x coordinates of the dynamical systems involved in the relation
      \param z user defined parameters (optional)
      \param[out] C the resulting matrix
  */
  void computeJachx(double time, const BlockVector &x, BlockVector &z,
                    SimpleMatrix &C);

  /** default function to compute \f$\nabla_z h\f$
      \param time current time value
      \param x coordinates of the dynamical systems involved in the relation
      \param z user defined parameters (optional)
      \param[out] F the resulting matrix
  */
  void computeJachz(double time, const BlockVector &x, BlockVector &z,
                    SimpleMatrix &F);

  /** to compute \f$ B = \nabla_{\lambda}g \f$
      \param time current time value
      \param x coordinates of the dynamical systems involved in the relation
      \param lambda interaction \f$\lambda\f$ vector
      \param z user defined parameters (optional)
      \param[out] B the resulting matrix
  */
  void computeJacglambda(double time, const SiconosVector &lambda,
                         BlockVector &z, SimpleMatrix &B);

  /** default function to compute y, using the data from the Interaction and DS
   *  \param time current time (not used)
   *  \param inter Interaction using this Relation
   *  \param level not used
   */
  void computeOutput(double time, Interaction &inter,
                             unsigned int level = 0) override;

  /** default function to compute r, using the data from the Interaction and DS
   *  \param time current time (not used)
   *  \param inter Interaction using this Relation
   *  \param level not used
   */
  void computeInput(double time, Interaction &inter,
                            unsigned int level = 0) override;

  void computeJach(double time, Interaction &inter) override;

  void computeJacg(double time, Interaction &inter) override;

  /** return true if the relation requires the computation of residu
      \return true if residu are required, false otherwise
   */
  bool requireResidu()  override { return true; }

  ACCEPT_STD_VISITORS();
};

TYPEDEF_SPTR(FirstOrderType1R)

#endif
