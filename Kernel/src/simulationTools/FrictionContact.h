/* Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
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
/*! \file
  Fricton-Contact Non-Smooth Problem
*/
#ifndef FrictionContact_H
#define FrictionContact_H

#include "LinearOSNS.h"

/** Pointer to function of the type used for drivers for FrictionContact problems in Numerics */
typedef int (*Driver)(FrictionContact_Problem*, double*, double*, Solver_Options*, Numerics_Options*);

/** Formalization and Resolution of a Friction-Contact Problem
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) Dec 15, 2005
 *
 * This class is devoted to the formalization and the resolution of
 * friction contact problems defined by :
 * \f{eqnarray*}
 * velocity =  q + M reaction \\
 * \\
 * velocity \geq 0, reaction \geq 0,  reaction^{T} velocity =0
 * \f}
 * and a Coulomb friction law.
 *
 * With:
 *    - \f$velocity \in R^{n} \f$  and \f$reaction \in R^{n} \f$ the unknowns,
 *    - \f$M \in R^{n \times n } \f$  and \f$q \in R^{n} \f$
 *
 * The dimension of the problem (2D or 3D) is given by the variable contactProblemDim and the right
 * Numerics driver will be called according to this value.
 *
 * \b Construction:
 *   - XML reading (inputs = xml node with tag "OneStepNSProblem" and a SP::Simulation)
 *   - Constructor from data (inputs = Simulations*, id, SP::NonSmoothSolver) - The solver is optional.
 * Main functions:
 *
 * \b Main functions:
 *  - formalization of the problem: computes M,q using the set of "active" UnitaryRelations from the simulation and \n
 *  the unitaryBlock-matrices saved in the field unitaryBlocks.\n
 *  Functions: initialize(), computeUnitaryBlock(), preCompute()
 *  - solving of the FrictionContact problem: function compute(), used to call solvers from Numerics through \n
 * the frictionContact2D_driver() or frictionContact3D_driver() interface of Numerics.
 *  - post-treatment of data: set values of y/lambda variables of the active UR (ie Interactions) using \n
 *  ouput results from the solver (velocity,reaction); function postCompute().
 *
 */
class FrictionContact : public LinearOSNS
{
protected:

  /** Type (dimension) of the contact problem (2D or 3D) */
  int _contactProblemDim;

  SP::MuStorage _mu;

  /** Pointer to the function used to call the Numerics driver to solve the problem */
  Driver _frictionContact_driver;

public:

  /** xml constructor
   *  \param SP::OneStepNSProblemXML : the XML linked-object
   */
  FrictionContact(SP::OneStepNSProblemXML);

  /** constructor from data
   *  \param int dim (2D or 3D) of the friction-contact problem
   *  \param Solver* pointer to object that contains solver algorithm and formulation \n
   *  (optional, default = NULL => read .opt file in Numerics)
   *  \param string id of the problem (optional)
   */
  FrictionContact(int dimPb, SP::NonSmoothSolver newSolver = SP::NonSmoothSolver(),
                  const std::string& newId = "unamed_friction_contact_problem"):
    LinearOSNS("FrictionContact", newSolver, newId), _contactProblemDim(dimPb) {};

  /** destructor
   */
  ~FrictionContact() {};

  // GETTERS/SETTERS

  /** get the type of FrictionContact problem (2D or 3D)
   *  \return an int (2 or 3)
   */
  inline int getFrictionContactDim() const
  {
    return _contactProblemDim;
  }

  // --- Mu ---
  /** get the vector mu, list of the friction coefficients
   *  \return a vector of double
   */
  inline const MuStorage getMu() const
  {
    return *_mu;
  }

  /** get a pointer to mu, the list of the friction coefficients
   *  \return pointer on a std::vector<double>
   */

  inline SP::MuStorage mu() const
  {
    return _mu;
  }

  /** get the value of the component number i of mu, the vector of the friction coefficients
   *  \return double value of mu
   */
  inline const double getMu(unsigned int i) const
  {
    return (*_mu)[i];
  }

  /** set the driver-function used to solve the problem
      \param a function of prototype Driver
  */
  inline void setNumericsDriver(Driver newFunction)
  {
    _frictionContact_driver = newFunction;
  };

  // --- Others functions ---

  /** initialize the FrictionContact problem(compute topology ...)
      \param the simulation, owner of this OSNSPB
   */
  void initialize(SP::Simulation);

  /** Compute the unknown reaction and velocity and update the Interaction (y and lambda )
   *  \param double current time
   *  \return int information about the solver convergence (0: ok, >0 problem, see Numerics documentation)
   */
  int compute(double time);

  /** print the data to the screen */
  void display() const;

  /** encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param OneStepNSProblem* : the one step problem which must be converted
   * \return a pointer on the problem if it is of the right type, NULL otherwise
   */
  static FrictionContact* convert(OneStepNSProblem* osnsp);
};

TYPEDEF_SPTR(FrictionContact);

#endif // FrictionContact_H
