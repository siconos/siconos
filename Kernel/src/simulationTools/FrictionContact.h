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

#include "OneStepNSProblem.h"
#include "SimpleVector.h"
#include "SimpleMatrix.h"

// Pointer to function of the type used for drivers for FrictionContact problems in Numerics
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
class FrictionContact : public OneStepNSProblem
{
protected:

  /** Type (dimension) of the contact problem (2D or 3D) */
  int contactProblemDim;

  /** contains the vector velocity of a FrictionContact system */
  SP::SimpleVector velocity;

  /** contains the vector reaction of a FrictionContact system */
  SP::SimpleVector reaction;

  /** contains the matrix M of a FrictionContact system */
  SP::OSNSMatrix M;

  /** contains the vector q of a FrictionContact system */
  SP::SimpleVector q;

  boost::shared_ptr< std::vector<double> > mu;

  /** Storage type for M - 0: SiconosMatrix (dense), 1: Sparse Storage (embedded into OSNSMatrix) */
  int MStorageType;

  /** Pointer to the function used to call the Numerics driver to solve the problem */
  Driver frictionContact_driver;

private:

  /** default constructor (private)
   */
  FrictionContact(const std::string& pbType = "FrictionContact2D");

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
  FrictionContact(int, SP::NonSmoothSolver = SP::NonSmoothSolver(), const std::string& = "unamed_friction_contact_problem");

  /** destructor
   */
  ~FrictionContact() {};

  // GETTERS/SETTERS

  /** get the type of FrictionContact problem (2D or 3D)
   *  \return an int (2 or 3)
   */
  inline int getFrictionContactDim() const
  {
    return contactProblemDim;
  }

  // --- Velocity ---
  /** get the value of velocity, the initial state of the DynamicalSystem
   *  \return SimpleVector
   *  \warning: SiconosVector is an abstract class => can not be an lvalue => return SimpleVector
   */
  inline const SimpleVector getVelocity() const
  {
    return *velocity;
  }

  /** get velocity, the initial state of the DynamicalSystem
   *  \return pointer on a SimpleVector
   */
  inline SP::SimpleVector getVelocityPtr() const
  {
    return velocity;
  }

  /** set the value of velocity to newValue
   *  \param SimpleVector newValue
   */
  void setVelocity(const SimpleVector&);

  /** set velocity to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  void setVelocityPtr(SP::SimpleVector);

  // --- Reaction ---
  /** get the value of reaction, the initial state of the DynamicalSystem
   *  \return SimpleVector
   *  \warning: SimpleVector is an abstract class => can not be an lvalue => return SimpleVector
   */
  inline const SimpleVector getReaction() const
  {
    return *reaction;
  }

  /** get reaction, the initial state of the DynamicalSystem
   *  \return pointer on a SimpleVector
   */
  inline SP::SimpleVector getReactionPtr() const
  {
    return reaction;
  }

  /** set the value of reaction to newValue
   *  \param SimpleVector newValue
   */
  void setReaction(const SimpleVector&);

  /** set reaction to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  void setReactionPtr(SP::SimpleVector) ;

  // --- M ---

  /** get M
   *  \return pointer on a OSNSMatrix
   */
  inline SP::OSNSMatrix getMPtr() const
  {
    return M;
  }

  /** set the value of M to newValue
   *  \param  newValue
   */
  void setM(const OSNSMatrix &);

  /** set M to pointer newPtr
   *  \param  OSNSMatrix* newPtr
   */
  void setMPtr(SP::OSNSMatrix);

  // --- Q ---
  /** get the value of q, the initial state of the DynamicalSystem
   *  \return SimpleVector
   *  \warning: SimpleVector is an abstract class => can not be an lvalue => return SimpleVector
   */
  inline const SimpleVector getQ() const
  {
    return *q;
  }

  /** get q, the initial state of the DynamicalSystem
   *  \return pointer on a SimpleVector
   */
  inline SP::SimpleVector getQPtr() const
  {
    return q;
  }

  /** set the value of q to newValue
   *  \param SimpleVector newValue
   */
  void setQ(const SimpleVector&);

  /** set q to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  void setQPtr(SP::SimpleVector);

  // --- Mu ---
  /** get the vector mu, list of the friction coefficients
   *  \return a vector of double
   */
  inline const std::vector<double> getMu() const
  {
    return *mu;
  }

  /** get a pointer to mu, the list of the friction coefficients
   *  \return pointer on a std::vector<double>
   */

  inline boost::shared_ptr< std::vector<double> > getMuPtr() const
  {
    return mu;
  }

  /** get the value of the component number i of mu, the vector of the friction coefficients
   *  \return SimpleVector
   *  \warning: SimpleVector is an abstract class => can not be an lvalue => return SimpleVector
   */
  inline const double getMu(unsigned int i) const
  {
    return (*mu)[i];
  }

  /** get the type of storage for M */
  inline const int getMStorageType() const
  {
    return MStorageType;
  };

  /** set which type of storage will be used for M - Note that this function does not
      allocate any memory for M, it just sets an indicator for future use */
  inline void setMStorageType(int i)
  {
    MStorageType = i;
  };

  /** set the driver-function used to solve the problem
      \param a function of prototype Driver
  */
  inline void setNumericsDriver(Driver newFunction)
  {
    frictionContact_driver = newFunction;
  };

  // --- Others functions ---

  /** initialize the FrictionContact problem(compute topology ...)
      \param the simulation, owner of this OSNSPB
   */
  void initialize(SP::Simulation);

  /** computes extra diagonal unitaryBlock-matrix that corresponds to UR1 and UR2
   *  Move this to Unitary Relation class?
   *  \param a pointer to UnitaryRelation
   *  \param a pointer to UnitaryRelation
   */
  void computeUnitaryBlock(SP::UnitaryRelation, SP::UnitaryRelation);

  /** To compute a part of the "q" vector of the OSNS
      \param SP::UnitaryRelation, the UR which corresponds to the considered block
       \param unsigned int, the position of the first element of yOut to be set
  */
  void computeQBlock(SP::UnitaryRelation, unsigned int);

  /** compute vector q
   *  \param double : current time
   */
  void computeQ(double time);

  /** pre-treatment for LCP
   *  \param double : current time
   */
  void preCompute(double time);

  /** Compute the unknown reaction and velocity and update the Interaction (y and lambda )
   *  \param double current time
   *  \return int information about the solver convergence (0: ok, >0 problem, see Numerics documentation)
   */
  int compute(double time);

  /** post-treatment of output from Numerics solver: \n
   *  set values of the unknowns of Interactions using (velocity,reaction)
   */
  void postCompute();

  /** copy the data of the OneStepNSProblem to the XML tree */
  void saveNSProblemToXML();

  /** print the data to the screen */
  void display() const;

  /** encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param OneStepNSProblem* : the one step problem which must be converted
   * \return a pointer on the problem if it is of the right type, NULL otherwise
   */
  static FrictionContact* convert(OneStepNSProblem* osnsp);
};

#endif // FrictionContact_H
