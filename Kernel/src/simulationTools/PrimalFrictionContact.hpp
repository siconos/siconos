/* Siconos-Kernel, Copyright INRIA 2005-2010.
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
/*! \file
  Primal Fricton-Contact Non-Smooth Problem
*/
#ifndef PrimalFrictionContact_H
#define PrimalFrictionContact_H

#include "LinearOSNS.hpp"
#include "SimpleVector.hpp"
#include "SimpleMatrix.hpp"

/** Pointer to function of the type used for drivers for PrimalFrictionContact problems in Numerics */
typedef int (*PFC_Driver)(PrimalFrictionContactProblem*, double*, double*, SolverOptions*, NumericsOptions*);

/** Formalization and Resolution of a Friction-Contact Problem
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) Dec 15, 2005
 *
 * This class is devoted to the formalization and the resolution of
 * primal friction contact problems defined by :
 * \f{eqnarray*}
 *  M velocity =  q +  reaction \\
 *  localVelocity = H^T velocity + tildeLocalVelocity\\
 *  reaction = H localReaction \\
 * \f}
 * and \f$localVelocity,localReaction\f$ belongs to the Coulomb friction law with unilateral contact.
 *
 * With:
 *    - \f$velocity \in R^{n} \f$  and \f$reaction \in R^{n} \f$ the unknowns,
 *    - \f$M \in R^{n \times n } \f$  and \f$q \in R^{n} \f$
 *    - \f$localVelocity \in R^{m} \f$  and \f$localReaction \in R^{m} \f$ the unknowns,
 *    - \f$tildeLocalVelocity \in R^{m} \f$ is the modified local velocity (\f$ e U_{N,k}\f$)
 *    - \f$M \in R^{n \times n } \f$  and \f$q \in R^{n} \f$
 *    - \f$H \in R^{n \times m } \f$
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
 *  - solving of the PrimalFrictionContact problem: function compute(), used to call solvers from Numerics through \n
 * the frictionContact2D_driver() or frictionContact3D_driver() interface of Numerics.
 *  - post-treatment of data: set values of y/lambda variables of the active UR (ie Interactions) using \n
 *  ouput results from the solver (velocity,reaction); function postCompute().
 *
 */
class PrimalFrictionContact : public LinearOSNS
{
protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(PrimalFrictionContact);


  /** Type (dimension) of the contact problem (2D or 3D) */
  int _contactProblemDim;

  /** size of the local problem to solve */
  unsigned int _sizeLocalOutput;

  /** contains the vector localVelocity of a PrimalFrictionContact system */
  SP::SiconosVector _localVelocity;

  /** contains the vector localReaction of a PrimalFrictionContact system */
  SP::SiconosVector _localReaction;

  /** contains the vector localVelocity of a PrimalFrictionContact system */
  SP::SiconosVector _tildeLocalVelocity;

  /** contains the matrix M of a PrimalFrictionContact system */
  SP::OSNSMatrix H;

  /** contains the vector q of a PrimalFrictionContact system */
  SP::MuStorage _mu;

  /** Pointer to the function used to call the Numerics driver to solve the problem */
  PFC_Driver primalFrictionContact_driver;

public:

  /** xml constructor
   *  \param SP::OneStepNSProblemXML : the XML linked-object
   */
  PrimalFrictionContact(SP::OneStepNSProblemXML);

  /** constructor from data
   *  \param int dim (2D or 3D) of the friction-contact problem
   *  \param string numericsSolvername
   *  \param string id of the problem (optional)
   */
  PrimalFrictionContact(int dimPb,
                        const int newNumericsSolverId =
                          SICONOS_FRICTION_3D_PRIMAL_NSGS_WR ,
                        const std::string& newId = "unamed_primal_friction_contact_problem");

  /** destructor
   */
  ~PrimalFrictionContact() {};

  // GETTERS/SETTERS

  /** get the type of PrimalFrictionContact problem (2D or 3D)
   *  \return an int (2 or 3)
   */
  inline int getPrimalFrictionContactDim() const
  {
    return _contactProblemDim;
  }

  /** get dimension of the problem
   *  \return an unsigned ing
   */
  inline unsigned int getLocalSizeOutput() const
  {
    return _sizeLocalOutput;
  }

  /** set the value of sizeOutput
   *  \param an unsigned int
   */
  inline void setLocalSizeOutput(const unsigned int newVal)
  {
    _sizeLocalOutput = newVal;
  }

  // --- LocalVelocity ---
  /** get the value of localVelocity, the initial state of the DynamicalSystem
   *  \return SimpleVector
   *  \warning: SiconosVector is an abstract class => can not be an lvalue => return SimpleVector
   */
  inline const SimpleVector getLocalVelocity() const
  {
    return *_localVelocity;
  }

  /** get localVelocity
   *  \return pointer on a SimpleVector
   */
  inline SP::SiconosVector localVelocity() const
  {
    return _localVelocity;
  }

  /** set the value of localVelocity to newValue
   *  \param SimpleVector newValue
   */
  void setLocalVelocity(const SiconosVector&);

  /** set localVelocity to pointer newPtr
   *  \param SP::SiconosVector  newPtr
   */
  inline void setLocalVelocityPtr(SP::SiconosVector newPtr)
  {
    _localVelocity = newPtr;
  }

  // --- Reaction ---
  /** get the value of reaction, the initial state of the DynamicalSystem
   *  \return SimpleVector
   *  \warning: SimpleVector is an abstract class => can not be an lvalue => return SimpleVector
   */
  inline const SimpleVector getLocalReaction() const
  {
    return *_localReaction;
  }

  /** get localReaction, the initial state of the DynamicalSystem
   *  \return pointer on a SimpleVector
   */
  inline SP::SiconosVector localReaction() const
  {
    return _localReaction;
  }

  /** set the value of localReaction to newValue
   *  \param SimpleVector newValue
   */
  void setLocalReaction(const SiconosVector&);

  /** set localReaction to pointer newPtr
   *  \param SP::SiconosVector  newPtr
   */
  inline void setLocalReactionPtr(SP::SiconosVector newPtr)
  {
    _localReaction = newPtr;
  }

  // --- TildeLocalVelocity ---

  /** get the value of tildeLocalVelocity, the initial state of the DynamicalSystem
   *  \return SimpleVector
   *  \warning: SiconosVector is an abstract class => can not be an lvalue => return SimpleVector
   */
  inline const SimpleVector getTildeLocalVelocity() const
  {
    return *_tildeLocalVelocity;
  }

  /** get tildeLocalVelocity
   *  \return pointer on a SimpleVector
   */
  inline SP::SiconosVector tildeLocalVelocity() const
  {
    return _tildeLocalVelocity;
  }

  /** set the value of tildeLocalVelocity to newValue
   *  \param SimpleVector newValue
   */
  void setTildeLocalVelocity(const SiconosVector&);

  /** set tildeLocalVelocity to pointer newPtr
   *  \param SP::SiconosVector  newPtr
   */
  inline void setTildeLocalVelocityPtr(SP::SiconosVector newPtr)
  {
    _tildeLocalVelocity = newPtr;
  }

  // --- H ---

  /** get H
   *  \return pointer on a OSNSMatrix
   */
  inline SP::OSNSMatrix h() const
  {
    return H;
  }

  /** set the value of H to newValue
   *  \param  newValue
   */
  void setH(const OSNSMatrix &);

  /** set H to pointer newPtr
   *  \param  SP::OSNSMatrix newPtr
   */
  inline void setHPtr(SP::OSNSMatrix newPtr)
  {
    H = newPtr;
  }

  // --- Mu ---
  /** get the vector mu, list of the friction coefficients
   *  \return a vector of double
   */
  inline const std::vector<double> getMu() const
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
   *  \return SimpleVector
   *  \warning: SimpleVector is an abstract class => can not be an lvalue => return SimpleVector
   */
  inline double getMu(unsigned int i) const
  {
    return (*_mu)[i];
  }

  /** set the driver-function used to solve the problem
      \param a function of prototype Driver
  */
  inline void setNumericsDriver(PFC_Driver newFunction)
  {
    primalFrictionContact_driver = newFunction;
  };

  // --- Others functions ---

  /** initialize the PrimalFrictionContact problem(compute topology ...)
    \param the simulation, owner of this OSNSPB
    */
  void initialize(SP::Simulation);

  /** computes extra diagonal unitaryBlock-matrix that corresponds to UR1 and UR2
   *  Move this to Unitary Relation class?
   *  \param a pointer to UnitaryRelation
   *  \param a pointer to UnitaryRelation
   */
  void computeUnitaryBlock(SP::UnitaryRelation, SP::UnitaryRelation);

  /** computes DSBlock-matrix that corresponds to DS1
   *  Move this to Unitary Relation class?
   *  \param a pointer to DynamicalSystem DS1
   */
  void computeDSBlock(SP::DynamicalSystem);

  /** computes  UnitaryDSBlock-matrix that corresponds to UR1 and DS2
   *  Move this to Unitary Relation class?
   *  \param a pointer to UnitaryRelation UR1
   *  \param a pointer to DynamicalSystems DS2
   */
  void computeUnitaryDSBlock(SP::UnitaryRelation, SP::DynamicalSystem);

  /** To compute a part of the "q" vector of the OSNS
      \param SP::UnitaryRelation, the UR which corresponds to the considered block
       \param unsigned int, the position of the first element of yOut to be set
  */
  void computeqBlock(SP::DynamicalSystem, unsigned int);

  /** compute vector q
   *  \param double : current time
   */
  void computeq(double time);

  /** To compute a part of the "tildeLovalVelocity" vector of the OSNS
      \param SP::UnitaryRelation, the UR which corresponds to the considered block
       \param unsigned int, the position of the first element of yOut to be set
  */
  void computeTildeLocalVelocityBlock(SP::UnitaryRelation, unsigned int);

  /** compute vector tildeLocalVelocity
   *  \param double : current time
   */
  void computeTildeLocalVelocity(double time);

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

  /** print the data to the screen */
  void display() const;

  /** encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param OneStepNSProblem* : the one step problem which must be converted
   * \return a pointer on the problem if it is of the right type, NULL otherwise
   */
  static PrimalFrictionContact* convert(OneStepNSProblem* osnsp);
};

#endif // PrimalFrictionContact_H
