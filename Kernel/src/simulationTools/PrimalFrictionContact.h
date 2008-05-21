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
  Primal Fricton-Contact Non-Smooth Problem
*/
#ifndef PrimalFrictionContact_H
#define PrimalFrictionContact_H

#include "OneStepNSProblem.h"
#include "SimpleVector.h"
#include "SimpleMatrix.h"

// Pointer to function of the type used for drivers for PrimalFrictionContact problems in Numerics
typedef int (*PFC_Driver)(PrimalFrictionContact_Problem*, double*, double*, Solver_Options*, Numerics_Options*);

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
 *   - XML reading (inputs = xml node with tag "OneStepNSProblem" and a Simulation*)
 *   - Constructor from data (inputs = Simulations*, id, NonSmoothSolver*) - The solver is optional.
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
class PrimalFrictionContact : public OneStepNSProblem
{
protected:

  /** Type (dimension) of the contact problem (2D or 3D) */
  int contactProblemDim;

  /** size of the local problem to solve */
  unsigned int sizeLocalOutput;

  /** contains the vector velocity of a PrimalFrictionContact system */
  SimpleVector *velocity;

  /** contains the vector reaction of a PrimalFrictionContact system */
  SimpleVector *reaction;

  /** contains the vector localVelocity of a PrimalFrictionContact system */
  SimpleVector *localVelocity;

  /** contains the vector localReaction of a PrimalFrictionContact system */
  SimpleVector *localReaction;

  /** contains the vector localVelocity of a PrimalFrictionContact system */
  SimpleVector *tildeLocalVelocity;

  /** contains the matrix M of a PrimalFrictionContact system */
  OSNSMatrix *M;

  /** contains the matrix M of a PrimalFrictionContact system */
  OSNSMatrix *H;

  /** contains the vector q of a PrimalFrictionContact system */
  SimpleVector *q;

  /** contains the vector q of a PrimalFrictionContact system */
  std::vector<double>* mu;

  /** Flags to check wheter pointers were allocated in class constructors or not */
  bool isVelocityAllocatedIn;
  bool isReactionAllocatedIn;
  bool isLocalVelocityAllocatedIn;
  bool isLocalReactionAllocatedIn;
  bool isTildeLocalVelocityAllocatedIn;
  bool isMAllocatedIn;
  bool isHAllocatedIn;
  bool isQAllocatedIn;

  /** Storage type for M - 0: SiconosMatrix (dense), 1: Sparse Storage (embedded into OSNSMatrix) */
  int MStorageType;

  /** Pointer to the function used to call the Numerics driver to solve the problem */
  PFC_Driver primalFrictionContact_driver;

private:

  /** default constructor (private)
   */
  PrimalFrictionContact(const std::string& pbType = "PrimalFrictionContact2D");

public:

  /** xml constructor
   *  \param OneStepNSProblemXML* : the XML linked-object
   *  \param Simulation * the simulation that owns the problem
   */
  PrimalFrictionContact(OneStepNSProblemXML*, Simulation*);

  /** constructor from data
   *  \param Simulation* the simulation that owns this problem
   *  \param int dim (2D or 3D) of the friction-contact problem
   *  \param Solver* pointer to object that contains solver algorithm and formulation \n
   *  (optional, default = NULL => read .opt file in Numerics)
   *  \param string id of the problem (optional)
   */
  PrimalFrictionContact(Simulation *, int, NonSmoothSolver* = NULL, const std::string& = "unamed_friction_contact_problem");

  /** destructor
   */
  ~PrimalFrictionContact();

  // GETTERS/SETTERS

  /** get the type of PrimalFrictionContact problem (2D or 3D)
   *  \return an int (2 or 3)
   */
  inline int getPrimalFrictionContactDim() const
  {
    return contactProblemDim;
  }

  /** get dimension of the problem
   *  \return an unsigned ing
   */
  inline const unsigned int getLocalSizeOutput() const
  {
    return sizeLocalOutput;
  }

  /** set the value of sizeOutput
   *  \param an unsigned int
   */
  inline void setLocalSizeOutput(const unsigned int newVal)
  {
    sizeLocalOutput = newVal;
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
  inline SimpleVector* getVelocityPtr() const
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
  void setVelocityPtr(SimpleVector*);

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
  inline SimpleVector* getReactionPtr() const
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
  void setReactionPtr(SimpleVector*) ;


  // --- LocalVelocity ---
  /** get the value of localVelocity, the initial state of the DynamicalSystem
   *  \return SimpleVector
   *  \warning: SiconosVector is an abstract class => can not be an lvalue => return SimpleVector
   */
  inline const SimpleVector getLocalVelocity() const
  {
    return *localVelocity;
  }

  /** get localVelocity
   *  \return pointer on a SimpleVector
   */
  inline SimpleVector* getLocalVelocityPtr() const
  {
    return localVelocity;
  }

  /** set the value of localVelocity to newValue
   *  \param SimpleVector newValue
   */
  void setLocalVelocity(const SimpleVector&);

  /** set localVelocity to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  void setLocalVelocityPtr(SimpleVector*);

  // --- Reaction ---
  /** get the value of reaction, the initial state of the DynamicalSystem
   *  \return SimpleVector
   *  \warning: SimpleVector is an abstract class => can not be an lvalue => return SimpleVector
   */
  inline const SimpleVector getLocalReaction() const
  {
    return *localReaction;
  }

  /** get localReaction, the initial state of the DynamicalSystem
   *  \return pointer on a SimpleVector
   */
  inline SimpleVector* getLocalReactionPtr() const
  {
    return localReaction;
  }

  /** set the value of localReaction to newValue
   *  \param SimpleVector newValue
   */
  void setLocalReaction(const SimpleVector&);

  /** set localReaction to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  void setLocalReactionPtr(SimpleVector*) ;
  // --- TildeLocalVelocity ---

  /** get the value of tildeLocalVelocity, the initial state of the DynamicalSystem
   *  \return SimpleVector
   *  \warning: SiconosVector is an abstract class => can not be an lvalue => return SimpleVector
   */
  inline const SimpleVector getTildeLocalVelocity() const
  {
    return *tildeLocalVelocity;
  }

  /** get tildeLocalVelocity
   *  \return pointer on a SimpleVector
   */
  inline SimpleVector* getTildeLocalVelocityPtr() const
  {
    return tildeLocalVelocity;
  }

  /** set the value of tildeLocalVelocity to newValue
   *  \param SimpleVector newValue
   */
  void setTildeLocalVelocity(const SimpleVector&);

  /** set tildeLocalVelocity to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  void setTildeLocalVelocityPtr(SimpleVector*);


  // --- M ---

  /** get M
   *  \return pointer on a OSNSMatrix
   */
  inline OSNSMatrix* getMPtr() const
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
  void setMPtr(OSNSMatrix*);


  // --- H ---

  /** get H
   *  \return pointer on a OSNSMatrix
   */
  inline OSNSMatrix* getHPtr() const
  {
    return H;
  }

  /** set the value of H to newValue
   *  \param  newValue
   */
  void setH(const OSNSMatrix &);

  /** set H to pointer newPtr
   *  \param  OSNSMatrix* newPtr
   */
  void setHPtr(OSNSMatrix*);

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
  inline SimpleVector* getQPtr() const
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
  void setQPtr(SimpleVector*);

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
  inline std::vector<double>* getMuPtr() const
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
  inline void setNumericsDriver(PFC_Driver newFunction)
  {
    primalFrictionContact_driver = newFunction;
  };

  // --- Others functions ---

  /** initialize the PrimalFrictionContact problem(compute topology ...)
   */
  void initialize();

  /** computes extra diagonal unitaryBlock-matrix that corresponds to UR1 and UR2
   *  Move this to Unitary Relation class?
   *  \param a pointer to UnitaryRelation
   *  \param a pointer to UnitaryRelation
   */
  void computeUnitaryBlock(UnitaryRelation*, UnitaryRelation*);

  /** computes DSBlock-matrix that corresponds to DS1
   *  Move this to Unitary Relation class?
   *  \param a pointer to DynamicalSystem DS1
   */
  void computeDSBlock(DynamicalSystem*);

  /** computes  UnitaryDSBlock-matrix that corresponds to UR1 and DS2
   *  Move this to Unitary Relation class?
   *  \param a pointer to UnitaryRelation UR1
   *  \param a pointer to DynamicalSystems DS2
   */
  void computeUnitaryDSBlock(UnitaryRelation*, DynamicalSystem*);

  /** To compute a part of the "q" vector of the OSNS
      \param UnitaryRelation*, the UR which corresponds to the considered block
       \param unsigned int, the position of the first element of yOut to be set
  */
  void computeQBlock(DynamicalSystem*, unsigned int);

  /** compute vector q
   *  \param double : current time
   */
  void computeQ(double time);

  /** To compute a part of the "tildeLovalVelocity" vector of the OSNS
      \param UnitaryRelation*, the UR which corresponds to the considered block
       \param unsigned int, the position of the first element of yOut to be set
  */
  void computeTildeLocalVelocityBlock(UnitaryRelation*, unsigned int);

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

  /** copy the data of the OneStepNSProblem to the XML tree */
  void saveNSProblemToXML();

  /** print the data to the screen */
  void display() const;

  /** encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param OneStepNSProblem* : the one step problem which must be converted
   * \return a pointer on the problem if it is of the right type, NULL otherwise
   */
  static PrimalFrictionContact* convert(OneStepNSProblem* osnsp);
};

#endif // PrimalFrictionContact_H
