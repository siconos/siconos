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
/*! \file OSNSMultipleImpact.hpp
 * \brief Linear Complementarity Problem formulation and solving
 */

#ifndef _OSNSMULTIPLEIMPACT_
#define _OSNSMULTIPLEIMPACT_

#include "LinearOSNS.hpp"
#include <string>

using namespace RELATION;

#define DEFAULT__tolImpact MACHINE_PREC
#define DEFAULT_TOL_VEL MACHINE_PREC
#define DEFAULT_TOL_ENER MACHINE_PREC

/** Formalization and Resolution of a Multiple Impact Non-Smooth problem.

\todo write a short introduction about OSNSMultipleImpact ...
 */
class OSNSMultipleImpact : public LinearOSNS
{
private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(OSNSMultipleImpact);

  //! Time-like variable (Impulse)
  double _impulseVariable;
  //! Time variable
  double _timeVariable;
  //! Number of contacts (only the active contacts)
  unsigned int _nContact;
  //! Maximal number of steps for each computation
  unsigned int _nStepMax;
  //! Tolerance to define zero
  double _tolImpact;
  //! Type of the compliance model
  std::string _typeCompLaw;
  //Velocity of bodies during impact
  //SP::SiconosVector VelAllBody;
  // Relative velocity at all Interactions (with or without contact)
  //SP::SiconosVector VelAllIteractions;
  //! Relative velocity during impact (at the end of each calculation step)
  SP::SiconosVector _velocityContact;
  //! Relative velocity during impact (at the beginning of each calculation step)
  SP::SiconosVector _oldVelocityContact;
  //! Potential energy during impact (at the end of each calculation step)
  SP::SiconosVector _energyContact;
  //! Work done during the last compression phase at contact
  SP::SiconosVector _WorkcContact;
  //! Distribution vector to distribute the incremental impulse at contact
  SP::SiconosVector _distributionVector;
  /** State of contacts at the beginning of impact
   if *_stateContact[i] = 0 => no impact at this contact (at contact with positive relative velocity and no potential energy, may be the impact has been terminated at this contact)
   if *_stateContact[i] = 1 => impact takes place at this contact without potential energy (beginning of impact or repeating impact)
   if *_stateContact[i] = 2 => impact takes place with not-zero potential energy */
  SP::IndexInt _stateContact;
  //!Stiffness at contacts
  SP::SiconosVector  _Kcontact;
  //! Restitution coefficient of contacts
  SP::SiconosVector _restitutionContact;
  //! Elasticity coefficient of contacts
  SP::SiconosVector _elasticyCoefficientcontact;
  //! Incremental impulse at contacts
  SP::SiconosVector _deltaImpulseContact;
  //! Total impulse at contacts
  SP::SiconosVector _tolImpulseContact;
  //! Impulse at contacts for each update time
  SP::SiconosVector _impulseContactUpdate;
  //! Force at contacts
  SP::SiconosVector _forceContact;
  /** Flag to select the primary contact based on the relative velocity or on the potential energy
   at contacts
   if _selectPrimaConInVel = true => select the primary contact according to the relative velocity
   at contact
   if _selectPrimaConInVel = false => select the primary contact according to the potential energy
   at contact */
  bool _selectPrimaConInVel;
  //! ID of the primary contact
  unsigned int _primaryContactId;
  /** Indicator about the selection of the primary contact
      true if primary contact is selected according to the potential energy
      false if primary contact is selected according to the relative velocity */
  bool _isPrimaryContactEnergy;
  //! Relative velocity at primary contact
  double _relativeVelocityPrimaryContact;
  //! Potential energy at primary contact
  double _energyPrimaryContact;
  //! Step size for the iterative calculation
  double _deltaP;
  //! Name of file into which the datat is writen
  std::string  _namefile;
  /** YesWriteData = true ==>save the data during impact
      YesWriteData = false ==> not save the data during impact */
  bool _saveData;
  /** bool variable to set the step size for multiple impact computation
      If IsNumberOfStepsEst = true ==> estimate the step size from the state of the dynamic system before impact and the number of step needed
      Number of steps after which the data is saved */
  unsigned int _nStepSave; //! If IsNumberOfStepsEst = false ==> user choose the step size
  //! Matrix on which the data during impact is saved
  SP::SiconosMatrix _DataMatrix;
  //! Number of points to be save during impacts
  unsigned int _sizeDataSave;
  /** indicator on the termination of the multiple impact process
      _IsImpactEnd = true: impact is terminated
      _IsImpactEnd = false: otherwise */
  bool _IsImpactEnd;
  //! Tolerance to define a negligeble value for a velocity grandeur
  double _Tol_Vel;
  //! Tolerance to define a negligeable value for a potential energy grandeur
  double _Tol_Ener;
  //! Epsilon to define a zero value for relative velocity in termination condition
  double _ZeroVel_EndIm;
  //! Epsilon to define a zero value for potential energy in termination condition
  double _ZeroEner_EndIm;
  //! we start to save data from _stepMinSave to _stepMaxSave
  unsigned int _stepMinSave, _stepMaxSave;
public:
  //!Default constructor
  OSNSMultipleImpact();
  /** Constructor from data (step size is required here)
   *  \param type  the type of the compliance law
   *  \param step step size estimated
   */
  OSNSMultipleImpact(std::string type, double step);

  //!Destructor
  ~OSNSMultipleImpact();
  /* To get the type of the compliance law at contact
   * \return std::string
   */
  inline std::string get_typeCompLaw() const
  {
    return _typeCompLaw;
  };

  /** To set the type of the compliance law
   * \param newTypeLaw
   */
  void set_typeCompLaw(std::string newTypeLaw);

  /** To set the tolerance to define zero
   * \param  newTolZero
   */
  void setTolImpact(double newTolZero);

  /** To get the tolerance to define zero
   * \return double
   */
  inline double getTolImpact()
  {
    return _tolImpact;
  };

  /** To set the flag to save the data during impact or not
   * \param var
   */
  void SetSaveData(bool var);

  /** To set the name for the output file
   * \param file_name
   */
  void SetNameOutput(std::string file_name);

  /** To get step size
   * \return double
   */
  inline double GetStepSize()
  {
    return _deltaP;
  };

  /* To get the duration of multiple impacts process
   * \return double
   */
  inline double DurationImpact()
  {
    return _timeVariable;
  };

  /** To set the variable _nStepSave
   * \param var
   */
  void SetNstepSave(unsigned int var);

  /** To set the maximal number of steps allowed for each computation
   * \param var
   */
  void SetNstepMax(unsigned int var);

  /** Set number of points to be saved during impact
   * \param var
   */
  void SetSizeDataSave(unsigned int var);

  /** Set tolerence to define whether or not a velocity is zero
   * \param var
   */
  void SetTolVel(double var);

  /** Set tolerence to define whether or not a potential energy is zero
   * \param var
   */
  void SetTolEner(double var);

  /** Set epsilon _ZeroVel_EndIm
   * \param var
   */
  void SetZeroVelEndImp(double var);

  /** Set epsilon _ZeroEner_EndIm
   * \param var
   */
  void SetZeroEnerEndImp(double var);

  /** Set the step number to start the data save and step number to stop save
   * \param min
   * \param max
   */
  void SetStepMinMaxSave(unsigned int min,  unsigned int max);

  /** To compare a double number with zero
   * \param var
   * \return bool   
   */
  bool isZero(double var);

  /** To compare a velocity value with zero
   * \param var 
   * \return bool   
   */
  bool isVelNegative(double var);

  /** To compare an energy value with zero
   * \param var
   * \return bool   
   */
  bool isEnerZero(double var);

  /** To select the pramary contact
   */
  void SelectPrimaContact();

  /** Calculate the vector of distributing rule */
  void Compute_distributionVector();

  /** Compute the normal imulse at contacts */
  void ComputeImpulseContact();

  /** Compute the relative velocity at contacts */
  void Compute_velocityContact();

  /** Compute the potential energy at contacts during each computation step */
  void Compute_energyContact();

  /** Compute the velocity of the bodies during impact */
  void UpdateDuringImpact();

  /** Run the iterative procedure to solve the multiple impact problem */
  void ComputeImpact();

  /** Post-compute for multiple impacts */
  void PostComputeImpact();

  /** Check if the multiple impacts process is terminated or not 
   * \return bool
   */
  bool IsMulImpactTerminate();

  /** To allocate the memory */
  void AllocateMemory();

  /** To build the vector of stiffnesses and restitution coefficient at contacts */
  void BuildParaContact();

  /** To get the velocity of bodies, relative velocity and potential energy at the beginning of impact */
  void InitializeInput();

  /** To check the state of contacts during impact */
  void Check_stateContact();

  /** Pre-compute for multiple impacs */
  void PreComputeImpact();

  /** To get the primary contact according to the relative velocity
      In this case, the primary contact correspond to the contact at which the relative velocity
      is minimum (the relative velocity for two approching bodies is negative so the magnitude of
      the relative velocity at the primary contact is maximum)
  */
  void PrimConVelocity();

  /** To get the primary contact according to the potential energy. In this case, the primary
      contact corresponds to the one at which the potential energy is maximum
  */
  void PrimConEnergy();

  /** To decide if the primary contact is selected according to the relative velocity or to the
   *   potential energy. The first case happens when there is no potential energy at any contact
   * \return bool
   */
  bool IsEnermaxZero();

  /** Verify if the minimum relative velocity at contacts is negative or not 
   * \return bool
   */
  bool IsVcminNegative();

  /** compute the unknown post-impact relative velocity and post-impact impulse
   * \param time
   * \return int
   */
  int compute(double time);

  /**initialize
   * \param sim
   */
  void initialize(SP::Simulation sim);

  /** print the data to the screen */
  void display() const;

  /** To write a SiconosVector into a matrix
   * \param v
   * \param row position starting to write
   * \param col position starting to write
   */

  void WriteVectorIntoMatrix(const SiconosVector& v, const unsigned int row, const unsigned int col);

  /** Save data for each step
   * \param i pointer to be save */
  void SaveDataOneStep(unsigned int i);

  /** Estimate size of data matrix
   * \return unsigned int
   */
  unsigned int EstimateNdataCols();

  ACCEPT_STD_VISITORS();
};

#endif
