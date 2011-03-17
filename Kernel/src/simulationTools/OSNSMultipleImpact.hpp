//This is the header file for the class OSNSMultipleImpact
//==========================================================================================
#ifndef _OSNSMULTIPLEIMPACT_
#define _OSNSMULTIPLEIMPACT_
//-----------------------------------------------------------------------------------------
#include "SimpleVector.hpp"
#include "SimpleMatrix.hpp"
#include "SiconosPointers.hpp"
#include "LinearOSNS.hpp"
#include<fstream>
#include <string>
using namespace std;
using namespace RELATION;

const double DEFAULT_TOL_IMPACT = 100 * MACHINE_PREC;
class OSNSMultipleImpact : public LinearOSNS
{
private:
  // Time-like variable (Impulse)
  double Impulse_variable;
  // Time variable
  double Time_variable;
  // Number of contacts (only the active contacts)
  unsigned int Ncontact;
  // Number of calculation steps estimated
  unsigned int NstepEst;
  // Maximal number of steps for each computation
  unsigned int NstepMax;
  // Tolerance to define zero
  double TOL_IMPACT;
  // Type of the compliance model
  std::string TypeCompLaw;
  //Velocity of bodies during impact
  //SP::SiconosVector VelAllBody;
  // Relative velocity at all Interactions (with or without contact)
  //SP::SiconosVector VelAllIteractions;
  // Relative velocity during impact (at the end of each calculation step)
  SP::SiconosVector VelContact;
  // Relative velocity during impact (at the beginning of each calculation step)
  SP::SiconosVector OldVelContact;
  // Potential energy during impact (at the end of each calculation step)
  SP::SiconosVector EnerContact;
  // Work done during the last compression phase at contact
  SP::SiconosVector WcContact;
  // Distribution vector to distribute the incremental impulse at contact
  SP::SiconosVector DistriVector;
  // State of contacts at the beginning of impact
  // if *StateContact[i] = 0 => no impact at this contact (at contact with positive relative velocity and no potential energy, may be the impact has been terminated at this contact)
  // if *StateContact[i] = 1 => impact takes place at this contact without potential energy (beginning of impact or repeating impact)
  // if *StateContact[i] = 2 => impact takes place with not-zero potential energy
  SP::IndexInt StateContact;
  //Stiffness at contacts
  SP::SiconosVector  Kcontact;
  // Restitution coefficient of contacts
  SP::SiconosVector ResContact;
  // Elasticity coefficient of contacts
  SP::SiconosVector ElasCoefContact;
  // Incremental impulse at contacts
  SP::SiconosVector DelImpulseContact;
  // Total impulse at contacts
  SP::SiconosVector TolImpulseContact;
  // Impulse at contacts for each update time
  SP::SiconosVector ImpulseContact_update;
  // Force at contacts
  SP::SiconosVector ForceContact;
  // Flag to select the primary contact based on the relative velocity or on the potential energy
  // at contacts
  // if SelectPrimaConInVel = true => select the primary contact according to the relative velocity
  // at contact
  // if SelectPrimaConInVel = false => select the primary contact according to the potential energy
  // at contact
  bool SelectPrimaConInVel;
  // ID of the primary contact
  unsigned int IdPrimaContact;
  // Indicator about the selection of the primary contact
  // true if primary contact is selected according to the potential energy
  // false if primary contact is selected according to the relative velocity
  bool IsPrimaConEnergy;
  // Relative velocity at primary contact
  double VelAtPrimaCon;
  // Potential energy at primary contact
  double EnerAtPrimaCon;
  // Step size for the iterative calculation
  double DeltaP;
  // ofstream objet to save the data during impact
  std::ofstream OutputFile;
  // Name of file into which the datat is writen
  std::string  NameFile;
  // YesWriteData = true ==>save the data during impact
  // YesWriteData = false ==> not save the data during impact
  bool YesSaveData;
  // bool variable to set the step size for multiple impact computation
  // If IsNumberOfStepsEst = true ==> estimate the step size from the state of the dynamic system before impact and the number of step needed
  // Number of steps after which the data is saved
  unsigned int NstepSave;
  // If IsNumberOfStepsEst = false ==> user choose the step size
  bool IsNumberOfStepsEst;
  // Matrix on which the data during impact is saved
  SP::SiconosMatrix _DataMatrix;
  // indicator to save the data
  // YesSaveByMatrix = true ==> we save data by a matrix allocated before impact computation
  // YesSaveByMatrix = false ==> we save data by writing in to a file.dat
  bool YesSaveByMatrix;
  // Number of points to be save during impacts
  unsigned int SizeDataSave;
  // indicator on the termination of the multiple impact process
  // _IsImpactEnd = true: impact is terminated
  // _IsImpactEnd = false: otherwise
  bool _IsImpactEnd;
public:
  //Default constructor
  OSNSMultipleImpact();
  //Constructor from data (number of steps is required here)
  //1st parameter: the type of the compliance law
  //3rd parameter: number of steps estimated
  OSNSMultipleImpact(string, unsigned int);
  //Constructor from data (step size is required here)
  //1st parameter: the type of the compliance law
  //3rd parameter: step size estimated
  OSNSMultipleImpact(string, double);
  //Destructor
  ~OSNSMultipleImpact();
  //To get the type of the compliance law at contact
  inline std::string getTypeCompLaw() const
  {
    return TypeCompLaw;
  };
  //To set the type of the compliance law
  void setTypeCompLaw(std::string newTypeLaw);
  // To set the tolerance to define zero
  void setTolImpact(double newTolZero);
  // To get the tolerance to define zero
  inline double getTolImpact()
  {
    return TOL_IMPACT;
  };
  // To set the flag to save the data during impact or not
  void SetYesSaveData(bool var);
  // To set the name for the output file
  void SetNameOutput(std::string file_name);
  // To get step size
  inline double GetStepSize()
  {
    return DeltaP;
  };
  // To get the duration of multiple impacts process
  inline double DurationImpact()
  {
    return Time_variable;
  };
  // To set the variable NstepSave
  void SetNstepSave(unsigned int var);
  // To set the maximal number of steps allowed for each computation
  void SetNstepMax(unsigned int var);
  // To set the indicator YesSaveByMatrix
  void SetYesSaveByMatrix(bool);
  // Set number of points to be saved during impact
  void SetSizeDataSave(unsigned int);
  // To compare a double number with zero
  bool isZero(const double);
  // To calculate the step size for the iterative procedure
  void ComputeStepSize();
  // To select the pramary contact
  void SelectPrimaContact();
  // Calculate the vector of distributing rule
  void ComputeDistriVector();
  // Compute the normal imulse at contacts
  void ComputeImpulseContact();
  // Compute the relative velocity at contacts
  void ComputeVelContact();
  // Compute the potential energy at contacts during each computation step
  void ComputeEnerContact();
  // Compute the velocity of the bodies during impact
  void UpdateDuringImpact();
  // Run the iterative procedure to solve the multiple impact problem
  void ComputeImpact();
  // Post-compute for multiple impacts
  void PostComputeImpact();
  // Check if the multiple impacts process is terminated or not
  //bool IsMulImpactTerminate();
  inline bool IsMulImpactTerminate()
  {
    return _IsImpactEnd;
  }
  // To allocate the memory
  void AllocateMemory();
  // To build the vector of stiffnesses and restitution coefficient at contacts
  void BuildParaContact();
  // To get the velocity of bodies, relative velocity and potential energy at the beginning of impact
  void InitializeInput();
  // To check the state of contacts during impact
  void CheckStateContact();
  // Pre-compute for multiple impacs
  void PreComputeImpact();
  // To get the primary contact according to the relative velocity
  // In this case, the primary contact correspond to the contact at which the relative velocity
  // is minimum (the relative velocity for two approching bodies is negative so the magnitude of
  // the relative velocity at the primary contact is maximum)
  void PrimConVelocity();
  // To get the primary contact according to the potential energy. In this case, the primary
  // contact corresponds to the one at which the potential energy is maximum
  void PrimConEnergy();
  // To decide if the primary contact is selected according to the relative velocity or to the
  // potential energy. The first case happens when there is no potential energy at any contact
  bool IsEnermaxZero();
  // Verify if the minimum relative velocity at contacts is negative or not
  bool IsVcminNegative();
  // compute the unknown post-impact relative velocity and post-impact impulse
  int compute(double);
  //initialize
  void initialize(SP::Simulation);
  // print the data to the screen
  void display() const;
  // To write a SiconosVector into an ouput file
  void WriteSiconosVector(const SiconosVector&);
  // To write a SiconosVector into a matrix
  // row and columns positions starting to write
  void WriteVectorIntoMatrix(const SimpleVector, const unsigned int, const unsigned int);
  // Save data for each step
  // parameter: ith pointer to be save
  void SaveDataOneStep(unsigned int);
  // Estimate size of data matrix
  unsigned int EstimateNdataCols();
  //==================== Friend methods ================================
  ACCEPT_STD_VISITORS();
};
DEFINE_SPTR(OSNSMultipleImpact);
#endif
