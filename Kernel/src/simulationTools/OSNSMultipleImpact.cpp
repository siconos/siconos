// This is the implementation for the class OSNSMultipleImpact
//========================================================================================
#include "OSNSMultipleImpact.hpp"
#include "LagrangianDS.hpp"
#include "MultipleImpactNSL.hpp"
#include "Simulation.hpp"
#include "ioMatrix.hpp"
#include <boost/numeric/ublas/io.hpp>
#include <boost/progress.hpp>
using namespace std;
//Default constructor
OSNSMultipleImpact::OSNSMultipleImpact(): LinearOSNS()
{
  TypeCompLaw = "BiStiffness";
  NstepEst = 10000;
  NstepSave = 100;
  TOL_IMPACT = DEFAULT_TOL_IMPACT;
  YesSaveData = false;
  IsNumberOfStepsEst = true;
  YesSaveByMatrix = true;
  SizeDataSave = 1000;
  NstepMax = 100000;
  NameFile = "DataMultipleImpact.dat";
}
//------------------------------ -------------------------------------------------------------
OSNSMultipleImpact::OSNSMultipleImpact(string newTypeLaw, unsigned int newNstepEst = 10000): LinearOSNS()
{
  TypeCompLaw = newTypeLaw;
  NstepEst = newNstepEst;
  NstepSave = 100;
  TOL_IMPACT = DEFAULT_TOL_IMPACT;
  YesSaveData = false;
  IsNumberOfStepsEst = true;
  NameFile = "DataMultipleImpact.dat";
  YesSaveByMatrix = true;
  SizeDataSave = 1000;
  NstepMax = 100000;
  if ((TypeCompLaw != "MonoStiffness") && (TypeCompLaw != "BiStiffness"))
    RuntimeException::selfThrow("OSNSMultipleImpact::TypeCompLaw type of the compliance model must be either MonoStiffness or BiStiffness!");
}
//------------------------------ -------------------------------------------------------------
OSNSMultipleImpact::OSNSMultipleImpact(string newTypeLaw, double newDelP = 1.0e-5): LinearOSNS()
{
  TypeCompLaw = newTypeLaw;
  DeltaP = newDelP;
  NstepSave = 100;
  TOL_IMPACT = DEFAULT_TOL_IMPACT;
  YesSaveData = false;
  IsNumberOfStepsEst = false;
  NameFile = "DataMultipleImpact.dat";
  YesSaveByMatrix = true;
  SizeDataSave = 1000;
  NstepMax = 100000;
  if ((TypeCompLaw != "MonoStiffness") && (TypeCompLaw != "BiStiffness"))
    RuntimeException::selfThrow("OSNSMultipleImpact::TypeCompLaw type of the compliance model must be either MonoStiffness or BiStiffness!");
}
//-------------------------------------------------------------------------------------------------
OSNSMultipleImpact::~OSNSMultipleImpact() {}
//------------------------------------------------------------------------------------------------

void OSNSMultipleImpact::setTolImpact(double newTolZero)
{
  TOL_IMPACT = newTolZero;
};

void OSNSMultipleImpact::SetYesSaveData(bool var)
{
  YesSaveData = var;
};

void OSNSMultipleImpact::SetNameOutput(std::string file_name)
{
  NameFile = file_name;
};

void OSNSMultipleImpact::SetNstepSave(unsigned int var)
{
  NstepSave = var;
};

void OSNSMultipleImpact::SetNstepMax(unsigned int var)
{
  NstepMax = var;
};

void OSNSMultipleImpact::setTypeCompLaw(std::string newTypeLaw)
{
  TypeCompLaw = newTypeLaw;
  if ((TypeCompLaw != "MonoStiffness") && (TypeCompLaw != "BiStiffness"))
    RuntimeException::selfThrow("OSNSMultipleImpact::TypeCompLaw type of the compliance model must be either MonoStiffness or BiStiffness!");
};

void OSNSMultipleImpact::SetYesSaveByMatrix(bool var)
{
  YesSaveByMatrix = var;
}

void OSNSMultipleImpact::SetSizeDataSave(unsigned int var)
{
  SizeDataSave = var;
}
//--------------------------------------------------------------------------------------------------
void OSNSMultipleImpact::WriteSiconosVector(const SiconosVector& m)
{
  DenseVect*  p = m.dense();
  std::copy(p->begin(), p->end(), std::ostream_iterator<double>(OutputFile, " "));
}
//---------------------------------------------------------------------------------------------------
void OSNSMultipleImpact::WriteVectorIntoMatrix(const SimpleVector m, const unsigned int pos_row, const unsigned int pos_col)
{
  for (unsigned int i = 0; i < m.size(); ++i)
  {
    (*_DataMatrix)(pos_row, pos_col + i) = m(i);
  }
}
//----------------------------------------------------------------------------------------------------
bool OSNSMultipleImpact::isZero(const double Var)
{
  if (std::abs(Var) <= TOL_IMPACT)
    return true;
  else
    return false;
}
//--------------------------------------------------------------------------------------------------
unsigned int OSNSMultipleImpact::EstimateNdataCols()
{
  unsigned int _numberCols = 2;
  // Number of columns for data at contacts
  SP::InteractionsGraph indexSet = simulation()->indexSet(0); // get indexSet[0]
  InteractionsGraph::VIterator ui, uiend;
  for (boost::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui)
  {
    _numberCols = _numberCols + 3 * (indexSet->bundle(*ui)->nonSmoothLaw()->size()) + 1;
  }
  // Number of columns for data at particles
  SP::DynamicalSystemsGraph DSG = simulation()->model()->nonSmoothDynamicalSystem()->dynamicalSystems();
  DynamicalSystemsGraph::VIterator dsi, dsiend;
  for (boost::tie(dsi, dsiend) = DSG->vertices(); dsi != dsiend; ++dsi)
  {
    _numberCols = _numberCols + (DSG->bundle(*dsi)->getDim());
  }
  return(_numberCols);
}
//-----------------------------------------------------------------------------------------------
void OSNSMultipleImpact::AllocateMemory()
{
  if (!VelContact)
    VelContact.reset(new SimpleVector(maxSize()));
  else
  {
    if (VelContact->size() != maxSize())
      VelContact->resize(maxSize());
  };
  //
  if (!OldVelContact)
    OldVelContact.reset(new SimpleVector(maxSize()));
  else
  {
    if (OldVelContact->size() != maxSize())
      OldVelContact->resize(maxSize());
  };
  //
  if (! EnerContact)
    EnerContact.reset(new SimpleVector(maxSize()));
  else
  {
    if (EnerContact->size() != maxSize())
      EnerContact->resize(maxSize());
  };
  //
  if (!WcContact)
    WcContact.reset(new SimpleVector(maxSize()));
  else
  {
    if (WcContact->size() != maxSize())
      WcContact->resize(maxSize());
  };
  //
  if (!DistriVector)
    DistriVector.reset(new SimpleVector(maxSize()));
  else
  {
    if (DistriVector->size() != maxSize())
      DistriVector->resize(maxSize());
  };
  //
  if (!StateContact)
    StateContact.reset(new IndexInt(maxSize()));
  else
  {
    if (StateContact->size() != maxSize())
      StateContact->resize(maxSize());
  };
  //
  if (!Kcontact)
    Kcontact.reset(new SimpleVector(maxSize()));
  else
  {
    if (Kcontact->size() != maxSize())
      Kcontact->resize(maxSize());
  };
  //
  if (!ResContact)
    ResContact.reset(new SimpleVector(maxSize()));
  else
  {
    if (ResContact->size() != maxSize())
      ResContact->resize(maxSize());
  };
  //
  if (!ElasCoefContact)
    ElasCoefContact.reset(new SimpleVector(maxSize()));
  else
  {
    if (ElasCoefContact->size() != maxSize())
      ElasCoefContact->resize(maxSize());
  };
  if (!TolImpulseContact)
    TolImpulseContact.reset(new SimpleVector(maxSize()));
  else
  {
    if (TolImpulseContact->size() != maxSize())
      TolImpulseContact->resize(maxSize());
  };
  //
  if (!DelImpulseContact)
    DelImpulseContact.reset(new SimpleVector(maxSize()));
  else
  {
    if (DelImpulseContact->size() != maxSize())
      DelImpulseContact->resize(maxSize());
  };
  //
  if (!ImpulseContact_update)
    ImpulseContact_update.reset(new SimpleVector(maxSize()));
  else
  {
    if (ImpulseContact_update->size() != maxSize())
      ImpulseContact_update->resize(maxSize());
  }
  //
  if (!ForceContact)
    ForceContact.reset(new SimpleVector(maxSize()));
  else
  {
    if (ForceContact->size() != maxSize())
      ForceContact->resize(maxSize());
  };
  // for the data matrix
  if (YesSaveByMatrix)
  {
    unsigned int _numberCols = EstimateNdataCols();
    if (!_DataMatrix)
      _DataMatrix.reset(new SimpleMatrix(SizeDataSave, _numberCols));
    else
    {
      if ((_DataMatrix->size(0) != SizeDataSave) || (_DataMatrix->size(1) != _numberCols))
        _DataMatrix->resize(SizeDataSave, _numberCols);
    }
  }
}
//=====================================================================================
void OSNSMultipleImpact::BuildParaContact()
{
  SP::InteractionsGraph indexSet = simulation()->indexSet(1); // get indexSet[1]
  //Loop over the Interactionof the indexSet(1)
  InteractionsGraph::VIterator ui, uiend;
  for (boost::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui)
  {
    SP::Interaction inter = indexSet->bundle(*ui);
    SP::NonSmoothLaw nslaw = inter->nslaw();
    SP::MultipleImpactNSL Mulnslaw = boost::dynamic_pointer_cast<MultipleImpactNSL>(nslaw);
    assert(Mulnslaw && "In OSNSMultipleImpact::BuildStiffResCofVec, non-smooth law used must be MultipleImpactNSL!!!");
    // Get the relative position of inter-interactionBlock in the vector VelContact
    unsigned int pos = _M->getPositionOfInteractionBlock(inter);
    (*ResContact)(pos) = Mulnslaw->ResCof();
    (*Kcontact)(pos) = Mulnslaw->Stiff();
    (*ElasCoefContact)(pos) = Mulnslaw->ElasCof();
  }
  /*
     cout << " Restitution coefficients: " << endl;
     ResContact->display();
     cout << "Stiffnesses: " << endl;
     Kcontact->display();
     cout << "Elasticity coeffients at contacts: " << endl;
     ElasCoefContact->display();
     */

}
//======================================================================================
void OSNSMultipleImpact::ComputeStepSize()
{
  SP::DynamicalSystemsGraph DSG = simulation()->model()->nonSmoothDynamicalSystem()->dynamicalSystems();
  SP::InteractionsGraph indexSet1 = simulation()->indexSet(1); // get indexSet[1]
  std::set<SP::DynamicalSystem> setOfDS;
  UpdateDuringImpact(); // Update the dynamical state
  //Loop over all Interactionto obtain the list of DS involved into impact
  unsigned int pos;
  SP::DynamicalSystem DS1;
  SP::DynamicalSystem DS2;
  InteractionsGraph::VIterator ui, uiend;
  for (boost::tie(ui, uiend) = indexSet1->vertices(); ui != uiend; ++ui)
  {
    SP::Interaction inter = indexSet1->bundle(*ui);
    // Get the relative position of inter-interactionBlock in the vector VelContact/TolImpulseContact
    pos = _M->getPositionOfInteractionBlock(inter);
    // Check if impact take place at this contact
    if ((*StateContact)[pos] != 0)
    {
      DS1 = indexSet1->properties(*ui).source;
      DS2 = indexSet1->properties(*ui).target;
      assert(DS1);
      assert(DS2);
      setOfDS.insert(DS1); // source DS on this Interaction
      setOfDS.insert(DS2); // target DS on this Interaction
    }
  }
  //
  std::set<SP::DynamicalSystem>::iterator it;
  double Pest = 0.0;
  for (it = setOfDS.begin(); it != setOfDS.end(); it++)
  {
    SP::LagrangianDS Lagds = boost::dynamic_pointer_cast<LagrangianDS>(*it);
    SP::SiconosMatrix mass_ds =  Lagds->mass();    // Mass matrix of DS
    SP::SiconosVector vel_ds = Lagds->velocity();  // Pre-impact velocity of DS

    SP::SiconosVector abs_vel_ds(new SimpleVector(vel_ds->size()));
    abs_wise(*vel_ds, *abs_vel_ds); // get the absolute (velocity vector of ds)
    SP::SimpleVector prod_mass_vel(new SimpleVector(mass_ds->size(0)));
    prod(*mass_ds, *abs_vel_ds, *prod_mass_vel);
    Pest = Pest + prod_mass_vel->sum();

    // cout << "DS involved into impact: " << (*it)->number() << endl;
    //cout << "mass of DS" << endl;
    // mass_ds->display();
    //cout << "velocity of DS" << endl;
    //vel_ds->display();
  }
  DeltaP = Pest / ((setOfDS.size()) * NstepEst);
  //cout << "Size of setOfDS" << setOfDS.size() << endl;
  //cout << "Step size is: " << DeltaP << endl;
  if (DeltaP <= 0.0)
    RuntimeException::selfThrow("OSNSMultipleImpact::ComputeStepSize the step size DeltaP must be positive  !!");


  // SP::DynamicalSystemsGraph DSG = simulation()->model()->nonSmoothDynamicalSystem()->dynamicalSystems();
  // SP::InteractionsGraph indexSet1 = simulation()->indexSet(1); // get indexSet[1]
  // //Loop over alls DS
  // DynamicalSystemsGraph::VIterator ui, uiend;
  // DynamicalSystemsGraph::OEIterator edi, ediend;
  // double Pest = 0.0;
  // for (boost::tie(ui,uiend) = DSG->vertices(); ui != uiend; ++ui)
  //   {
  //     // Check if this DS is involved in the impact or not
  //     bool found = false;
  //     // Loop over all edges comming from this DS
  //     for(boost::tie(edi,ediend) = DSG->out_edges(*ui); edi != ediend; ++edi)
  //       {
  //         SP::Interaction interp = DSG->bundle(*edi);
  //         if (indexSet1->is_vertex(inter)) // there exist at least one edge of this DS in the indexSet[1]
  //           {
  //             found = true;
  //             break;
  //           }
  //       }
  //     //
  //     if (found) // if this DS involved in the impact
  //       {
  //         SP::DynamicalSystem ds = DSG->bundle(*ui); // DS
  //         SP::LagrangianDS Lagds = boost::dynamic_pointer_cast<LagrangianDS>(ds);
  //         SP::SiconosMatrix mass_ds =  Lagds->mass();    // Mass matrix of DS
  //         SP::SiconosVector vel_ds = Lagds->velocityMemory()->getSiconosVector(1);  // Pre-impact velocity of DS
  //         SP::SiconosVector abs_vel_ds(new SimpleVector(vel_ds->size()));
  //         abs_wise(*vel_ds,*abs_vel_ds);// get the absolute (velocity vector of ds)
  //         SP::SimpleVector prod_mass_vel(new SimpleVector(mass_ds->size(0)));
  //         prod(*mass_ds,*abs_vel_ds,*prod_mass_vel);
  //         Pest = Pest + prod_mass_vel->sum();
  //       }
  //   };
  // DeltaP = Pest/NstepEst;
  // if (DeltaP <= 0.0)
  //   RuntimeException::selfThrow("OSNSMultipleImpact::ComputeStepSize the step size DeltaP must be positive  !!");
  // //
  //cout << "Step size:" << DeltaP << endl;
}
//========================================================================================
void OSNSMultipleImpact::PreComputeImpact()
{
  //1. Get the number of contacts and bodies involved in the impact
  if (levelMin() != 1)
    RuntimeException::selfThrow("OSNSMultipleImpact::PreComputeImpact==> the levelMin must be equal to 1 in the multiple impact model !!");
  SP::InteractionsGraph indexSet = simulation()->indexSet(levelMin()); // get indexSet[1]
  Ncontact = indexSet->size();
  //2. Compute matrix _M
  SP::Topology topology = simulation()->model()->nonSmoothDynamicalSystem()->topology();
  bool isLinear = simulation()->model()->nonSmoothDynamicalSystem()->isLinear();
  if (!_hasBeenUpdated || !isLinear)
  {
    // Computes new _interactionBlocks if required
    updateInteractionBlocks();
    // Updates matrix M
    _M->fill(indexSet, !_hasBeenUpdated);
    _sizeOutput = _M->size();
  }
  if (Ncontact != _sizeOutput)
    RuntimeException::selfThrow("OSNSMultipleImpact::ComputeWMinvWtrans: number of contacts different from the size of output--> this case is not yet implemented!");
  //3. Checks size of vectors
  if (VelContact->size() != _sizeOutput)
  {
    VelContact->resize(_sizeOutput);
  }
  VelContact->zero();
  //
  if (OldVelContact->size() != _sizeOutput)
  {
    OldVelContact->resize(_sizeOutput);
  }
  OldVelContact->zero();
  //
  if (EnerContact->size() != _sizeOutput)
  {
    EnerContact->resize(_sizeOutput);
  }
  EnerContact->zero();
  //
  if (WcContact->size() != _sizeOutput)
  {
    WcContact->resize(_sizeOutput);
  }
  WcContact->zero();
  //
  if (DistriVector->size() != _sizeOutput)
  {
    DistriVector->resize(_sizeOutput);
  }
  DistriVector->zero();
  //
  if (StateContact->size() != _sizeOutput)
  {
    StateContact->resize(_sizeOutput);
  }
  //
  if (Kcontact->size() != _sizeOutput)
  {
    Kcontact->resize(_sizeOutput);
  }
  Kcontact->zero();
  //
  if (ResContact->size() != _sizeOutput)
  {
    ResContact->resize(_sizeOutput);
  }
  ResContact->zero();
  //
  if (ElasCoefContact->size() != _sizeOutput)
  {
    ElasCoefContact->resize(_sizeOutput);
  }
  ElasCoefContact->zero();
  //
  if (TolImpulseContact->size() != _sizeOutput)
  {
    TolImpulseContact->resize(_sizeOutput);
  }
  TolImpulseContact->zero();
  //
  if (DelImpulseContact->size() != _sizeOutput)
  {
    DelImpulseContact->resize(_sizeOutput);
  }
  DelImpulseContact->zero();
  //
  if (ImpulseContact_update->size() != _sizeOutput)
  {
    ImpulseContact_update->resize(_sizeOutput);
  }
  ImpulseContact_update->zero();
  //
  if (ForceContact->size() != _sizeOutput)
  {
    ForceContact->resize(_sizeOutput);
  }
  ForceContact->zero();
  //4. Initialize the relative velocity, potential energy, impulse at contacts
  InitializeInput();
  // //5. Compute the step size
  // if (IsNumberOfStepsEst) // We need to estimate the step size from the initial date of the dynamic system
  //   {
  //     ComputeStepSize();
  //   };
  //6. Build the vectors of stifnesseses, of restitution coefficients, and of elaticity coefficients
  BuildParaContact();
  /*
     cout<<"Matrix OSNS problem:  " << endl;
     _M->display();
     */
}
//=======================================================================================
void OSNSMultipleImpact::InitializeInput()
{
  //Loop over alls Interactioninvolved in the indexSet[1]
  SP::InteractionsGraph indexSet = simulation()->indexSet(levelMin()); // get indexSet[1]
  InteractionsGraph::VIterator ui, uiend;
  for (boost::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui)
  {
    SP::Interaction inter = indexSet->bundle(*ui);
    //SP::SiconosVector Vc0 = inter->y(1); // Relative velocity at beginning of impact
    SP::SiconosVector Vc0 = inter->yOld(1); // Relative velocity at beginning of impact
    unsigned int pos_inter = _M->getPositionOfInteractionBlock(inter);
    setBlock(*Vc0, VelContact, Vc0->size(), 0, pos_inter);
    SP::SiconosVector ener0(new SimpleVector(Vc0->size()));
    ener0->zero(); // We suppose that the initial potential energy before impact is equal to zero at any contact
    // at the beginning of impact
    setBlock(*ener0, EnerContact, ener0->size(), 0, pos_inter);
    //SP::SiconosVector impulse0= (inter)->lambda(1))->vector(inter->number());
    SP::SiconosVector impulse0(new SimpleVector(Vc0->size()));
    impulse0->zero(); // We suppose that the impulse before impact is equal to zero at any contact
    // at the beginning of impact
    setBlock(*impulse0, TolImpulseContact, impulse0->size(), 0, pos_inter);
  };
  /*
     cout << "Initial relative velocity at contacts" << endl;
     VelContact->display();
     cout<< "Initial energy at contacts" << endl;
     EnerContact->display();
     cout << "Impulse at contact" << endl;
     TolImpulseContact->display();
     */
}
//=========================================================================================
void OSNSMultipleImpact::initialize(SP::Simulation sim)
{

  // General initialize for OneStepNSProblem
  OneStepNSProblem::initialize(sim);
  // Allocate the memory
  AllocateMemory();
  // get topology
  SP::Topology topology = simulation()->model()->nonSmoothDynamicalSystem()->topology();
  // Note that _interactionBlocks is up to date since updateInteractionBlocks
  // has been called during OneStepNSProblem::initialize()

  if (! _M)
  {
    if (_MStorageType == 0)
      _M.reset(new OSNSMatrix(maxSize(), 0));

    else // if(_MStorageType == 1) size = number of _interactionBlocks
      // = number of Interactionin the largest considered indexSet
      _M.reset(new OSNSMatrix(simulation()->indexSet(levelMin())->size(), 1));
  }

};
//========================================================================================
void OSNSMultipleImpact::PrimConVelocity()
{
  getMin(*VelContact, VelAtPrimaCon, IdPrimaContact);
  EnerAtPrimaCon = (*EnerContact)(IdPrimaContact);
  if (VelAtPrimaCon >= 0.0)
    RuntimeException::selfThrow("OSNSMultipleImpact::PrimConVelocity the velocity at the primary contact must be negative !!");
  /*
     cout << "Primary contact according to relative velocity: " << IdPrimaContact << endl;
     cout << "Relative velocity at the primary contact: " << VelAtPrimaCon << endl;
     cout << "Potential energy at the primary contact: " << EnerAtPrimaCon << endl;
     */
}
//=======================================================================================
void OSNSMultipleImpact::PrimConEnergy()
{
  getMax(*EnerContact, EnerAtPrimaCon, IdPrimaContact);
  VelAtPrimaCon = (*VelContact)(IdPrimaContact);
  if ((EnerAtPrimaCon < 0.0) || (isZero(EnerAtPrimaCon)))
    RuntimeException::selfThrow("OSNSMultipleImpact::PrimConEnergy the potential energy at the primary contact must be different from zero !!");
  /*
     cout << "Primary contact according to potenial energy: " << IdPrimaContact << endl;
     cout << "Relative velocity at the primary contact: " << VelAtPrimaCon << endl;
     cout << "Potential energy at the primary contact: " << EnerAtPrimaCon << endl;
     */
}
//======================================================================================
bool OSNSMultipleImpact::IsEnermaxZero()
{
  double MaxEner;
  unsigned int IdMax;
  getMax(*EnerContact, MaxEner, IdMax);
  if (isZero(MaxEner))
    return true;
  else
    return false;
}
//======================================================================================
bool OSNSMultipleImpact::IsVcminNegative()
{
  double MinVelCon;
  unsigned int IdConVmin;
  getMin(*VelContact, MinVelCon, IdConVmin);
  if (MinVelCon < 0.0)
    return true;
  else
    return false;
}
//=======================================================================================
void OSNSMultipleImpact::CheckStateContact()
{
  _IsImpactEnd = true;
  for (unsigned int i = 0; i < Ncontact; ++i)
  {
    if (isZero((*EnerContact)(i)))
    {
      if ((*VelContact)(i) >= 0.0) // no impact at this contact
        (*StateContact)[i] = 0;
      else  // impact happens without potential energy
      {
        (*StateContact)[i] = 1;
        _IsImpactEnd = false;
      }
    }
    else // impact happens with not zero potential energy
    {
      _IsImpactEnd = false;
      if ((*StateContact)[i] != 2)
      {
        (*StateContact)[i] = 2;
      }
    }
  };
  //
  /*
     cout << "State at contacts: ";
     for (unsigned int i = 0; i < Ncontact;++i)
     {
     cout << (*StateContact)[i] << " ";
     };
     cout << endl;
     */
  //
}
//=======================================================================================
// bool OSNSMultipleImpact::IsMulImpactTerminate(){
//   bool var = true;
//   for(unsigned int i = 0; i < Ncontact;++i)
//     {
//       if ((*StateContact)[i] != 0)
//         {
//           var = false;
//           break;
//         };
//     };
//   return var;
//   /*
//   cout << "Is the multiple impacts is terminated: " << var << endl;
//   */
// }
//=======================================================================================
void OSNSMultipleImpact::SelectPrimaContact()
{
  if (IsEnermaxZero()) // case of no potential energy at any contact
  {
    PrimConVelocity(); // Select the primary contact according to the relative velocity at contact
    IsPrimaConEnergy = false;
  }
  else
  {
    PrimConEnergy(); // Select the primary contact according to the potential energy at contacts
    IsPrimaConEnergy = true;
  }
  //
  // cout << "The primary contact is :" << IdPrimaContact << endl;
  // cout << "Is the primary contact is selected according to the potential energy: " << IsPrimaConEnergy << endl;
}
//=======================================================================================
void OSNSMultipleImpact::ComputeDistriVector()
{
  //Case 1: if no potential energy at any contact
  double _ratio_mu, ratio_stiff, ratio_ener;
  double mu_prima = (*ElasCoefContact)(IdPrimaContact); // Elasticity coefficient at the primary contact
  double stiff_prima = (*Kcontact)(IdPrimaContact);     // Stiffness at the primary contact
  double _mu, _stiff, _vel, _energy;
  if (!IsPrimaConEnergy) // case of primary contact selected according to the relative velocity
  {
    double ratio_vel;
    for (unsigned int i = 0; i < Ncontact; ++i)
    {
      if ((*StateContact)[i] != 0) // the impact can takes place at this contact
      {
        _mu = (*ElasCoefContact)(i); // Elasticity coefficient at the current contact
        _stiff = (*Kcontact)(i);     // Stiffness at the current contact
        _vel = (*VelContact)(i);     // Relative velocity at the current contact
        _ratio_mu = (std::pow(_mu + 1, (_mu / (_mu + 1)))) / (std::pow(mu_prima + 1, (mu_prima / (mu_prima + 1))));
        ratio_stiff = (std::pow(_stiff, (1 / (1 + _mu)))) / (std::pow(stiff_prima, (1 / (1 + mu_prima))));
        if (_vel > 0.0)
        {
          RuntimeException::selfThrow("OSNSMultipleImpact::ComputeDistriVector, the relative velocity when particle starts to impact must be negative!!");
        }
        ratio_vel = (std::pow(std::fabs(_vel), (_mu / (_mu + 1)))) / (std::pow(std::fabs(VelAtPrimaCon), (mu_prima / (1 + mu_prima))));
        (*DistriVector)(i) = std::pow((_ratio_mu * ratio_stiff * ratio_vel), (1 + _mu)) * std::pow(DeltaP, ((_mu - mu_prima) / (1 + mu_prima)));
      }
      else
      {
        (*DistriVector)(i) = 0.0;
      }
      if ((*DistriVector)(i) < 0.0)
        RuntimeException::selfThrow("OSNSMultipleImpact::ComputeDistriVector the component of DistriVector must be positive !!");
    };
  }
  //Case 2: case of primary contact selected according to the potential energy
  else
  {
    for (unsigned int i = 0; i < Ncontact; ++i)
    {
      _mu = (*ElasCoefContact)(i);
      _stiff = (*Kcontact)(i);
      _ratio_mu = (std::pow(_mu + 1, (_mu / (_mu + 1)))) / (std::pow(mu_prima + 1, (mu_prima / (mu_prima + 1))));
      ratio_stiff = (std::pow(_stiff, (1 / (1 + _mu)))) / (std::pow(stiff_prima, (1 / (1 + mu_prima))));
      if ((*StateContact)[i] == 1) // no potential energy at this contact, including the contacts at which impact repeats
      {
        if ((*VelContact)(i) > 0)
        {
          RuntimeException::selfThrow("OSNSMultipleImpact::ComputeDistriVector, the pre-impact velocity must be negative!!");
        }
        else
        {
          _vel = (*VelContact)(i);
          ratio_ener = (std::pow(std::fabs(_vel), (_mu / (_mu + 1)))) / (std::pow(EnerAtPrimaCon, (mu_prima / (mu_prima + 1))));
          (*DistriVector)(i) = std::pow((_ratio_mu * ratio_stiff * ratio_ener), (1 + _mu)) * std::pow(DeltaP, _mu);
        }
      }
      else if ((*StateContact)[i] == 2) // potential is not zero at this contact
      {
        _energy = (*EnerContact)(i); // Potential energy at the current contact
        ratio_ener = (std::pow(_energy, (_mu / (_mu + 1)))) / (std::pow(EnerAtPrimaCon, (mu_prima / (mu_prima + 1))));
        (*DistriVector)(i) = _ratio_mu * ratio_stiff * ratio_ener;
      }
      else // no impact at this contact
      {
        (*DistriVector)(i) = 0.0;
      };
      if ((*DistriVector)(i) < 0.0)
        RuntimeException::selfThrow("OSNSMultipleImpact::ComputeDistriVector the component of DistriVector must be positive !!");
    };
  };
  // //Case 1: if no potential energy at any contact
  // if (!IsPrimaConEnergy) // case of primary contact selected according to the relative velocity
  //   {
  //     for(unsigned int i = 0; i< Ncontact; ++i)
  //       {
  //         if ((*StateContact)[i] != 0) // the impact can takes place at this contact
  //           {
  //             (*DistriVector)(i) = ((*Kcontact)(i)/(*Kcontact)(IdPrimaContact))*std::pow(((*VelContact)(i)/VelAtPrimaCon),PowCompLaw);
  //           }
  //         else
  //           {
  //             (*DistriVector)(i) = 0.0;
  //           }
  //         if ((*DistriVector)(i) < 0.0)
  //           RuntimeException::selfThrow("OSNSMultipleImpact::ComputeDistriVector the component of DistriVector must be positive !!");
  //       };
  //   }
  // //Case 2: case of primary contact selected according to the potential energy
  // else
  //   {
  //     for (unsigned int i = 0; i < Ncontact; ++i)
  //       {
  //         if ((*StateContact)[i]==1) // no potential energy at this contact, including the contacts at which impact repeats
  //           {
  //             if ((*VelContact)(i) > 0)
  //               {
  //                 RuntimeException::selfThrow("OSNSMultipleImpact::ComputeDistriVector, the pre-impact velocity must be negative!!");
  //               }
  //             else
  //               {
  //                 (*DistriVector)(i) = ((*Kcontact)(i)/(*Kcontact)(IdPrimaContact))*std::pow(((std::fabs((*VelContact)(i))*DeltaP)/EnerAtPrimaCon),PowCompLaw);
  //               }
  //           }
  //         else if ((*StateContact)[i]==2) // potential is not zero at this contact
  //           {
  //             (*DistriVector)(i) = std::pow((*Kcontact)(i)/(*Kcontact)(IdPrimaContact),1.0/(1.0 + PowCompLaw))*std::pow((*EnerContact)(i)/EnerAtPrimaCon,PowCompLaw/(PowCompLaw + 1.0));
  //           }
  //         else // no impact at this contact
  //           {
  //             (*DistriVector)(i) = 0.0;
  //           };
  //         if ((*DistriVector)(i) < 0.0)
  //           RuntimeException::selfThrow("OSNSMultipleImpact::ComputeDistriVector the component of DistriVector must be positive !!");
  //       };
  //   };
  // /*
  //cout << "Tensor of distributing rule:" << endl;
  //DistriVector->display();
  // */

}
//=======================================================================================
void OSNSMultipleImpact::ComputeImpulseContact()
{
  (*DelImpulseContact) = (*DistriVector) * DeltaP;
  (*TolImpulseContact) = (*TolImpulseContact) + (*DelImpulseContact);
  (*ImpulseContact_update) = (*ImpulseContact_update) + (*DelImpulseContact);
  // Compute the contact force
  double PowCompLaw;
  for (unsigned i = 0; i < Ncontact; ++i)
  {
    PowCompLaw = (*ElasCoefContact)(i);
    if (isZero((*EnerContact)(i))) // if potential energy at this contact is zero
    {
      if ((*VelContact)(i) < 0.0) // if the relative velocity at contact is negative
      {
        (*ForceContact)(i) = std::pow((1.0 + PowCompLaw), PowCompLaw / (1.0 + PowCompLaw)) * std::pow((*Kcontact)(i), 1.0 / (1.0 + PowCompLaw)) * std::pow((std::fabs((*VelContact)(i)) * (*DelImpulseContact)(i)), PowCompLaw / (1.0 + PowCompLaw));
      }
      else
      {
        (*ForceContact)(i) = 0.0;
      };
    }
    else
    {
      (*ForceContact)(i) = std::pow((1.0 + PowCompLaw), PowCompLaw / (1.0 + PowCompLaw)) * std::pow((*Kcontact)(i), 1.0 / (1.0 + PowCompLaw)) * std::pow((*EnerContact)(i), PowCompLaw / (1.0 + PowCompLaw));
    }
  };
  //
  /*
     cout << "Increment of normal impulse at contacts:" << endl;
     DelImpulseContact->display();
     cout << "Total normal impulse at contacts:" << endl;
     TolImpulseContact->display();
     cout << "Force at contacts: " << endl;
     ForceContact->display();
     */
  //

}
//=======================================================================================
void OSNSMultipleImpact::ComputeVelContact()
{
  (*OldVelContact) = (*VelContact); //save the relative velocity at the beginning of the step
  (*VelContact) = (*VelContact) + prod(*(_M->defaultMatrix()), *DelImpulseContact); // compute the relative velocity at the end of the step
  //
  /*
     cout << "Relative velocity at contacts at the beginning of step:" << endl;
     OldVelContact->display();
     cout << "Relative velocity at contacts at the end of step:" << endl;
     VelContact->display();
     */
  //
}
//=======================================================================================
void OSNSMultipleImpact::ComputeEnerContact()
{
  if (TypeCompLaw == "BiStiffness")
    // For Bistiffness model
  {
    for (unsigned i = 0; i < Ncontact; ++i)
    {
      if ((*OldVelContact)(i) < 0.0) // Contact located in the compression phase
      {
        (*EnerContact)(i) = (*EnerContact)(i) - 0.5 * ((*OldVelContact)(i) + (*VelContact)(i)) * ((*DelImpulseContact)(i));
      }
      else                       // Contact located in the expansion phase
      {
        if (!isZero((*ResContact)(i)))
          (*EnerContact)(i) = (*EnerContact)(i) - (1.0 / std::pow((*ResContact)(i), 2)) * 0.5 * ((*OldVelContact)(i) + (*VelContact)(i)) * ((*DelImpulseContact)(i));
        else
          (*EnerContact)(i) = 0.0; // In this case, no potential energy at contacts when the contact is located in the compression phase

      };
      if ((*EnerContact)(i) < 0.0)
      {
        (*EnerContact)(i) = 0.0;
      };
    };
  }
  else
    // For the mono-stiffness model
  {
    for (unsigned i = 0; i < Ncontact; ++i)
    {
      //1: dertermine the work done by the last compression phase at contacts (nessessary for the mono-stiffness compliance model)
      // if Vc(k) < 0 and Vc(k +1) >= 0 ==> transition from the compression phase to the expansion phase, Wc = E(k)
      if (((*OldVelContact)(i) < 0.0) && ((*VelContact)(i) >= 0.0))
      {
        (*WcContact)(i) = (*EnerContact)(i);
      };
      //2: Calculate the potential energy at the end of stap
      (*EnerContact)(i) = (*EnerContact)(i) - 0.5 * ((*OldVelContact)(i) + (*VelContact)(i)) * ((*DelImpulseContact)(i));
      //3: Check if the termination condition is verified or not (if Vc(k+1) > 0.0 and E(k+1) <= (1-e^2)*Wc). If yes, discard the potential energy
      // in order to respect the energetic constraint
      if (((*StateContact)[i] == 2) && (((*VelContact)(i) > 0.0) && ((*EnerContact)(i) <= ((1.0 - std::pow((*ResContact)(i), 2)) * (*WcContact)(i)))))
      {
        (*EnerContact)(i) = 0.0; // potential energy at this contact is completely dissipated before the compression phase finishes
      };
    };
  }

  //
  /*
     cout << "Potential energy at contacts at the end of step:" << endl;
     EnerContact->display();
     cout << "Work done during the compression phase at contacts" << endl;
     WcContact->display();
     */
  //
}
//======================================================================================
void OSNSMultipleImpact::UpdateDuringImpact()
{
  //1. Copy VelContact/DelImpulseContact into the vector y/lambda for Interactions
  SP::InteractionsGraph indexSet = simulation()->indexSet(levelMin());
  // y and lambda vectors
  SP::SiconosVector lambda;
  SP::SiconosVector y;
  // === Loop through "active" Interactions (ie present in indexSets[1]) ===
  unsigned int pos;
  InteractionsGraph::VIterator ui, uiend;
  for (boost::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui)
  {
    SP::Interaction inter = indexSet->bundle(*ui);
    // Get the relative position of inter-interactionBlock in the vector VelContact/TolImpulseContact
    pos = _M->getPositionOfInteractionBlock(inter);
    // Get Y and Lambda for the current Interaction
    y = inter->y(levelMin());
    lambda = inter->lambda(levelMin());
    // Copy VelContact/TolImpulseContact, starting from index pos into y/lambda
    // save into y !!
    setBlock(*VelContact, y, y->size(), pos, 0);
    // saved into lambda[1] !!
    setBlock(*ImpulseContact_update, lambda, lambda->size(), pos, 0);
    //setBlock(*DelImpulseContact, lambda, lambda->size(), pos, 0);
  };
  //2. Update the Input[1], state of DS systems, Output[1]
  simulation()->update(levelMin());
  ImpulseContact_update->zero(); // reset input[1] to zero after each update
}
//--------------------------------------------------------------------------------------
void OSNSMultipleImpact::SaveDataOneStep(unsigned int _ithPoint)
{
  // Save the total impulse at the primary contacts (time-like independent variable) and the time evolution during impact
  if (YesSaveByMatrix)
  {
    if (_ithPoint >= _DataMatrix->size(0))
      RuntimeException::selfThrow("In OSNSMultipleImpact::ComputeImpact, number of points saved exceeds the size of matrix allocated!!!");
    (*_DataMatrix)(_ithPoint, 0) =  Impulse_variable;
    (*_DataMatrix)(_ithPoint, 1) = Time_variable;
  }
  else
  {
    OutputFile << Impulse_variable << " ";
    OutputFile << Time_variable << " ";
  }
  // Save the data related to Interactions
  SP::InteractionsGraph indexSet0 = simulation()->indexSet(0);
  SP::InteractionsGraph indexSet1 = simulation()->indexSet(levelMin());
  unsigned int pos;
  int Status_inter;
  InteractionsGraph::VIterator ui, uiend;
  unsigned int col_pos = 2;
  for (boost::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui)
  {
    SP::Interaction inter = indexSet0->bundle(*ui); // Interaction
    unsigned int id_contact = inter->number();
    SP::SiconosVector ydot = inter->y(1);
    SP::SiconosVector P_inter(new SimpleVector(inter->nonSmoothLaw()->size()));
    SP::SiconosVector F_inter(new SimpleVector(inter->nonSmoothLaw()->size()));
    SP::SiconosVector E_inter(new SimpleVector(1));

    if (indexSet1->is_vertex(inter)) // if Interactionbelongs to the IndexSet[1]
    {
      pos = _M->getPositionOfInteractionBlock(inter);
      Status_inter = (*StateContact)[pos]; //status of Interaction
      setBlock(*TolImpulseContact, P_inter, P_inter->size(), pos, 0);
      setBlock(*ForceContact, F_inter, F_inter->size(), pos, 0);
      setBlock(*EnerContact, E_inter, E_inter->size(), pos, 0);
    }
    else
    {
      Status_inter = -1; // no contact at this Interaction
      P_inter->zero();   // no impulse at this Interaction
      F_inter->zero();   // no force at this Interaction
      E_inter->zero();   // no potential at this Interaction
    };
    //
    if (YesSaveByMatrix)
    {
      // Write the relative velocity at the Interaction
      WriteVectorIntoMatrix(*ydot, _ithPoint, col_pos);
      col_pos = col_pos + ydot->size();
      // Write the total impulse at the Interaction
      WriteVectorIntoMatrix(*P_inter, _ithPoint, col_pos);
      col_pos = col_pos + P_inter->size();
      // Write the force at the Interaction
      WriteVectorIntoMatrix(*F_inter, _ithPoint, col_pos);
      col_pos = col_pos + F_inter->size();
      // Write the potential energy at the Interaction
      WriteVectorIntoMatrix(*E_inter, _ithPoint, col_pos);
      col_pos = col_pos + E_inter->size();
    }
    else
    {
      OutputFile << id_contact  << " "; // Write the Id of the its parent interaction
      OutputFile << Status_inter << " "; // Write the status of the Interaction
      WriteSiconosVector(*ydot); // Write the relative velocity at the Interaction
      WriteSiconosVector(*P_inter); // Write the total impulse at the Interaction
      WriteSiconosVector(*F_inter); // Write the force at the Interaction
      WriteSiconosVector(*E_inter); // Write the potential energy at the Interaction
    }
  }
  // Save the data related to DS
  SP::DynamicalSystemsGraph DSG = simulation()->model()->nonSmoothDynamicalSystem()->dynamicalSystems();
  DynamicalSystemsGraph::VIterator dsi, dsiend;
  for (boost::tie(dsi, dsiend) = DSG->vertices(); dsi != dsiend; ++dsi)
  {
    SP::DynamicalSystem ds = DSG->bundle(*dsi); // DS
    unsigned int id_ds = ds->number();
    SP::LagrangianDS Lagds = boost::dynamic_pointer_cast<LagrangianDS>(ds);
    SP::SiconosVector qdot = Lagds->velocity();
    // Write
    if (YesSaveByMatrix)
    {
      WriteVectorIntoMatrix(*qdot, _ithPoint, col_pos);
      col_pos = col_pos + qdot->size();
    }
    else
    {
      OutputFile << id_ds  << " "; // number of DS
      WriteSiconosVector(*qdot);
    }
  }
  if (!YesSaveByMatrix)
  {
    // Terminate the line
    OutputFile << "\n";
  }
}
//=======================================================================================
void OSNSMultipleImpact::ComputeImpact()
{
  // open the stream file
  if ((YesSaveData) && (!YesSaveByMatrix))
  {
    OutputFile.open(NameFile.c_str(), ios::out | ios::trunc);
    if (!OutputFile.is_open())
      RuntimeException::selfThrow("OSNSMultipleImpact::ComputeImpact==> write error : Fail to open \"" + NameFile + "\"");
    OutputFile.precision(15);
  }
  //
  Impulse_variable = 0.0;
  Time_variable = 0.0;
  unsigned int number_step = 1;
  unsigned int point_save = 0;
  unsigned int _counterstepsave = 1;
  unsigned int Old_prima_contact;
  bool IsPrimaConChange = true;
  // Show computation progress
  cout << "*********** Impact computation progress *************" << endl;
  boost::progress_display show_progress(NstepMax);
  //
  /*
     cout << "----------Before multiple impacts computation---------------" << endl;
     cout << "Velocity at contacts: ";
     VelContact->display();
     cout << "Impulse at contact: ";
     TolImpulseContact->display();
     */
  //cout << "-------------------Multiple impacts computation starts:-----------------------" << endl;
  // First save at the beginning of impact computation
  if (YesSaveData)
  {
    SaveDataOneStep(point_save); // Save the data
    point_save++;
  }
  //
  while (1 != 0)
  {
    //
    /*
       cout << "==================Step==================:  " << number_step << endl;
       cout << "Impulse variable: " << Impulse_variable << endl;
       cout << "Time_variable: " << Time_variable << endl;
       */
    //Step 1: check the state at contacts
    CheckStateContact();
    //Step 2: check if the multiple impact is terminated or not
    if (IsMulImpactTerminate()) // multiple impact terminated
    {
      if (YesSaveData) // Save the date at the end of impact
      {
        UpdateDuringImpact(); // Update state of dynamical system
        SaveDataOneStep(point_save); // Save the data
      }
      break;
    }
    // Select the primary contact
    SelectPrimaContact();
    // If the primary contact change, we have to re-estimate the step size
    if ((IsPrimaConChange) && (IsNumberOfStepsEst))
    {
      ComputeStepSize();
    }
    //Step 3: compute the vector of distributing law
    ComputeDistriVector();
    //Check if the primary contact change during computation or not
    if (number_step == 1) // for first step
    {
      Old_prima_contact = IdPrimaContact;
      IsPrimaConChange = false;
    }
    else
    {
      if (Old_prima_contact == IdPrimaContact)
      {
        IsPrimaConChange = false;
      }
      else
      {
        IsPrimaConChange = true;
        Old_prima_contact = IdPrimaContact;
      }
    }
    //cout << "Id of the primary contact: " << IdPrimaContact << endl;
    //cout << "Is the primary contact change: " << IsPrimaConChange << endl;

    //Step 4: compute the increment of normal impulse and the total impulse at contacts
    ComputeImpulseContact();
    // Step 5: compute the relative velocity at contacts
    ComputeVelContact();
    // Step 6: compute the potential energy at contacts
    ComputeEnerContact();
    //Step 7: Update the time-like variable
    ++number_step;
    ++_counterstepsave;
    ++show_progress;
    Impulse_variable = Impulse_variable + DeltaP;
    Time_variable = Time_variable + DeltaP / (*ForceContact)(IdPrimaContact);
    // Step 8: update the state of DS and output during impact and write data into output file at the beginning of each step
    if ((YesSaveData) & (_counterstepsave >= NstepSave))
    {
      UpdateDuringImpact(); // Update state of dynamical system
      SaveDataOneStep(point_save); // Save the data
      point_save++;
      _counterstepsave = 0; // reset the counter to 0
    }
    //
    if (number_step > NstepMax)
      // RuntimeException::selfThrow("In OSNSMultipleImpact::ComputeImpact, number of integration steps perfomed exceeds the maximal number of steps allowed!!!");
      break;
  };
  //
  cout << "*****************Impact computation is terminated******************" << endl;
  cout << "Number of integration steps: " << number_step << endl;
  /*
     cout << "Velocity at contacts: ";
     VelContact->display();
     cout << "Impulse at contact: ";
     TolImpulseContact->display();
     cout << "Duration of the multiple impacts process: " << Time_variable << " s" << endl;
     */
  // Close the stream file
  if (YesSaveData)
  {
    if (YesSaveByMatrix)
    {
      ioMatrix io(NameFile.c_str(), "ascii");
      io.write(*_DataMatrix, "noDim");
    }
    else
      OutputFile.close();
  }
}
//=======================================================================================
void OSNSMultipleImpact::PostComputeImpact()
{
  // === Get index set from Topology ===
  SP::InteractionsGraph indexSet = simulation()->indexSet(levelMin());
  // y and lambda vectors
  SP::SiconosVector lambda;
  SP::SiconosVector y;
  // === Loop through "active" Interactions (ie present in indexSets[1]) ===
  unsigned int pos;
  InteractionsGraph::VIterator ui, uiend;
  for (boost::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui)
  {
    SP::Interaction inter = indexSet->bundle(*ui);
    // Get the relative position of inter-interactionBlock in the vector VelContact/TolImpulseContact
    pos = _M->getPositionOfInteractionBlock(inter);
    // Get Y and Lambda for the current Interaction
    y = inter->y(levelMin());
    lambda = inter->lambda(levelMin());
    // Copy VelContact/TolImpulseContact, starting from index pos into y/lambda
    // save into y !!
    setBlock(*VelContact, y, y->size(), pos, 0);// Warning: yEquivalent is
    // saved into lambda[1] !!
    setBlock(*ImpulseContact_update, lambda, lambda->size(), pos, 0);
    // If the update is performed at the end of the impact process, we update the total normal impulse at contacts
    // from the beginning to the end of impact (vector TolImpulseContact). Otherwise, we must reset the lambda[1] to zero because
    // the post-impact velocity has been calculated during impact
    // if (!YesSaveData) // we update the impact state at the end of impact
    //   {
    //     // Copy VelContact/TolImpulseContact, starting from index pos into y/lambda
    //     // save into y !!
    //     setBlock(*VelContact, y, y->size(), pos, 0);// Warning: yEquivalent is
    //     // saved into lambda[1] !!
    //     setBlock(*TolImpulseContact, lambda, lambda->size(), pos, 0);
    //   }
    // else //
    //   lambda->zero();
  }
}
//========================================================================================
int OSNSMultipleImpact::compute(double time)
{
  // Pre-compute for impact
  PreComputeImpact();
  // solve the multiple impacts
  if ((Ncontact != 0) && IsVcminNegative()) // if there is at least one contact and the vilocity before impact is negative
  {
    ComputeImpact();
  };
  // Post-compute for multiple impacts
  PostComputeImpact();
  return  0;
}

//========================================================================================
void OSNSMultipleImpact::display() const
{
  cout << "<<<<<<<<<<<<<<<<< Information about the multiple impact >>>>>>>>>>>>>>>>>>>>>" << endl;
  cout << "Type of the contact compliance law: " << TypeCompLaw << endl;
  cout << "Number of contacts involved into impacts: " << Ncontact << endl;
  cout << "Step size used: " << DeltaP << endl;
  cout << "Primary impulse at the end of impact: " << Impulse_variable << endl;
  cout << "Duration of the multiple impacs process: " << Time_variable << endl;
};
