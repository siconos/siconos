// This is the implementation for the class OSNSMultipleImpact
//========================================================================================
#include "OSNSMultipleImpact.hpp"
#include "LagrangianDS.hpp"
#include "MultipleImpactNSL.hpp"
#include "Simulation.hpp"
#include <boost/numeric/ublas/io.hpp>
using namespace std;
//Default constructor
OSNSMultipleImpact::OSNSMultipleImpact(): LinearOSNS()
{
  TypeCompLaw = "BiStiffness";
  PowCompLaw = 1.0;
  NstepEst = 10000;
  NstepSave = 100;
  TOL_IMPACT = DEFAULT_TOL_IMPACT;
  UpdateEndImpact = true;
  YesSaveData = false;
  IsNumberOfStepsEst = true;
  NameFile = "DataMultipleImpact.dat";
}
//------------------------------ -------------------------------------------------------------
OSNSMultipleImpact::OSNSMultipleImpact(string newTypeLaw, double newPowerLaw, unsigned int newNstepEst = 10000): LinearOSNS()
{
  TypeCompLaw = newTypeLaw;
  PowCompLaw = newPowerLaw;
  NstepEst = newNstepEst;
  NstepSave = 100;
  TOL_IMPACT = DEFAULT_TOL_IMPACT;
  UpdateEndImpact = true;
  YesSaveData = false;
  IsNumberOfStepsEst = true;
  NameFile = "DataMultipleImpact.dat";
  if ((TypeCompLaw != "MonoStiffness") && (TypeCompLaw != "BiStiffness"))
    RuntimeException::selfThrow("OSNSMultipleImpact::TypeCompLaw type of the compliance model must be either MonoStiffness or BiStiffness!");
}
//------------------------------ -------------------------------------------------------------
OSNSMultipleImpact::OSNSMultipleImpact(string newTypeLaw, double newPowerLaw, double newDelP = 1.0e-5): LinearOSNS()
{
  TypeCompLaw = newTypeLaw;
  PowCompLaw = newPowerLaw;
  DeltaP = newDelP;
  NstepSave = 100;
  TOL_IMPACT = DEFAULT_TOL_IMPACT;
  UpdateEndImpact = true;
  YesSaveData = false;
  IsNumberOfStepsEst = false;
  NameFile = "DataMultipleImpact.dat";
  if ((TypeCompLaw != "MonoStiffness") && (TypeCompLaw != "BiStiffness"))
    RuntimeException::selfThrow("OSNSMultipleImpact::TypeCompLaw type of the compliance model must be either MonoStiffness or BiStiffness!");
}
//-------------------------------------------------------------------------------------------------
OSNSMultipleImpact::~OSNSMultipleImpact() {}
//--------------------------------------------------------------------------------------------------
void OSNSMultipleImpact::WriteSiconosVector(const SiconosVector& m)
{
  DenseVect*  p = m.dense();
  std::copy(p->begin(), p->end(), std::ostream_iterator<double>(OutputFile, " "));
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
void OSNSMultipleImpact::AllocateMemory()
{
  if (! VelContact)
    VelContact.reset(new SimpleVector(maxSize()));
  else
  {
    if (VelContact->size() != maxSize())
      VelContact->resize(maxSize());
  };
  //
  if (! OldVelContact)
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
  if (! WcContact)
    WcContact.reset(new SimpleVector(maxSize()));
  else
  {
    if (WcContact->size() != maxSize())
      WcContact->resize(maxSize());
  };
  //
  if (! DistriVector)
    DistriVector.reset(new SimpleVector(maxSize()));
  else
  {
    if (DistriVector->size() != maxSize())
      DistriVector->resize(maxSize());
  };
  //
  if (! StateContact)
    StateContact.reset(new IndexInt(maxSize()));
  else
  {
    if (StateContact->size() != maxSize())
      StateContact->resize(maxSize());
  };
  //
  if (! Kcontact)
    Kcontact.reset(new SimpleVector(maxSize()));
  else
  {
    if (Kcontact->size() != maxSize())
      Kcontact->resize(maxSize());
  };
  //
  if (! ResContact)
    ResContact.reset(new SimpleVector(maxSize()));
  else
  {
    if (ResContact->size() != maxSize())
      ResContact->resize(maxSize());
  };
  //
  if (! TolImpulseContact)
    TolImpulseContact.reset(new SimpleVector(maxSize()));
  else
  {
    if (TolImpulseContact->size() != maxSize())
      TolImpulseContact->resize(maxSize());
  };
  //
  if (! DelImpulseContact)
    DelImpulseContact.reset(new SimpleVector(maxSize()));
  else
  {
    if (DelImpulseContact->size() != maxSize())
      DelImpulseContact->resize(maxSize());
  };
  //
  if (! ForceContact)
    ForceContact.reset(new SimpleVector(maxSize()));
  else
  {
    if (ForceContact->size() != maxSize())
      ForceContact->resize(maxSize());
  };
}
//=====================================================================================
void OSNSMultipleImpact::BuildStiffResCofVec()
{
  SP::UnitaryRelationsGraph indexSet = simulation()->indexSet(1); // get indexSet[1]
  //Loop over the UR of the indexSet(1)
  UnitaryRelationsGraph::VIterator ui, uiend;
  for (boost::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui)
  {
    SP::UnitaryRelation ur = indexSet->bundle(*ui);
    SP::NonSmoothLaw nslaw = ur->interaction()->nslaw();
    SP::MultipleImpactNSL Mulnslaw = boost::dynamic_pointer_cast<MultipleImpactNSL>(nslaw);
    // Get the relative position of UR-unitaryBlock in the vector VelContact
    unsigned int pos = _M->getPositionOfUnitaryBlock(ur);
    (*ResContact)(pos) = Mulnslaw->ResCof();
    (*Kcontact)(pos) = Mulnslaw->Stiff();
  }
  /*
  cout << " Restitution coefficients: " << endl;
  ResContact->display();
  cout << "Stiffnesses: " << endl;
  Kcontact->display();
  */
}
//======================================================================================
void OSNSMultipleImpact::ComputeStepSize()
{
  SP::DynamicalSystemsGraph DSG = simulation()->model()->nonSmoothDynamicalSystem()->dynamicalSystems();
  SP::UnitaryRelationsGraph indexSet1 = simulation()->indexSet(1); // get indexSet[1]
  //Loop over alls DS
  DynamicalSystemsGraph::VIterator ui, uiend;
  DynamicalSystemsGraph::OEIterator edi, ediend;
  double Pest = 0.0;
  for (boost::tie(ui, uiend) = DSG->vertices(); ui != uiend; ++ui)
  {
    // Check if this DS is involved in the impact or not
    bool found = false;
    // Loop over all edges comming from this DS
    for (boost::tie(edi, ediend) = DSG->out_edges(*ui); edi != ediend; ++edi)
    {
      SP::UnitaryRelation urp = DSG->bundle(*edi);
      if (indexSet1->is_vertex(urp)) // there exist at least one edge of this DS in the indexSet[1]
      {
        found = true;
        break;
      }
    }
    //
    if (found) // if this DS involved in the impact
    {
      SP::DynamicalSystem ds = DSG->bundle(*ui); // DS
      SP::LagrangianDS Lagds = boost::dynamic_pointer_cast<LagrangianDS>(ds);
      SP::SiconosMatrix mass_ds =  Lagds->mass();    // Mass matrix of DS
      SP::SiconosVector vel_ds = Lagds->velocityMemory()->getSiconosVector(1);  // Pre-impact velocity of DS
      SP::SiconosVector abs_vel_ds(new SimpleVector(vel_ds->size()));
      abs_wise(*vel_ds, *abs_vel_ds); // get the absolute (velocity vector of ds)
      SP::SimpleVector prod_mass_vel(new SimpleVector(mass_ds->size(0)));
      prod(*mass_ds, *abs_vel_ds, *prod_mass_vel);
      Pest = Pest + prod_mass_vel->sum();
    }
  };
  DeltaP = Pest / NstepEst;
  if (DeltaP <= 0.0)
    RuntimeException::selfThrow("OSNSMultipleImpact::ComputeStepSize the step size DeltaP must be positive  !!");
  //
  //cout << "Step size:" << DeltaP << endl;
}
//========================================================================================
void OSNSMultipleImpact::PreComputeImpact()
{
  //1. Get the number of contacts and bodies involved in the impact
  if (levelMin() != 1)
    RuntimeException::selfThrow("OSNSMultipleImpact::PreComputeImpact==> the levelMin must be equal to 1 in the multiple impact model !!");
  SP::UnitaryRelationsGraph indexSet = simulation()->indexSet(levelMin()); // get indexSet[1]
  Ncontact = indexSet->size();
  //2. Compute matrix _M
  SP::Topology topology = simulation()->model()->nonSmoothDynamicalSystem()->topology();
  bool b = topology->isTimeInvariant();
  bool isLinear = simulation()->model()->nonSmoothDynamicalSystem()->isLinear();
  if (!b || !isLinear)
  {
    // Computes new _unitaryBlocks if required
    updateUnitaryBlocks();
    // Updates matrix M
    _M->fill(indexSet);
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
  //
  if (ForceContact->size() != _sizeOutput)
  {
    ForceContact->resize(_sizeOutput);
  }
  ForceContact->zero();
  //4. Initialize the relative velocity, potential energy, impulse at contacts
  InitializeInput();
  //5. Compute the step size
  if (IsNumberOfStepsEst) // We need to estimate the step size from the initial date of the dynamic system
  {
    ComputeStepSize();
  };
  //6. Build the vectors of stifnesseses and of restitution coefficients
  BuildStiffResCofVec();
  /*
  cout<<"Matrix OSNS problem:  " << endl;
   _M->display();
  */
}
//=======================================================================================
void OSNSMultipleImpact::InitializeInput()
{
  //Loop over alls UR involved in the indexSet[1]
  SP::UnitaryRelationsGraph indexSet = simulation()->indexSet(levelMin()); // get indexSet[1]
  UnitaryRelationsGraph::VIterator ui, uiend;
  for (boost::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui)
  {
    SP::UnitaryRelation ur = indexSet->bundle(*ui);
    //SP::SiconosVector Vc0 = ur->y(1); // Relative velocity at beginning of impact
    SP::SiconosVector Vc0 = ur->yOld(1); // Relative velocity at beginning of impact
    unsigned int pos_ur = _M->getPositionOfUnitaryBlock(ur);
    setBlock(*Vc0, VelContact, Vc0->size(), 0, pos_ur);
    SP::SiconosVector ener0(new SimpleVector(Vc0->size()));
    ener0->zero(); // We suppose that the initial potential energy before impact is equal to zero at any contact
    // at the beginning of impact
    setBlock(*ener0, EnerContact, ener0->size(), 0, pos_ur);
    //SP::SiconosVector impulse0= (ur->interaction())->lambda(1))->vector(ur->number());
    SP::SiconosVector impulse0(new SimpleVector(Vc0->size()));
    impulse0->zero(); // We suppose that the impulse before impact is equal to zero at any contact
    // at the beginning of impact
    setBlock(*impulse0, TolImpulseContact, impulse0->size(), 0, pos_ur);
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
  // Note that _unitaryBlocks is up to date since updateUnitaryBlocks
  // has been called during OneStepNSProblem::initialize()

  // If the topology is TimeInvariant ie if M structure does not
  // change during simulation:
  bool b = topology->isTimeInvariant();
  bool isLinear = simulation()->model()->nonSmoothDynamicalSystem()->isLinear();
  if ((b || !isLinear) &&   !interactions()->isEmpty())
  {
    updateM();
  }
  else // in that case, M will be updated during preCompute
  {
    // Default size for M = maxSize()
    if (! _M)
    {
      if (_MStorageType == 0)
        _M.reset(new OSNSMatrix(maxSize(), 0));

      else // if(_MStorageType == 1) size = number of _unitaryBlocks
        // = number of UR in the largest considered indexSet
        _M.reset(new OSNSMatrix(simulation()->indexSet(levelMin())->size(), 1));
    }
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
  for (unsigned int i = 0; i < Ncontact; ++i)
  {
    if (isZero((*EnerContact)(i)))
    {
      if ((*VelContact)(i) >= 0.0) // no impact at this contact
        (*StateContact)[i] = 0;
      else  // impact happens without potential energy
        (*StateContact)[i] = 1;
    }
    else if ((*StateContact)[i] != 2) // impact happens with not zero potential energy
    {
      (*StateContact)[i] = 2;
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
bool OSNSMultipleImpact::IsMulImpactTerminate()
{
  bool var = true;
  for (unsigned int i = 0; i < Ncontact; ++i)
  {
    if ((*StateContact)[i] != 0)
    {
      var = false;
      break;
    };
  };
  return var;
  /*
  cout << "Is the multiple impacts is terminated: " << var << endl;
  */
}
//=======================================================================================
void OSNSMultipleImpact::ComputeDistriVector()
{
  //Case 1: if no potential energy at any contact
  if (IsEnermaxZero()) // case of no potential energy at any contact
  {
    PrimConVelocity(); // Select the primary contact according to the relative velocity at contact
    for (unsigned int i = 0; i < Ncontact; ++i)
    {
      if ((*StateContact)[i] != 0) // the impact can takes place at this contact
      {
        (*DistriVector)(i) = ((*Kcontact)(i) / (*Kcontact)(IdPrimaContact)) * std::pow(((*VelContact)(i) / VelAtPrimaCon), PowCompLaw);
      }
      else
      {
        (*DistriVector)(i) = 0.0;
      }
      if ((*DistriVector)(i) < 0.0)
        RuntimeException::selfThrow("OSNSMultipleImpact::ComputeDistriVector the component of DistriVector must be positive !!");
    };
  }
  //Case 2: if the potential energy is not zero at contacts
  else
  {
    PrimConEnergy(); // Select the primary contact according to the potential energy at contacts
    for (unsigned int i = 0; i < Ncontact; ++i)
    {
      if ((*StateContact)[i] == 1) // no potential energy at this contact, including the contacts at which impact repeats
      {
        if ((*VelContact)(i) > 0)
        {
          RuntimeException::selfThrow("OSNSMultipleImpact::ComputeDistriVector, the pre-impact velocity must be negative!!");
        }
        else
        {
          (*DistriVector)(i) = ((*Kcontact)(i) / (*Kcontact)(IdPrimaContact)) * std::pow(((std::fabs((*VelContact)(i)) * DeltaP) / EnerAtPrimaCon), PowCompLaw);
        }
      }
      else if ((*StateContact)[i] == 2) // potential is not zero at this contact
      {
        (*DistriVector)(i) = std::pow((*Kcontact)(i) / (*Kcontact)(IdPrimaContact), 1.0 / (1.0 + PowCompLaw)) * std::pow((*EnerContact)(i) / EnerAtPrimaCon, PowCompLaw / (PowCompLaw + 1.0));
      }
      else // no impact at this contact
      {
        (*DistriVector)(i) = 0.0;
      };
      if ((*DistriVector)(i) < 0.0)
        RuntimeException::selfThrow("OSNSMultipleImpact::ComputeDistriVector the component of DistriVector must be positive !!");
    };
  };
  /*
  cout << "Tensor of distributing rule:" << endl;
  DistriVector->display();
  */

}
//=======================================================================================
void OSNSMultipleImpact::ComputeImpulseContact()
{
  (*DelImpulseContact) = (*DistriVector) * DeltaP;
  (*TolImpulseContact) = (*TolImpulseContact) + (*DelImpulseContact);
  // Compute the contact force
  for (unsigned i = 0; i < Ncontact; ++i)
  {
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
  //1. Copy VelContact/DelImpulseContact into the vector y/lambda for unitary relations
  SP::UnitaryRelationsGraph indexSet = simulation()->indexSet(levelMin());
  // y and lambda vectors
  SP::SiconosVector lambda;
  SP::SiconosVector y;
  // === Loop through "active" Unitary Relations (ie present in indexSets[1]) ===
  unsigned int pos;
  UnitaryRelationsGraph::VIterator ui, uiend;
  for (boost::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui)
  {
    SP::UnitaryRelation ur = indexSet->bundle(*ui);
    // Get the relative position of UR-unitaryBlock in the vector VelContact/TolImpulseContact
    pos = _M->getPositionOfUnitaryBlock(ur);
    // Get Y and Lambda for the current Unitary Relation
    y = ur->y(levelMin());
    lambda = ur->lambda(levelMin());
    // Copy VelContact/TolImpulseContact, starting from index pos into y/lambda
    // save into y !!
    setBlock(*VelContact, y, y->size(), pos, 0);// Warning: yEquivalent is
    // saved into lambda[1] !!
    setBlock(*DelImpulseContact, lambda, lambda->size(), pos, 0);
  };
  //2. Update the Input[1], state of DS systems, Output[1]
  simulation()->update(levelMin());
}
//--------------------------------------------------------------------------------------
void OSNSMultipleImpact::SaveDataOneStep()
{
  // Save the total impulse at the primary contacts (time-like independent variable) and the time evolution during impact
  OutputFile << Impulse_variable << " ";
  OutputFile << Time_variable << " ";
  // Save the data related to UnitaryRelations
  SP::UnitaryRelationsGraph indexSet0 = simulation()->indexSet(0);
  SP::UnitaryRelationsGraph indexSet1 = simulation()->indexSet(levelMin());
  unsigned int pos;
  int Status_ur;
  UnitaryRelationsGraph::VIterator ui, uiend;
  for (boost::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui)
  {
    SP::UnitaryRelation ur = indexSet0->bundle(*ui); // UR
    SP::SiconosVector ydot = ur->y(1);
    SP::SiconosVector P_ur(new SimpleVector(ur->interaction()->nonSmoothLaw()->size()));
    SP::SiconosVector F_ur(new SimpleVector(ur->interaction()->nonSmoothLaw()->size()));
    SP::SiconosVector E_ur(new SimpleVector(ur->interaction()->nonSmoothLaw()->size()));

    if (indexSet1->is_vertex(ur)) // if UR belongs to the IndexSet[1]
    {
      pos = _M->getPositionOfUnitaryBlock(ur);
      Status_ur = (*StateContact)[pos]; //status of UR
      setBlock(*TolImpulseContact, P_ur, P_ur->size(), pos, 0);
      setBlock(*ForceContact, F_ur, F_ur->size(), pos, 0);
      setBlock(*EnerContact, E_ur, E_ur->size(), pos, 0);
    }
    else
    {
      Status_ur = -1; // no contact at this UR
      P_ur->zero();   // no impulse at this UR
      F_ur->zero();   // no force at this UR
      E_ur->zero();   // no potential at this UR
    };
    //
    OutputFile << ur->interaction()->number() << " "; // Write the Id of the its parent interaction
    OutputFile << Status_ur << " "; // Write the status of UR
    WriteSiconosVector(*ydot); // Write the relative velocity at the UR
    WriteSiconosVector(*P_ur); // Write the total impulse at the UR
    WriteSiconosVector(*F_ur); // Write the force at the UR
    WriteSiconosVector(*E_ur); // Write the potential energy at the UR
  }
  // Save the data related to DS
  SP::DynamicalSystemsGraph DSG = simulation()->model()->nonSmoothDynamicalSystem()->dynamicalSystems();
  DynamicalSystemsGraph::VIterator dsi, dsiend;
  for (boost::tie(dsi, dsiend) = DSG->vertices(); dsi != dsiend; ++dsi)
  {
    SP::DynamicalSystem ds = DSG->bundle(*dsi); // DS
    SP::LagrangianDS Lagds = boost::dynamic_pointer_cast<LagrangianDS>(ds);
    SP::SiconosVector qdot = Lagds->velocity();
    // Write
    OutputFile << ds->number() << " "; // number of DS
    WriteSiconosVector(*qdot);
  }
  // Terminate the line
  OutputFile << "\n";
}
//=======================================================================================
void OSNSMultipleImpact::ComputeImpact()
{
  // open the stream file
  if (YesSaveData)
  {
    OutputFile.open(NameFile.c_str(), ios::out | ios::trunc);
    if (!OutputFile.is_open())
      SiconosVectorException::selfThrow("OSNSMultipleImpact::ComputeImpact==> write error : Fail to open \"" + NameFile + "\"");
    OutputFile.precision(15);
  }
  //
  Impulse_variable = 0.0;
  Time_variable = 0.0;
  unsigned int number_step = 1;
  unsigned int _counterstepsave = NstepSave - 1;
  //
  /*
  cout << "----------Before multiple impacts computation---------------" << endl;
  cout << "Velocity at contacts: ";
  VelContact->display();
  cout << "Impulse at contact: ";
  TolImpulseContact->display();
  */
  //cout << "-------------------Multiple impacts computation starts:-----------------------" << endl;

  while (1 != 0)
  {
    //
    /*
    cout << "==================Step==================:  " << number_step << endl;
    cout << "Impulse variable: " << Impulse_variable << endl;
    cout << "Time_variable: " << Time_variable << endl;
    */
    //
    //Step 1: Write data into output file at the beginning of each step
    if ((YesSaveData) && (!UpdateEndImpact))
    {
      ++_counterstepsave;
      if (_counterstepsave >= NstepSave)
      {
        OutputFile << number_step << " "; // number of calculation step
        SaveDataOneStep(); // Save the date
        _counterstepsave = 0; // reset the counter to 0
      };
    }
    //Step 2: check the state at contacts
    CheckStateContact();
    //Step 3: check if the multiple impact is terminated or not
    if (IsMulImpactTerminate()) // multiple impact terminated
      break;
    //Step 4: compute the vector of distributing law
    ComputeDistriVector();
    //Step 5: compute the increment of normal impulse and the total impulse at contacts
    ComputeImpulseContact();
    // Step 6: compute the relative velocity at contacts
    ComputeVelContact();
    // Step 7: compute the potential energy at contacts
    ComputeEnerContact();
    // Step 8: update the state of DS and output during impact
    if (!UpdateEndImpact)
    {
      UpdateDuringImpact();
    };
    //Step 9: Update the time-like variable
    Impulse_variable = Impulse_variable + DeltaP;
    Time_variable = Time_variable + DeltaP / (*ForceContact)(IdPrimaContact);
    ++number_step;
  };
  //
  /*
  cout << "----------After multiple impacts computation---------------" << endl;
  cout << "Number of iterations: " << number_step << endl;
  cout << "Velocity at contacts: ";
  VelContact->display();
  cout << "Impulse at contact: ";
  TolImpulseContact->display();
  cout << "Duration of the multiple impacts process: " << Time_variable << " s" << endl;
  */
  // Close the stream file
  if (YesSaveData)
  {
    OutputFile.close();
  }
}
//=======================================================================================
void OSNSMultipleImpact::PostComputeImpact()
{
  // === Get index set from Topology ===
  SP::UnitaryRelationsGraph indexSet = simulation()->indexSet(levelMin());
  // y and lambda vectors
  SP::SiconosVector lambda;
  SP::SiconosVector y;
  // === Loop through "active" Unitary Relations (ie present in indexSets[1]) ===
  unsigned int pos;
  UnitaryRelationsGraph::VIterator ui, uiend;
  for (boost::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui)
  {
    SP::UnitaryRelation ur = indexSet->bundle(*ui);
    // Get the relative position of UR-unitaryBlock in the vector VelContact/TolImpulseContact
    pos = _M->getPositionOfUnitaryBlock(ur);
    // Get Y and Lambda for the current Unitary Relation
    y = ur->y(levelMin());
    lambda = ur->lambda(levelMin());
    // If the update is performed at the end of the impact process, we update the total normal impulse at contacts
    // from the beginning to the end of impact (vector TolImpulseContact). Otherwise, we must reset the lambda[1] to zero because
    // the post-impact velocity has been calculated during impact
    if (UpdateEndImpact)
    {
      // Copy VelContact/TolImpulseContact, starting from index pos into y/lambda
      // save into y !!
      setBlock(*VelContact, y, y->size(), pos, 0);// Warning: yEquivalent is
      // saved into lambda[1] !!
      setBlock(*TolImpulseContact, lambda, lambda->size(), pos, 0);
    }
    else //
      lambda->zero();
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
  cout << "Power of the compliance law: " << PowCompLaw << endl;
  cout << "Number of contacts involved into impacts: " << Ncontact << endl;
  cout << "Step size used: " << DeltaP << endl;
  cout << "Primary impulse at the end of impact: " << Impulse_variable << endl;
  cout << "Duration of the multiple impacs process: " << Time_variable << endl;
};
