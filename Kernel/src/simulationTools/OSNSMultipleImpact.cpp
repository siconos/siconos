
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
  _Tol_Vel = DEFAULT_TOL_VEL;
  _Tol_Ener = DEFAULT_TOL_ENER;
  _ZeroVel_EndIm = DEFAULT_TOL_VEL;
  _ZeroEner_EndIm = DEFAULT_TOL_ENER;
  YesSaveData = false;
  IsNumberOfStepsEst = true;
  SizeDataSave = 1000;
  NstepMax = 100000;
  Step_min_save = 1;
  Step_max_save = NstepMax;
  NameFile = "DataMultipleImpact.dat";
}
//------------------------------ -------------------------------------------------------------
OSNSMultipleImpact::OSNSMultipleImpact(string newTypeLaw, unsigned int newNstepEst = 10000): LinearOSNS()
{
  TypeCompLaw = newTypeLaw;
  NstepEst = newNstepEst;
  NstepSave = 100;
  TOL_IMPACT = DEFAULT_TOL_IMPACT;
  _Tol_Vel = DEFAULT_TOL_VEL;
  _Tol_Ener = DEFAULT_TOL_ENER;
  _ZeroVel_EndIm = DEFAULT_TOL_VEL;
  _ZeroEner_EndIm = DEFAULT_TOL_ENER;
  YesSaveData = false;
  IsNumberOfStepsEst = true;
  NameFile = "DataMultipleImpact.dat";
  SizeDataSave = 1000;
  NstepMax = 100000;
  Step_min_save = 1;
  Step_max_save = NstepMax;
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
  _Tol_Vel = DEFAULT_TOL_VEL;
  _Tol_Ener = DEFAULT_TOL_ENER;
  _ZeroVel_EndIm = DEFAULT_TOL_VEL;
  _ZeroEner_EndIm = DEFAULT_TOL_ENER;
  YesSaveData = false;
  IsNumberOfStepsEst = false;
  NameFile = "DataMultipleImpact.dat";
  SizeDataSave = 1000;
  NstepMax = 100000;
  Step_min_save = 1;
  Step_max_save = NstepMax;
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

void OSNSMultipleImpact::SetTolVel(double _var)
{
  _Tol_Vel = _var;
};

void OSNSMultipleImpact::SetTolEner(double _var)
{
  _Tol_Ener = _var;
};

void OSNSMultipleImpact::SetZeroVelEndImp(double _var)
{
  _ZeroVel_EndIm = _var;
};

void OSNSMultipleImpact::SetZeroEnerEndImp(double _var)
{
  _ZeroEner_EndIm = _var;
};

void OSNSMultipleImpact::SetNstepSave(unsigned int var)
{
  NstepSave = var;
};

void OSNSMultipleImpact::SetNstepMax(unsigned int var)
{
  NstepMax = var;
};

void OSNSMultipleImpact::SetStepMinMaxSave(unsigned int var1, unsigned int var2)
{
  Step_min_save = var1;
  Step_max_save = var2;
}

void OSNSMultipleImpact::setTypeCompLaw(std::string newTypeLaw)
{
  TypeCompLaw = newTypeLaw;
  if ((TypeCompLaw != "MonoStiffness") && (TypeCompLaw != "BiStiffness"))
    RuntimeException::selfThrow("OSNSMultipleImpact::TypeCompLaw type of the compliance model must be either MonoStiffness or BiStiffness!");
};

void OSNSMultipleImpact::SetSizeDataSave(unsigned int var)
{
  SizeDataSave = var;
}
//---------------------------------------------------------------------------------------------------
void OSNSMultipleImpact::WriteVectorIntoMatrix(const SiconosVector m, const unsigned int pos_row, const unsigned int pos_col)
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
//------------------------------------------------------------------------------------------------
bool OSNSMultipleImpact::isVelNegative(const double Var)
{
  if (Var < - _Tol_Vel)
    return true;
  else
    return false;
}
//-------------------------------------------------------------------------------------------------

bool OSNSMultipleImpact::isEnerZero(const double Var)
{
  if (std::abs(Var) <= _Tol_Ener)
    return true;
  else
    return false;
}
//--------------------------------------------------------------------------------------------------
unsigned int OSNSMultipleImpact::EstimateNdataCols()
{
  unsigned int _numberCols = 1;
  // Number of columns for data at contacts
  SP::InteractionsGraph indexSet = simulation()->indexSet(0); // get indexSet[0]
  InteractionsGraph::VIterator ui, uiend;
  for (std11::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui)
  {
    //_numberCols = _numberCols + 3*(indexSet->bundle(*ui)->interaction()->nonSmoothLaw()->size()) + 1;
    _numberCols = _numberCols + (indexSet->bundle(*ui)->getSizeOfY());
  }
  // Number of columns for data at particles
  SP::DynamicalSystemsGraph DSG = simulation()->model()->nonSmoothDynamicalSystem()->dynamicalSystems();
  DynamicalSystemsGraph::VIterator dsi, dsiend;
  for (std11::tie(dsi, dsiend) = DSG->vertices(); dsi != dsiend; ++dsi)
  {
    _numberCols = _numberCols + (DSG->bundle(*dsi)->getDim());
  }
  return(_numberCols);
}
//-----------------------------------------------------------------------------------------------
void OSNSMultipleImpact::AllocateMemory()
{
  if (!VelContact)
    VelContact.reset(new SiconosVector(maxSize()));
  else
  {
    if (VelContact->size() != maxSize())
      VelContact->resize(maxSize());
  };
  //
  if (!OldVelContact)
    OldVelContact.reset(new SiconosVector(maxSize()));
  else
  {
    if (OldVelContact->size() != maxSize())
      OldVelContact->resize(maxSize());
  };
  //
  if (! EnerContact)
    EnerContact.reset(new SiconosVector(maxSize()));
  else
  {
    if (EnerContact->size() != maxSize())
      EnerContact->resize(maxSize());
  };
  //
  if (!WcContact)
    WcContact.reset(new SiconosVector(maxSize()));
  else
  {
    if (WcContact->size() != maxSize())
      WcContact->resize(maxSize());
  };
  //
  if (!DistriVector)
    DistriVector.reset(new SiconosVector(maxSize()));
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
    Kcontact.reset(new SiconosVector(maxSize()));
  else
  {
    if (Kcontact->size() != maxSize())
      Kcontact->resize(maxSize());
  };
  //
  if (!ResContact)
    ResContact.reset(new SiconosVector(maxSize()));
  else
  {
    if (ResContact->size() != maxSize())
      ResContact->resize(maxSize());
  };
  //
  if (!ElasCoefContact)
    ElasCoefContact.reset(new SiconosVector(maxSize()));
  else
  {
    if (ElasCoefContact->size() != maxSize())
      ElasCoefContact->resize(maxSize());
  };
  if (!TolImpulseContact)
    TolImpulseContact.reset(new SiconosVector(maxSize()));
  else
  {
    if (TolImpulseContact->size() != maxSize())
      TolImpulseContact->resize(maxSize());
  };
  //
  if (!DelImpulseContact)
    DelImpulseContact.reset(new SiconosVector(maxSize()));
  else
  {
    if (DelImpulseContact->size() != maxSize())
      DelImpulseContact->resize(maxSize());
  };
  //
  if (!ImpulseContact_update)
    ImpulseContact_update.reset(new SiconosVector(maxSize()));
  else
  {
    if (ImpulseContact_update->size() != maxSize())
      ImpulseContact_update->resize(maxSize());
  }
  //
  if (!ForceContact)
    ForceContact.reset(new SiconosVector(maxSize()));
  else
  {
    if (ForceContact->size() != maxSize())
      ForceContact->resize(maxSize());
  };
  // for the data matrix
  unsigned int _numberCols = EstimateNdataCols();
  if (!_DataMatrix)
    _DataMatrix.reset(new SimpleMatrix(SizeDataSave, _numberCols));
  else
  {
    if ((_DataMatrix->size(0) != SizeDataSave) || (_DataMatrix->size(1) != _numberCols))
      _DataMatrix->resize(SizeDataSave, _numberCols);
  }
}
//=====================================================================================
void OSNSMultipleImpact::BuildParaContact()
{
  SP::InteractionsGraph indexSet = simulation()->indexSet(1); // get indexSet[1]
  //Loop over the Interactionof the indexSet(1)
  InteractionsGraph::VIterator ui, uiend;
  for (std11::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui)
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
  //Loop over alls DS
  DynamicalSystemsGraph::VIterator ui, uiend;
  DynamicalSystemsGraph::OEIterator edi, ediend;
  double Pest = 0.0;
  for (std11::tie(ui, uiend) = DSG->vertices(); ui != uiend; ++ui)
  {
    // Check if this DS is involved in the impact or not
    bool found = false;
    // Loop over all edges comming from this DS
    for (boost::tie(edi, ediend) = DSG->out_edges(*ui); edi != ediend; ++edi)
    {
      SP::Interaction urp = DSG->bundle(*edi);
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
      SP::SiconosVector abs_vel_ds(new SiconosVector(vel_ds->size()));
      abs_wise(*vel_ds, *abs_vel_ds); // get the absolute (velocity vector of ds)
      SP::SiconosVector prod_mass_vel(new SiconosVector(mass_ds->size(0)));
      prod(*mass_ds, *abs_vel_ds, *prod_mass_vel);
      Pest = Pest + prod_mass_vel->sum();
    }
  };
  DeltaP = Pest / NstepEst;
  //cout << "Size of setOfDS" << setOfDS.size() << endl;
  //cout << "Step size is: " << DeltaP << endl;
  if (DeltaP <= 0.0)
    RuntimeException::selfThrow("OSNSMultipleImpact::ComputeStepSize the step size DeltaP must be positive  !!");
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
    // Computes new _unitaryBlocks if required
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
  //5. Compute the step size
  if (IsNumberOfStepsEst) // We need to estimate the step size from the initial date of the dynamic system
  {
    ComputeStepSize();
  };
  //6. Build the vectors of stifnesseses, of restitution coefficients, and of elaticity coefficients
  BuildParaContact();
}
//=======================================================================================
void OSNSMultipleImpact::InitializeInput()
{
  //Loop over alls Interactioninvolved in the indexSet[1]
  SP::InteractionsGraph indexSet = simulation()->indexSet(levelMin()); // get indexSet[1]
  InteractionsGraph::VIterator ui, uiend;
  for (std11::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui)
  {
    SP::Interaction inter = indexSet->bundle(*ui);
    //SP::SiconosVector Vc0 = inter->y(1); // Relative velocity at beginning of impact
    SP::SiconosVector Vc0 = inter->yOld(1); // Relative velocity at beginning of impact
    unsigned int pos_inter = _M->getPositionOfInteractionBlock(inter);
    setBlock(*Vc0, VelContact, Vc0->size(), 0, pos_inter);
    SP::SiconosVector ener0(new SiconosVector(Vc0->size()));
    ener0->zero(); // We suppose that the initial potential energy before impact is equal to zero at any contact
    // at the beginning of impact
    setBlock(*ener0, EnerContact, ener0->size(), 0, pos_inter);
    //SP::SiconosVector impulse0= (inter)->lambda(1))->vector(inter->number());
    SP::SiconosVector impulse0(new SiconosVector(Vc0->size()));
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
  if (!isVelNegative(VelAtPrimaCon))
  {
    RuntimeException::selfThrow("OSNSMultipleImpact::PrimConVelocity, the velocity at the primary contact must be negative !!");
  }
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
  if (EnerAtPrimaCon < 0.0)
  {
    RuntimeException::selfThrow("OSNSMultipleImpact::PrimConEnergy the potential energy at the primary contact must be positive !!");
  }
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
  if (isEnerZero(MaxEner))
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
  if (isVelNegative(MinVelCon))
    return true;
  else
    return false;
}
//=======================================================================================
void OSNSMultipleImpact::CheckStateContact()
{
  for (unsigned int i = 0; i < Ncontact; ++i)
  {
    if (isEnerZero((*EnerContact)(i))) // potential energy is zero
    {
      if (!isVelNegative((*VelContact)(i))) // relative velocity is positive or equal to zero
        (*StateContact)[i] = 0; // no impact at this contact
      else  // impact happens without potential energy
      {
        (*StateContact)[i] = 1;
      }
    }
    else // impact happens with not zero potential energy
    {
      if ((*StateContact)[i] != 2)
      {
        (*StateContact)[i] = 2;
      }
    }
  }
}
//=======================================================================================
bool OSNSMultipleImpact::IsMulImpactTerminate()
{
  _IsImpactEnd = true;
  for (unsigned int i = 0; i < Ncontact; ++i)
  {
    if (((*EnerContact)(i) > _ZeroEner_EndIm) || ((*VelContact)(i) < -_ZeroVel_EndIm)) // if potential energy is not equal to zero or the relative velocity is negative
    {
      _IsImpactEnd = false;
    }
  }
  return _IsImpactEnd;
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
  //
  //cout << "Is the multiple impacts is terminated: " << _IsImpactEnd << endl;
  //
}
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
        _ratio_mu = (std::pow(_mu + 1.0, (_mu / (_mu + 1.0)))) / (std::pow(mu_prima + 1.0, (mu_prima / (mu_prima + 1.0))));
        ratio_stiff = (std::pow(_stiff, (1.0 / (1.0 + _mu)))) / (std::pow(stiff_prima, (1.0 / (1.0 + mu_prima))));
        if (!isVelNegative(_vel))
        {
          RuntimeException::selfThrow("OSNSMultipleImpact::ComputeDistriVector, the relative velocity when particle starts to impact must be negative!!");
        }

        ratio_vel = (std::pow(std::fabs(_vel), (_mu / (_mu + 1.0)))) / (std::pow(std::fabs(VelAtPrimaCon), (mu_prima / (1.0 + mu_prima))));
        (*DistriVector)(i) = std::pow((_ratio_mu * ratio_stiff * ratio_vel), (1.0 + _mu)) * std::pow(DeltaP, ((_mu - mu_prima) / (1.0 + mu_prima)));
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
      //
      _mu = (*ElasCoefContact)(i);
      _stiff = (*Kcontact)(i);
      _ratio_mu = (std::pow(_mu + 1.0, (_mu / (_mu + 1.0)))) / (std::pow(mu_prima + 1.0, (mu_prima / (mu_prima + 1.0))));
      ratio_stiff = (std::pow(_stiff, (1.0 / (1.0 + _mu)))) / (std::pow(stiff_prima, (1.0 / (1.0 + mu_prima))));
      if ((*StateContact)[i] == 1) // no potential energy at this contact, including the contacts at which impact repeats
      {
        if (!isVelNegative((*VelContact)(i)))
        {
          RuntimeException::selfThrow("OSNSMultipleImpact::ComputeDistriVector, the pre-impact velocity must be negative!!");
        }
        else
        {
          _vel = (*VelContact)(i);
          ratio_ener = (std::pow(std::fabs(_vel * DeltaP), (_mu / (_mu + 1.0)))) / (std::pow(EnerAtPrimaCon, (mu_prima / (mu_prima + 1.0))));
          //
          // cout << "_ratio_m: " << _ratio_mu << endl;
          // cout << "Stiff: " << _stiff << endl;
          // cout << "ratio_stiff: " << ratio_stiff << endl;
          // cout << "energy ratio: " << ratio_ener << endl;

          //
          (*DistriVector)(i) = std::pow((_ratio_mu * ratio_stiff * ratio_ener), (1.0 + _mu));
        }
      }
      else if ((*StateContact)[i] == 2) // potential is not zero at this contact
      {
        _energy = (*EnerContact)(i); // Potential energy at the current contact
        ratio_ener = (std::pow(_energy, (_mu / (_mu + 1.0)))) / (std::pow(EnerAtPrimaCon, (mu_prima / (mu_prima + 1.0))));
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
}
//=======================================================================================
void OSNSMultipleImpact::ComputeImpulseContact()
{
  (*DelImpulseContact) = (*DistriVector) * DeltaP;
  (*TolImpulseContact) = (*TolImpulseContact) + (*DelImpulseContact);
  (*ImpulseContact_update) = (*ImpulseContact_update) + (*DelImpulseContact);
  // Compute the contact force
  double PowCompLaw;
  for (unsigned int i = 0; i < Ncontact; ++i)
  {
    PowCompLaw = (*ElasCoefContact)(i);
    if (isEnerZero((*EnerContact)(i))) // if potential energy at this contact is zero
    {
      if (isVelNegative((*VelContact)(i))) // if the relative velocity at contact is negative
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
    if ((*ForceContact)(i) < 0.0)
    {
      RuntimeException::selfThrow("OSNSMultipleImpact::ComputeImpulseContact, the contact force must be positive or equal to zero!!!");
    }
  };
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
    for (unsigned int i = 0; i < Ncontact; ++i)
    {
      if ((0.5 * ((*OldVelContact)(i) + (*VelContact)(i))) <= 0.0) // Contact located in the compression phase
      {
        (*EnerContact)(i) = (*EnerContact)(i) - 0.5 * ((*OldVelContact)(i) + (*VelContact)(i)) * ((*DelImpulseContact)(i));
      }
      else                       // Contact located in the expansion phase
      {
        if (!isZero((*ResContact)(i)))
        {
          (*EnerContact)(i) = (*EnerContact)(i) - (1.0 / std::pow((*ResContact)(i), 2)) * 0.5 * ((*OldVelContact)(i) + (*VelContact)(i)) * ((*DelImpulseContact)(i));
          //
          if ((*EnerContact)(i) <  0.0)
          {
            (*EnerContact)(i) = 0.0;
          };
        }
        else // restitution coefficient equal to zero
        {
          (*EnerContact)(i) = 0.0; // In this case, no potential energy at contacts when the contact is located in the compression phase
        }
      };
      //
      if ((*EnerContact)(i) <  0.0)
      {
        RuntimeException::selfThrow("OSNSMultipleImpact::ComputeEnerContact, the potential energy during compression phase must be positive!!!");
      };
    };
  }
  else
    // For the mono-stiffness model
  {
    for (unsigned int i = 0; i < Ncontact; ++i)
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
  /*

    cout << "Potential energy at contacts at the end of step:" << endl;
    EnerContact->display();
    cout << "Work done during the compression phase at contacts" << endl;
    WcContact->display();

  */
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
  for (std11::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui)
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
  if (_ithPoint >= _DataMatrix->size(0))
    RuntimeException::selfThrow("In OSNSMultipleImpact::ComputeImpact, number of points saved exceeds the size of matrix allocated!!!");
  (*_DataMatrix)(_ithPoint, 0) = Time_variable;
  //(*_DataMatrix)(_ithPoint,0) = Impulse_variable;
  // Save the data related to UnitaryRelations
  SP::InteractionsGraph indexSet0 = simulation()->indexSet(0);
  SP::InteractionsGraph indexSet1 = simulation()->indexSet(levelMin());
  unsigned int pos;
  InteractionsGraph::VIterator ui, uiend;
  unsigned int col_pos = 1;
  for (std11::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui)
  {
    SP::Interaction ur = indexSet0->bundle(*ui); // UR
    SP::SiconosVector ydot = ur->y(1);
    SP::SiconosVector P_ur(new SiconosVector(ur->getSizeOfY()));
    SP::SiconosVector F_ur(new SiconosVector(ur->getSizeOfY()));
    SP::SiconosVector E_ur(new SiconosVector(1));
    if (indexSet1->is_vertex(ur)) // if UR belongs to the IndexSet[1]
    {
      pos = _M->getPositionOfInteractionBlock(ur);
      setBlock(*TolImpulseContact, P_ur, P_ur->size(), pos, 0);
      setBlock(*ForceContact, F_ur, F_ur->size(), pos, 0);
      setBlock(*EnerContact, E_ur, E_ur->size(), pos, 0);
    }
    else
    {
      P_ur->zero();   // no impulse at this UR
      F_ur->zero();   // no force at this UR
      E_ur->zero();   // no potential at this UR
    };
    // Write the force at the UR
    WriteVectorIntoMatrix(*F_ur, _ithPoint, col_pos);
    col_pos = col_pos + F_ur->size();
  } // Save the data related to DS
  SP::DynamicalSystemsGraph DSG = simulation()->model()->nonSmoothDynamicalSystem()->dynamicalSystems();
  DynamicalSystemsGraph::VIterator dsi, dsiend;
  for (std11::tie(dsi, dsiend) = DSG->vertices(); dsi != dsiend; ++dsi)
  {
    SP::DynamicalSystem ds = DSG->bundle(*dsi); // DS
    SP::LagrangianDS Lagds = boost::dynamic_pointer_cast<LagrangianDS>(ds);
    SP::SiconosVector qdot = Lagds->velocity();
    // Write

    WriteVectorIntoMatrix(*qdot, _ithPoint, col_pos);
    col_pos = col_pos + qdot->size();
  }
}
//=======================================================================================
void OSNSMultipleImpact::ComputeImpact()
{
  Impulse_variable = 0.0;
  Time_variable = 0.0;
  unsigned int number_step = 1;
  unsigned int point_save = 0;
  unsigned int _counterstepsave = 0;
  // Show computation progress
  //cout << "*********** Impact computation progress *************" << endl;
  //boost::progress_display show_progress(NstepMax);
  /*
     cout << "----------Before multiple impacts computation---------------" << endl;
     cout << "Velocity at contacts: ";
     VelContact->display();
     cout << "Impulse at contact: ";
     TolImpulseContact->display();
  */
  //cout << "-------------------Multiple impacts computation starts:-----------------------" << endl;
  // First save at the beginning of impact computation
  if ((YesSaveData) && (Step_min_save == 1))
  {
    SaveDataOneStep(point_save); // Save the data
    point_save++;
  }
  //
  while (1 != 0)
  {
    // cout << "==================Step==================:  " << number_step << endl;
    // cout << "Impulse variable: " << Impulse_variable << endl;
    // cout << "Time_variable: " << Time_variable << endl;

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
    //Step 3: compute the vector of distributing law
    ComputeDistriVector();
    //Step 4: compute the increment of normal impulse and the total impulse at contacts
    ComputeImpulseContact();
    // Step 5: compute the relative velocity at contacts
    ComputeVelContact();
    // Step 6: compute the potential energy at contacts
    ComputeEnerContact();
    //Step 7: Update the time-like variable
    ++number_step;
    ++_counterstepsave;
    //++show_progress;
    Impulse_variable = Impulse_variable + DeltaP;
    Time_variable = Time_variable + DeltaP / (*ForceContact)(IdPrimaContact);
    // Step 8: update the state of DS and output during impact and write data into output file at the beginning of each step
    if ((YesSaveData) & (_counterstepsave >= NstepSave))
    {
      if ((number_step >= Step_min_save) && (number_step <= Step_max_save))
      {
        UpdateDuringImpact(); // Update state of dynamical system
        SaveDataOneStep(point_save); // Save the data
        point_save++;
        _counterstepsave = 0; // reset the counter to 0
      }
    }
    //
    if (number_step > NstepMax)
    {
      RuntimeException::selfThrow("In OSNSMultipleImpact::ComputeImpact, number of integration steps perfomed exceeds the maximal number of steps allowed!!!");
      //cout << "Causion: so long computation, the computation is stopped even when the impact is not yet terminated!!! " << endl;
      // break;
    }
    // cout << "Distribution vector: ";
    // DistriVector->display();
    // cout << "Incremental Impulse: ";
    // DelImpulseContact->display();
    // cout << "Impulse at contact: ";
    // TolImpulseContact->display();
    // cout << "Velocity at contacts: ";
    // VelContact->display();
    // cout << "Potential energy at contacts: ";
    // EnerContact->display();

  }

  //
  // cout << "*****************Impact computation is terminated******************" << endl;
  // cout << "Number of integration steps: " << number_step << endl;

  // cout << "Velocity at contacts: ";
  // VelContact->display();
  // cout << "Impulse at contact: ";
  // TolImpulseContact->display();
  // cout << "Duration of the multiple impacts process: " << Time_variable << " s" << endl;

  // Close the stream file
  if (YesSaveData)
  {
    ioMatrix::write(NameFile.c_str(), "ascii", *_DataMatrix, "noDim");
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
  for (std11::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui)
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
