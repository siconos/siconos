/* Siconos-Kernel, Copyright INRIA 2005-2012.
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
#include "OSNSMultipleImpact.hpp"
#include "LagrangianDS.hpp"
#include "MultipleImpactNSL.hpp"
#include "Simulation.hpp"
#include "ioMatrix.hpp"
#include <boost/numeric/ublas/io.hpp>
#include <boost/progress.hpp>
#include "OSNSMatrix.hpp"
#include "Model.hpp"
#include "NonSmoothDynamicalSystem.hpp"

//Default constructor
OSNSMultipleImpact::OSNSMultipleImpact(): LinearOSNS()
{
  TypeCompLaw = "BiStiffness";
  NstepSave = 100;
  TOL_IMPACT = DEFAULT_TOL_IMPACT;
  _Tol_Vel = DEFAULT_TOL_VEL;
  _Tol_Ener = DEFAULT_TOL_ENER;
  _ZeroVel_EndIm = DEFAULT_TOL_VEL;
  _ZeroEner_EndIm = DEFAULT_TOL_ENER;
  YesSaveData = false;
  SizeDataSave = 1000;
  NstepMax = 100000;
  Step_min_save = 1;
  Step_max_save = NstepMax;
  NameFile = "DataMultipleImpact.dat";
}
//------------------------------ -------------------------------------------------------------
OSNSMultipleImpact::OSNSMultipleImpact(std::string newTypeLaw, double newDelP = 1.0e-5): LinearOSNS()
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
bool OSNSMultipleImpact::isZero(double Var)
{
  if (std::abs(Var) <= TOL_IMPACT)
    return true;
  else
    return false;
}
//------------------------------------------------------------------------------------------------
bool OSNSMultipleImpact::isVelNegative(double Var)
{
  if (Var < - _Tol_Vel)
    return true;
  else
    return false;
}
//-------------------------------------------------------------------------------------------------

bool OSNSMultipleImpact::isEnerZero(double Var)
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
  SP::DynamicalSystemsGraph DSG = simulation()->nonSmoothDynamicalSystem()->dynamicalSystems();
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
    SP::MultipleImpactNSL Mulnslaw = std11::dynamic_pointer_cast<MultipleImpactNSL>(nslaw);
    assert(Mulnslaw && "In OSNSMultipleImpact::BuildStiffResCofVec, non-smooth law used must be MultipleImpactNSL!!!");
    // Get the relative position of inter-interactionBlock in the vector VelContact
    unsigned int pos = _M->getPositionOfInteractionBlock(*inter);
    (*ResContact)(pos) = Mulnslaw->ResCof();
    (*Kcontact)(pos) = Mulnslaw->Stiff();
    (*ElasCoefContact)(pos) = Mulnslaw->ElasCof();
  }
  /*
    std::cout << " Restitution coefficients: " <<std::endl;
    ResContact->display();
    std::cout << "Stiffnesses: " <<std::endl;
    Kcontact->display();
    std::cout << "Elasticity coeffients at contacts: " <<std::endl;
    ElasCoefContact->display();
  */

}
//========================================================================================
void OSNSMultipleImpact::PreComputeImpact()
{
  //1. Get the number of contacts and bodies involved in the impact
  if (indexSetLevel() != 1)
    RuntimeException::selfThrow("OSNSMultipleImpact::PreComputeImpact==> the levelMin must be equal to 1 in the multiple impact model !!");
  SP::InteractionsGraph indexSet = simulation()->indexSet(indexSetLevel()); // get indexSet[1]
  Ncontact = indexSet->size();
  //2. Compute matrix _M
  SP::Topology topology = simulation()->nonSmoothDynamicalSystem()->topology();
  bool isLinear = simulation()->nonSmoothDynamicalSystem()->isLinear();
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
  //5. Build the vectors of stifnesseses, of restitution coefficients, and of elaticity coefficients
  BuildParaContact();
}
//=======================================================================================
void OSNSMultipleImpact::InitializeInput()
{
  //Loop over alls Interactioninvolved in the indexSet[1]
  SP::InteractionsGraph indexSet = simulation()->indexSet(indexSetLevel()); // get indexSet[1]
  InteractionsGraph::VIterator ui, uiend;
  for (std11::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui)
  {
    SP::Interaction inter = indexSet->bundle(*ui);
    //SP::SiconosVector Vc0 = inter->y(1); // Relative velocity at beginning of impact
    SP::SiconosVector Vc0 = inter->yOld(1); // Relative velocity at beginning of impact
    unsigned int pos_inter = _M->getPositionOfInteractionBlock(*inter);
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
    std::cout << "Initial relative velocity at contacts" <<std::endl;
    VelContact->display();
    std::cout<< "Initial energy at contacts" <<std::endl;
    EnerContact->display();
    std::cout << "Impulse at contact" <<std::endl;
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
  SP::Topology topology = simulation()->nonSmoothDynamicalSystem()->topology();
  // Note that _interactionBlocks is up to date since updateInteractionBlocks
  // has been called during OneStepNSProblem::initialize()

  if (! _M)
  {
    if (_MStorageType == 0)
      _M.reset(new OSNSMatrix(maxSize(), 0));

    else // if(_MStorageType == 1) size = number of _interactionBlocks
      // = number of Interactionin the largest considered indexSet
      _M.reset(new OSNSMatrix(simulation()->indexSet(indexSetLevel())->size(), 1));
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
    std::cout << "Primary contact according to relative velocity: " << IdPrimaContact <<std::endl;
    std::cout << "Relative velocity at the primary contact: " << VelAtPrimaCon <<std::endl;
    std::cout << "Potential energy at the primary contact: " << EnerAtPrimaCon <<std::endl;
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
    std::cout << "Primary contact according to potenial energy: " << IdPrimaContact <<std::endl;
    std::cout << "Relative velocity at the primary contact: " << VelAtPrimaCon <<std::endl;
    std::cout << "Potential energy at the primary contact: " << EnerAtPrimaCon <<std::endl;
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
  //cout << "Is the multiple impacts is terminated: " << _IsImpactEnd <<std::endl;
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
  // std::cout << "The primary contact is :" << IdPrimaContact <<std::endl;
  // std::cout << "Is the primary contact is selected according to the potential energy: " << IsPrimaConEnergy <<std::endl;
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
          // std::cout << "_ratio_m: " << _ratio_mu <<std::endl;
          // std::cout << "Stiff: " << _stiff <<std::endl;
          // std::cout << "ratio_stiff: " << ratio_stiff <<std::endl;
          // std::cout << "energy ratio: " << ratio_ener <<std::endl;

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
    std::cout << "Relative velocity at contacts at the beginning of step:" <<std::endl;
    OldVelContact->display();
    std::cout << "Relative velocity at contacts at the end of step:" <<std::endl;
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

    std::cout << "Potential energy at contacts at the end of step:" <<std::endl;
    EnerContact->display();
    std::cout << "Work done during the compression phase at contacts" <<std::endl;
    WcContact->display();

  */
}
//======================================================================================
void OSNSMultipleImpact::UpdateDuringImpact()
{
  //1. Copy VelContact/DelImpulseContact into the vector y/lambda for Interactions
  SP::InteractionsGraph indexSet = simulation()->indexSet(indexSetLevel());
  // y and lambda vectors
  SP::SiconosVector lambda;
  SP::SiconosVector y;
  // === Loop through "active" Interactions (ie present in indexSets[1]) ===
  unsigned int pos;
  InteractionsGraph::VIterator ui, uiend;
  for (std11::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui)
  {
    Interaction& inter = *indexSet->bundle(*ui);
    // Get the relative position of inter-interactionBlock in the vector VelContact/TolImpulseContact
    pos = _M->getPositionOfInteractionBlock(inter);
    // Get Y and Lambda for the current Interaction
    y = inter.y(inputOutputLevel());
    lambda = inter.lambda(inputOutputLevel());
    // Copy VelContact/TolImpulseContact, starting from index pos into y/lambda
    // save into y !!
    setBlock(*VelContact, y, y->size(), pos, 0);
    // saved into lambda[1] !!
    setBlock(*ImpulseContact_update, lambda, lambda->size(), pos, 0);
    //setBlock(*DelImpulseContact, lambda, lambda->size(), pos, 0);
  };
  //2. Update the Input[1], state of DS systems, Output[1]
  simulation()->update(inputOutputLevel());
  ImpulseContact_update->zero(); // reset input[1] to zero after each update
}
//--------------------------------------------------------------------------------------
void OSNSMultipleImpact::SaveDataOneStep(unsigned int _ithPoint)
{
  // Save the total impulse at the primary contacts (time-like independent variable) and the time evolution during impact
  if (_ithPoint >= _DataMatrix->size(0))
    RuntimeException::selfThrow("In OSNSMultipleImpact::ComputeImpact, number of points saved exceeds the size of matrix allocated!!!");
  //(*_DataMatrix)(_ithPoint,0) = Time_variable;
  (*_DataMatrix)(_ithPoint, 0) = Impulse_variable;
  // Save the data related to UnitaryRelations
  SP::InteractionsGraph indexSet0 = simulation()->indexSet(0);
  SP::InteractionsGraph indexSet1 = simulation()->indexSet(indexSetLevel());
  unsigned int pos;
  InteractionsGraph::VIterator ui, uiend;
  unsigned int col_pos = 1;
  for (std11::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui)
  {
    SP::Interaction inter = indexSet0->bundle(*ui);
    SP::SiconosVector ydot = inter->y(1);
    SP::SiconosVector P_inter(new SiconosVector(inter->getSizeOfY()));
    SP::SiconosVector F_inter(new SiconosVector(inter->getSizeOfY()));
    SP::SiconosVector E_inter(new SiconosVector(1));
    if (indexSet1->is_vertex(inter)) // if Interaction belongs to the IndexSet[1]
    {
      pos = _M->getPositionOfInteractionBlock(*inter);
      setBlock(*TolImpulseContact, P_inter, P_inter->size(), pos, 0);
      setBlock(*ForceContact, F_inter, F_inter->size(), pos, 0);
      setBlock(*EnerContact, E_inter, E_inter->size(), pos, 0);
    }
    else
    {
      P_inter->zero();   // no impulse at this Interaction
      F_inter->zero();   // no force at this Interaction
      E_inter->zero();   // no potential at this Interaction
    };
    // Write the force at the Interaction
    WriteVectorIntoMatrix(*F_inter, _ithPoint, col_pos);
    //WriteVectorIntoMatrix(*P_inter, _ithPoint, col_pos);
    //WriteVectorIntoMatrix(*E_inter, _ithPoint, col_pos);
    col_pos = col_pos + F_inter->size();
  } // Save the data related to DS
  SP::DynamicalSystemsGraph DSG = simulation()->nonSmoothDynamicalSystem()->dynamicalSystems();
  DynamicalSystemsGraph::VIterator dsi, dsiend;
  for (std11::tie(dsi, dsiend) = DSG->vertices(); dsi != dsiend; ++dsi)
  {
    SP::DynamicalSystem ds = DSG->bundle(*dsi); // DS
    SP::LagrangianDS Lagds = std11::dynamic_pointer_cast<LagrangianDS>(ds);
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
  //cout << "*********** Impact computation progress *************" <<std::endl;
  //boost::progress_display show_progress(NstepMax);
  /*
     std::cout << "----------Before multiple impacts computation---------------" <<std::endl;
     std::cout << "Velocity at contacts: ";
     VelContact->display();
     std::cout << "Impulse at contact: ";
     TolImpulseContact->display();
  */
  //cout << "-------------------Multiple impacts computation starts:-----------------------" <<std::endl;
  // First save at the beginning of impact computation
  if ((YesSaveData) && (Step_min_save == 1))
  {
    SaveDataOneStep(point_save); // Save the data
    point_save++;
  }
  //
  while (1 != 0)
  {
    // std::cout << "==================Step==================:  " << number_step <<std::endl;
    // std::cout << "Impulse variable: " << Impulse_variable <<std::endl;
    // std::cout << "Time_variable: " << Time_variable <<std::endl;

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
      //cout << "Causion: so long computation, the computation is stopped even when the impact is not yet terminated!!! " <<std::endl;
      break;
    }
    // std::cout << "Distribution vector: ";
    // DistriVector->display();
    // std::cout << "Incremental Impulse: ";
    // DelImpulseContact->display();
    // std::cout << "Impulse at contact: ";
    // TolImpulseContact->display();
    // std::cout << "Velocity at contacts: ";
    // VelContact->display();
    // std::cout << "Potential energy at contacts: ";
    // EnerContact->display();

  }

  //
  // std::cout << "*****************Impact computation is terminated******************" <<std::endl;
  // std::cout << "Number of integration steps: " << number_step <<std::endl;
  // std::cout << "Velocity at contacts: ";
  // VelContact->display();
  // std::cout << "Impulse at contact: ";
  // TolImpulseContact->display();
  // std::cout << "Duration of the multiple impacts process: " << Time_variable << " s" <<std::endl;

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
  SP::InteractionsGraph indexSet = simulation()->indexSet(indexSetLevel());
  // y and lambda vectors
  SP::SiconosVector lambda;
  SP::SiconosVector y;
  // === Loop through "active" Interactions (ie present in indexSets[1]) ===
  unsigned int pos;
  InteractionsGraph::VIterator ui, uiend;
  for (std11::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui)
  {
    Interaction& inter = *indexSet->bundle(*ui);
    // Get the relative position of inter-interactionBlock in the vector VelContact/TolImpulseContact
    pos = _M->getPositionOfInteractionBlock(inter);
    // Get Y and Lambda for the current Interaction
    y = inter.y(inputOutputLevel());
    lambda = inter.lambda(inputOutputLevel());
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
  std::cout << "<<<<<<<<<<<<<<<<< Information about the multiple impact >>>>>>>>>>>>>>>>>>>>>" <<std::endl;
  std::cout << "Type of the contact compliance law: " << TypeCompLaw <<std::endl;
  std::cout << "Number of contacts involved into impacts: " << Ncontact <<std::endl;
  std::cout << "Step size used: " << DeltaP <<std::endl;
  std::cout << "Primary impulse at the end of impact: " << Impulse_variable <<std::endl;
  std::cout << "Duration of the multiple impacs process: " << Time_variable <<std::endl;
  // Display post-impact velocities
  SP::DynamicalSystemsGraph DSG0 = simulation()->nonSmoothDynamicalSystem()->topology()->dSG(0);
  DynamicalSystemsGraph::VIterator ui, uiend;
  for (std11::tie(ui, uiend) = DSG0->vertices(); ui != uiend; ++ui)
  {
    SP::DynamicalSystem ds = DSG0->bundle(*ui);
    SP::LagrangianDS lag_ds = std11::dynamic_pointer_cast<LagrangianDS>(ds);
    std::cout << "DS number: " << ds->number() <<std::endl;
    std::cout << "Pre-impact velocity: ";
    (lag_ds->velocityMemory()->getSiconosVector(1))->display();
    std::cout << "Post-impact velocity: ";
    (lag_ds->velocity())->display();
  }
  // Display impulses at contact points
  SP::InteractionsGraph IndexSet0 = simulation()->nonSmoothDynamicalSystem()->topology()->indexSet(0);
  InteractionsGraph::VIterator vi, viend;
  for (std11::tie(vi, viend) = IndexSet0->vertices(); vi != viend; ++vi)
  {
    SP::Interaction inter = IndexSet0->bundle(*vi);
    std::cout << "Impulse at contact point " << inter->number() << ":";
    (inter->lambda(1))->display();
  }
};
