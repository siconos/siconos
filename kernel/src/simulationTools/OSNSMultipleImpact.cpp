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
OSNSMultipleImpact::OSNSMultipleImpact(): LinearOSNS(), _typeCompLaw("BiStiffness")
{
  _nStepSave = 100;
  _tolImpact = DEFAULT__tolImpact;
  _Tol_Vel = DEFAULT_TOL_VEL;
  _Tol_Ener = DEFAULT_TOL_ENER;
  _ZeroVel_EndIm = DEFAULT_TOL_VEL;
  _ZeroEner_EndIm = DEFAULT_TOL_ENER;
  _saveData = false;
  _sizeDataSave = 1000;
  _nStepMax = 100000;
  _stepMinSave = 1;
  _stepMaxSave = _nStepMax;
  _namefile = "DataMultipleImpact.dat";
}
//------------------------------ -------------------------------------------------------------
OSNSMultipleImpact::OSNSMultipleImpact(std::string newTypeLaw, double newDelP = 1.0e-5): LinearOSNS()
{
  _typeCompLaw = newTypeLaw;
  _deltaP = newDelP;
  _nStepSave = 100;
  _tolImpact = DEFAULT__tolImpact;
  _Tol_Vel = DEFAULT_TOL_VEL;
  _Tol_Ener = DEFAULT_TOL_ENER;
  _ZeroVel_EndIm = DEFAULT_TOL_VEL;
  _ZeroEner_EndIm = DEFAULT_TOL_ENER;
  _saveData = false;
  _namefile = "DataMultipleImpact.dat";
  _sizeDataSave = 1000;
  _nStepMax = 100000;
  _stepMinSave = 1;
  _stepMaxSave = _nStepMax;
  if ((_typeCompLaw != "MonoStiffness") && (_typeCompLaw != "BiStiffness"))
    RuntimeException::selfThrow("OSNSMultipleImpact::_typeCompLaw type of the compliance model must be either MonoStiffness or BiStiffness!");
}
//-------------------------------------------------------------------------------------------------
OSNSMultipleImpact::~OSNSMultipleImpact() {}
//------------------------------------------------------------------------------------------------

void OSNSMultipleImpact::setTolImpact(double newTolZero)
{
  _tolImpact = newTolZero;
};

void OSNSMultipleImpact::SetSaveData(bool var)
{
  _saveData = var;
};

void OSNSMultipleImpact::SetNameOutput(std::string file_name)
{
  _namefile = file_name;
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
  _nStepSave = var;
};

void OSNSMultipleImpact::SetNstepMax(unsigned int var)
{
  _nStepMax = var;
};

void OSNSMultipleImpact::SetStepMinMaxSave(unsigned int var1, unsigned int var2)
{
  _stepMinSave = var1;
  _stepMaxSave = var2;
}

void OSNSMultipleImpact::set_typeCompLaw(std::string newTypeLaw)
{
  _typeCompLaw = newTypeLaw;
  if ((_typeCompLaw != "MonoStiffness") && (_typeCompLaw != "BiStiffness"))
    RuntimeException::selfThrow("OSNSMultipleImpact::_typeCompLaw type of the compliance model must be either MonoStiffness or BiStiffness!");
};

void OSNSMultipleImpact::SetSizeDataSave(unsigned int var)
{
  _sizeDataSave = var;
}
//---------------------------------------------------------------------------------------------------
void OSNSMultipleImpact::WriteVectorIntoMatrix(const SiconosVector& m, const unsigned int pos_row, const unsigned int pos_col)
{
  for (unsigned int i = 0; i < m.size(); ++i)
  {
    (*_DataMatrix)(pos_row, pos_col + i) = m(i);
  }
}
//----------------------------------------------------------------------------------------------------
bool OSNSMultipleImpact::isZero(double Var)
{
  if (std::abs(Var) <= _tolImpact)
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
    _numberCols = _numberCols + (DSG->bundle(*dsi)->dimension());
  }
  return(_numberCols);
}
//-----------------------------------------------------------------------------------------------
void OSNSMultipleImpact::AllocateMemory()
{
  if (!_velocityContact)
    _velocityContact.reset(new SiconosVector(maxSize()));
  else
  {
    if (_velocityContact->size() != maxSize())
      _velocityContact->resize(maxSize());
  };
  //
  if (!_oldVelocityContact)
    _oldVelocityContact.reset(new SiconosVector(maxSize()));
  else
  {
    if (_oldVelocityContact->size() != maxSize())
      _oldVelocityContact->resize(maxSize());
  };
  //
  if (! _energyContact)
    _energyContact.reset(new SiconosVector(maxSize()));
  else
  {
    if (_energyContact->size() != maxSize())
      _energyContact->resize(maxSize());
  };
  //
  if (!_WorkcContact)
    _WorkcContact.reset(new SiconosVector(maxSize()));
  else
  {
    if (_WorkcContact->size() != maxSize())
      _WorkcContact->resize(maxSize());
  };
  //
  if (!_distributionVector)
    _distributionVector.reset(new SiconosVector(maxSize()));
  else
  {
    if (_distributionVector->size() != maxSize())
      _distributionVector->resize(maxSize());
  };
  //
  if (!_stateContact)
    _stateContact.reset(new IndexInt(maxSize()));
  else
  {
    if (_stateContact->size() != maxSize())
      _stateContact->resize(maxSize());
  };
  //
  if (!_Kcontact)
    _Kcontact.reset(new SiconosVector(maxSize()));
  else
  {
    if (_Kcontact->size() != maxSize())
      _Kcontact->resize(maxSize());
  };
  //
  if (!_restitutionContact)
    _restitutionContact.reset(new SiconosVector(maxSize()));
  else
  {
    if (_restitutionContact->size() != maxSize())
      _restitutionContact->resize(maxSize());
  };
  //
  if (!_elasticyCoefficientcontact)
    _elasticyCoefficientcontact.reset(new SiconosVector(maxSize()));
  else
  {
    if (_elasticyCoefficientcontact->size() != maxSize())
      _elasticyCoefficientcontact->resize(maxSize());
  };
  if (!_tolImpulseContact)
    _tolImpulseContact.reset(new SiconosVector(maxSize()));
  else
  {
    if (_tolImpulseContact->size() != maxSize())
      _tolImpulseContact->resize(maxSize());
  };
  //
  if (!_deltaImpulseContact)
    _deltaImpulseContact.reset(new SiconosVector(maxSize()));
  else
  {
    if (_deltaImpulseContact->size() != maxSize())
      _deltaImpulseContact->resize(maxSize());
  };
  //
  if (!_impulseContactUpdate)
    _impulseContactUpdate.reset(new SiconosVector(maxSize()));
  else
  {
    if (_impulseContactUpdate->size() != maxSize())
      _impulseContactUpdate->resize(maxSize());
  }
  //
  if (!_forceContact)
    _forceContact.reset(new SiconosVector(maxSize()));
  else
  {
    if (_forceContact->size() != maxSize())
      _forceContact->resize(maxSize());
  };
  // for the data matrix
  unsigned int _numberCols = EstimateNdataCols();
  if (!_DataMatrix)
    _DataMatrix.reset(new SimpleMatrix(_sizeDataSave, _numberCols));
  else
  {
    if ((_DataMatrix->size(0) != _sizeDataSave) || (_DataMatrix->size(1) != _numberCols))
      _DataMatrix->resize(_sizeDataSave, _numberCols);
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
    SP::NonSmoothLaw nslaw = inter->nonSmoothLaw();
    SP::MultipleImpactNSL Mulnslaw = std11::dynamic_pointer_cast<MultipleImpactNSL>(nslaw);
    assert(Mulnslaw && "In OSNSMultipleImpact::BuildStiffResCofVec, non-smooth law used must be MultipleImpactNSL!!!");
    // Get the relative position of inter-interactionBlock in the vector _velocityContact
    unsigned int pos = _M->getPositionOfInteractionBlock(*inter);
    (*_restitutionContact)(pos) = Mulnslaw->ResCof();
    (*_Kcontact)(pos) = Mulnslaw->Stiff();
    (*_elasticyCoefficientcontact)(pos) = Mulnslaw->ElasCof();
  }
  /*
    std::cout << " Restitution coefficients: " <<std::endl;
    _restitutionContact->display();
    std::cout << "Stiffnesses: " <<std::endl;
    _Kcontact->display();
    std::cout << "Elasticity coeffients at contacts: " <<std::endl;
    _elasticyCoefficientcontact->display();
  */

}
//========================================================================================
void OSNSMultipleImpact::PreComputeImpact()
{
  //1. Get the number of contacts and bodies involved in the impact
  if (indexSetLevel() != 1)
    RuntimeException::selfThrow("OSNSMultipleImpact::PreComputeImpact==> the levelMin must be equal to 1 in the multiple impact model !!");
  SP::InteractionsGraph indexSet = simulation()->indexSet(indexSetLevel()); // get indexSet[1]
  _nContact = indexSet->size();
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
  if (_nContact != _sizeOutput)
    RuntimeException::selfThrow("OSNSMultipleImpact::ComputeWMinvWtrans: number of contacts different from the size of output--> this case is not yet implemented!");
  //3. Checks size of vectors
  if (_velocityContact->size() != _sizeOutput)
  {
    _velocityContact->resize(_sizeOutput);
  }
  _velocityContact->zero();
  //
  if (_oldVelocityContact->size() != _sizeOutput)
  {
    _oldVelocityContact->resize(_sizeOutput);
  }
  _oldVelocityContact->zero();
  //
  if (_energyContact->size() != _sizeOutput)
  {
    _energyContact->resize(_sizeOutput);
  }
  _energyContact->zero();
  //
  if (_WorkcContact->size() != _sizeOutput)
  {
    _WorkcContact->resize(_sizeOutput);
  }
  _WorkcContact->zero();
  //
  if (_distributionVector->size() != _sizeOutput)
  {
    _distributionVector->resize(_sizeOutput);
  }
  _distributionVector->zero();
  //
  if (_stateContact->size() != _sizeOutput)
  {
    _stateContact->resize(_sizeOutput);
  }
  //
  if (_Kcontact->size() != _sizeOutput)
  {
    _Kcontact->resize(_sizeOutput);
  }
  _Kcontact->zero();
  //
  if (_restitutionContact->size() != _sizeOutput)
  {
    _restitutionContact->resize(_sizeOutput);
  }
  _restitutionContact->zero();
  //
  if (_elasticyCoefficientcontact->size() != _sizeOutput)
  {
    _elasticyCoefficientcontact->resize(_sizeOutput);
  }
  _elasticyCoefficientcontact->zero();
  //
  if (_tolImpulseContact->size() != _sizeOutput)
  {
    _tolImpulseContact->resize(_sizeOutput);
  }
  _tolImpulseContact->zero();
  //
  if (_deltaImpulseContact->size() != _sizeOutput)
  {
    _deltaImpulseContact->resize(_sizeOutput);
  }
  _deltaImpulseContact->zero();
  //
  if (_impulseContactUpdate->size() != _sizeOutput)
  {
    _impulseContactUpdate->resize(_sizeOutput);
  }
  _impulseContactUpdate->zero();
  //
  if (_forceContact->size() != _sizeOutput)
  {
    _forceContact->resize(_sizeOutput);
  }
  _forceContact->zero();
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
    setBlock(*Vc0, _velocityContact, Vc0->size(), 0, pos_inter);
    SP::SiconosVector ener0(new SiconosVector(Vc0->size()));
    ener0->zero(); // We suppose that the initial potential energy before impact is equal to zero at any contact
    // at the beginning of impact
    setBlock(*ener0, _energyContact, ener0->size(), 0, pos_inter);
    //SP::SiconosVector impulse0= (inter)->lambda(1))->vector(inter->number());
    SP::SiconosVector impulse0(new SiconosVector(Vc0->size()));
    impulse0->zero(); // We suppose that the impulse before impact is equal to zero at any contact
    // at the beginning of impact
    setBlock(*impulse0, _tolImpulseContact, impulse0->size(), 0, pos_inter);
  };
  /*
    std::cout << "Initial relative velocity at contacts" <<std::endl;
    _velocityContact->display();
    std::cout<< "Initial energy at contacts" <<std::endl;
    _energyContact->display();
    std::cout << "Impulse at contact" <<std::endl;
    _tolImpulseContact->display();
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
  getMin(*_velocityContact, _relativeVelocityPrimaryContact, _primaryContactId);
  _energyPrimaryContact = (*_energyContact)(_primaryContactId);
  if (!isVelNegative(_relativeVelocityPrimaryContact))
  {
    RuntimeException::selfThrow("OSNSMultipleImpact::PrimConVelocity, the velocity at the primary contact must be negative !!");
  }
  /*
    std::cout << "Primary contact according to relative velocity: " << _primaryContactId <<std::endl;
    std::cout << "Relative velocity at the primary contact: " << _relativeVelocityPrimaryContact <<std::endl;
    std::cout << "Potential energy at the primary contact: " << _energyPrimaryContact <<std::endl;
  */
}
//=======================================================================================
void OSNSMultipleImpact::PrimConEnergy()
{
  getMax(*_energyContact, _energyPrimaryContact, _primaryContactId);
  _relativeVelocityPrimaryContact = (*_velocityContact)(_primaryContactId);
  if (_energyPrimaryContact < 0.0)
  {
    RuntimeException::selfThrow("OSNSMultipleImpact::PrimConEnergy the potential energy at the primary contact must be positive !!");
  }
  /*
    std::cout << "Primary contact according to potenial energy: " << _primaryContactId <<std::endl;
    std::cout << "Relative velocity at the primary contact: " << _relativeVelocityPrimaryContact <<std::endl;
    std::cout << "Potential energy at the primary contact: " << _energyPrimaryContact <<std::endl;
  */

}
//======================================================================================
bool OSNSMultipleImpact::IsEnermaxZero()
{
  double MaxEner;
  unsigned int IdMax;
  getMax(*_energyContact, MaxEner, IdMax);
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
  getMin(*_velocityContact, MinVelCon, IdConVmin);
  if (isVelNegative(MinVelCon))
    return true;
  else
    return false;
}
//=======================================================================================
void OSNSMultipleImpact::Check_stateContact()
{
  for (unsigned int i = 0; i < _nContact; ++i)
  {
    if (isEnerZero((*_energyContact)(i))) // potential energy is zero
    {
      if (!isVelNegative((*_velocityContact)(i))) // relative velocity is positive or equal to zero
        (*_stateContact)[i] = 0; // no impact at this contact
      else  // impact happens without potential energy
      {
        (*_stateContact)[i] = 1;
      }
    }
    else // impact happens with not zero potential energy
    {
      if ((*_stateContact)[i] != 2)
      {
        (*_stateContact)[i] = 2;
      }
    }
  }
}
//=======================================================================================
bool OSNSMultipleImpact::IsMulImpactTerminate()
{
  _IsImpactEnd = true;
  for (unsigned int i = 0; i < _nContact; ++i)
  {
    if (((*_energyContact)(i) > _ZeroEner_EndIm) || ((*_velocityContact)(i) < -_ZeroVel_EndIm)) // if potential energy is not equal to zero or the relative velocity is negative
    {
      _IsImpactEnd = false;
    }
  }
  return _IsImpactEnd;
  //   bool var = true;
  //   for(unsigned int i = 0; i < _nContact;++i)
  //     {
  //       if ((*_stateContact)[i] != 0)
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
    _isPrimaryContactEnergy = false;
  }
  else
  {
    PrimConEnergy(); // Select the primary contact according to the potential energy at contacts
    _isPrimaryContactEnergy = true;
  }
  //
  // std::cout << "The primary contact is :" << _primaryContactId <<std::endl;
  // std::cout << "Is the primary contact is selected according to the potential energy: " << _isPrimaryContactEnergy <<std::endl;
}
//=======================================================================================
void OSNSMultipleImpact::Compute_distributionVector()
{
  //Case 1: if no potential energy at any contact
  double _ratio_mu, ratio_stiff, ratio_ener;
  double mu_prima = (*_elasticyCoefficientcontact)(_primaryContactId); // Elasticity coefficient at the primary contact
  double stiff_prima = (*_Kcontact)(_primaryContactId);     // Stiffness at the primary contact
  double _mu, _stiff, _vel, _energy;
  if (!_isPrimaryContactEnergy) // case of primary contact selected according to the relative velocity
  {
    double ratio_vel;
    for (unsigned int i = 0; i < _nContact; ++i)
    {
      if ((*_stateContact)[i] != 0) // the impact can takes place at this contact
      {
        _mu = (*_elasticyCoefficientcontact)(i); // Elasticity coefficient at the current contact
        _stiff = (*_Kcontact)(i);     // Stiffness at the current contact
        _vel = (*_velocityContact)(i);     // Relative velocity at the current contact
        _ratio_mu = (std::pow(_mu + 1.0, (_mu / (_mu + 1.0)))) / (std::pow(mu_prima + 1.0, (mu_prima / (mu_prima + 1.0))));
        ratio_stiff = (std::pow(_stiff, (1.0 / (1.0 + _mu)))) / (std::pow(stiff_prima, (1.0 / (1.0 + mu_prima))));
        if (!isVelNegative(_vel))
        {
          RuntimeException::selfThrow("OSNSMultipleImpact::Compute_distributionVector, the relative velocity when particle starts to impact must be negative!!");
        }

        ratio_vel = (std::pow(std::fabs(_vel), (_mu / (_mu + 1.0)))) / (std::pow(std::fabs(_relativeVelocityPrimaryContact), (mu_prima / (1.0 + mu_prima))));
        (*_distributionVector)(i) = std::pow((_ratio_mu * ratio_stiff * ratio_vel), (1.0 + _mu)) * std::pow(_deltaP, ((_mu - mu_prima) / (1.0 + mu_prima)));
      }
      else
      {
        (*_distributionVector)(i) = 0.0;
      }
      if ((*_distributionVector)(i) < 0.0)
        RuntimeException::selfThrow("OSNSMultipleImpact::Compute_distributionVector the component of _distributionVector must be positive !!");
    };
  }
  //Case 2: case of primary contact selected according to the potential energy
  else
  {
    for (unsigned int i = 0; i < _nContact; ++i)
    {
      //
      _mu = (*_elasticyCoefficientcontact)(i);
      _stiff = (*_Kcontact)(i);
      _ratio_mu = (std::pow(_mu + 1.0, (_mu / (_mu + 1.0)))) / (std::pow(mu_prima + 1.0, (mu_prima / (mu_prima + 1.0))));
      ratio_stiff = (std::pow(_stiff, (1.0 / (1.0 + _mu)))) / (std::pow(stiff_prima, (1.0 / (1.0 + mu_prima))));
      if ((*_stateContact)[i] == 1) // no potential energy at this contact, including the contacts at which impact repeats
      {
        if (!isVelNegative((*_velocityContact)(i)))
        {
          RuntimeException::selfThrow("OSNSMultipleImpact::Compute_distributionVector, the pre-impact velocity must be negative!!");
        }
        else
        {
          _vel = (*_velocityContact)(i);
          ratio_ener = (std::pow(std::fabs(_vel * _deltaP), (_mu / (_mu + 1.0)))) / (std::pow(_energyPrimaryContact, (mu_prima / (mu_prima + 1.0))));
          //
          // std::cout << "_ratio_m: " << _ratio_mu <<std::endl;
          // std::cout << "Stiff: " << _stiff <<std::endl;
          // std::cout << "ratio_stiff: " << ratio_stiff <<std::endl;
          // std::cout << "energy ratio: " << ratio_ener <<std::endl;

          //
          (*_distributionVector)(i) = std::pow((_ratio_mu * ratio_stiff * ratio_ener), (1.0 + _mu));
        }
      }
      else if ((*_stateContact)[i] == 2) // potential is not zero at this contact
      {
        _energy = (*_energyContact)(i); // Potential energy at the current contact
        ratio_ener = (std::pow(_energy, (_mu / (_mu + 1.0)))) / (std::pow(_energyPrimaryContact, (mu_prima / (mu_prima + 1.0))));
        (*_distributionVector)(i) = _ratio_mu * ratio_stiff * ratio_ener;
      }
      else // no impact at this contact
      {
        (*_distributionVector)(i) = 0.0;
      };
      if ((*_distributionVector)(i) < 0.0)
        RuntimeException::selfThrow("OSNSMultipleImpact::Compute_distributionVector the component of _distributionVector must be positive !!");
    };
  };
}
//=======================================================================================
void OSNSMultipleImpact::ComputeImpulseContact()
{
  (*_deltaImpulseContact) = (*_distributionVector) * _deltaP;
  (*_tolImpulseContact) = (*_tolImpulseContact) + (*_deltaImpulseContact);
  (*_impulseContactUpdate) = (*_impulseContactUpdate) + (*_deltaImpulseContact);
  // Compute the contact force
  double PowCompLaw;
  for (unsigned int i = 0; i < _nContact; ++i)
  {
    PowCompLaw = (*_elasticyCoefficientcontact)(i);
    if (isEnerZero((*_energyContact)(i))) // if potential energy at this contact is zero
    {
      if (isVelNegative((*_velocityContact)(i))) // if the relative velocity at contact is negative
      {
        (*_forceContact)(i) = std::pow((1.0 + PowCompLaw), PowCompLaw / (1.0 + PowCompLaw)) * std::pow((*_Kcontact)(i), 1.0 / (1.0 + PowCompLaw)) * std::pow((std::fabs((*_velocityContact)(i)) * (*_deltaImpulseContact)(i)), PowCompLaw / (1.0 + PowCompLaw));
      }
      else
      {
        (*_forceContact)(i) = 0.0;
      };
    }
    else
    {
      (*_forceContact)(i) = std::pow((1.0 + PowCompLaw), PowCompLaw / (1.0 + PowCompLaw)) * std::pow((*_Kcontact)(i), 1.0 / (1.0 + PowCompLaw)) * std::pow((*_energyContact)(i), PowCompLaw / (1.0 + PowCompLaw));
    }
    if ((*_forceContact)(i) < 0.0)
    {
      RuntimeException::selfThrow("OSNSMultipleImpact::ComputeImpulseContact, the contact force must be positive or equal to zero!!!");
    }
  };
}
//=======================================================================================
void OSNSMultipleImpact::Compute_velocityContact()
{
  (*_oldVelocityContact) = (*_velocityContact); //save the relative velocity at the beginning of the step
  (*_velocityContact) = (*_velocityContact) + prod(*(_M->defaultMatrix()), *_deltaImpulseContact); // compute the relative velocity at the end of the step
  //
  /*
    std::cout << "Relative velocity at contacts at the beginning of step:" <<std::endl;
    _oldVelocityContact->display();
    std::cout << "Relative velocity at contacts at the end of step:" <<std::endl;
    _velocityContact->display();
  */
  //
}
//=======================================================================================
void OSNSMultipleImpact::Compute_energyContact()
{
  if (_typeCompLaw == "BiStiffness")
    // For Bistiffness model
  {
    for (unsigned int i = 0; i < _nContact; ++i)
    {
      if ((0.5 * ((*_oldVelocityContact)(i) + (*_velocityContact)(i))) <= 0.0) // Contact located in the compression phase
      {
        (*_energyContact)(i) = (*_energyContact)(i) - 0.5 * ((*_oldVelocityContact)(i) + (*_velocityContact)(i)) * ((*_deltaImpulseContact)(i));
      }
      else                       // Contact located in the expansion phase
      {
        if (!isZero((*_restitutionContact)(i)))
        {
          (*_energyContact)(i) = (*_energyContact)(i) - (1.0 / std::pow((*_restitutionContact)(i), 2)) * 0.5 * ((*_oldVelocityContact)(i) + (*_velocityContact)(i)) * ((*_deltaImpulseContact)(i));
          //
          if ((*_energyContact)(i) <  0.0)
          {
            (*_energyContact)(i) = 0.0;
          };
        }
        else // restitution coefficient equal to zero
        {
          (*_energyContact)(i) = 0.0; // In this case, no potential energy at contacts when the contact is located in the compression phase
        }
      };
      //
      if ((*_energyContact)(i) <  0.0)
      {
        RuntimeException::selfThrow("OSNSMultipleImpact::Compute_energyContact, the potential energy during compression phase must be positive!!!");
      };
    };
  }
  else
    // For the mono-stiffness model
  {
    for (unsigned int i = 0; i < _nContact; ++i)
    {
      //1: dertermine the work done by the last compression phase at contacts (nessessary for the mono-stiffness compliance model)
      // if Vc(k) < 0 and Vc(k +1) >= 0 ==> transition from the compression phase to the expansion phase, Wc = E(k)
      if (((*_oldVelocityContact)(i) < 0.0) && ((*_velocityContact)(i) >= 0.0))
      {
        (*_WorkcContact)(i) = (*_energyContact)(i);
      };
      //2: Calculate the potential energy at the end of stap
      (*_energyContact)(i) = (*_energyContact)(i) - 0.5 * ((*_oldVelocityContact)(i) + (*_velocityContact)(i)) * ((*_deltaImpulseContact)(i));
      //3: Check if the termination condition is verified or not (if Vc(k+1) > 0.0 and E(k+1) <= (1-e^2)*Wc). If yes, discard the potential energy
      // in order to respect the energetic constraint
      if (((*_stateContact)[i] == 2) && (((*_velocityContact)(i) > 0.0) && ((*_energyContact)(i) <= ((1.0 - std::pow((*_restitutionContact)(i), 2)) * (*_WorkcContact)(i)))))
      {
        (*_energyContact)(i) = 0.0; // potential energy at this contact is completely dissipated before the compression phase finishes
      };
    };
  }
  /*

    std::cout << "Potential energy at contacts at the end of step:" <<std::endl;
    _energyContact->display();
    std::cout << "Work done during the compression phase at contacts" <<std::endl;
    _WorkcContact->display();

  */
}
//======================================================================================
void OSNSMultipleImpact::UpdateDuringImpact()
{
  //1. Copy _velocityContact/_deltaImpulseContact into the vector y/lambda for Interactions
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
    // Get the relative position of inter-interactionBlock in the vector _velocityContact/_tolImpulseContact
    pos = _M->getPositionOfInteractionBlock(inter);
    // Get Y and Lambda for the current Interaction
    y = inter.y(inputOutputLevel());
    lambda = inter.lambda(inputOutputLevel());
    // Copy _velocityContact/_tolImpulseContact, starting from index pos into y/lambda
    // save into y !!
    setBlock(*_velocityContact, y, y->size(), pos, 0);
    // saved into lambda[1] !!
    setBlock(*_impulseContactUpdate, lambda, lambda->size(), pos, 0);
    //setBlock(*_deltaImpulseContact, lambda, lambda->size(), pos, 0);
  };
  //2. Update the Input[1], state of DS systems, Output[1]
  simulation()->update(inputOutputLevel());
  _impulseContactUpdate->zero(); // reset input[1] to zero after each update
}
//--------------------------------------------------------------------------------------
void OSNSMultipleImpact::SaveDataOneStep(unsigned int _ithPoint)
{
  // Save the total impulse at the primary contacts (time-like independent variable) and the time evolution during impact
  if (_ithPoint >= _DataMatrix->size(0))
    RuntimeException::selfThrow("In OSNSMultipleImpact::ComputeImpact, number of points saved exceeds the size of matrix allocated!!!");
  //(*_DataMatrix)(_ithPoint,0) = _timeVariable;
  (*_DataMatrix)(_ithPoint, 0) = _impulseVariable;
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
      setBlock(*_tolImpulseContact, P_inter, P_inter->size(), pos, 0);
      setBlock(*_forceContact, F_inter, F_inter->size(), pos, 0);
      setBlock(*_energyContact, E_inter, E_inter->size(), pos, 0);
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
  _impulseVariable = 0.0;
  _timeVariable = 0.0;
  unsigned int number_step = 1;
  unsigned int point_save = 0;
  unsigned int _counterstepsave = 0;
  // Show computation progress
  //cout << "*********** Impact computation progress *************" <<std::endl;
  //boost::progress_display show_progress(_nStepMax);
  /*
     std::cout << "----------Before multiple impacts computation---------------" <<std::endl;
     std::cout << "Velocity at contacts: ";
     _velocityContact->display();
     std::cout << "Impulse at contact: ";
     _tolImpulseContact->display();
  */
  //cout << "-------------------Multiple impacts computation starts:-----------------------" <<std::endl;
  // First save at the beginning of impact computation
  if ((_saveData) && (_stepMinSave == 1))
  {
    SaveDataOneStep(point_save); // Save the data
    point_save++;
  }
  //
  while (1 != 0)
  {
    // std::cout << "==================Step==================:  " << number_step <<std::endl;
    // std::cout << "Impulse variable: " << _impulseVariable <<std::endl;
    // std::cout << "_timeVariable: " << _timeVariable <<std::endl;

    //Step 1: check the state at contacts
    Check_stateContact();
    //Step 2: check if the multiple impact is terminated or not
    if (IsMulImpactTerminate()) // multiple impact terminated
    {
      if (_saveData) // Save the date at the end of impact
      {
        UpdateDuringImpact(); // Update state of dynamical system
        SaveDataOneStep(point_save); // Save the data
      }
      break;
    }
    // Select the primary contact
    SelectPrimaContact();
    //Step 3: compute the vector of distributing law
    Compute_distributionVector();
    //Step 4: compute the increment of normal impulse and the total impulse at contacts
    ComputeImpulseContact();
    // Step 5: compute the relative velocity at contacts
    Compute_velocityContact();
    // Step 6: compute the potential energy at contacts
    Compute_energyContact();
    //Step 7: Update the time-like variable
    ++number_step;
    ++_counterstepsave;
    //++show_progress;
    _impulseVariable = _impulseVariable + _deltaP;
    _timeVariable = _timeVariable + _deltaP / (*_forceContact)(_primaryContactId);
    // Step 8: update the state of DS and output during impact and write data into output file at the beginning of each step
    if ((_saveData) & (_counterstepsave >= _nStepSave))
    {
      if ((number_step >= _stepMinSave) && (number_step <= _stepMaxSave))
      {
        UpdateDuringImpact(); // Update state of dynamical system
        SaveDataOneStep(point_save); // Save the data
        point_save++;
        _counterstepsave = 0; // reset the counter to 0
      }
    }
    //
    if (number_step > _nStepMax)
    {
      RuntimeException::selfThrow("In OSNSMultipleImpact::ComputeImpact, number of integration steps perfomed exceeds the maximal number of steps allowed!!!");
      //cout << "Causion: so long computation, the computation is stopped even when the impact is not yet terminated!!! " <<std::endl;
      break;
    }
    // std::cout << "Distribution vector: ";
    // _distributionVector->display();
    // std::cout << "Incremental Impulse: ";
    // _deltaImpulseContact->display();
    // std::cout << "Impulse at contact: ";
    // _tolImpulseContact->display();
    // std::cout << "Velocity at contacts: ";
    // _velocityContact->display();
    // std::cout << "Potential energy at contacts: ";
    // _energyContact->display();

  }

  //
  // std::cout << "*****************Impact computation is terminated******************" <<std::endl;
  // std::cout << "Number of integration steps: " << number_step <<std::endl;
  // std::cout << "Velocity at contacts: ";
  // _velocityContact->display();
  // std::cout << "Impulse at contact: ";
  // _tolImpulseContact->display();
  // std::cout << "Duration of the multiple impacts process: " << _timeVariable << " s" <<std::endl;

  // Close the stream file
  if (_saveData)
  {
    ioMatrix::write(_namefile.c_str(), "ascii", *_DataMatrix, "noDim");
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
    // Get the relative position of inter-interactionBlock in the vector _velocityContact/_tolImpulseContact
    pos = _M->getPositionOfInteractionBlock(inter);
    // Get Y and Lambda for the current Interaction
    y = inter.y(inputOutputLevel());
    lambda = inter.lambda(inputOutputLevel());
    // Copy _velocityContact/_tolImpulseContact, starting from index pos into y/lambda
    // save into y !!
    setBlock(*_velocityContact, y, y->size(), pos, 0);// Warning: yEquivalent is
    // saved into lambda[1] !!
    setBlock(*_impulseContactUpdate, lambda, lambda->size(), pos, 0);
    // If the update is performed at the end of the impact process, we update the total normal impulse at contacts
    // from the beginning to the end of impact (vector _tolImpulseContact). Otherwise, we must reset the lambda[1] to zero because
    // the post-impact velocity has been calculated during impact
    // if (!_saveData) // we update the impact state at the end of impact
    //   {
    //     // Copy _velocityContact/_tolImpulseContact, starting from index pos into y/lambda
    //     // save into y !!
    //     setBlock(*_velocityContact, y, y->size(), pos, 0);// Warning: yEquivalent is
    //     // saved into lambda[1] !!
    //     setBlock(*_tolImpulseContact, lambda, lambda->size(), pos, 0);
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
  if ((_nContact != 0) && IsVcminNegative()) // if there is at least one contact and the vilocity before impact is negative
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
  std::cout << "Type of the contact compliance law: " << _typeCompLaw <<std::endl;
  std::cout << "Number of contacts involved into impacts: " << _nContact <<std::endl;
  std::cout << "Step size used: " << _deltaP <<std::endl;
  std::cout << "Primary impulse at the end of impact: " << _impulseVariable <<std::endl;
  std::cout << "Duration of the multiple impacs process: " << _timeVariable <<std::endl;
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
