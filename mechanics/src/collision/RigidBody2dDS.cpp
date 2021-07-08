#include "MechanicsFwd.hpp"
#include "NewtonEulerDS.hpp"
#include "RigidBody2dDS.hpp"
#include "SiconosContactor.hpp"
RigidBody2dDS::RigidBody2dDS(SP::SiconosVector position,
                             SP::SiconosVector velocity,
                             SP::SiconosMatrix  mass)
  : LagrangianLinearTIDS(position, velocity, mass)
  , _contactors(std::make_shared<SiconosContactorSet>())
  , _useContactorInertia(true)
  , _allowSelfCollide(true)
{
  // Check size of positions, velocities and mass matrix
  if((position->size() !=3) or (velocity->size() !=3))
  {
    THROW_EXCEPTION("RigidBody2dDS::RigidBody2dDS(...). The size of position and velocity must of size 3");
  }

  _scalarMass = mass->getValue(0,0);
}

RigidBody2dDS::RigidBody2dDS(SP::SiconosVector position,
                             SP::SiconosVector velocity,
                             double mass,
                             double inertia)
  : LagrangianLinearTIDS(position, velocity, SP::SimpleMatrix(new SimpleMatrix(3,3)))
  ,_scalarMass(mass)
  , _contactors(std::make_shared<SiconosContactorSet>())
  , _useContactorInertia(true)
  , _allowSelfCollide(true)
{
  _mass->setValue(0,0,mass);
  _mass->setValue(1,1,mass);
  _mass->setValue(2,2,inertia);

  // Check size of positions, velocities and mass matrix
    if((position->size() !=3) or (velocity->size() !=3))
  {
    THROW_EXCEPTION("RigidBody2dDS::RigidBody2dDS(...). The size of position and velocity must of size 3");
  }


}

RigidBody2dDS::~RigidBody2dDS()
{
}
