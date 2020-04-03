#include "MechanicsFwd.hpp"
#include "NewtonEulerDS.hpp"
#include "RigidBody2dDS.hpp"
#include "SiconosContactor.hpp"
RigidBody2dDS::RigidBody2dDS(SP::SiconosVector position,
                             SP::SiconosVector velocity,
                             SP::SiconosMatrix  mass)
  : LagrangianDS(position, velocity, mass)
  , _contactors(std::make_shared<SiconosContactorSet>())
  , _useContactorInertia(true)
  , _allowSelfCollide(true)
{
  // Check size of positions, velocities and mass matrix
  _scalarMass = mass->getValue(0,0);
}

RigidBody2dDS::RigidBody2dDS(SP::SiconosVector position,
                             SP::SiconosVector velocity,
                             double mass,
                             double inertia)
  : LagrangianDS(position, velocity)
  ,_scalarMass(mass)
  , _contactors(std::make_shared<SiconosContactorSet>())
  , _useContactorInertia(true)
  , _allowSelfCollide(true)
{


  SP::SiconosMatrix mass_matrix(new SimpleMatrix(3,3));
  mass_matrix->setValue(0,0,mass);
  mass_matrix->setValue(1,1,mass);
  mass_matrix->setValue(2,2,inertia);

  _mass = mass_matrix;


  // Check size of positions, velocities and mass matrix
}

RigidBody2dDS::~RigidBody2dDS()
{
}
