
#include "MechanicsFwd.hpp"
#include "NewtonEulerDS.hpp"
#include "BodyDS.hpp"
#include "SiconosContactor.hpp"

// TODO
SP::SiconosMatrix identity_inertia()
{
  SP::SiconosMatrix i(new SimpleMatrix(3,3));
  i->zero();
  (*i)(0,0) = 1.0;
  (*i)(1,1) = 1.0;
  (*i)(2,2) = 1.0;
  return i;
}

BodyDS::BodyDS(SP::SiconosVector position,
               SP::SiconosVector velocity,
               double mass)
  : NewtonEulerDS(position, velocity, mass, identity_inertia())
{
}

BodyDS::~BodyDS()
{
}
