#include "OccBody.hpp"
#include "OccBody_impl.hpp"

OccBody::OccBody(SP::SiconosVector position,
                 SP::SiconosVector velocity,
                 double mass ,
                 SP::SiconosMatrix inertia) :
  NewtonEulerDS(position, velocity, mass, inertia),
  _contactShapes(new ContactShapes())
  {};


void OccBody::addContactShape(SP::OccContactShape shape)
{
  this->_contactShapes->push_back(shape);
};

void OccBody::updateContactShapes()
{
  for (ContactShapes::iterator csi = _contactShapes->begin();
       csi != _contactShapes->end(); ++ csi)
  {
    (**csi).move(*_q);
  }
}

const OccContactShape& OccBody::contactShape(unsigned int id) const
{
  return *(*this->_contactShapes)[id];
};
