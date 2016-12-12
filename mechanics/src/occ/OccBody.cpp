#include "OccBody.hpp"
#include "OccBody_impl.hpp"
#include "OccUtils.hpp"

OccBody::OccBody(SP::SiconosVector position,
                 SP::SiconosVector velocity,
                 double mass ,
                 SP::SiconosMatrix inertia) :
  NewtonEulerDS(position, velocity, mass, inertia),
  _contactShapes(new ContactShapes()),
  _shapes(new TopoDS_Shapes())
{}


void OccBody::addContactShape(SP::OccContactShape shape,
                              SP::SiconosVector position,
                              SP::SiconosVector orientation,
                              unsigned int group)
{
  this->_contactShapes->push_back(shape);
  this->updateContactShapes();
  shape->computeUVBounds();
}

void OccBody::addShape(SP::TopoDS_Shape shape,
                       SP::SiconosVector position,
                       SP::SiconosVector orientation)
{
  this->_shapes->push_back(shape);
  this->updateShapes();
}

void OccBody::updateContactShapes()
{
  for (ContactShapes::iterator csi = _contactShapes->begin();
       csi != _contactShapes->end(); ++ csi)
  {
    occ_move((**csi).data(), *_q);
  }
}

void OccBody::updateShapes()
{
  for (TopoDS_Shapes::iterator csi = _shapes->begin();
       csi != _shapes->end(); ++ csi)
  {
    occ_move(**csi, *_q);
  }
}

const OccContactShape& OccBody::contactShape(unsigned int id) const
{
  return *(*this->_contactShapes)[id];
}

const TopoDS_Shape& OccBody::shape(unsigned int id) const
{
  return *(*this->_shapes)[id];
}

