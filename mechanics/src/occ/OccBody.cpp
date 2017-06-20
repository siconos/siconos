#include "OccBody.hpp"
#include "OccBody_impl.hpp"
#include "OccUtils.hpp"

#include <boost/math/quaternion.hpp>

OccBody::OccBody(SP::SiconosVector position,
                 SP::SiconosVector velocity,
                 double mass ,
                 SP::SiconosMatrix inertia) :
  NewtonEulerDS(position, velocity, mass, inertia),
  _contactShapes(new ContactShapes()),
  _shapes(new TopoDS_Shapes())
{}


void OccBody::addContactShape(SP::OccContactShape shape,
                              SP::SiconosVector pos,
                              SP::SiconosVector ori,
                              unsigned int group)
{
  OffSet offset = {0, 0, 0, 1, 0, 0 ,0};
  if (pos) {
    offset[0] = (*pos)(0);
    offset[1] = (*pos)(1);
    offset[2] = (*pos)(2);
  }
  if (ori) {
    offset[3] = (*ori)(0);
    offset[4] = (*ori)(1);
    offset[5] = (*ori)(2);
    offset[6] = (*ori)(3);
  }

  this->_contactShapes->push_back(
    boost::tuple<SP::OccContactShape, OffSet, int>(shape, offset, group));

  this->updateContactShapes();
  shape->computeUVBounds();
}

void OccBody::addShape(SP::TopoDS_Shape shape,
                       SP::SiconosVector pos,
                       SP::SiconosVector ori)
{
  OffSet offset = {0, 0, 0, 1, 0, 0 ,0};
  if (pos) {
    offset[0] = (*pos)(0);
    offset[1] = (*pos)(1);
    offset[2] = (*pos)(2);
  }
  if (ori) {
    offset[3] = (*ori)(0);
    offset[4] = (*ori)(1);
    offset[5] = (*ori)(2);
    offset[6] = (*ori)(3);
  }

  this->_shapes->push_back(
    boost::tuple<SP::TopoDS_Shape, OffSet>(shape, offset));

  this->updateShapes();
}

void OccBody::updateContactShapes()
{
  boost::math::quaternion<double> q((*_q)(3), (*_q)(4), (*_q)(5), (*_q)(6));

  for (ContactShapes::iterator csi = _contactShapes->begin();
       csi != _contactShapes->end(); ++ csi)
  {
    OffSet offset = boost::get<1>(*csi);

    boost::math::quaternion<double> pv = boost::math::quaternion<double>(0, offset[0], offset[1], offset[2]);

    boost::math::quaternion<double> rv = q * pv * boost::math::conj(q);

    boost::math::quaternion<double> r = q * boost::math::quaternion<double>(offset[3], offset[4], offset[5], offset[6]);

    SiconosVector fp = SiconosVector(7);
    fp(0) = (*_q)(0)+rv.R_component_2();
    fp(1) = (*_q)(1)+rv.R_component_3();
    fp(2) = (*_q)(2)+rv.R_component_4();
    fp(3) = r.R_component_1();
    fp(4) = r.R_component_2();
    fp(5) = r.R_component_3();
    fp(6) = r.R_component_4();

    occ_move(boost::get<0>(*csi)->data(), fp);
  }
}

void OccBody::updateShapes()
{

   boost::math::quaternion<double> q((*_q)(3), (*_q)(4), (*_q)(5), (*_q)(6));

  for (TopoDS_Shapes::iterator csi = _shapes->begin();
       csi != _shapes->end(); ++ csi)
  {

    OffSet offset = boost::get<1>(*csi);

    boost::math::quaternion<double> pv = boost::math::quaternion<double>(0, offset[0], offset[1], offset[2]);

    boost::math::quaternion<double> rv = q * pv * boost::math::conj(q);

    boost::math::quaternion<double> r = q * boost::math::quaternion<double>(offset[3], offset[4], offset[5], offset[6]);

    SiconosVector fp = SiconosVector(7);
    fp(0) = (*_q)(0)+rv.R_component_2();
    fp(1) = (*_q)(1)+rv.R_component_3();
    fp(2) = (*_q)(2)+rv.R_component_4();
    fp(3) = r.R_component_1();
    fp(4) = r.R_component_2();
    fp(5) = r.R_component_3();
    fp(6) = r.R_component_4();

    occ_move(*(boost::get<0>(*csi)), fp);
  }
}

const OccContactShape& OccBody::contactShape(unsigned int id) const
{
  return *boost::get<0>((*this->_contactShapes)[id]);
}

const TopoDS_Shape& OccBody::shape(unsigned int id) const
{
  return *boost::get<0>((*this->_shapes)[id]);
}
