#include "OccR.hpp"
#include "OccR_impl.hpp"
#include "OccContactShape.hpp"

#include <limits>

OccR::OccR(const ContactPoint& contact1,
           const ContactPoint& contact2) :
  NewtonEulerFrom3DLocalFrameR(),
  _contact1(contact1), _contact2(contact2),
  _normalFromFace1(true),
  _offsetp1(false),
  _offset(0.)
{
}

void OccR::computeh(double time, BlockVector& q0, SiconosVector& y)
{
  double X1,X2,Y1,Y2,Z1,Z2,n1x,n1y,n1z;

  double dist = std::numeric_limits<double>::infinity();

  _contact1.shape.distance(_contact2.shape,
                           X1,Y1,Z1,
                           X2,Y2,Z2,
                           n1x,n1y,n1z,
                           _normalFromFace1, dist);

  //_Pc1->setValue(0,X1); _Pc1->setValue(1,Y1); _Pc1->setValue(2,Z1);
  if(_offsetp1)
  {
    _Pc1->setValue(0,X1+_offset*n1x);
    _Pc1->setValue(1,Y1+_offset*n1y);
    _Pc1->setValue(2,Z1+_offset*n1z);
    _Pc2->setValue(0,X2);
    _Pc2->setValue(1,Y2);
    _Pc2->setValue(2,Z2);
  }
  else
  {
    _Pc1->setValue(0,X1);
    _Pc1->setValue(1,Y1);
    _Pc1->setValue(2,Z1);
    _Pc2->setValue(0,X2-_offset*n1x);
    _Pc2->setValue(1,Y2-_offset*n1y);
    _Pc2->setValue(2,Z2-_offset*n1z);
  }

  /*Because in CAD model, the normal is going outside of the body.*/
  _Nc->setValue(0,-n1x);
  _Nc->setValue(1,-n1y);
  _Nc->setValue(2,-n1z);


  dist -= _offset;

  y.setValue(0, dist);

}
