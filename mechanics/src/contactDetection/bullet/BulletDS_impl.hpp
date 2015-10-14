#ifndef BulletDS_impl_hpp
#define BulletDS_impl_hpp

#include <map>

#include <boost/tuple/tuple.hpp>

#include "BulletDS.hpp"
typedef boost::array<double, 7> OffSet;

class CollisionObjects :
  public std::map<const btCollisionObject*,
                  boost::tuple<SP::btCollisionObject, OffSet , long unsigned int> >
{
private:
  ACCEPT_SERIALIZATION(CollisionObjects);
};
#endif
