#ifndef BulletSpaceFilter_impl_hpp
#define BulletSpaceFilter_impl_hpp


#include <map>
#include <Question.hpp>

struct StaticObjects :
  public std::map<const btCollisionObject*,
                  std::pair<SP::btCollisionObject, long unsigned int> >
{
  ACCEPT_SERIALIZATION(StaticObjects);
};


struct ForStaticObjects : public Question< SP::StaticObjects >
{
  ANSWER(BulletSpaceFilter, staticObjects());
};


#endif
