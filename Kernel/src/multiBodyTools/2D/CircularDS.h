#ifndef CircularDS_h
#define CircularDS_h

#include "LagrangianDS.h"

class CircularDS : public LagrangianDS
{
protected:
  double radius;
  double massValue;

public:

  CircularDS(double, double, const SiconosVector&, const SiconosVector&);

  ~CircularDS();

  inline double getQ(unsigned int pos)
  {
    return (*q[0])(pos);
  };
  inline double getVelocity(unsigned int pos)
  {
    return (*q[1])(pos);
  };

  inline double getMassValue()
  {
    return massValue;
  };

  inline double getRadius()
  {
    return radius;
  };
};

TYPEDEF_SPTR(CircularDS);

#endif /* CircularDS_h */
