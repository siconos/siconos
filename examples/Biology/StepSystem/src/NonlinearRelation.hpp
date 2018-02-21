#ifndef NONLINEARRELATION_H
#define NONLINEARRELATION_H

#include "SiconosKernel.hpp"

class NonlinearRelation : public FirstOrderType2R
{
protected:
public:
  NonlinearRelation();
  virtual ~NonlinearRelation() {};

  /** default function to compute h */
  virtual void computeh(double t, SiconosVector& x, SiconosVector& lambda, SiconosVector& y);

  /** default function to compute g */
  virtual void computeg(double t, SiconosVector& lambda, SiconosVector& r);

  /** default function to compute jacobian of h w.r.t x and lambda */
  virtual void computeJachx(double t, SiconosVector& x, SiconosVector& lambda, SimpleMatrix& C);
  virtual void computeJachlambda(double t, SiconosVector& x, SiconosVector& lambda, SimpleMatrix& D);

  /** default function to compute jacobian of g  w.r.t  lambda  */
  virtual void computeJacglambda(double t, SiconosVector& lambda, SimpleMatrix& B);

};

TYPEDEF_SPTR(NonlinearRelation);

#endif
