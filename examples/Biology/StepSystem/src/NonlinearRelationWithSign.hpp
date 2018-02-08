#ifndef NONLINEARRELATIONWITHSIGN_H
#define NONLINEARRELATIONWITHSIGN_H

#include "SiconosKernel.hpp"

class NonlinearRelationWithSign : public FirstOrderType2R
{
protected:
public:
  NonlinearRelationWithSign();
  virtual ~NonlinearRelationWithSign() {};

  /** default function to compute h
   *  \param double : current time
   */
  virtual void computeh(double t, SiconosVector& x, SiconosVector& lambda, SiconosVector& y);

  /** default function to compute g
   *  \param double : current time
   */
  virtual void computeg(double t, SiconosVector& lambda, SiconosVector& r);

  /** default function to compute jacobianH
   *  \param double : current time
   *  \param index for jacobian (0: jacobian according to x, 1 according to lambda)
   */
  virtual void computeJachx(double t, SiconosVector& x, SiconosVector& lambda, SimpleMatrix& C);
  virtual void computeJachlambda(double t, SiconosVector& x, SiconosVector& lambda, SimpleMatrix& D);

  /** default function to compute jacobianG according to lambda
   *  \param double : current time
   *  \param index for jacobian: at the time only one possible jacobian => i = 0 is the default value .
   */
  virtual void computeJacglambda(double t, SiconosVector& lambda, SimpleMatrix& B);

};

TYPEDEF_SPTR(NonlinearRelationWithSign);

#endif
