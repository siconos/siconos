#ifndef ADJOINTINPUT_H
#define ADJOINTINPUT_H

#include "SiconosKernel.hpp"

class adjointInput : public FirstOrderNonLinearR
{
protected:
  SP::SimpleMatrix  K2;

public:
  adjointInput();
  virtual ~adjointInput() {};


  virtual void initialize(Interaction& inter);


  /** default function to compute h
   *  \param double : current time
   */
  virtual void computeh(double time, SiconosVector& x, SiconosVector& lambda, SiconosVector& z, SiconosVector& y);

  /** default function to compute g
   *  \param double time, Interaction& inter : current time
   */
  virtual void computeg(double time, SiconosVector& x, SiconosVector& lambda, SiconosVector& z, SiconosVector& r);

   /** default function to compute jacobianH
   *  \param double time, Interaction& inter : current time
   *  \param index for jacobian (0: jacobian according to x, 1 according to lambda)
   */
  virtual void computeJachx(double time, SiconosVector& x, SiconosVector& lambda, SiconosVector& z, SimpleMatrix& C);
  virtual void computeJachlambda(double time, SiconosVector& x, SiconosVector& lambda, SiconosVector& z, SimpleMatrix& D);

  /** default function to compute jacobianG according to lambda
   *  \param double time, Interaction& inter : current time
   *  \param index for jacobian: at the time only one possible jacobian => i = 0 is the default value .
   */
  virtual void computeJacgx(double time, SiconosVector& x, SiconosVector& lambda, SiconosVector& z, SimpleMatrix& K);
  virtual void computeJacglambda(double time, SiconosVector& x, SiconosVector& lambda, SiconosVector& z, SimpleMatrix& B);


  double source(double t);

  void beta(double t, SiconosVector& xvalue, SP::SiconosVector alpha);

  void JacobianXbeta(double t, SiconosVector& xvalue, SP::SimpleMatrix JacbetaX);

};

TYPEDEF_SPTR(adjointInput);

#endif
