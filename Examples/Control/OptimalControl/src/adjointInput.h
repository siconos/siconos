#ifndef ADJOINTINPUT_H
#define ADJOINTINPUT_H

#include "SiconosKernel.hpp"

class adjointInput : public FirstOrderType2R
{
protected:
  SimpleMatrix   * K2 ;

public:
  adjointInput();
  virtual ~adjointInput() {};


  virtual void initialize(SP::Interaction inter);


  /** default function to compute h
   *  \param double : current time
   */
  virtual void computeh(double) ;

  /** default function to compute g
   *  \param double : current time
   */
  virtual void computeg(double) ;

  /** default function to compute jacobianH
   *  \param double : current time
   *  \param index for jacobian (0: jacobian according to x, 1 according to lambda)
   */
  virtual void computeJachx(double);
  virtual void computeJachlambda(double);

  /** default function to compute jacobianG according to lambda
   *  \param double : current time
   *  \param index for jacobian: at the time only one possible jacobian => i = 0 is the default value .
   */
  virtual void computeJacgx(double);
  virtual void computeJacglambda(double);


  double source(double t);

  void beta(double t, SP::SiconosVector xvalue, SP::SiconosVector alpha);

  void JacobianXbeta(double t, SP::SiconosVector xvalue, SP::SiconosMatrix JacbetaX);

};

TYPEDEF_SPTR(adjointInput);

#endif
