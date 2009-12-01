#ifndef ADJOINTINPUT_H
#define ADJOINTINPUT_H

#include "SiconosKernel.hpp"

class adjointInput : public FirstOrderType2R
{
protected:
public:
  adjointInput();
  virtual ~adjointInput() {};


  virtual void initialize(SP::Interaction inter);


  /** default function to compute h
   *  \param double : current time
   */
  virtual void computeH(double) ;

  /** default function to compute g
   *  \param double : current time
   */
  virtual void computeG(double) ;

  /** default function to compute jacobianH
   *  \param double : current time
   *  \param index for jacobian (0: jacobian according to x, 1 according to lambda)
   */
  virtual void computeJacXH(double);
  virtual void computeJacLH(double);

  /** default function to compute jacobianG according to lambda
   *  \param double : current time
   *  \param index for jacobian: at the time only one possible jacobian => i = 0 is the default value .
   */
  virtual void computeJacXG(double);
  virtual void computeJacLG(double);


  double source(double t);

};

TYPEDEF_SPTR(adjointInput);

#endif
