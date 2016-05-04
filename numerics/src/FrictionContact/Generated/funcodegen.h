#ifndef FUNCODEGEN
#define FUNCODEGEN

/* Definitions are incomplete in frama-c Magnesium, so we use our
 * specs for math functions */
#ifndef __FRAMAC__
#include <math.h>
#endif

#include <float.h>
#include <assert.h>

/*@
  axiomatic Max {
  logic real fmax(real x, real y);
  axiom fmax_1: \forall real x, real y; fmax(x, y) >= x && fmax(y, y) >= y;
  axiom fmax_2: \forall real x, real y; fmax(x, y) == x || fmax(y, y) == y;
  }
  axiomatic sqrt {
  logic real sqrt(real x);
  axiom sqrt_1: \forall real x; x >= 0 <==> x == sqrt(x) * sqrt(x);
  axiom sqrt_2: \forall real x; x > 0 <==> sqrt(x) > 0;
  }
  axiomatic Heaviside {
  logic real Heaviside(real x);
  axiom Heaviside_1: \forall real x; (x < 0 ==> Heaviside(x) == 0);
  axiom Heaviside_2: \forall real x; (x > 0 ==> Heaviside(x) == 1);
  axiom Heaviside_3: \forall real x; (x == 0 ==> 0 < Heaviside(x) < 1);
  }
  axiomatic pow {
  logic real pow(real x, real y);
  axiom pow_1: \forall real x, real y; x >=0 ==> pow(x, y) >= 0;
  }
  axiomatic general {
  axiom sq: \forall real x; x*x >= 0.;
  }
*/

/*@ lemma one: \forall real mu; 2*mu*mu - 2*mu + 1 >= 0.5; */

#define Sign(x) ((double)(x>0) - (double)(x<0))
#define Max fmax
#define Abs(x) (x < 0 ? -x : x)
#define Heaviside(x) (x < 0 ? 0 : ((x > 0) ? 1 : .5))

/*@ requires \is_finite((double) x);
    requires x >= 0.;
    assigns \nothing;
    ensures \is_finite((double) \result);
    ensures 0. <= \result <= 1. + x;
    ensures x <= \result <= 1. || 1. <= \result <= x;
    ensures \result * \result == x;
    ensures \result > 2.22044604925e-16 ==> x > 4.930380657631323783823303533017413935457540219431394e-32;
*/
extern double sqrt(double x);

/*@ requires \is_finite((double) x);
    requires \is_finite((double) y);
    assigns \nothing;
    ensures \is_finite(\result);
    ensures (\result == x || \result == y);
    ensures \result >= x;
    ensures \result >= y;
 */
extern double fmax(double x, double y);

#ifndef epsilon
#define epsilon 0x1.0000000000000p-52
#endif

#ifdef __FRAMAC__
#include <__fc_builtin.h>
#endif

#endif
