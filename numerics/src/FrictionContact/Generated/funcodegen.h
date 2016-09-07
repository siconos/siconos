#ifndef FUNCODEGEN
#define FUNCODEGEN

/* Definitions are incomplete in frama-c Magnesium, so we use our
 * specs for math functions */
#ifndef __FRAMAC__
#include <math.h>
#endif

#include <float.h>
#include <assert.h>

#ifndef epsilon
#define epsilon 0x1.0000000000000p-52
#endif

/*@ lemma one: \forall real mu; 2*mu*mu - 2*mu + 1 >= 0.5; */

#define Sign(x) ((double)(x>0) - (double)(x<0))
#define Max fmax
#define Abs(x) (x < 0 ? -x : x)
#define Heaviside(x) (x < 0 ? 0 : ((x > 0) ? 1 : 0.))

/*@ requires \is_finite((double) x);
    requires x >= 0.;
    assigns \nothing;
    ensures \is_finite((double) \result);
    ensures 0. <= \result <= 1. + x;
    ensures x <= \result <= 1. || 1. <= \result <= x;
    ensures \result * \result == x;
    ensures \result > epsilon ==> x > epsilon*epsilon;
    ensures x > epsilon ==> \result > epsilon;
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



#ifdef __FRAMAC__
#include <__fc_builtin.h>
#endif

#endif
