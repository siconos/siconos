#ifndef FUNCODEGEN
#define FUNCODEGEN

/* Definitions are incomplete in frama-c Magnesium, so we use our
 * specs for math functions */
#ifndef __FRAMAC__
#include <math.h>
#endif

#include <float.h>
#include <assert.h>

#define Sign(x) ((double)(x>0) - (double)(x<0))
#define Max fmax
#define Heaviside(x) (x < 0 ? 0 : ((x > 0) ? 1 : .5))
#define Rand(x) ((double) rand()/ (double) RAND_MAX)

/*@ requires \is_finite((double) x) && x >= 0.;
    assigns \nothing;
    ensures \is_finite((double) \result);
    ensures 0. <= \result <= 1. + x;
    ensures \result * \result == x; */
extern double sqrt(double x);

/*@ requires \is_finite((double) x) && \is_finite((double) y);
    assigns \nothing;
    ensures \is_finite(\result) && (\result == x || \result == y) && (\result >= x && \result >= y);
 */
extern double fmax(double x, double y);

#ifdef __FRAMAC__
#include <__fc_builtin.h>
#endif

#endif
