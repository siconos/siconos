#ifndef UTILITIES_H
#define UTILITIES_H

#define PI 3.1415926535

/* Sum of beads by floor */
unsigned int sum(unsigned int);

/* Sum of beads of all floors */
unsigned int SUM(unsigned int);

unsigned int qlq(unsigned int, unsigned int);

unsigned int L(unsigned int, unsigned);

/* Function to compute the 3 beads under each bead */
void Q(unsigned int, unsigned int, unsigned int, double *, double, double, double, double);

#endif
