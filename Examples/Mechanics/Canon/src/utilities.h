#ifndef UTILITIES_H
#define UTILITIES_H

#define PI 3.1415926535

/* Number of beads in floor N (floor 1 = top floor) */
unsigned int sum(unsigned int N);

/* Number of beads between floor 1 and floor N (both included) */
unsigned int SUM(unsigned int N);

unsigned int qlq(unsigned int, unsigned int);

unsigned int L(unsigned int, unsigned);

/* Function to compute the 3 beads under each bead */
void Q(unsigned int, unsigned int, unsigned int, double *, double, double, double, double);

#endif
