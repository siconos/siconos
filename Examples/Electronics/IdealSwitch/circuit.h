#ifndef CIRCUIT_H
#define CIRCUIT_H

#define CLSC_CIRCUIT

//#define SICONOS_DEBUG

#ifdef CLSC_CIRCUIT
/*
 *An example from Paolo Maffezzoni nov 2006.
 *
 *
 */

static double sR = 1;
static double sR1s = 0.001;
static double sR1d = 0.001;
static double sR2 = 1000;
static double sL = 200e-6;
static double sT = 20e-6;
static double sStep = 0.1e-6;
static double sE_plus = 7.5;
static double sE_moins = -2.5;
static double sTf = 10 * sT;
static int sNSLawSize = 9;
static int sN = 4;
static int sM = 5;
static double sAmpli = 100.0;

#else
/*    V1     /      V2
 * ___|_____/  _____|________
 * |         ^               |
 * |         |               |
 * |         V1              |
 * |                         |
 *_____                     ---
 *|u=  |                       C
 *|S(t)|                    ___
 *-----                      |
 * |                         |
 * |                         |
 * |                         |
 * |__________________________
 *
 *
 * x=V2
 * L=(V1,i,L3,L4)^t=(L1,L2,L3,L4)^t=
 *
 * C dx/dt = 0.x +r
 * r = g(L)=L2
 *
 *           |-L1+s(t)
 * Y=h(x,L)= |L1-X+(L4+R1)L2
 *           |R2-L4-R1
 *           |-L1+L3
 *
 *
 *
 */
static int sNSLawSize = 4;
static int sN = 2;
static int sM = 2;

static double sR1 = 1;
static double sR2 = 1000;
static double sC = 1e-2;
static double sW = 1000;
static double sTf = 0.07;
static double sStep = 0.00001;



#endif


#endif
