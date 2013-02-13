/*
 *
 */
/**
 *\file ace.h
  \brief This file contains declations of constants and tools used in the ACEF module.


*/

#ifndef ACE_H
#define ACE_H

/**
 *
 * define PRE_COMPUTE_ADAPTIVE to use the precompute adaptive time stepping.
 *
 *
 */

//#define PRE_COMPUTE_ADAPTIVE


#include <vector>
#include <fstream>





#include "acetime.h"


#define ACE_TIMER_MAIN 0
#define ACE_TIMER_GRAPHIC 1
#define ACE_TIMER_DIST 2
#define ACE_TIMER_UPDATE_POS 3
#define ACE_TIMER_SICONOS 4
#define ACE_TIMER_UV_TO_3D 5
#define ACE_TIMER_UV_POLYNOME 6
#define ACE_TIMER_UV_CLASS 7
#define ACE_TIMER_UV_GRAD 8
#define ACE_TIMER_3D_PROJ 9
#define ACE_TIMER_CAD_VALUE 10
#define ACE_TIMER_CAD_VALUES 11
#define ACE_TIMER_CAD_1 12
#define ACE_TIMER_CAD_12 13
#define ACE_TIMER_CAD_13 14
#define ACE_TIMER_CAD_14 15
#define ACE_TIMER_CAD_15 16
#define ACE_TIMER_CAD_16 17
#define ACE_TIMER_CAD_17 18
#define ACE_TIMER_CAD_OK 19
#define ACE_TIMER_LAST 20




extern aceTime ACE_times[];



//TIME FUNCTION
void ACE_INIT_TIME();
void ACE_PRINT_TIME();


#endif //ACE_H
