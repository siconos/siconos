/* Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
 */

/*! \file DrawSiconos.h
  \brief Driver for Siconos Simulation, using drawstuff for visualisation
*/

#ifndef DrawSiconos_H
#define DrawSiconos_H

#include "CanonBallsModel.h"

/** DrawSiconos

    \author F. Perignon
    \version 3.0.0.
    \date (Creation) May 2008



*/
/* class DrawSiconos */
/* { */
/*  private:  */

/*   /\** *\/ */
/*   bool GLOB_COMPUTE; */
/*   bool GLOB_STEP; */

/*   /\** The object used to define the Siconos Model *\/ */
/*   BeadsModel* canon; */

/*  public: */

/*   /\** Basic constructor *\/ */
/*  DrawSiconos() : GLOB_COMPUTE(false), GLOB_STEP(false), canon(NULL){}; */

/*   /\** Destructor *\/ */
/*   ~DrawSiconos(); */

/* } */


/** Function called before the simulation: build model, initialize ... */
void Start();

/** Draw environment (move this elsewhere?) */
void DrawBox(float alpha);

/** Run Siconos simulation */
void computeSiconos();

/** call computeSiconos + draw */
void SimuLoop(int pause);

/** Command list display */
void Command(int cmd);

/** Actions performs when the simulation is over */
void End();

/** Start, run, draw and end ... */
void run(int argc, char* argv[]);

#endif
