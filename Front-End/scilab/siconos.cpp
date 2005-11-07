/* Siconos version 1.0, Copyright INRIA 2005.
 * Siconos is a program dedicated to the modeling, the simulation and the control
 * of non smooth dynamical systems
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
#include <iostream>
#include "Model.h"
#include "check.h"

#include "LagrangianDS.h"
#include "LinearDS.h"
//#include "LinearSystemDS.h"

#include "LagrangianLinearTIDS.h"
//#include "LagrangianNonLinearR.h"
#include "LagrangianR.h"

#include <libxml/parser.h>
#include "SiconosVector.h"
//#include "NewSiconosVector.h"
#include "SimpleVector.h"
#include "SiconosMatrix.h"
#include "SiconosDOMTreeTools.h"
#include <math.h>
#include <stdio.h>
#include "LCP.h"

#include <sys/time.h>




using namespace std;


void cartouche()
{
  cout << endl;
  cout << endl;
  cout << "+---------------------------------------------------------------+\n";
  cout << "+                                                               +\n";
  cout << "+                        SICONOS / WP2                          +\n";
  cout << "+                       INRIA - 2005 (c)                        +\n";
  cout << "+                                                               +\n";
  cout << "+---------------------------------------------------------------+\n";
  cout << endl;
  cout << endl;
}

//
// TODO: add a structure with array and manage index
//
Model *GLOB_MODEL;
Strategy *GLOB_STRATEGIE;
FILE *GLOB_FD;

extern "C" void sicLoadModel(int *ret, char ModelXmlFile[])
{

  GLOB_FD = fopen("result.dat", "w");

  if (GLOB_FD == NULL)
  {
    printf("error:: result.dat write\n");
  }

  try
  {

    GLOB_MODEL = new Model(ModelXmlFile);
    *ret = 0;

  }
  catch (SiconosException e)
  {
    cout << e.report() << endl;
    *ret = -1;
  }
  catch (...)
  {
    cout << "Exception caught in sicLoadModel" << endl;
  }

}

extern "C" void sicInitStrategy()
{
  try
  {
    if (GLOB_MODEL != NULL)
    {
      GLOB_STRATEGIE = GLOB_MODEL->getStrategyPtr();
      GLOB_STRATEGIE->initialize();
    }
    else
    {
      RuntimeException::selfThrow("siconos/C:: sicInitStrategy failed");
    }

  }
  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught in sicInitStrategy" << endl;
  }
}

extern "C" void sicTimeGetH(double *H)
{
  *H = GLOB_STRATEGIE->getTimeDiscretisationPtr()->getH();
}

extern "C" void sicTimeGetN(int *N)
{
  *N = GLOB_STRATEGIE->getTimeDiscretisationPtr()->getNSteps();
}

extern "C" void sicTimeGetK(int *K)
{
  *K = GLOB_STRATEGIE->getTimeDiscretisationPtr()->getK();
}

extern "C" void  sicSTNextStep(int *ret)
{
  GLOB_STRATEGIE->nextStep();

  *ret = 0;
}

extern "C" void sicSTComputeFreeState(int *ret)
{
  GLOB_STRATEGIE->computeFreeState();

  *ret = 0;
}

// extern "C" void sicSTformalisePb(int *ret)
// {

// try{

//   GLOB_STRATEGIE->formaliseOneStepNSProblem();
//   *ret = 0;

// }
//    catch(SiconosException e)
//     {cout << e.report() << endl;}
//    catch(...)
//     {cout << "Exception caught in sicSTformalisePb" << endl;}


// }

extern "C" void sicSTcomputePb(int *ret)
{
  GLOB_STRATEGIE->computeOneStepNSProblem();

  *ret = 0;
}

extern "C" void sicSTupdateState(int *ret)
{
  GLOB_STRATEGIE->update();

  *ret = 0;
}


extern "C" void sicDebug(int *ret)
{
  // GLOB_STRATEGIE->update();
  cout << "--- Debug" << endl;

  // LagrangianDS* ball1 = static_cast<LagrangianDS*> (GLOB_MODEL->getNonSmoothDynamicalSystemPtr()->getDynamicalSystemPtr(0));

  //   vector<Interaction*> vIS= (GLOB_MODEL->getNonSmoothDynamicalSystemPtr()->getInteractions());

  //   for (int index = 0; index < vIS.size(); index++)
  //     {
  //       vIS[index]->display();
  //     }

  fclose(GLOB_FD);

  *ret = 0;
}

extern "C" void sicModelgetQ(double *value, int *index)
{


  try
  {

    LagrangianDS* system = static_cast<LagrangianDS*>(GLOB_MODEL->getNonSmoothDynamicalSystemPtr()->getDynamicalSystemPtr(*index));

    if (system != NULL)
    {
      int size = (system->getQ()).size();
      // cout <<"SIZE "<<size<<endl;
      // for(int i=0;i<size;i++)
      //  value[i]=system->getQ()(i);
      *value = system->getQ()(0);
    }

  }
  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught in sicSTformalisePb" << endl;
  }

}
extern "C" void simul(int *ret)
{
  FILE *fd;
  double plot[4], H;
  int index;


  try
  {

    // --- Model loading from xml file ---
    //Model uBeads("./ThreeBeadsColumn.xml");
    //cout << "\n *** ThreeBeadsColumn.xml loaded ***" << endl;

    Model *uBeads = GLOB_MODEL;
    Strategy *s = GLOB_STRATEGIE;


    if ((uBeads == NULL) || (s == NULL))
    {
      *ret = -1;
      RuntimeException::selfThrow("siconos/C:: simul failed, GLOB_MODEL not created");
    }
    // --- Get and initialize the strategy ---

    cout << "\n **** the strategy is ready ****" << endl;

    // --- Get the time discretisation scheme ---
    TimeDiscretisation* t = s->getTimeDiscretisationPtr();
    int k = t->getK(); // Current step
    int N = t->getNSteps(); // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    //SiconosMatrix dataPlot(N+1, 11);


    cout << " Computation 1... " << endl;
    // while(k < N)
    // {
    // transfer of state i+1 into state i and time incrementation
    s->nextStep();

    // get current time step
    k = t->getK();

    // solve ...
    s->computeFreeState();
    s->computeOneStepNSProblem();

    // update
    s->update();

    sicTimeGetH(&H);

    plot[0] = k * H;
    index = 0;
    sicModelgetQ(&plot[1], &index);
    index = 1;
    sicModelgetQ(&plot[2], &index);
    index = 2;
    sicModelgetQ(&plot[3], &index);

    fprintf(GLOB_FD, "%lf %lf %lf %lf \n", plot[0], plot[1], plot[2], plot[3]);
    //  }
    cout << "End of computation - Number of iterations  done: " << k << endl;
  }

  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught in \'sample/ThreeBeadsColumn\'" << endl;
  }

  *ret = 0;
}

extern "C" void simul1(int *ret)
{

  try
  {

    // --- Model loading from xml file ---
    Model uBeads("./ThreeBeadsColumn.xml");
    cout << "\n *** ThreeBeadsColumn.xml loaded ***" << endl;

    // --- Get and initialize the strategy ---
    Strategy* s = uBeads.getStrategyPtr();
    s->initialize();
    cout << "\n **** the strategy is ready ****" << endl;

    // --- Get the time discretisation scheme ---
    TimeDiscretisation* t = s->getTimeDiscretisationPtr();
    int k = t->getK(); // Current step
    int N = t->getNSteps(); // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    SiconosMatrix dataPlot(N + 1, 11);

    cout << "Prepare data for plotting ... " << endl;
    // For the initial time step:
    // time
    dataPlot(k, 0) = k * t->getH();

    // state q and velocity for the first dynamical system
    LagrangianLinearTIDS* bead = static_cast<LagrangianLinearTIDS*>(uBeads.getNonSmoothDynamicalSystemPtr()->getDynamicalSystemPtr(0));
    dataPlot(k, 1) = (bead->getQ())(0);
    dataPlot(k, 2) = (bead->getVelocity())(0);

    // state q and velocity for the second dynamical system
    LagrangianLinearTIDS* bead2 = static_cast<LagrangianLinearTIDS*>(uBeads.getNonSmoothDynamicalSystemPtr()->getDynamicalSystemPtr(1));
    dataPlot(k, 3) = (bead2->getQ())(0);
    dataPlot(k, 4) = (bead2->getVelocity())(0);

    // state q and velocity for the third dynamical system
    LagrangianLinearTIDS* bead3 = static_cast<LagrangianLinearTIDS*>(uBeads.getNonSmoothDynamicalSystemPtr()->getDynamicalSystemPtr(2));
    dataPlot(k, 5) = (bead3->getQ())(0);
    dataPlot(k, 6) = (bead3->getVelocity())(0);

    // state q and velocity for the fourth dynamical system (ground)
    LagrangianLinearTIDS* ground = static_cast<LagrangianLinearTIDS*>(uBeads.getNonSmoothDynamicalSystemPtr()->getDynamicalSystemPtr(3));
    dataPlot(k, 7) = (ground->getQ())(0);
    dataPlot(k, 8) = (ground->getVelocity())(0);

    // state q and velocity for the fourth dynamical system (ceiling)
    LagrangianLinearTIDS* ceiling = static_cast<LagrangianLinearTIDS*>(uBeads.getNonSmoothDynamicalSystemPtr()->getDynamicalSystemPtr(4));
    dataPlot(k, 9) = (ceiling->getQ())(0);
    dataPlot(k, 10) = (ceiling->getVelocity())(0);

    cout << " Computation ... " << endl;
    while (k < N)
    {
      // transfer of state i+1 into state i and time incrementation
      s->nextStep();

      // get current time step
      k = t->getK();

      // solve ...
      s->computeFreeState();
      s->computeOneStepNSProblem();

      // update
      s->update();

      // --- Get values to be plotted ---

      dataPlot(k, 0) = k * t->getH();
      dataPlot(k, 1) = (bead->getQ())(0);
      dataPlot(k, 2) = (bead->getVelocity())(0);
      dataPlot(k, 3) = (bead2->getQ())(0);
      dataPlot(k, 4) = (bead2->getVelocity())(0);
      dataPlot(k, 5) = (bead3->getQ())(0);
      dataPlot(k, 6) = (bead3->getVelocity())(0);
      dataPlot(k, 7) = (ground->getQ())(0);
      dataPlot(k, 8) = (ground->getVelocity())(0);
      dataPlot(k, 9) = (ceiling->getQ())(0);
      dataPlot(k, 10) = (ceiling->getVelocity())(0);
      //  dataPlot(k, 5) = (uBeads.getNonSmoothDynamicalSystemPtr()->getInteraction(0)->getLambda())(0);
    }
    cout << "End of computation - Number of iterations  done: " << k << endl;

    dataPlot.rawWrite("result.dat", "ascii");

    //    uBeads.saveToXMLFile("./ThreeBeadsColumn.xml.output");

  }

  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught in \'sample/ThreeBeadsColumn\'" << endl;
  }

  *ret = 0;
}
