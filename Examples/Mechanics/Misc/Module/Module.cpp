/* Siconos version 1.0, Copyright INRIA 2005.
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
#include "SiconosKernel.h"

using namespace std;


char* itoa(char* s, int e)
{
  int tmp;
  int i = 0;

  tmp = e;
  while (tmp != 0)
  {
    tmp = tmp / 10;
    i = i + 1;
  }

  s = (char*)malloc(i + 1);
  s[i] = 0;
  i = i - 1;

  tmp = e;
  while (tmp != 0)
  {
    s[i] = tmp % 10 + '0';
    i = i - 1;
    tmp = tmp / 10;
  }

  return(s);
}



int main(int argc, char* argv[])
{

  FILE *f1, *f2;
  int i, n_ds, sizeQ;
  double dataQ[6];

  try
  {

    // --- Model loading from xml file ---
    Model module3D("./module_siconos.xml");
    cout << "\n *** module_siconos.xml file loaded ***" << endl;

    // Get number of DS in model
    int dsNumber = module3D.getNonSmoothDynamicalSystemPtr()->getNumberOfDS();
    unsigned int nDof;

    // Get ds of the model and save them into vectorDS
    vector<LagrangianDS *> vectorDS; // the list of DS
    vectorDS.resize(dsNumber, NULL);
    vector<SimpleVector *> F0;
    F0.resize(dsNumber, NULL);


    for (i = 0; i < dsNumber - 1; i++)
    {

      char *ch1 = NULL;
      char ch2[9]  = "fext" ;

      vectorDS[i] = static_cast<LagrangianDS*>(module3D.getNonSmoothDynamicalSystemPtr()->getDynamicalSystemPtr(i));
      nDof = vectorDS[i]->getNdof();

      F0[i] = new SimpleVector(nDof);
      ch1 = itoa(ch1, i + 1);
      strcat(ch2, ch1);

      F0[i] = new SimpleVector(strcat(ch2, ".dat"), 1);

      printf("taille FO %d\n", F0[i]->size());

      for (unsigned int j = 0; j < F0[i]->size(); j++)
        printf("DSPlugin ... %f\n", (*(F0[i]))(j));

      // 2 - Assign this param to the function FExt
      vectorDS[i]->setZPtr(F0[i]);
    }

    // --- Get and initialize the simulation ---
    TimeStepping* s = static_cast<TimeStepping*>(module3D.getSimulationPtr());

    cout << "simulation initialization" << endl;

    s->initialize();

    cout << "\n **** the simulation is ready ****" << endl;

    // --- Get the time discretisation scheme ---
    TimeDiscretisation* t = s->getTimeDiscretisationPtr();
    int k = 0;
    int N = t->getNSteps(); // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    SimpleMatrix dataPlot(N + 1, 9); // h + ddl bar bloquee

    cout << "Prepare data for plotting ... " << endl;
    // For the initial time step:
    // state q for the first dynamical system (bar2D)

    f1 = fopen("module_time.dat", "w+") ;
    f2 = fopen("module_deplacement.dat", "w+");


    n_ds = 4;

    fprintf(f1, "%d\n", N + 1);
    fprintf(f2, "%d %d\n", N + 1, n_ds);

    //time
    dataPlot(k, 0) = k * t->getH();

    for (i = 0; i < 8; i++)
    {

      dataPlot(k, i + 1) = (module3D.getNonSmoothDynamicalSystemPtr()->getInteractionPtr(i)->getLambda(1))(0);

    }


    fprintf(f1, "%10.4e\n", k * t->getH());



    for (i = 0; i < 4; i++)
    {

      LagrangianDS* module = static_cast<LagrangianDS*>(module3D.getNonSmoothDynamicalSystemPtr()->getDynamicalSystemPtr(i));

      sizeQ = 3;


      if (sizeQ == 3)
      {

        dataQ[0] = 0.;
        dataQ[1] = 0.;
        dataQ[2] = 0.;
        dataQ[3] = (module->getQ())(0);
        dataQ[4] = (module->getQ())(1);
        dataQ[5] = (module->getQ())(2);


      }
      else if (sizeQ == 6)
      {

        dataQ[0] = (module->getQ())(0);
        dataQ[1] = (module->getQ())(1);
        dataQ[1] = (module->getQ())(2);
        dataQ[3] = (module->getQ())(3);
        dataQ[4] = (module->getQ())(4);
        dataQ[5] = (module->getQ())(5);


      }
      else
      {

        printf(" Not available \n");
        break;

      }

      fprintf(f2, "%10.4e %10.4e %10.4e %10.4e %10.4e %10.4e\n", dataQ[0], dataQ[1], dataQ[2], dataQ[3], dataQ[4], dataQ[5]);


    }

    // --- Compute elapsed time ---
    double t1, t2, elapsed;
    struct timeval tp;
    int rtn;
    clock_t start, end;
    double elapsed2;
    start = clock();
    rtn = gettimeofday(&tp, NULL);
    t1 = (double)tp.tv_sec + (1.e-6) * tp.tv_usec;

    n_ds = 4;

    cout << "Computation ... " << endl;

    // --- Time loop  ---
    while (k < N)
    {
      k++;
      // solve ...
      s->computeOneStep();

      // --- Get values to be plotted ---
      //time
      dataPlot(k, 0) = k * t->getH();
      fprintf(f1, "%10.4e\n", k * t->getH());




      for (i = 0; i < 8; i++)
      {

        dataPlot(k, i + 1) = (module3D.getNonSmoothDynamicalSystemPtr()->getInteractionPtr(i)->getLambda(1))(0);

      }


      for (i = 0; i < 4; i++)
      {

        LagrangianDS* module = static_cast<LagrangianDS*>(module3D.getNonSmoothDynamicalSystemPtr()->getDynamicalSystemPtr(i));


        if (sizeQ == 3)
        {

          dataQ[0] = 0.;
          dataQ[1] = 0.;
          dataQ[2] = 0.;
          dataQ[3] = (module->getQ())(0);
          dataQ[4] = (module->getQ())(1);
          dataQ[5] = (module->getQ())(2);


        }
        else if (sizeQ == 6)
        {

          dataQ[0] = (module->getQ())(0);
          dataQ[1] = (module->getQ())(1);
          dataQ[1] = (module->getQ())(2);
          dataQ[3] = (module->getQ())(3);
          dataQ[4] = (module->getQ())(4);
          dataQ[5] = (module->getQ())(5);


        }
        else
        {

          printf(" Not available \n");
          break;

        }

        fprintf(f2, "%10.4e %10.4e %10.4e %10.4e %10.4e %10.4e\n", dataQ[0], dataQ[1], dataQ[2], dataQ[3], dataQ[4], dataQ[5]);

      }

      // transfer of state i+1 into state i and time incrementation
      s->nextStep();


    }



    fclose(f1);
    fclose(f2);

    // --- elapsed time computing ---
    end = clock();

    rtn = gettimeofday(&tp, NULL);
    t2 = (double)tp.tv_sec + (1.e-6) * tp.tv_usec;
    elapsed = t2 - t1;
    elapsed2 = (end - start) / (double)CLOCKS_PER_SEC;
    cout << "time = " << elapsed << " --- cpu time " << elapsed2 << endl;


    // Number of time iterations
    cout << "Number of iterations done: " << k << endl;

    // dataPlot (ascii) output
    ioMatrix io("result.dat", "ascii");
    io.write(dataPlot, "noDim");

    for (i = 0; i < dsNumber; i++)
      delete F0[i];
  }


  // --- Exceptions handling ---
  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)

  {
    cout << "Exception caught in \'sample/MODULE\'" << endl;
  }


}
