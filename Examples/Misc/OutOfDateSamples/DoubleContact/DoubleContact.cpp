/* Siconos-sample version 2.0.1, Copyright INRIA 2005-2006.
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
#include <iostream>
#include "Model.h"
#include "check.h"

#include "LagrangianDS.h"
#include "LinearDS.h"
#include "LagrangianLinearTIDS.h"
#include "LagrangianNonLinearR.h"

#include <libxml/parser.h>
//#include "SiconosVector.h"
#include "NewSiconosVector.h"
#include "SimpleVector.h"
#include "SiconosMatrix.h"
#include "SiconosDOMTreeTools.h"
#include <math.h>
#include <stdio.h>
#include "LCP.h"



using namespace std;


void cartouche()
{
  cout << endl;
  cout << endl;
  cout << "+---------------------------------------------------------------+\n";
  cout << "+                                                               +\n";
  cout << "+                        SICONOS / WP2                          +\n";
  cout << "+                       INRIA - 2004 (c)                        +\n";
  cout << "+                                                               +\n";
  cout << "+---------------------------------------------------------------+\n";
  cout << endl;
  cout << endl;
}



int main(int argc, char* argv[])
{
  cartouche();

  try
  {

    Model doubleContact("./DoubleContact.xml");

    cout << "\n *** DoubleContact.xml loaded ***" << endl;
    //getchar();
    Strategy* s = doubleContact.getStrategy();

    cout << "the strategy will be initialized" << endl;
    s->initialize();
    cout << "\n **** the strategy is ready ****" << endl;
    //cout<<"K "<<k<<endl;

    //        static_cast<LCP*>(s->getOneStepNSProblem())->display();
    //        cout<<" @ LCP = "<<s->getOneStepNSProblem()<<endl;
    //        cout<<"Press <<enter>> to continue"<<endl;
    //        getchar();


    TimeDiscretisation* t = s->getTimeDiscretisation();
    int k = t->getK();
    int N = t->getN();

    // Trace Values
    SiconosMatrix m(N + 1, 7);

    //time
    m(k, 0) = k * t->getH();

    // position
    LagrangianLinearTIDS* ball = static_cast<LagrangianLinearTIDS*>(doubleContact.getNonSmoothDynamicalSystem()->getDynamicalSystem(0));
    m(k, 1) = (ball->getQ())(1);
    m(k, 2) = (ball->getVelocity())(1);

    // position
    LagrangianLinearTIDS* ball2 = static_cast<LagrangianLinearTIDS*>(doubleContact.getNonSmoothDynamicalSystem()->getDynamicalSystem(1));
    m(k, 3) = (ball2->getQ())(1);
    m(k, 4) = (ball2->getVelocity())(1);

    // position
    LagrangianLinearTIDS* ball3 = static_cast<LagrangianLinearTIDS*>(doubleContact.getNonSmoothDynamicalSystem()->getDynamicalSystem(2));
    m(k, 5) = (ball3->getQ())(1);
    m(k, 6) = (ball3->getVelocity())(1);

    /*SiconosVector*/
    SimpleVector vv;

    Interaction* I =  doubleContact.getNonSmoothDynamicalSystem()->getInteraction(0);
    I->display();
    //cout << "vv " << vv << endl;

    //    m(k,7) =0.0;


    while (k < N)
    {
      s->nextStep();
      cout << "NextStep done" << endl;

      k = t->getK();
      cout << "iteration : " << k << endl;

      s->computeFreeState();

      s->formaliseOneStepNSProblem();

      s->computeOneStepNSProblem();

      s->updateState();


      // Trace Values
      //time
      m(k, 0) = k * t->getH();

      // position
      LagrangianDS* ball = static_cast<LagrangianDS*>(doubleContact.getNonSmoothDynamicalSystem()->getDynamicalSystem(0));
      m(k, 1) = (ball->getQ())(1);
      m(k, 2) = (ball->getVelocity())(1);

      m(k, 3) = (ball2->getQ())(1);
      m(k, 4) = (ball2->getVelocity())(1);

      m(k, 5) = (ball3->getQ())(1);
      m(k, 6) = (ball3->getVelocity())(1);

      //  m(k, 5) = (doubleContact.getNonSmoothDynamicalSystem()->getInteraction(0)->getLambda())(0);

      //ball->display();
      //ground->display();

      //      cout<<"Press a key to begin the computations   <<enter>>"<<endl;
      //getchar();

    }
    cout << "iterations  done: " << k << endl;

    m.write("resultDC.dat", "ascii");

    doubleContact.saveToXMLFile("./DoubleContact.xml.output");

  }

  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught in \'sample/DoubleContact\'" << endl;
  }
}
