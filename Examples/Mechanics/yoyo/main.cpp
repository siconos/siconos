#include "SiconosKernel.hpp"
#include "donnee.h"

using namespace std;


int main(int argc, char* argv[])
{
  try
  {

    unsigned int nDof = 3;  // nombre de degrés de liberté du yoyo
    double t0 = 0;      //  instants initial et final de la simulation
    double T = 50;
    double h = 0.001;       // pas de discrétisation du temps
    const double theta = 0.5;   // coefficient pour le générateur de simulation
    double e = 0; // coefficient de restitution
    int N = (int)((T - t0) / h);


    // déclarations

    DynamicalSystemsSet allDS;
    SiconosMatrix *M = new SimpleMatrix(nDof, nDof);
    SiconosVector * q0 = new SimpleVector(nDof);
    SiconosVector* v0 = new SimpleVector(nDof);
    LagrangianDS* yoyo;
    InteractionsSet allInteractions;

    NonSmoothLaw * loi = new NewtonImpactNSL(e);
    Relation * relation = new LagrangianRheonomousR("YoyoPlugin:h1", "YoyoPlugin:G11", "YoyoPlugin:G10");
    Interaction * inter ;
    NonSmoothDynamicalSystem * system;
    TimeDiscretisation * t;
    TimeStepping* s ;
    Moreau * OSI ;
    Model * jeu ;
    OneStepNSProblem * osnspb;

    SiconosMatrix *H = new SimpleMatrix(1, nDof);
    (*H)(0, 0) = 0;
    (*H)(0, 1) = 0;
    (*H)(0, 2) = 0;
    NonSmoothLaw * loi0 = new NewtonImpactNSL(e);
    Relation * relation0 = new LagrangianLinearTIR(*H);


    // pramètres du solveur siconos

    IntParameters iparam(5);
    iparam[0] = 1000; // Max number of iteration
    DoubleParameters dparam(5);
    dparam[0] = 1e-15; // Tolerance
    string solverName = "Lemke" ;
    NonSmoothSolver * mySolver = new NonSmoothSolver(solverName, iparam, dparam);


    unsigned int outputSize = 9;
    SimpleMatrix dataPlot(N, outputSize);


    SiconosVector * q ;
    SiconosVector * v ;
    SiconosVector * p ;
    SiconosVector * lambda;


    boost::progress_display show_progress(N);
    boost::timer time;
    time.restart();


    (*q0)(0) = L / (2 * r); //  vlaeur de  teta( 0 )
    (*q0)(1) = -L + r * (*q0)(0); // y ( 0 )
    (*q0)(2) = 0;  // valeur initiale de h
    // vitesses initiels de teta , y et  h
    (*v0)(0) = 0;
    (*v0)(1) = r * (*v0)(0);
    (*v0)(2) = 0;

    // Objectifs du joueurs avec les instants correspendants
    SimpleMatrix Controle(G + 1, 2);
    Controle(0, 1) = (*q0)(1) - (*q0)(2);
    for (int i = 0; i < G ; i++)
    {
      Controle(i + 1, 0) = temps[i];
      Controle(i + 1, 1) = Som[i] - L;
    }

    int k = 0;

    while (k < N)
    {


      ///////////////////////////////////////Phase contrainte //////////////////////////////////////

      M->eye();
      (*M)(0, 0) = I + m * r * r;
      (*M)(0, 2) = m * r;
      (*M)(1, 0) = -r;
      (*M)(1, 1) = 1;
      (*M)(1, 2) = -1;

      if (k != 0)
      {
        t0 = dataPlot(k - 1, 0);
        (*q0)(0) = (*q)(0);
        (*q0)(2) = (*q)(2);
        (*q0)(1) = (*q0)(2) - L + r * (*q0)(0);
        (*v0)(0) = (*v)(0);
        (*v0)(2) = (*v)(2);
        (*v0)(1) = r * (*v0)(0) + (*v0)(2);


        allDS.clear();
        allInteractions.clear();
        delete inter;
        delete osnspb;
        delete t;
        delete s;
        delete OSI;
        delete jeu;
        delete system;
        delete yoyo;
      }

      // création et insertion  du système dynamique représentant la yoyo dans le récipient allDS
      yoyo = new LagrangianDS(0, *q0, *v0, *M);

      yoyo->setComputeFExtFunction("YoyoPlugin.so", "force_ext");
      yoyo->setComputeFIntFunction("YoyoPlugin.so", "F_int");
      yoyo->setComputeJacobianFIntFunction(1, "YoyoPlugin.so", "jacobianVFInt");
      yoyo->setComputeJacobianFIntFunction(0, "YoyoPlugin.so", "jacobianQFInt");

      allDS.insert(yoyo);


      ////////////////  loi d'impact et relations /////////////////////////////////

      inter = new Interaction("impact", allDS, 0, 1, loi0, relation0);
      allInteractions.insert(inter);
      system = new NonSmoothDynamicalSystem(allDS, allInteractions);


      /////////////////////////  MODEL //////////////////////////////////////////////////

      jeu = new Model(t0, T);
      jeu->setNonSmoothDynamicalSystemPtr(system); // set NonSmoothDynamicalSystem of this model

      ///////////////////// SIMULATION /////////////////////////////////


      // déscrétisation du temps
      t = new TimeDiscretisation(h, jeu);
      s = new TimeStepping(t);
      // -- OneStepIntegrators --
      OSI = new Moreau(yoyo, theta, s);
      // -- OneStepNsProblem --
      osnspb = new LCP(s, mySolver);
      // --- Simulation initialization ---
      s->initialize();


      q = yoyo->getQPtr();
      v = yoyo->getVelocityPtr();
      p = yoyo->getPPtr();
      lambda = inter->getLambdaPtr(1);


      // --- sauver les valeurs dans une matrice dataPlot

      dataPlot(k, 0) = jeu->getT0();
      dataPlot(k, 1) = (*q)(0);
      dataPlot(k, 2) = (*v)(0);
      dataPlot(k, 3) = (*p)(0);
      dataPlot(k, 4) = (*lambda)(0);
      dataPlot(k, 5) = (*q)(1);
      dataPlot(k, 6) = (*v)(1);
      dataPlot(k, 7) =  L + (*q)(1) - r * (*q)(0) - (*q)(2);
      dataPlot(k, 8) = (*q)(1) - (*q)(2);
      k++;


      while (k < N  && (*q)(0) > 0.0)
      {
        s->computeOneStep();
        // --- Get values to be plotted ---
        dataPlot(k, 0) =  s->getNextTime();
        dataPlot(k, 1) = (*q)(0);
        dataPlot(k, 2) = (*v)(0);
        dataPlot(k, 3) = (*p)(0);
        dataPlot(k, 4) = (*lambda)(0); // multiplicateur de lagrange
        dataPlot(k, 5) = (*q)(1);
        dataPlot(k, 6) = (*v)(1);
        dataPlot(k, 7) =  L + (*q)(1) - r * (*q)(0) - (*q)(2); // contrainte géométrique
        dataPlot(k, 8) = (*q)(1) - (*q)(2);
        if (abs((*v)(0)) <= 0.05) cout << "valeur max de teta est : "  << (*q)(0) << endl;
        k++;
        s->nextStep();
        ++show_progress;
      }


      //////////////////////////////////////////Phase libre//////////////////////////////////////////////////////


      M->eye();
      (*M)(1, 1) = m;
      (*M)(0, 0) = I;

      t0 = dataPlot(k - 1, 0) + h;
      (*q0)(0) = 0;
      (*q0)(2) = (*q)(2);
      (*q0)(1) = (*q0)(2) - L;
      (*v0)(0) = -(*v)(0);
      (*v0)(2) = (*v)(2);
      (*v0)(1) = -r * (*v0)(0) + (*v0)(2);

      // libération de mémoire
      allDS.clear();
      allInteractions.clear();
      delete inter;
      delete osnspb;
      delete t;
      delete s;
      delete OSI;
      delete jeu;
      delete system;
      delete yoyo;



      yoyo = new LagrangianDS(0, *q0, *v0, *M);
      yoyo->setComputeFExtFunction("YoyoPlugin.so", "force_extf");
      yoyo->setComputeFIntFunction("YoyoPlugin.so", "F_intf");
      yoyo->setComputeJacobianFIntFunction(1, "YoyoPlugin.so", "jacobianVFIntf");
      yoyo->setComputeJacobianFIntFunction(0, "YoyoPlugin.so", "jacobianQFIntf");

      allDS.insert(yoyo);

      inter = new Interaction("impact", allDS, 0, 1, loi, relation);
      allInteractions.insert(inter);

      system = new NonSmoothDynamicalSystem(allDS, allInteractions);
      jeu = new Model(t0, T);
      jeu->setNonSmoothDynamicalSystemPtr(system);


      t = new TimeDiscretisation(h, jeu);
      s = new TimeStepping(t);
      OSI = new Moreau(yoyo, theta, s);
      osnspb = new LCP(s, mySolver);
      s->initialize();


      q = yoyo->getQPtr();
      v = yoyo->getVelocityPtr();
      p = yoyo->getPPtr();
      lambda = inter->getLambdaPtr(1);

      while (k < N)
      {
        s->computeOneStep();
        dataPlot(k, 0) =  s->getNextTime();
        dataPlot(k, 1) = (*q)(0);
        dataPlot(k, 2) = (*v)(0);
        dataPlot(k, 3) = (*p)(0);
        dataPlot(k, 4) = (*lambda)(0);
        dataPlot(k, 5) = (*q)(1);
        dataPlot(k, 6) = (*v)(1);
        dataPlot(k, 7) = L + (*q)(1) - r * (*q)(0) - (*q)(2);
        dataPlot(k, 8) = (*q)(1) - (*q)(2);
        k++;
        if ((*lambda)(0) > 0 && (-r * (*v)(0) + (*v)(1) - (*v)(2)) < 10e-14) break;
        s->nextStep();
        ++show_progress;
      }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    cout << endl << "End of computation - Number of iterations done: " << k - 1 << endl;
    cout << "Computation Time " << time.elapsed()  << endl;

    // écrire les valeur de la matrice dataPlot dans un fichier result.dat

    cout << "====> Output file writing ..." << endl;
    ioMatrix io1("result.dat", "ascii");
    io1.write(dataPlot, "noDim");
    ioMatrix io2("fichier.dat", "ascii");
    io2.write(Controle, "noDim");

    // --- Libérer de la mémoire

    delete mySolver;
    delete relation;
    delete relation0;
    delete loi;
    delete loi0;
    delete q0;
    delete v0;
    delete M;
    delete H;
  }

  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught in main.cpp" << endl;
  }

}
