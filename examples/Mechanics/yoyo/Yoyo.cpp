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
    int N = ceil((T - t0) / h) + 1;


    // déclarations

    SP::SiconosMatrix M(new SimpleMatrix(nDof, nDof));
    SP::SiconosVector q0(new SiconosVector(nDof));
    SP::SiconosVector v0(new SiconosVector(nDof));

    SP::NonSmoothLaw loi(new NewtonImpactNSL(e));
    SP::Relation relation(new LagrangianRheonomousR("YoyoPlugin:h1", "YoyoPlugin:G10", "YoyoPlugin:G11"));

    SP::SimpleMatrix H(new SimpleMatrix(1, nDof));
    (*H)(0, 0) = 0;
    (*H)(0, 1) = 0;
    (*H)(0, 2) = 0;
    SP::NonSmoothLaw loi0(new NewtonImpactNSL(e));
    SP::Relation relation0(new LagrangianLinearTIR(H));


    // pramètres du solveur siconos


    unsigned int outputSize = 9;
    SimpleMatrix dataPlot(N, outputSize);


    SP::SiconosVector q ;
    SP::SiconosVector v ;
    SP::SiconosVector p ;
    SP::SiconosVector lambda;


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
      }

      // création et insertion  du système dynamique représentant la yoyo dans le récipient allDS
      SP::LagrangianDS yoyo(new LagrangianDS(q0, v0, M));

      yoyo->setComputeFExtFunction("YoyoPlugin", "force_ext");
      yoyo->setComputeFIntFunction("YoyoPlugin", "F_int");
      yoyo->setComputeJacobianFIntqDotFunction("YoyoPlugin", "jacobianVFInt");
      yoyo->setComputeJacobianFIntqFunction("YoyoPlugin", "jacobianFIntq");


      ////////////////  loi d'impact et relations /////////////////////////////////

      SP::Interaction inter(new Interaction(loi0, relation0));

      /////////////////////////  MODEL //////////////////////////////////////////////////
      SP::Model jeu(new Model(t0, T));
      jeu->nonSmoothDynamicalSystem()->insertDynamicalSystem(yoyo);
      jeu->nonSmoothDynamicalSystem()->link(inter, yoyo);
      ///////////////////// SIMULATION /////////////////////////////////


      // déscrétisation du temps
      SP::TimeDiscretisation t(new TimeDiscretisation(t0, h));

      // -- OneStepIntegrators --
      SP::MoreauJeanOSI OSI(new MoreauJeanOSI(theta));

      // -- OneStepNsProblem --
      SP::OneStepNSProblem osnspb(new LCP());

      SP::TimeStepping s(new TimeStepping(t, OSI, osnspb));

      jeu->setSimulation(s);

      // --- Model initialization ---
      jeu->initialize();

      q = yoyo->q();
      v = yoyo->velocity();
      p = yoyo->p(1);
      lambda = inter->lambda(1);


      // --- sauver les valeurs dans une matrice dataPlot

      dataPlot(k, 0) = jeu->t0();
      dataPlot(k, 1) = (*q)(0);
      dataPlot(k, 2) = (*v)(0);
      dataPlot(k, 3) = (*p)(0);
      dataPlot(k, 4) = (*lambda)(0);
      dataPlot(k, 5) = (*q)(1);
      dataPlot(k, 6) = (*v)(1);
      dataPlot(k, 7) =  L + (*q)(1) - r * (*q)(0) - (*q)(2);
      dataPlot(k, 8) = (*q)(1) - (*q)(2);
      k++;


      while (s->hasNextEvent() && (*q)(0) > 0.0)
      {
        s->computeOneStep();
        // --- Get values to be plotted ---
        dataPlot(k, 0) =  s->nextTime();
        dataPlot(k, 1) = (*q)(0);
        dataPlot(k, 2) = (*v)(0);
        dataPlot(k, 3) = (*p)(0);
        dataPlot(k, 4) = (*lambda)(0); // multiplicateur de lagrange
        dataPlot(k, 5) = (*q)(1);
        dataPlot(k, 6) = (*v)(1);
        dataPlot(k, 7) =  L + (*q)(1) - r * (*q)(0) - (*q)(2); // contrainte géométrique
        dataPlot(k, 8) = (*q)(1) - (*q)(2);
        if (abs((*v)(0)) <= 0.05) cout << "valeur max de theta est : "  << (*q)(0) << endl;
        k++;
        s->nextStep();
        ++show_progress;
      }


      //////////////////////////////////////////Phase libre//////////////////////////////////////////////////////


      M->eye();
      (*M)(1, 1) = m;
      (*M)(0, 0) = I;

      t0 = dataPlot(k - 1, 0) + h;
      if (t0 + h < T)
      {
        (*q0)(0) = 0;
        (*q0)(2) = (*q)(2);
        (*q0)(1) = (*q0)(2) - L;
        (*v0)(0) = -(*v)(0);
        (*v0)(2) = (*v)(2);
        (*v0)(1) = -r * (*v0)(0) + (*v0)(2);

        yoyo.reset(new LagrangianDS(q0, v0, M));
        yoyo->setComputeFExtFunction("YoyoPlugin", "force_extf");
        yoyo->setComputeFIntFunction("YoyoPlugin", "F_intf");
        yoyo->setComputeJacobianFIntqDotFunction("YoyoPlugin", "jacobianVFIntf");
        yoyo->setComputeJacobianFIntqFunction("YoyoPlugin", "jacobianFIntqf");

        inter.reset(new Interaction(loi, relation));

        jeu.reset(new Model(t0, T));
        jeu->nonSmoothDynamicalSystem()->insertDynamicalSystem(yoyo);
        jeu->nonSmoothDynamicalSystem()->link(inter, yoyo);

        t.reset(new TimeDiscretisation(t0, h));
        OSI.reset(new MoreauJeanOSI(theta));
        osnspb.reset(new LCP());
        s.reset(new TimeStepping(t, OSI, osnspb));
	jeu->setSimulation(s);
        jeu->initialize();

        q = yoyo->q();
        v = yoyo->velocity();
        p = yoyo->p(1);
        lambda = inter->lambda(1);

        while (s->hasNextEvent())
        {
          s->computeOneStep();
          dataPlot(k, 0) =  s->nextTime();
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
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    cout << endl << "End of computation - Number of iterations done: " << k - 1 << endl;
    cout << "Computation Time " << time.elapsed()  << endl;

    // écrire les valeur de la matrice dataPlot dans un fichier result.dat

    cout << "====> Output file writing ..." << endl;
    ioMatrix::write("result.dat", "ascii", dataPlot, "noDim");
    ioMatrix::write("fichier.dat", "ascii", Controle, "noDim");

    // --- Libérer de la mémoire
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
