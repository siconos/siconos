/* Siconos-sample , Copyright INRIA 2005-2011.
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

/*!\file BouncingBallNETS.cpp
  \brief \ref EMBouncingBall - C++ input file, Time-Stepping version -
  V. Acary, O. Bonnefon.

  A Ball bouncing on the ground.
  Direct description of the model without XML input.
  Simulation with a Time-Stepping scheme.
*/
#include "SphereNEDSPlanR.hpp"
#include "SiconosKernel.hpp"
static double sBallRadius = 0.1;
#define WITH_PROJ
#define WITH_FC3D
using namespace std;
#ifdef WITH_FC3D
#define R_CLASS NewtonEulerFrom3DLocalFrameR
#else
#define R_CLASS NewtonEulerFrom1DLocalFrameR
#endif
class my_NERFC3D : public R_CLASS
{
public:
  my_NERFC3D(): R_CLASS() {};
  virtual void computeh(double t)
  {
    SP::SiconosVector y = interaction()->y(0);
    SP::NewtonEulerR ner = (boost::static_pointer_cast<NewtonEulerR>(interaction()->relation()));
    double hpc = fabs(data[q0]->getValue(0)) - sBallRadius;
    y->setValue(0, hpc);
    ner->yProj()->setValue(0, hpc);
    _Nc->setValue(0, 1);
    _Nc->setValue(1, 0);
    _Nc->setValue(2, 0);
    _Pc1->setValue(0, hpc);
    _Pc1->setValue(1, data[q0]->getValue(1));
    _Pc1->setValue(2, data[q0]->getValue(2));
    _Pc2->setValue(0, hpc);
    _Pc2->setValue(1, data[q0]->getValue(1));
    _Pc2->setValue(2, data[q0]->getValue(2));
    //printf("my_NERFC3D N, Pc\n");
    //_Nc->display();
    //_Pc1->display();
  }
  //ACCEPT_VISITORS();
};
TYPEDEF_SPTR(my_NERFC3D);





int main(int argc, char* argv[])
{
  try
  {


    // ================= Creation of the model =======================

    // User-defined main parameters
    unsigned int nDof = 3;           // degrees of freedom for the ball
    unsigned int qDim = 7;           // degrees of freedom for the ball
    unsigned int nDim = 6;           // degrees of freedom for the ball
    double t0 = 0;                   // initial computation time
    double T = 10.0;                  // final computation time
    double h = 0.005;                // time step
    double position_init = 1.0;      // initial position for lowest bead.
    double velocity_init = 2.0;      // initial velocity for lowest bead.
    double omega_initx = 0.0;
    double omega_initz = 0.0;// initial velocity for lowest bead.
    double theta = 0.5;              // theta for Moreau integrator
    double m = 1; // Ball mass
    double g = 9.81; // Gravity

    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    cout << "====> Model loading ..." << endl << endl;

    // -- Initial positions and velocities --
    SP::SimpleVector q0(new SimpleVector(qDim));
    SP::SimpleVector v0(new SimpleVector(nDim));
    SP::SimpleMatrix I(new SimpleMatrix(3, 3));
    v0->zero();
    q0->zero();
    I->eye();
    (*q0)(0) = position_init;
    /*initial quaternion equal to (1,0,0,0)*/
    (*q0)(3) = 1.0;

    (*v0)(0) = velocity_init;
    (*v0)(3) = omega_initx;
    (*v0)(5) = omega_initz;
    // -- The dynamical system --
    SP::NewtonEulerDS ball(new NewtonEulerDS(q0, v0, m, I));

    // -- Set external forces (weight) --
    SP::SimpleVector weight(new SimpleVector(nDof));
    (*weight)(0) = -m * g;
    ball->setFExtPtr(weight);

    // --------------------
    // --- Interactions ---
    // --------------------

    // -- nslaw --
    double e = 0.9;

    // Interaction ball-floor
    //

    //     vector<SP::SiconosMatrix> vecMatrix1;
    //     vecMatrix1.push_back(H);
    //     SP::BlockMatrix H_block(new BlockMatrix(vecMatrix1,1,1));

    //     SP::SiconosMatrix HT(new SimpleMatrix(1,nDim));
    //     vector<SP::SiconosMatrix> vecMatrix2;
    //     vecMatrix2.push_back(HT);
    //     SP::BlockMatrix HT_block(new BlockMatrix(vecMatrix2,1,1));

#ifdef WITH_FC3D
    int nslawsize = 3;
    SP::NonSmoothLaw nslaw0(new NewtonImpactFrictionNSL(e, e, 0.6, 3));
#else
    int nslawsize = 1;
    SP::NonSmoothLaw nslaw0(new NewtonImpactNSL(e));
#endif
    SP::SimpleMatrix H(new SimpleMatrix(nslawsize, qDim));
    H->zero();
    (*H)(0, 0) = 1.0;
#ifdef WITH_FC3D
    (*H)(1, 1) = 1.0;
    (*H)(2, 2) = 1.0;
#endif
    SP::NewtonEulerR relation0(new my_NERFC3D());
    //SP::NewtonEulerR relation0(new SphereNEDSPlanR(0.1,1.0,0.0,0.0,0.0));
    //SP::NewtonEulerR relation0(new NewtonEulerR());
    relation0->setJachq(H);
    //    relation0->setJacQH(H_block);
    //    relation0->setJacQHT(HT_block);
    cout << "main jacQH" << endl;
    relation0->jachq()->display();
    SP::Interaction inter(new Interaction(nslawsize, nslaw0, relation0));

    // -------------
    // --- Model ---
    // -------------
    SP::Model bouncingBall(new Model(t0, T));
    // add the dynamical system in the non smooth dynamical system
    bouncingBall->nonSmoothDynamicalSystem()->insertDynamicalSystem(ball);

    // link the interaction and the dynamical system
    bouncingBall->nonSmoothDynamicalSystem()->link(inter, ball);

    // ------------------
    // --- Simulation ---
    // ------------------

    // -- (1) OneStepIntegrators --
    SP::Moreau OSI(new Moreau(ball, theta));

    // -- (2) Time discretisation --
    SP::TimeDiscretisation t(new TimeDiscretisation(t0, h));

    // -- (3) one step non smooth problem
    SP::OneStepNSProblem osnspb(new GenericMechanical());
#ifdef WITH_PROJ
    SP::OneStepNSProblem osnspb_pos(new MLCPProjectOnConstraints(SICONOS_MLCP_ENUM));
#endif
    // -- (4) Simulation setup with (1) (2) (3)
#ifdef WITH_PROJ
    SP::TimeStepping s(new TimeSteppingProjectOnConstraints(t, OSI, osnspb, osnspb_pos));
#else
    SP::TimeStepping s(new TimeStepping(t, OSI, osnspb));
#endif
    s->setNewtonTolerance(1e-4);
    s->setNewtonMaxIteration(10);

    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    // --- Simulation initialization ---

    cout << "====> Initialisation ..." << endl << endl;
    bouncingBall->initialize(s);
    int N = ceil((T - t0) / h); // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 16;
    SimpleMatrix dataPlot(N + 1, outputSize);

    SP::SiconosVector q = ball->q();
    SP::SiconosVector v = ball->velocity();
    SP::SiconosVector p = ball->p(1);
    SP::SiconosVector lambda = inter->lambda(1);

    dataPlot(0, 0) = bouncingBall->t0();
    dataPlot(0, 1) = (*q)(0);
    dataPlot(0, 2) = (*v)(0);
    dataPlot(0, 3) = (*p)(0);
    dataPlot(0, 4) = (*lambda)(0);
    dataPlot(0, 5) = acos((*q)(3));
    dataPlot(0, 6) = relation0->contactForce()->norm2();
    dataPlot(0, 7) = (*q)(0);
    dataPlot(0, 8) = (*q)(1);
    dataPlot(0, 9) = (*q)(2);
    dataPlot(0, 10) = (*q)(3);
    dataPlot(0, 11) = (*q)(4);
    dataPlot(0, 12) = (*q)(5);
    dataPlot(0, 13) = (*q)(6);
    dataPlot(0, 14) = (*v)(1);
    dataPlot(0, 15) = (*v)(2);

    // --- Time loop ---
    cout << "====> Start computation ... " << endl << endl;
    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 1;
    boost::progress_display show_progress(N);

    boost::timer time;
    time.restart();
    dataPlot(k, 6) = relation0->contactForce()->norm2();
    while (s->nextTime() < T)
    {
      //      s->computeOneStep();
      s->advanceToEvent();
      // --- Get values to be plotted ---
      dataPlot(k, 0) =  s->nextTime();
      dataPlot(k, 1) = (*q)(0);
      dataPlot(k, 2) = (*v)(0);
      dataPlot(k, 3) = (*p)(0);
      dataPlot(k, 4) = (*lambda)(0);
      dataPlot(k, 5) = (*q)(3);
      dataPlot(k, 6) = relation0->contactForce()->norm2();
      dataPlot(k, 7) = (*q)(0);
      dataPlot(k, 8) = (*q)(1);
      dataPlot(k, 9) = (*q)(2);
      dataPlot(k, 10) = (*q)(3);
      dataPlot(k, 11) = (*q)(4);
      dataPlot(k, 12) = (*q)(5);
      dataPlot(k, 13) = (*q)(6);
      dataPlot(k, 14) = (*v)(1);
      dataPlot(k, 15) = (*v)(2);
      s->nextStep();
      ++show_progress;
      k++;
    }
    cout << endl << "End of computation - Number of iterations done: " << k - 1 << endl;
    cout << "Computation Time " << time.elapsed()  << endl;

    // --- Output files ---
    cout << "====> Output file writing ..." << endl;
    ioMatrix io("result.dat", "ascii");
    dataPlot.resize(k, outputSize);
    io.write(dataPlot, "noDim");

    // Comparison with a reference file
    SimpleMatrix dataPlotRef(dataPlot);
    dataPlotRef.zero();
#ifdef WITH_PROJ
    ioMatrix ref("resultNETS-WITHPROJ.ref", "ascii");
#else
    ioMatrix ref("resultNETS.ref", "ascii");
#endif
    ref.read(dataPlotRef);
    //std::cout << (dataPlot-dataPlotRef).normInf() <<std::endl;
    if ((dataPlot - dataPlotRef).normInf() > 1e-10)
    {
      std::cout << "Warning. The results is rather different from the reference file. err = " << (dataPlot - dataPlotRef).normInf() << std::endl;
      return 1;
    }
  }

  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught in BouncingBallTS.cpp" << endl;
  }

}
