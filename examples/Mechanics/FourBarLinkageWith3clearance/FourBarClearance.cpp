#include "SiconosKernel.hpp"
#include <stdlib.h>
using namespace std;

#include <boost/progress.hpp>
#define PI 3.14159265
// parameters according to Table 1
// geometrical characteristics
double l1 = 1.0;
double l2 = 4.0;
double l3 = 2.5;
double l0 = 3.0;
double r1 = 0.0;
double r3 = 0.0;
double r5 = 0.0;
double Kp = 0.0;
double lmd = 0.0;
// force elements
double gravity = 9.81;
double m1 = 1.0;
double m2 = 1.0;
double m3 = 1.0;
double I1 = m1*l1*l1/3.0;
double J2 = m2*l2*l2/12.0;
double I3 = m3*l3*l3/3.0;
//double PI = 22.0/7.0;


int main(int argc, char* argv[])
{
  try
  {

    // ================= Creation of the model =======================

    // User-defined main parameters
    unsigned int nDof = 7;           // degrees of freedom for robot arm
    double t0 = 0;                   // initial computation time
    double T = 1.0;                   // final computation time
    double h = 1e-4;                // time step
    //double criterion = 1e-6;
    //unsigned int maxIter = 2000;
    char  filename[50] = "simu_";
    double eN = 0.0;
    //eN1 = 0.1;
    double eT = 0.0;
    double mu = 0.1;

    cout << "argc :" << argc << endl;
    if (argc < 2)
    {
      cout << "Using default arguments" << endl;

      ::r1 = 0.055;
      ::r3 = 0.055;
      ::r5 = 0.055;
      ::Kp = 200;
      ::lmd = 10;
      int sizeofargv1 = strlen("0.055");
      int sizeofargv2 = strlen("0.055");
      int sizeofargv3 = strlen("0.055");
      int sizeofargv4 = strlen("200");
      int sizeofargv5 = strlen("10");
      strncpy(&filename[5],"0.055",sizeofargv1);
      strncpy(&filename[11],"0.055",sizeofargv2);
      strncpy(&filename[17],"0.055",sizeofargv3);
      strncpy(&filename[23],"200",sizeofargv4);
      strncpy(&filename[26],"10",sizeofargv5);

    }
    else if (argc==6)
    {

      ::r1 = atof(argv[1]);
      ::r3 = atof(argv[2]);
      ::r5 = atof(argv[3]);
      ::Kp = atof(argv[4]);
      ::lmd = atof(argv[5]);
      int sizeofargv1 = strlen(argv[1]);
      int sizeofargv2 = strlen(argv[2]);
      int sizeofargv3 = strlen(argv[3]);
      int sizeofargv4 = strlen(argv[4]);
      int sizeofargv5 = strlen(argv[5]);
      strncpy(&filename[5],argv[1],sizeofargv1);
      strncpy(&filename[11],argv[2],sizeofargv2);
      strncpy(&filename[17],argv[3],sizeofargv3);
      strncpy(&filename[23],argv[4],sizeofargv4);
      strncpy(&filename[26],argv[5],sizeofargv5);
    }
    else
    {
      cout << "Wrong number of arguments" << endl;
      return 1;
    }

    // -> mind to set the initial conditions below.

    // -------------------------
    // --- Dynamical systems ---
    // -------------------------


    // --- DS: slidercrank ---

    // Initial position (angles in radian)
    SP::SiconosVector q0(new SiconosVector(nDof));
    SP::SiconosVector v0(new SiconosVector(nDof));
    q0->zero();
    v0->zero();
    (*q0)(0) = 1.570823772407980;//1.5708;
    (*q0)(1) = 0.3532842020624460;//0.3533;
    (*q0)(2) = 1.264872058968431;//1.2649;
    (*q0)(3) = 1.876454585097650;//1.87647;
    (*q0)(4) = 1.691962091335582;//1.69199;
    (*q0)(5) = 0.3764686197082958;//0.3764+3.5e-5;
    (*q0)(6) = 1.191962183453451;//1.19197;
    (*v0)(0) = 0.0;
    (*v0)(1) = 0.0;
    (*v0)(2) = 0.0;

    SP::SiconosVector z(new SiconosVector(3));

    SP::LagrangianDS fourbar(new LagrangianDS(q0, v0, "FourBarClearancePlugin:mass"));

    // external plug-in
    fourbar->setComputeFGyrFunction("FourBarClearancePlugin.so", "FGyr");
    fourbar->setComputeJacobianFGyrqFunction("FourBarClearancePlugin.so", "jacobianFGyrq");
    fourbar->setComputeJacobianFGyrqDotFunction("FourBarClearancePlugin.so", "jacobianFGyrVelocity");
    fourbar->setComputeFIntFunction("FourBarClearancePlugin.so", "FInt");
    fourbar->setComputeJacobianFIntqDotFunction("FourBarClearancePlugin.so", "jacobianFIntqDot");
    fourbar->setComputeJacobianFIntqFunction("FourBarClearancePlugin.so", "jacobianFIntq");


    // -------------------
    // --- Interactions---
    // -------------------


    //InteractionsSet allInteractions;

    // -- relations --

    SP::NonSmoothLaw nslaw(new NewtonImpactFrictionNSL(eN, eT, mu, 2)); //EqualityConditionNSL NewtonImpactNSL
    SP::Relation relation(new LagrangianScleronomousR("FourBarClearancePlugin:g1", "FourBarClearancePlugin:W1"));
    SP::Interaction inter(new Interaction(nslaw, relation));

    SP::NonSmoothLaw nslaw1(new NewtonImpactFrictionNSL(eN, eT, mu, 2)); //EqualityConditionNSL NewtonImpactNSL
    SP::Relation relation1(new LagrangianScleronomousR("FourBarClearancePlugin:g2", "FourBarClearancePlugin:W2"));
    SP::Interaction inter1(new Interaction(nslaw1, relation1));

    SP::NonSmoothLaw nslaw2(new NewtonImpactFrictionNSL(eN, eT, mu, 2)); //EqualityConditionNSL NewtonImpactNSL
    SP::Relation relation2(new LagrangianScleronomousR("FourBarClearancePlugin:g3", "FourBarClearancePlugin:W3"));
    SP::Interaction inter2(new Interaction(nslaw2, relation2));


    //SP::NonSmoothLaw nslaw2(new EqualityConditionNSL(e1)); //EqualityConditionNSL NewtonImpactNSL
    //SP::Relation relation2(new LagrangianScleronomousR("FourBarClearancePlugin:g3", "FourBarClearancePlugin:W3"));
    //SP::Interaction inter2(new Interaction(nslaw2, relation2));

    //SP::NonSmoothLaw nslaw3(new EqualityConditionNSL(e1)); //EqualityConditionNSL NewtonImpactNSL
    //SP::Relation relation3(new LagrangianScleronomousR("FourBarClearancePlugin:g4", "FourBarClearancePlugin:W4"));
    //SP::Interaction inter3(new Interaction(nslaw3, relation3));

    //allInteractions.insert(inter);
    //allInteractions.insert(inter1);
    //allInteractions.insert(inter2);
    //allInteractions.insert(inter3);

    // -------------
    // --- Model ---
    // -------------


    SP::Model FourBarIdeal(new Model(t0, T));
    FourBarIdeal->nonSmoothDynamicalSystem()->insertDynamicalSystem(fourbar);
    FourBarIdeal->nonSmoothDynamicalSystem()->link(inter, fourbar);
    FourBarIdeal->nonSmoothDynamicalSystem()->link(inter1, fourbar);
    FourBarIdeal->nonSmoothDynamicalSystem()->link(inter2, fourbar);
    //FourBarIdeal->nonSmoothDynamicalSystem()->link(inter3, fourbar);
    // ----------------
    // --- Simulation ---
    // ----------------
    fourbar-> setzPtr(z);

    // -- Time discretisation --
    SP::TimeDiscretisation t(new TimeDiscretisation(t0, h));

    SP::MoreauJeanCombinedProjectionOSI OSI(new MoreauJeanCombinedProjectionOSI(0.5));
    // -- set the integrator for the four bar linkage --

    SP::OneStepNSProblem impact(new FrictionContact(2,SICONOS_FRICTION_2D_PGS)); /*,SICONOS_FRICTION_2D_ENUM
                                                                                   SICONOS_FRICTION_2D_LEMKE notworking //
                                                                                   SICONOS_FRICTION_2D_PGS=>working 0.75PI */
    impact->numericsSolverOptions()->dparam[0] = 1e-6;
    impact->numericsSolverOptions()->iparam[0] = 2000;
    impact->numericsSolverOptions()->iparam[2] = 1; // random
    SP::OneStepNSProblem position(new MLCPProjectOnConstraints(SICONOS_MLCP_ENUM));

    SP::TimeSteppingCombinedProjection s(new TimeSteppingCombinedProjection(t, OSI, impact, position,2));
    s->setProjectionMaxIteration(2000);
    s->setConstraintTolUnilateral(1e-6);
    s->setConstraintTol(1e-6);

    cout << "=== End of model loading === " << endl;
/////////////////////////////////////////////////////////////////////////
    // -- Time discretisation --
    //SP::TimeDiscretisation t(new TimeDiscretisation(t0, h));

    //SP::TimeStepping s(new TimeStepping(t));


    // -- OneStepIntegrators --//


   // double theta = 0.500001;

    //SP::Moreau OSI(new Moreau(fourbar, theta));
    //s->insertIntegrator(OSI);

    // -- OneStepNsProblem --//

    //SP::OneStepNSProblem osnspb(new FrictionContact(2));

    //s->insertNonSmoothProblem(osnspb);

    //cout << "=== End of model loading === " << endl;
    FourBarIdeal->setSimulation(s);

////////////////////////////////////////////////////////////////////////

    // =========================== End of model definition ===========================  dataPlot(k,7) = (*inter->y(0))(0);


    // ================================= Computation =================================

    // --- Simulation initialization ---
    FourBarIdeal->initialize();
    cout << "End of simulation initialisation" << endl;

    int k = 1;
    int N = ceil((T - t0) / h) + 1;
    cout << "Number of time step   " << N << endl;
    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 42;
    SimpleMatrix dataPlot((N/500)+1 , outputSize);
    SimpleMatrix beam1Plot(2,3*((N/500)+1));
    SimpleMatrix beam2Plot(2,3*((N/500)+1));
    SimpleMatrix beam3Plot(2,3*((N/500)+1));
    SimpleMatrix beam4Plot(2,3*((N/500)+1));
    SimpleMatrix beam5Plot(1,2*((N/500)+1));
    SimpleMatrix beam6Plot(1,1*((N/500)+1));
    SimpleMatrix beam7Plot(1,2*((N/500)+1));
    SimpleMatrix beam8Plot(1,2*((N/500)+1));
    SimpleMatrix beam9Plot(1,2*((N/500)+1));
    SimpleMatrix beam10Plot(1,2*((N/500)+1));
    SimpleMatrix beam11Plot(1,2*((N/500)+1));
    SimpleMatrix beam12Plot(1,2*((N/500)+1));
    SimpleMatrix beam13Plot(1,2*((N/500)+1));
	//cout << "size here " << (N/20) + 1 << endl;
    // For the initial time step:
    // time
    SP::SiconosVector q = fourbar->q();
    SP::SiconosVector v = fourbar->velocity();
    dataPlot(0, 0) = FourBarIdeal->t0();
    dataPlot(0, 1) = (*q)(0); // crank revolution
    dataPlot(0, 2) = (*q)(1);
    dataPlot(0, 3) = (*q)(2);
    dataPlot(0, 4) = (*q)(3);
    dataPlot(0, 5) = (*q)(4);
    dataPlot(0, 6) = (*q)(5);
    dataPlot(0, 7) = (*q)(6);
    dataPlot(0, 8) = (*v)(0);
    dataPlot(0, 9) = (*v)(1);
    dataPlot(0, 10) = (*v)(2);
    dataPlot(0, 11) = (*inter1->y(1))(0);
    dataPlot(0, 12) = (*inter1->y(1))(1);
    dataPlot(0, 13) = (*q)(3)-l1 * cos((*q)(0)) - 0.5*l2 * cos((*q)(1));
    dataPlot(0, 14) = (*q)(4) -l1 * sin((*q)(0)) - 0.5*l2 * sin((*q)(1));
    dataPlot(0, 15) = l0 + (*q)(5) + 0.5*l3*cos((*q)(2)) - 0.5*l2 * cos((*q)(1))-(*q)(3);
    dataPlot(0, 16) = (*q)(6) + 0.5*l3 * sin((*q)(2)) - 0.5*l2 * sin((*q)(1))-(*q)(4);
    dataPlot(0, 17) = (*q)(5)-0.5*l3 * cos((*q)(2));
    dataPlot(0, 18) = (*q)(6) -0.5*l3 * sin((*q)(2));
    dataPlot(0, 19) = (*inter->y(1))(0);
    dataPlot(0, 20) = (*inter->y(1))(1);
    dataPlot(0, 21) = (*inter->lambda(1))(0) ; // lambda1
    dataPlot(0, 22) = (*inter1->lambda(1))(0) ; // lambda2
    dataPlot(0, 23) = (*inter2->y(1))(0);
    dataPlot(0, 24) = (*inter2->y(1))(1);
    dataPlot(0, 25) = (*inter2->lambda(1))(0) ; // lambda2
    dataPlot(0, 26) = (*inter2->lambda(1))(1) ; // lambda2
    dataPlot(0, 27) = (*inter1->lambda(1))(1) ; // lambda2
    dataPlot(0, 28) = (*inter->lambda(1))(1) ; // lambda2
    dataPlot(0, 29) = sqrt(pow(((*q)(3) - 0.5 * l2 * cos((*q)(1)) - l1 * cos((*q)(0))),2) + pow(((*q)(4) - 0.5 * l2 * sin((*q)(1)) - l1 * sin((*q)(0))),2));
    dataPlot(0, 30) = sqrt(pow((l0 + (*q)(5) + 0.5*l3 * cos((*q)(2)) - 0.5*l2 * cos((*q)(1))-(*q)(3)),2) + pow((-(*q)(4) + (*q)(6)+0.5*l3 * sin((*q)(2)) - 0.5*l2 * sin((*q)(1))),2));
    dataPlot(0, 31) = sqrt((pow(((*q)(5) - 0.5 * l3 * cos((*q)(2))),2) + pow(((*q)(6) - 0.5 * l3 * sin((*q)(2))),2)));
    dataPlot(0, 32) = (*inter->y(0))(0);
    dataPlot(0, 33) = (*inter1->y(0))(0);
    dataPlot(0, 34) = (*inter2->y(0))(0);

    dataPlot(0, 35) = (*fourbar->fInt())(0) - (0.5 * m1) * gravity * l1 * cos((*q)(0));
    dataPlot(0, 36) =0.5*((*v)(0)-6.0*0.75*PI*cos(0.75*PI*h))*((*v)(0)-6.0*0.75*PI*cos(0.75*PI*h)) + 0.5*3000.0*((*q)(0)-6.0*sin(0.75*PI*h))*((*q)(0)-6.0*sin(0.75*PI*h)) + 0.5*10.0*((*v)(0)-6.0*0.75*PI*cos(0.75*PI*h))*((*q)(0)-6.0*sin(0.75*PI*h));
    dataPlot(0, 37) = 6.0*sin(0.75*PI*h);
    dataPlot(0, 38) = 6.0*0.75*PI*cos(0.75*PI*h);
    dataPlot(0, 39) = (*z)(2);
    dataPlot(0, 40) =0.5*(((*v)(0)-6.0*0.75*PI*cos(0.75*PI*h))+ 10.0*((*q)(0)-6.0*sin(0.75*PI*h)))*(((*v)(0)-6.0*0.75*PI*cos(0.75*PI*h))+ 10.0*((*q)(0)-6.0*sin(0.75*PI*h)));
    dataPlot(0, 41) =0.5*((*inter->y(1))(0)*(*inter->y(1))(0));
    boost::timer time;
    time.restart();

    // --- Time loop ---
    cout << "Start computation ... " << endl;

    boost::progress_display show_progress(N);
	double tt = 0;
	int kk = 1;
    while (s->hasNextEvent())
    {
      //k++;
      // ++show_progress;
      // if (!(div(k,1000).rem))  cout <<"Step number "<< k << "\n";
      s->advanceToEvent();
      // Solve problem
      //s->newtonSolve(criterion, maxIter); //2000000
      //std::cout << "jachq:" <<std::endl;
      //std11::static_pointer_cast<LagrangianScleronomousR>(relation)->jachq()->display();
      //std::cout << "=================================" << std::endl;
      //std::cout << "jachq:" <<std::endl;
      //std11::static_pointer_cast<LagrangianScleronomousR>(relation1)->jachq()->display();
      //std::cout << "*********************************" << std::endl;
      // Data Output

      tt = s->nextTime();
      if(((k%500)==0) && (k != 0)){
        // cout << "size here " << kk << endl;
      dataPlot(kk, 0) = tt;
      dataPlot(kk, 1) = (*q)(0); // crank revolution
      dataPlot(kk, 2) = (*q)(1);
      dataPlot(kk, 3) = (*q)(2);
      dataPlot(kk, 4) = (*q)(3);
      dataPlot(kk, 5) = (*q)(4);
      dataPlot(kk, 6) = (*q)(5);
      dataPlot(kk, 7) = (*q)(6);
      dataPlot(kk, 8) = (*v)(0);
      dataPlot(kk, 9) = (*v)(1);
      dataPlot(kk, 10) = (*v)(2);
      dataPlot(kk, 11) = (*inter1->y(1))(0);
      dataPlot(kk, 12) = (*inter1->y(1))(1);
      dataPlot(kk, 13) = (*q)(3)-l1 * cos((*q)(0)) - 0.5*l2 * cos((*q)(1));
      dataPlot(kk, 14) = (*q)(4) -l1 * sin((*q)(0)) - 0.5*l2 * sin((*q)(1));
      dataPlot(kk, 15) = l0 + (*q)(5) + 0.5*l3*cos((*q)(2)) - 0.5*l2 * cos((*q)(1))-(*q)(3);
      dataPlot(kk, 16) = (*q)(6) + 0.5*l3 * sin((*q)(2)) - 0.5*l2 * sin((*q)(1))-(*q)(4);
      dataPlot(kk, 17) = (*q)(5)-0.5*l3 * cos((*q)(2));
      dataPlot(kk, 18) = (*q)(6) -0.5*l3 * sin((*q)(2));
      dataPlot(kk, 19) = (*inter->y(1))(0);
      dataPlot(kk, 20) = (*inter->y(1))(1);
      dataPlot(kk, 21) = (*inter->lambda(1))(0) ; // lambda1
      dataPlot(kk, 22) = (*inter1->lambda(1))(0) ; // lambda2
      dataPlot(kk, 23) = (*inter2->y(1))(0);
      dataPlot(kk, 24) = (*inter2->y(1))(1);
      dataPlot(kk, 25) = (*inter2->lambda(1))(0) ; // lambda2
      dataPlot(kk, 26) = (*inter2->lambda(1))(1) ; // lambda2
      dataPlot(kk, 27) = (*inter1->lambda(1))(1) ; // lambda2
      dataPlot(kk, 28) = (*inter->lambda(1))(1) ; // lambda2
      dataPlot(kk, 29) = sqrt(pow(((*q)(3) - 0.5 * l2 * cos((*q)(1)) - l1 * cos((*q)(0))),2) + pow(((*q)(4) - 0.5 * l2 * sin((*q)(1)) - l1 * sin((*q)(0))),2));
      dataPlot(kk, 30) = sqrt(pow((l0 + (*q)(5) + 0.5*l3 * cos((*q)(2)) - 0.5*l2 * cos((*q)(1))-(*q)(3)),2) + pow((-(*q)(4) + (*q)(6)+0.5*l3 * sin((*q)(2)) - 0.5*l2 * sin((*q)(1))),2));
      dataPlot(kk, 31) = sqrt((pow(((*q)(5) - 0.5 * l3 * cos((*q)(2))),2) + pow(((*q)(6) - 0.5 * l3 * sin((*q)(2))),2)));
      dataPlot(kk, 32) = (*inter->y(0))(0);
      dataPlot(kk, 33) = (*inter1->y(0))(0);
      dataPlot(kk, 34) = (*inter2->y(0))(0);
      dataPlot(kk, 35) = (*fourbar->fInt())(0) - (0.5 * m1) * gravity * l1 * cos((*q)(0));
      dataPlot(kk, 36) = 0.5*((*v)(0)-6.0*0.75*PI*cos(0.75*PI*tt))*((*v)(0)-6.0*0.75*PI*cos(0.75*PI*tt)) + 0.5*3000.0*((*q)(0)-6.0*sin(0.75*PI*tt))*((*q)(0)-6.0*sin(0.75*PI*tt)) + 0.5*10.0*((*v)(0)-6.0*0.75*PI*cos(0.75*PI*tt))*((*q)(0)-6.0*sin(0.75*PI*tt));
      dataPlot(kk, 37) = 6.0*sin(0.75*PI*tt);
      dataPlot(kk, 38) = 6.0*0.75*PI*cos(0.75*PI*tt);
      dataPlot(kk, 39) = (*z)(2);
      dataPlot(kk, 40) = 0.5*(((*v)(0)-6.0*0.75*PI*cos(0.75*PI*tt))+ 10.0*((*q)(0)-6.0*sin(0.75*PI*tt)))*(((*v)(0)-6.0*0.75*PI*cos(0.75*PI*tt))+ 10.0*((*q)(0)-6.0*sin(0.75*PI*tt)));
      dataPlot(kk, 41) = 0.5*((*inter->y(1))(0)*(*inter->y(1))(0));

      beam1Plot(0,3*kk) = 0.0;
      beam1Plot(0,3*kk+1) = 0.0;
      beam1Plot(0,3*kk+2) = 0.0;
      beam1Plot(1,3*kk) = l1*cos((*q)(0));
      beam1Plot(1,3*kk+1) =l1*sin((*q)(0));
      beam1Plot(1,3*kk+2) = 0.0;

      beam2Plot(0,3*kk) = ((*q)(3) - 0.5*l2*cos((*q)(1))) + (-(*q)(3) + 0.5 * l2 * cos((*q)(1)) + l1 * cos((*q)(0)));
      beam2Plot(0,3*kk+1) = ((*q)(4) - 0.5*l2*sin((*q)(1))) + (-(*q)(4) + 0.5 * l2 * sin((*q)(1)) + l1 * sin((*q)(0)));
      beam2Plot(0,3*kk+2) = 0.0;
      beam2Plot(1,3*kk) = (*q)(3) + 0.5*l2 * cos((*q)(1));
      beam2Plot(1,3*kk+1) =(*q)(4) + 0.5*l2 * sin((*q)(1));
      beam2Plot(1,3*kk+2) = 0.0;

      beam3Plot(0,3*kk) = l0 + l3*cos((*q)(2));
      beam3Plot(0,3*kk+1) =l3*sin((*q)(2));
      beam3Plot(0,3*kk+2) = 0.0;
      beam3Plot(1,3*kk) = l0;
      beam3Plot(1,3*kk+1) =0.0;
      beam3Plot(1,3*kk+2) = 0.0;

      beam4Plot(0,3*kk) = 0.0;
      beam4Plot(0,3*kk+1) =0.0;
      beam4Plot(0,3*kk+2) = 0.0;
      beam4Plot(1,3*kk) = l0;
      beam4Plot(1,3*kk+1) =0.0;
      beam4Plot(1,3*kk+2) = 0.0;

      beam5Plot(0,2*kk) = l1*cos((*q)(0));
      beam5Plot(0,2*kk+1) = l1*sin((*q)(0));
     // beam5Plot(0,3*kk+2) = 0.0;
      //beam5Plot(1,3*kk) = 0.0;
      //beam5Plot(1,3*kk+1) =0.0;
      //beam5Plot(1,3*kk+2) = 0.0;

      beam6Plot(0,1*kk) = tt;

      beam7Plot(0,2*kk) = ((*q)(3) - 0.5*l2*cos((*q)(1))) - (-(*q)(3) + 0.5 * l2 * cos((*q)(1)) + l1 * cos((*q)(0)));
      beam7Plot(0,2*kk+1) = ((*q)(4) - 0.5*l2*sin((*q)(1))) - (-(*q)(4) + 0.5 * l2 * sin((*q)(1)) + l1 * sin((*q)(0)));

      beam8Plot(0,2*kk) = -2.0- 10*(-(*q)(3) + 0.5 * l2 * cos((*q)(1)) + l1 * cos((*q)(0)));
      beam8Plot(0,2*kk+1) = 2.5 - 10*(-(*q)(4) + 0.5 * l2 * sin((*q)(1)) + l1 * sin((*q)(0)));

      beam9Plot(0,2*kk) = -2.0;
      beam9Plot(0,2*kk+1) = 2.5;

      beam10Plot(0,2*kk) = (l0 + l3*cos((*q)(2)));
      beam10Plot(0,2*kk+1) = (l3*sin((*q)(2)));

      beam11Plot(0,2*kk) = ((*q)(3) + 0.5*l2*cos((*q)(1))) - (-l0-l3 * cos((*q)(2)) + 0.5*l2 * cos((*q)(1))+(*q)(3));
      beam11Plot(0,2*kk+1) = ((*q)(4) + 0.5*l2*sin((*q)(1))) - ((*q)(4) - l3 * sin((*q)(2)) + 0.5*l2 * sin((*q)(1)));

      beam12Plot(0,2*kk) = -2.0+ 10*(-l0-l3 * cos((*q)(2)) + 0.5*l2 * cos((*q)(1))+(*q)(3));
      beam12Plot(0,2*kk+1) = -1.5 + 10*((*q)(4) - l3 * sin((*q)(2)) + 0.5*l2 * sin((*q)(1)));

      beam13Plot(0,2*kk) = -2.0;
      beam13Plot(0,2*kk+1) = -1.5;

       kk++;
       }
      //dataPlot(k, 17) = s->getNewtonNbSteps();
      //dataPlot(k, 18) = s->nbProjectionIteration();
      //dataPlot(k, 19) = s->maxViolationUnilateral();
      //dataPlot(k, 20) = s->nbIndexSetsIteration();
      //dataPlot(k, 21) = s->cumulatedNewtonNbSteps();
      //dataPlot(k, 22) = s->nbCumulatedProjectionIteration();
      s->processEvents();
      ++show_progress;
      k++;
    }

    cout << "\nEnd of computation - Number of iterations done: " << k << endl;
    cout << "Computation Time " << time.elapsed()  << endl;

    // --- Output files ---
    dataPlot.resize(kk, outputSize);
    ioMatrix::write(filename, "ascii", dataPlot, "noDim");
    ioMatrix::write("Link1.dat", "ascii", beam1Plot, "noDim");
    ioMatrix::write("Link2.dat", "ascii", beam2Plot, "noDim");
    ioMatrix::write("Link3.dat", "ascii", beam3Plot, "noDim");
    ioMatrix::write("Link4.dat", "ascii", beam4Plot, "noDim");
    ioMatrix::write("ex_ey.dat", "ascii", beam5Plot, "noDim");
    ioMatrix::write("time.dat", "ascii", beam6Plot, "noDim");
    ioMatrix::write("ex1_ey1.dat", "ascii", beam7Plot, "noDim");
    ioMatrix::write("ex1x_ey1x.dat", "ascii", beam8Plot, "noDim");
    ioMatrix::write("ex1xy_ey1xy.dat", "ascii", beam9Plot, "noDim");
    ioMatrix::write("ex1xy1_ey1xy1.dat", "ascii", beam10Plot, "noDim");
    ioMatrix::write("ex1xy2_ey1xy2.dat", "ascii", beam11Plot, "noDim");
    ioMatrix::write("ex1xy3_ey1xy3.dat", "ascii", beam12Plot, "noDim");
    ioMatrix::write("ex1xy4_ey1xy4.dat", "ascii", beam13Plot, "noDim");
  }


  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught in FourBarD1MinusLinear.cpp" << endl;
  }
}
