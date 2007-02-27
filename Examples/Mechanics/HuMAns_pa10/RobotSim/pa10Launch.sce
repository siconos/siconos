global GLOBAL_U;
global GLOBAL_DT;
global GLOBAL_HIST;


getf 'lib/modules.sci'
getf 'lib/simulation.sci'
getf 'lib/plot.sci'


getf 'lib/pa10/modele.sci'
getf 'lib/pa10/pa10Jac.sci'
getf 'pa10Control.sci'
getf 'pa10Plot.sci'


// vector x=(q,qdot)
// Initial Values
//

pa10ModelInit();

GLOBAL_DT=0.001;

q0=[0 0.4 0 0 0 0 0];

T=10;

simulationLoop(q0,T,pa10Free, pa10Model,pa10Spy,%F,pa10Ndof());
//simulationLoop(q0,T,pa10JointControl, pa10Model,pa10Spy,%F,pa10Ndof());
//simulationLoop(q0,T,pa10CartesianControl, pa10Model,pa10Spy,%F,pa10Ndof());
 
//pa10Plot2D(2);

pa10Plot3D();

