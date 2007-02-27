global GLOBAL_U;
global GLOBAL_DT;
global GLOBAL_HIST;

getf 'lib/arm2/modele.sci'
getf 'lib/modules.sci'
getf 'lib/simulation.sci'
getf 'lib/plot.sci'

getf 'arm2Plot.sci'
getf 'arm2Control.sci'

// vector x=(q,qdot)
// Initial Values
//

arm2ModelInit();

GLOBAL_DT=0.1;

q0=zeros(1,arm2Ndof());

T=5.0;

simulationLoop(q0,T,arm2JointControl, arm2Model,arm2Ndof());
 
//simulationPlot2D(1,1);

//simulationPlot2D(1,2);

arm2Plot3D();
