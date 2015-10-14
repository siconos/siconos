// Siconos-sample version 3.0.0, Copyright INRIA 2005-2008.
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.	
// Siconos is a free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
// Siconos is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Siconos; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
//
// Contact: Vincent ACARY vincent.acary@inrialpes.fr 
//	
//
//


exec '../../Front-End/scilab/loaderSiconos.sce';

function Simul()

// Load Siconos

[t,Q]=Robot(%F);

endfunction


function [Time,Q] = Robot(bRealtime)

// Bug dgetri ?
errcatch(998,'continue');

// --- DS: robot arm ---
    
    // The dof are angles between ground and arm and between differents parts of the arm. (See corresponding .pdf for more details)

 // Initial position (angles in radian)
q0= [0.05,0.05,0];
v0= [0,0,0];
nDof=3;
t0 = 0;                   // initial computation time
T = 10;                   // final computation time 
h = 0.002;                // time step
criterion = 0.001;
maxIter = 50;
e = 0.9;                  // nslaw 
e2 = 0.0; 


nId=sicLagrangianDS(nDof,q0,v0);

// external plug-in
sicSetMass(nId,"RobotPlugin","mass");
sicSetNNL(nId,"RobotPlugin", "NNL");
sicSetJacQNNL(nId,"RobotPlugin", "jacobianQNNL");
sicSetJacVelNNL(nId,"RobotPlugin", "jacobianVNNL");
sicSetFInt(nId,"RobotPlugin", "FInt");
sicSetJacQFInt(nId,"RobotPlugin", "jacobianQFInt");
sicSetJacVelFInt(nId,"RobotPlugin", "jacobianQFInt");
sicSetFExt(nId,"RobotPlugin", "FExt");

// -------------------
// --- Interactions---
// -------------------

// Two interactions:
//  - one with Lagrangian non linear relation to define contact with ground
//  - the other to define angles limitations (articular stops), with lagrangian linear relation
//  Both with newton impact nslaw. 

DS=[0];
idInter=sicInteraction("floor-arm",1,DS,2);
sicLagrangianR(idInter,"scleronomic", "RobotPlugin:h2","RobotPlugin:G2");
sicNewtonImpactLawNSL(idInter,e);

H=[[-1,1,-1,1,-1,1],[1.7,1.7,0.3,0.3,3.14,3.14]];
b=[1.7,1.7,0.3,0.3,3.14,3.14];
idInter=sicInteraction("floor-arm2",1,DS,6);
sicLagrangianLinearR(idInter,H,b);
sicNewtonImpactLawNSL(idInter,e2);


// Construct NSDS 
sicNonSmoothDynamicalSystem(0);
// Construct Model 
sicModel(t0,T);

// Simulation Model 
sicSimulationTimeStepping(h);
Theta=[0.5];
sicOneStepIntegratorMoreau(Theta);
sicOneStepNSProblemLCP("NLGS",101,0.001);


// Simulation 

sicInitSimulation();

k = sicTimeGetK();
N = sicTimeGetN();

winId=waitbar('Siconos Computation');


dixpc=0;

while k < N do 

  // transfer of state i+1 into state i and time incrementation
  sicSTNextStep();

  // get current time step
  k = sicTimeGetK();
  
  // solve ... 
  sicSTnewtonSolve(criterion,maxIter);

  // update
  sicSTupdateState();
	
 // --- Get values to be plotted ---

 dixpc=dixpc+1/N;
 
 if (dixpc>=0.01) then
   waitbar(k/N,winId);
   dixpc=0;
 end
 
 // get Time and axis values
 Time(k+1)=k*sicTimeGetH();
 Q(1,k+1)=sicModelgetQ(0,0); 
 Q(2,k+1)=sicModelgetQ(0,1); 
 Q(3,k+1)=sicModelgetQ(0,2); 
   
end   //while

  winclose(winId);

endfunction


