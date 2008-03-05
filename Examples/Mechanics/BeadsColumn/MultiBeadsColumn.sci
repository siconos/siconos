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
getf("../../Front-End/scilab/visuballs.sci");

function Simul()

// Load Siconos

[t,qx,qy]=MultiBeadsColumn(5,%F);

VisuBalls(qx,qy,3,3,0.1);

endfunction


function [Time,Qx,Qy] = MultiBeadsColumn(nBalls,realtime)

// Bug dgetri ?
errcatch(998,'continue');

// Dynamical System (Beads) parameters
Mass=[1,0,0,0,1,0,0,0,1];
C=[0,0,0,0,0,0,0,0,0];
K=[0,0,0,0,1,0,0,0,0];
// Dynamical System initial conditions 
q0= [1,0,0];
v0= [1,0,0];
nDof=3;
inc_pos=0.5;
inc_vel=0;
// Simulation parameters 
Theta=[0,0,0];
// Interaction parameters
DS=[0,0];
H=[0,0,0,0,0,0];
b=[0];

// Creation of Dynamical Systems (Beads) 
for i = 1:nBalls
  sicLagrangianLinearTIDS(nDof,q0,v0,Mass,K,C,"BeadsPlugin.so","beadsFExt");  
  q0(1)=q0(1)+inc_pos;
  v0(1)=v0(1)+inc_vel;
  Theta(i)=0.500001;
end //for 

// Creation of Interactions 
// Bead and floor
DS=[0,0];
H(1)=1;
b(1)=-0.1;
idInter=sicInteraction("floor",1,DS,1);
sicLagrangianLinearR(idInter,H,b);
sicNewtonImpactLawNSL(idInter,0.9);
// Last Bead and ceiling 
DS=[nBalls-1,0]; 
H(1)=-1;
b(1)=1.5;
idInter=sicInteraction("ceiling",1,DS,1);
sicLagrangianLinearR(idInter,H,b);
sicNewtonImpactLawNSL(idInter,0.9);
// Between beads 
H(1)=-1; H(4)=1;
b(1)=-0;1;
for i = 1:nBalls-1
  DS=[i-1,i];
  nameInter=sprintf("inter%d\0",i);
  idInter=sicInteraction(nameInter,2,DS,1);
  sicLagrangianLinearR(idInter,H,b);
  sicNewtonImpactLawNSL(idInter,0.9);
end //for

// Construct NSDS 
sicNonSmoothDynamicalSystem(0);
// Construct Model 
sicModel(0.0,3.0);

// Simulation Model 
sicSimulationTimeStepping(0.001);
sicOneStepIntegratorMoreau(Theta);
sicOneStepNSProblemLCP("NSQP",101,0.0001);


// Simulation 

sicInitSimulation();

k = sicTimeGetK();
N = sicTimeGetN();

winId=waitbar('Siconos Computation');

if realtime then 
  arc=InitDrawBall(nBalls,3,0.1); 
end;
dixpc=0;

while k < N do 

  // transfer of state i+1 into state i and time incrementation
  sicSTNextStep(),

  // get current time step
  k = sicTimeGetK(),
  
  // solve ... 
  sicSTComputeFreeState(),
  sicSTcomputePb(),

  // update
  sicSTupdateState(),
	
 // --- Get values to be plotted ---

 dixpc=dixpc+1/N;
 
 if (dixpc>=0.01) then
   waitbar(k/N,winId);
   dixpc=0;
 end
 
 for i = 1:nBalls,
   Qx(i,k+1)=sicModelgetQ(i-1,0); 
   Qy(i,k+1)=sicModelgetQ(i-1,1);
   if realtime then 
     DrawBall(arc,i,Qx(i,k+1),Qy(i,k+1));
   end
 end //for
 
 Time(k+1)=k*sicTimeGetH();

end   

  winclose(winId);

endfunction


