// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.
//
// Copyright 2018 INRIA.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//	

exec '../../Front-End/scilab/loaderSiconos.sce';

getf("../../Front-End/scilab/visuballs.sci");



function Simul()

[t,qx,qy]=BouncingBall();

VisuBalls(qx,qy,1,1,0.1);


endfunction


function [Time,Qx,Qy] = BouncingBall()

// Bug dgetri ?
errcatch(998,'continue');

sicLoadModel('./BouncingBall_TIDS_sci.xml');
sicInitSimulation();

k = sicTimeGetK(); 
N = sicTimeGetN();

arc=InitDrawBall(1,1,0.1);

winId=waitbar('Siconos Computation');

dixpc=0;

while k < N do 

  // transfer of state i+1 into state i and time incrementation
  sicSTNextStep();
  
  // get current time step
  k = sicTimeGetK();

  // solve ... 
  sicSTComputeFreeState();
  sicSTcomputePb();

  // update
  sicSTupdateState();

  // --- Get values to be plotted ---
  dixpc=dixpc+1/N;
 if (dixpc>0.01) then
   waitbar(k/N,winId);
   dixpc=0;
 end
 
  Qx(1,k+1)=sicModelgetQ(0,0);
  Qy(1,k+1)=0;
  Time(k+1)=k*sicTimeGetH();
  DrawBall(arc,1,Qx(1,k+1),Qy(1,k+1));
  
end   

endfunction




