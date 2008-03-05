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




