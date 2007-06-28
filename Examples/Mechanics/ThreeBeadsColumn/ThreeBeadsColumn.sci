// Siconos-sample version 2.1.0, Copyright INRIA 2005-2006.
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


exec '../../../Front-End/scilab/loaderSiconos.sce';
getf("../../../Front-End/scilab/visuballs.sci");

function Simul()

[t,qx,qy]=ThreeBeadsColumn(%F);

VisuBalls(qx,qy,5,4,0.1);

endfunction


function [Time,Qx,Qy] = ThreeBeadsColumn(realtimeFlag)

nBalls=3;

sicLoadModel('./ThreeBeadsColumn_sci.xml');

sicInitSimulation();

k = 1;
N = sicTimeGetN();
h = sicTimeGetH();

winId=waitbar('Siconos Computation');

if realtimeFlag then 
  arc=InitDrawBall(nBalls,4,0.1); 
end;

dixpc=0;

while k < N do 

  // solve ... 
  sicSTComputeOneStep();
 // transfer of state i+1 into state i and time incrementation
  sicSTNextStep();
	
 // --- Get values to be plotted ---

 dixpc=dixpc+1/N;
 
 if (dixpc>=0.01) then
   waitbar(k/N,winId);
   dixpc=0;
 end
 
 for i = 1:nBalls,
   Qx(i,k+1)=sicModelgetQ(i-1,0); 
   Qy(i,k+1)=sicModelgetQ(i-1,1);
   if realtimeFlag then 
     DrawBall(arc,i,Qx(i,k+1),Qy(i,k+1));
   end
 end //for
 
 Time(k+1)=k*h;
 
 k=k+1;

end // while

  winclose(winId);

endfunction


