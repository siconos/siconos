// Siconos version 1.0, Copyright INRIA 2005.  
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

getf("./siconos.sci");
getf("./visuballs.sci");

cd "../../sample/ThreeBeadsColumn";


function Simul()

sicLink();

[t,q]=ThreeBeadsColumn(%F);

VisuBalls(q,10,40,1);

endfunction


function [Time,Q] = ThreeBeadsColumn(realtime)

// Bug dgetri ?
errcatch(998,'continue');

nBalls=3;

sicLink();
sicLoadModel('./ThreeBeadsColumn.xml');

sicInitStrategy();

k = sicTimeGetK();
N = sicTimeGetN();

winId=waitbar('Siconos Computation');

if realtime then 
  arc=InitDrawBall(nBalls,40,1); 
end;

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

  waitbar(k/N,winId);

  for i = 1:nBalls,
    Q(i,k+1)=sicModelgetQ(i-1);
    if realtime then 
      DrawBall(arc,i,0,Q(i,k+1));

    end;
     if realtime then 
      xset("wshow");
    end;
  end; //for
  
  Time(k+1)=k*sicTimeGetH(),
 
 // sicDebug,

end   

  winclose(winId);

endfunction


