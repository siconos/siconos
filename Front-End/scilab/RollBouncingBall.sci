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
cd  "../../sample/RollBouncingBall"

function Simul()

sicLink();

[t,q]=RollBouncingBall();

xset("pixmap",0);
xset("window",3);
plot2d(t,q);

endfunction

function [arc]=InitDrawBall()

set("figure_style","new");
xset("pixmap",1);
plot2d([-1;1],[0;0],2,"031"," ",[-1,-1,1,1]);
xarc(0,0,0.1,0.1,0,360*64);
arc=get("hdl"); //get handle on current entity (here the arc entity)
arc.fill_mode="on";
arc.foreground=2;
arc.visible="off";  
xset("wshow");

endfunction

function DrawBall(arc,x,y)

arc.visible="on";  
arc.data(1)=x;
arc.data(2)=y+0.1;
xset("wshow");

endfunction


function [Time,Q] = RollBouncingBall()

// Bug dgetri ?
errcatch(998,'continue');


sicLink();
sicLoadModel('./Ball.xml');
sicInitStrategy();

k = sicTimeGetK();
N = sicTimeGetN();

arc=InitDrawBall();

while k < N do 

  Q(k+1)=sicModelgetQ(0),
  
  Time(k+1)=k*sicTimeGetH(),
  
  DrawBall(arc,0,Q(k+1));

  sicSTNextStep(),

  k = sicTimeGetK(),

  sicSTComputeFreeState(),

//  sicSTformalisePb(),

  sicSTcomputePb(),

  sicSTupdateState(),

end   

endfunction


