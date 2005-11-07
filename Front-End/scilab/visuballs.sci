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

function [arc]=InitDrawBall(nBalls,high,sizeball)

set("figure_style","new");
xset("pixmap",1);
plot2d([-1;1],[0;0],2,"031"," ",[-1,-1,1,high]);

arc=get("hdl"); 

for i = 1:nBalls,
  xarc(0,0,sizeball,sizeball,0,360*64);
  //get handle on current entity (here the arc entity)
  arc(i)=get("hdl");
  arc(i).fill_mode="on";
  arc(i).foreground=i;
  arc(i).visible="on";  
end; //for

xset("wshow");
 

endfunction


function DrawBall(arc,index,x,y)


if ((index>0)&(index<=size(arc,1))) then
 arc(index).data(1)=x;
 arc(index).data(2)=y;
end; //if
xset("wshow");

endfunction

function VisuBalls(q,step,high,sizeball)
 
 nBalls=size(q,1);
 N =size(q,2);

 arc=InitDrawBall(nBalls,high,sizeball);

 xset("pixmap",0);
 xset("window",3);
 
 for i = 1:nBalls,
    plot2d(t,q(i,:),style=i);
 end; //for
   
k=0;
while k <= N do 
 
  for i = 1:nBalls,
    DrawBall(arc,i,0,q(i,k+1));
  end
  k=k+step;
  
 end;
 
endfunction
