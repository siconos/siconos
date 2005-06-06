//
// getf("./sample.sci");
//

getf("./siconos.sci");

function Simul()

sicLink();

[t,q]=BouncingBall();

xset("pixmap",0);
xset("window",3);
plot2d(t,q);

endfunction

function [arc]=InitDrawBall()

set("figure_style","new");
xset("pixmap",1);
plot2d([-1;1],[0;0],2,"031"," ",[-1,0,1,1]);
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

function [Time,Q] = BouncingBall()

// Bug dgetri ?
errcatch(998,'continue');


sicLink();
sicLoadModel('./BouncingBall_TIDS.xml');
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

  sicSTformalisePb(),

  sicSTcomputePb(),

  sicSTupdateState(),

end   

endfunction

function BouncingBallC()
sicLink();
call("simul");
endfunction

function plotBallC()
// Think to cut the first line of result.dat
Balls=fscanfMat('result.dat');
plot2d(Balls(:,2));
endfunction
