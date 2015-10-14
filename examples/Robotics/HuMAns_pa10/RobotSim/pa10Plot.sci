
function[]=pa10Plot2D(index)

// 1=q, 2=qdot, 3=U, 4=User (qd)

// q, qd axis index
subplot(2,1,1);
simulationPlot2D(1,index);
//simulationPlot2D(4,index);

// Error xe
simulationPlot2D(4,index);

// U
subplot(2,1,2);
simulationPlot2D(3,index);

endfunction

function[]=pa10Plot3D()

   pa10IndexTags=[1,7,9,15,17,23];
   simulationPlot3D(pa10Tags,pa10IndexTags,1);

endfunction

function[]=pa10PlotQ(Q)

   pa10IndexTags=[1,7,9,15,17,23];
   simulationPlotQ(Q,pa10Tags,pa10IndexTags,1)

endfunction
