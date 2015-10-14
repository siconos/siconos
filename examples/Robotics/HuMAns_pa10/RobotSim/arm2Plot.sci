
function[]=arm2Plot2D()


// Compute time scale absisse
time=(0:GLOBAL_DT:(length(GLOBAL_HIST(1)(2,:))-1)*GLOBAL_DT);

//
// plot Q
//
subplot(2,1,1)
simulationPlot2D(1,1);
simulationPlot2D(1,2);

//
// plot QD
//
xset('color',1)
simulationPlot2D(3,1);
simulationPlot2D(3,2);

//
// plot Torque
//
subplot(2,1,2)
simulationPlot2D(2,1);
simulationPlot2D(2,2);

endfunction

function[]=arm2Plot3D()

// arm2IndexTags=[1,3,3,5,5,8];
// simulationPlot3D(arm2Tags,arm2IndexTags,1);
   printf("Tags not generated\n");

endfunction

function[]=arm2PlotQ(Q)

//   arm2IndexTags=[1,3,3,5,5,8];
//   simulationPlotQ(Q,arm2Tags,arm2IndexTags,1);
     printf("Tags not generated\n");

endfunction
