function []=simulationPlot2D(indextype,axis)

 xset("window",0);

// Compute time scale absisse
time=(0:GLOBAL_DT:(length(GLOBAL_HIST(1)(2,:))-1)*GLOBAL_DT);


plot2d(time,GLOBAL_HIST(indextype)(axis,:)',indextype);


endfunction


function []=simulationPlot3D(funcTags,tabTags,sw)

 xset("window",1)
 xset('pixmap', 1);

for k = 1:size(GLOBAL_HIST(1), 2),
    tag = funcTags(GLOBAL_HIST(1)(:,k)); 
    xset("wwpc");
    param3d1(tag(tabTags(1):tabTags(2)),tag(tabTags(3):tabTags(4)),tag(tabTags(5):tabTags(6)),75, 75,"X@Y@Z",[1,4],[-sw,sw,-sw,sw,-sw,sw]);
    xset("wshow");
    xbasc();
end;

param3d1(tag(tabTags(1):tabTags(2)),tag(tabTags(3):tabTags(4)),tag(tabTags(5):tabTags(6)),75, 75,"X@Y@Z",[1,4],[-sw,sw,-sw,sw,-sw,sw]);

xset('pixmap', 0);

endfunction



function []=simulationPlotQ(Q,funcTags,tabTags,sw)

tag = funcTags(Q);

param3d1(tag(1:7),tag(9:15),tag(17:23),75, 75,"X@Y@Z",[1,4],[-sw,sw,-sw,sw,-sw,sw]);

endfunction
