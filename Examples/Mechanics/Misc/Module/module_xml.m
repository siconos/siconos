
%%%%%%%%%%%            %%%%%%%%%%%%%%
%            Données                %
%%%%%%%%%%             %%%%%%%%%%%%%%

t0 = 0.0;
T  = 0.30;
title = 'ModuleLevogyre';
author = 'SN_DD';
description='description';
date1 =date;
nb_ds = 1;
id = 'Bar2D';
step_mem= 1;
q0_size = 6;
q0 = [0.0 0.0 0.1 0.0 0.0 0.1];
V0_size = 6;
V0 =[0.0  0.0 -0.5 0.0 0.0 -0.5 ];
M_row = 6;
M_col = 6;
M = [0.728 0.0 0.0  0.364 0.0 0.0
     0.0 0.728 0.0 0.0 0.364 0.0
     0.0 0.0 0.728 0.0 0.0 0.328
     0.364 0.0 0.0 0.728 0.0 0.0
     0.0 0.364 0.0 0.0 0.728 0.0
     0.0 0.0 0.364 0.0 0.0 0.728];
		
ndof = 6;
K_row = 6;
K_col = 6;
K = [0.0 0.0 0.0 0.0 0.0 0.0
     0.0 5.6e7 0.0 0.0 -5.6e7 0.0
     0.0 0.0 0.0 0.0 0.0 0.0
     0.0 0.0 0.0 0.0 0.0 0.0
     0.0 -5.6e7 0.0 0.0 5.6e7 0.0
     0.0 0.0 0.0 0.0 0.0 0.0];


C_row = 6;
C_col = 6;
C = zeros(6,6);
nb_inter = 1;
nInter =1;
DSList_size =2;
DSList = [1 2];
H_row = 1;
H_col = 7;
H=[-0.7071 0.0 0.7071 0.0 0.0 0.0 0.0];
b_size = 1;
b = 0.0;
e = 0.5;
simulation_type='TimeStepping';
true = 'true';
h = 0.001;
DS_concerned_size =1;
DS_number = 2;
theta = 0.5;
n = 1;
Inter_num =1;
maxiter = 10;
tol = 0.0001;
%%%%%%%%%%%            %%%%%%%%%%%%%%
%           Fin données             %
%%%%%%%%%%             %%%%%%%%%%%%%%


fid=fopen('fichier_xml.xml','w');
fprintf(fid,'%s\n','<SiconosModel>');
fprintf(fid,'  %s%s%s\n','<Title>',title,'</Title>');
fprintf(fid,'  %s%s%s\n','<Author>',author,'</Author>');
fprintf(fid,'  %s%s%s\n','<Description>',description,'</Description>');
fprintf(fid,'  %s%s%s\n','<Date>',date1,'</Date>');
fprintf(fid,'  %s\n',' <SchemaXML>/share/SICONOS/SiconosModelSchema-V1.2.xsd</SchemaXML>');
fprintf(fid,'     %s\n','<Time>');
fprintf(fid,'        %s %g %s\n','<t0>',t0,'</t0>');
fprintf(fid,'        %s %g %s\n','<T>',T,'</T>');
fprintf(fid,'     %s\n','</Time>');
fprintf(fid,'     %s%s%s\n','<NSDS bvp=''','false','''>');
fprintf(fid,'        %s\n','<DS_Definition>');

%%%%%%%%%%%            %%%%%%%%%%%%%%
%         Boucle sur les DS         %
%%%%%%%%%%             %%%%%%%%%%%%%%
for i_ds = 1: nb_ds

fprintf(fid,'           %s%i%s\n','<LagrangianLinearTIDS number=''',i_ds,'''>');
fprintf(fid,'           %s %s %s\n','<Id>',id,'</Id>');
fprintf(fid,'           %s%i%s\n','<StepsInMemory>',step_mem,'</StepsInMemory>');
fprintf(fid,'           %s%i%s\n','<q0 vectorSize=''',q0_size,'''>');
fprintf(fid,'              ');
fprintf(fid,'%g ',q0);
%for i = 1: q0_size
%fprintf(fid,'%i ',q0(i));
%end
fprintf(fid,'\n            %s \n','</q0>');
fprintf(fid,'           %s%i%s\n','<Velocity0 vectorSize=''',V0_size,'''>');
fprintf(fid,'              ');
fprintf(fid,'%g ',V0);
fprintf(fid,'\n            %s \n','</Velocity0>');
fprintf(fid,'              %s \n','<Fext vectorPlugin="Bar2DPlugin:bar2DFExt"/>');
fprintf(fid,'              %s%i%s%i%s\n','<M matrixRowSize=''',M_row,''' matrixColSize=''',M_col,'''>');
for i = 1:M_row
fprintf(fid,'                 %s','<row>');
fprintf(fid,'%g ',M(i,:));
fprintf(fid,'%s','</row>');
fprintf(fid,'\n');
end
fprintf(fid,'              %s\n','</M>');
fprintf(fid,'              %s%i%s\n','<ndof>',ndof,'</ndof>');
fprintf(fid,'              %s%i%s%i%s\n','<K matrixRowSize=''',K_row,''' matrixColSize=''',K_col,'''>');
for i = 1:K_row
fprintf(fid,'                 %s','<row>');
fprintf(fid,'%g ',K(i,:));
fprintf(fid,'%s','</row>');
fprintf(fid,'\n');
end
fprintf(fid,'              %s\n','</K>');
fprintf(fid,'              %s%i%s%i%s\n','<C matrixRowSize=''',C_row,''' matrixColSize=''',C_col,'''>');
for i = 1:C_row
fprintf(fid,'                 %s','<row>');
fprintf(fid,'%g ',C(i,:));
fprintf(fid,'%s','</row>');
fprintf(fid,'\n');
end
fprintf(fid,'              %s\n','</C>');
fprintf(fid,'           %s\n','</LagrangianLinearTIDS>');

end %for i_ds
%%%%%%%%%%%            %%%%%%%%%%%%%%
%       Fin boucle sur les DS       %
%%%%%%%%%%             %%%%%%%%%%%%%%


nb_ds =nb_ds +1;
i_ds = nb_ds;
q0_size = 1;
q0 = 0.0;
V0_size = 1;
V0 = 0.0;
M_row =1;
M_col=1;
M= 1.0;
ndof =1;
K_row= 1;
K_col=1;
K= 0.0;
C_row=1;
C_col=1;
C=0.0;




fprintf(fid,'           %s%i%s\n','<LagrangianLinearTIDS number=''',i_ds,'''>');
fprintf(fid,'           %s %s %s\n','<Id>','Ground','</Id>');
fprintf(fid,'           %s %i %s\n','<StepsInMemory>',step_mem,'</StepsInMemory>');
fprintf(fid,'           %s%i%s\n','<q0 vectorSize=''',q0_size,'''>');
fprintf(fid,'              ');
fprintf(fid,'%g ',q0);
%for i = 1: q0_size
%fprintf(fid,'%i ',q0(i));
%end
fprintf(fid,'\n            %s \n','</q0>');
fprintf(fid,'           %s%i%s\n','<Velocity0 vectorSize=''',V0_size,'''>');
fprintf(fid,'              ');
fprintf(fid,'%g ',V0);
fprintf(fid,'\n            %s \n','</Velocity0>');
fprintf(fid,'              %s \n','<Fext vectorPlugin="Bar2DPlugin:groundFExt"/>');
fprintf(fid,'              %s%i%s%i%s\n','<M matrixRowSize=''',M_row,''' matrixColSize=''',M_col,'''>');
for i = 1:M_row
fprintf(fid,'                 %s','<row>');
fprintf(fid,'%g ',M(i,:));
fprintf(fid,'%s','</row>');
fprintf(fid,'\n');
end
fprintf(fid,'              %s\n','</M>');
fprintf(fid,'              %s%i%s\n','<ndof>',ndof,'</ndof>');
fprintf(fid,'              %s%i%s%i%s\n','<K matrixRowSize=''',K_row,''' matrixColSize=''',K_col,'''>');
for i = 1:K_row
fprintf(fid,'                 %s','<row>');
fprintf(fid,'%g ',K(i,:));
fprintf(fid,'%s','</row>');
fprintf(fid,'\n');
end
fprintf(fid,'              %s\n','</K>');
fprintf(fid,'              %s%i%s%i%s\n','<C matrixRowSize=''',C_row,''' matrixColSize=''',C_col,'''>');
for i = 1:C_row
fprintf(fid,'                 %s','<row>');
fprintf(fid,'%g ',C(i,:));
fprintf(fid,'%s','</row>');
fprintf(fid,'\n');
end
fprintf(fid,'              %s\n','</C>');
fprintf(fid,'           %s\n','</LagrangianLinearTIDS>');



fprintf(fid,'        %s\n','</DS_Definition>');

fprintf(fid,'\n\n        %s\n','<Interaction_Definition>');

%%%%%%%%%%%            %%%%%%%%%%%%%%
%     Boucle sur les interactions   %
%%%%%%%%%%             %%%%%%%%%%%%%%

for i_inter =1:nb_inter

fprintf(fid,'            %s%i%s\n','<Interaction number=''',i_inter,'''>');
fprintf(fid,'               %s%s%s%s%s\n','<Id>',id,'-','Ground','</Id>');
fprintf(fid,'               %s %i %s\n','<nInter>',nInter,'</nInter>');
fprintf(fid,'               %s\n','<DS_Concerned>');

fprintf(fid,'                  %s%i%s','<DSList vectorSize=''',DSList_size,'''>');
%fprintf(fid,'%i ',DSList);
i=1;
fprintf(fid,'%i',DSList(i));
for i = 2: DSList_size
fprintf(fid,' %i',DSList(i));
end
fprintf(fid,'%s\n','</DSList>');
fprintf(fid,'               %s\n','</DS_Concerned>');
fprintf(fid,'               %s\n','<Interaction_Content>');
fprintf(fid,'                  %s\n','<LagrangianLinearRelation>');
fprintf(fid,'                     %s%i%s%i%s\n','<H matrixRowSize=''',H_row,''' matrixColSize=''',H_col,'''>');
fprintf(fid,'                     %s','<row>');
fprintf(fid,'%g ',H);
fprintf(fid,'%s\n ','</row>');
fprintf(fid,'                     %s\n','</H>');
fprintf(fid,'                     %s%i%s\n','<b vectorSize=''', b_size,'''>');
fprintf(fid,'                        ');
fprintf(fid,'%g ',b);
fprintf(fid,'\n                     %s\n','</b>');
fprintf(fid,'                  %s\n','</LagrangianLinearRelation>');
fprintf(fid,'                  %s\n','<NewtonImpactLaw>');
fprintf(fid,'                  %s %g %s\n','<e>',e,'</e>');
fprintf(fid,'                  %s\n','</NewtonImpactLaw>');
fprintf(fid,'               %s\n','</Interaction_Content>');
fprintf(fid,'            %s\n','</Interaction>');

end % i_inter

%%%%%%%%%%%            %%%%%%%%%%%%%%
%  Fin boucle sur les interactions  %
%%%%%%%%%%             %%%%%%%%%%%%%%


fprintf(fid,'        %s\n\n','</Interaction_Definition>');

fprintf(fid,'     %s\n\n','</NSDS>');

fprintf(fid,'     %s%s%s\n\n','<Simulation type=''',simulation_type,'''>');
fprintf(fid,'        %s%s%s\n','<TimeDiscretisation isConstant=''',true,'''>');
fprintf(fid,'           %s %g %s\n','<h>',h,'</h>');
fprintf(fid,'        %s\n','</TimeDiscretisation>');

fprintf(fid,'        %s\n','<OneStepIntegrator_Definition>');

%%%%%%%%%%%            %%%%%%%%%%%%%%
%        Boucle sur les DS          %
%%%%%%%%%%             %%%%%%%%%%%%%%

for DS_number = 1:nb_ds 

fprintf(fid,'           %s\n','<Moreau>');
fprintf(fid,'              %s%i%s\n','<DS_Concerned size=''',DS_concerned_size,'''>');
fprintf(fid,'                 %s%i%s\n','<DS number=''',DS_number,'''/>');
fprintf(fid,'              %s\n','</DS_Concerned>');
fprintf(fid,'              %s %g %s\n','<Theta>',theta,'</Theta>');
fprintf(fid,'           %s\n','</Moreau>');

end % DS_number
%%%%%%%%%%%            %%%%%%%%%%%%%%
%       Fin  boucle sur les DS      %
%%%%%%%%%%             %%%%%%%%%%%%%%


fprintf(fid,'        %s\n','</OneStepIntegrator_Definition>');

fprintf(fid,'        %s\n','<OneStepNSProblem>');
fprintf(fid,'           %s\n','<LCP>');

%%%%%%%%%%%            %%%%%%%%%%%%%%
%    Boucle sur les interactions    %
%%%%%%%%%%             %%%%%%%%%%%%%%

for Inter_num = 1: nb_inter

fprintf(fid,'              %s%i%s\n','<n>',n,'</n>');
fprintf(fid,'              %s\n','<Interaction_Concerned>');
fprintf(fid,'                 %s%i%s\n','<Interaction number=''',Inter_num,'''/>');
fprintf(fid,'              %s\n','</Interaction_Concerned>');

end %Inter_num
%%%%%%%%%%%            %%%%%%%%%%%%%%
%   Fin boucle sur les interactions %
%%%%%%%%%%             %%%%%%%%%%%%%%

fprintf(fid,'           %s\n','</LCP>');
fprintf(fid,'           %s\n','<Solver>');


fprintf(fid,'              %s\n','<LcpSolving>');
fprintf(fid,'                 %s%i%s%g%s\n','<NSQP maxIter="', maxiter,'" tolerance="', tol,'"/>');
fprintf(fid,'              %s\n','</LcpSolving>');
fprintf(fid,'           %s\n','</Solver>');


fprintf(fid,'        %s\n','</OneStepNSProblem>');
fprintf(fid,'     %s\n','</Simulation>');
fprintf(fid,'%s\n','</SiconosModel>');
fclose(fid);
