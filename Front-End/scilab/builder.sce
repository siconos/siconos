// This is the builder.sce 
// must be run from this directory 

lines(0);

ilib_name  = 'libSiconosApiSci' 		// interface library name 


// objects files 

files = ['SiconosSci.o'];

libs  = [] 				// other libs needed for linking

// table of (scilab_name,interface-name) 
// for fortran coded interface use 'C2F(name)'

table =['sicLoadModel',	'sicLoadModelInterface';
        'sicInitSimulation', 'sicInitSimulationInterface';
	'sicTimeGetH', 'sicTimeGetHInterface';
	'sicTimeGetN', 'sicTimeGetNInterface';
	'sicSTNextStep','sicSTNextStepInterface';
	'sicSTSaveInMemory','sicSTSaveInMemoryInterface';
	'sicSTComputeOneStep','sicSTComputeOneStepInterface';
	'sicSTnewtonSolve','sicSTnewtonSolveInterface';
	'sicSTupdateState', 'sicSTupdateStateInterface';
	'sicAdvanceToEvent','sicAdvanceToEventInterface';
	'sicProcessEvents','sicProcessEventsInterface';
	'sicHasNextEvent','sicHasNextEventInterface';
	'sicGetTypeEvent','sicGetTypeEventInterface';
	'sicModelgetQ','sicModelgetQInterface';
	'sicLagrangianLinearTIDS','sicLagrangianLinearTIDSInterface';
	'sicLagrangianDS','sicLagrangianDSInterface';
	'sicSetMass','sicSetMassInterface';
	'sicSetNNL','sicSetNNLInterface';
	'sicSetJacQNNL','sicSetJacQNNLInterface';
	'sicSetJacVelNNL','sicSetJacVelNNLInterface';
	'sicSetFInt', 'sicSetFIntInterface';
	'sicSetJacQFInt', 'sicSetJacQFIntInterface';
	'sicSetJacVelFInt', 'sicSetJacVelFIntInterface';
	'sicSetFExt','sicSetFExtInterface';
	'sicInteraction', 'sicInteractionInterface';
	'sicLagrangianLinearR','sicLagrangianLinearRInterface';
	'sicNewtonImpactLawNSL','sicNewtonImpactLawNSLInterface';
	'sicNonSmoothDynamicalSystem','sicNonSmoothDynamicalSystemInterface';
	'sicModel', 'sicModelInterface';
	'sicSimulationTimeStepping','sicSimulationTimeSteppingInterface';
	'sicOneStepIntegratorMoreau','sicOneStepIntegratorMoreauInterface';
	'sicOneStepNSProblemLCP','sicOneStepNSProblemLCPInterface';
	'sicClean','sicCleanInterface'];
    

// extra parameters can be transmited to the linker 
// and to the C and Fortran compilers with 
// ldflags,cflags,fflags 
// for example to link a set of routines using the 
// ImageMagick library 
//  ldflags = "`Magick-config --ldflags --libs`"; 
//  cflags  = "`Magick-config --cppflags`"; 
//  fflags   = ""; 

ldflags = "";
cflags ="";
fflags ="";

// do not modify below 
// ----------------------------------------------
ilib_build(ilib_name,table,files,libs,'Makelib',ldflags,cflags,fflags);

quit;








