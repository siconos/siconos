//
// getf("./siconos.sci");
//
// LinkSiconos();
// sicLoadModel('./BouncingBall_TIDS.xml');
//

//
//
// configure may generate paths :
//   for standard lib 
//   for siconos and numerics lib
//

function LinkSiconos()
	// dynamic lib
	link("/usr/lib/libdl.so");
	// xml
	link("/usr/lib/libxml2.so");
	// lapack++
	link("//usr/lib/liblapack.so");
	link("/usr/local/lib/liblamatrix++.so");
	link("/usr/local/lib/libblas++.so");
	link("/usr/local/lib/liblapack++.so");
	// Numerics and Siconos
	link("/local_home/pissard/Workspace/siconos/trunk/Numerics/lib/libNumerics.so");
	link("/local_home/pissard/Workspace/siconos-user/lib/libSiconosKernel.so");
	link("/local_home/pissard/Workspace/siconos/trunk/Front-End/scilab/siconos.so",['simul','sicLoadModel'],'C');
endfunction



function BouncingBall()
        call("simul");
endfunction

function plotBall()
	Balls=fscanfMat('result.dat');
	plot2d(Balls(:,2));
endfunction


function ret_val=sicLoadModel(ModelXml)
        ret_val=call("sicLoadModel",ModelXml,2,"c","out",[1,1],1,"i");
endfunction

// Test
//
// LinkSiconos();
// sicLoadModel('./BouncingBall_TIDS.xml');
// 
