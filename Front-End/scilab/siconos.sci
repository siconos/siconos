//
// getf("./siconos.sci");
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
	link("/home/sed/pissard/Local/Workspace/siconos-user/lib/libSiconosKernel.so");
	link("/home/sed/pissard/Local/Workspace/siconos/trunk/Front-End/scilab/siconos.so",['simul'],'C');
endfunction



function BouncingBall()
        call("simul");
endfunction
