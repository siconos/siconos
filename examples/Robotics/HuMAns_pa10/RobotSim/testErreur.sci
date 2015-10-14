getf 'lib/modules.sci'
getf 'lib/pa10/modele.sci'
getf 'lib/pa10/pa10Jac.sci'



function [err] = xerr(x1,x2)

m1=f_Hmat(x1);
m2=f_Hmat(x2);

err= f_xerror(m1,m2);

endfunction

xerr([0 0 0 0 0 0],[0 0 0 0 0 0])

xerr([0 0 0 0 0 0],[0 0 0 0.1 0 0])

xerr([0 0 0 0 0 0],[0 0 0 0 0.1 0])

xerr([0 0 0 0 0 0.2],[0 0 0 0 0.0 0.1])

xerr([0 0 0 0 0.01 0.02],[0 0 0  0.01 0 0.01])
