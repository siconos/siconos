//
// Matrice de passage DH modifié
//
function [M]= DHmat(theta,r,alpha,lambda)

M=[cos(theta) -sin(theta) 0 r;
   cos(alpha)*sin(theta) cos(alpha)*cos(theta) -sin(alpha) -sin(alpha)*lambda;
   sin(alpha)*sin(theta) sin(alpha)*cos(theta) cos(alpha)  cos(alpha)*lambda;
   0 0 0 1];
endfunction


//
// Composition des matrices de DH
//
function[M]=robotKin(theta,r,alpha,lambda,ndof)

M=list();

M(1)=DHmat(theta(1),r(1),alpha(1),lambda(1));

for i = 2:ndof,
   M(i)=M(i-1)*DHmat(theta(i),r(i),alpha(i),lambda(i));
end;

endfunction

//
// Modele geometrique du PA10
//
function[M]=pa10Kin(q)

r=[0 0 0 0 0 0 0 0];
lambda=[0 0 0.45 0 0.48 0 0];
alpha=[0 -%pi/2 %pi/2 -%pi/2 %pi/2 -%pi/2 %pi/2];

M=robotKin(q,r,alpha,lambda,7);

endfunction;


//
// Jacobienne du PA10
//
//
function[J]=pa10Jac(q)

M=pa10Kin(q);
ndof=7;

J=zeros(6,7);

for i = 1:ndof,
 x=M(ndof)(1,4) - M(i)(1,4);
 y=M(ndof)(2,4) - M(i)(2,4);
 z=M(ndof)(3,4) - M(i)(3,4);

 J(1,i)= M(i)(2,3)*z - M(i)(3,3)*y;
 J(2,i)= M(i)(3,3)*x - M(i)(1,3)*z;
 J(3,i)= M(i)(1,3)*y - M(i)(2,3)*x;

 J(4,i)=M(i)(1,3);
 J(5,i)=M(i)(2,3);
 J(6,i)=M(i)(3,3);

end

endfunction;
