
//
// arm2
// 

function [xd] = f_ctraj(t)
 xd=[1.0,0]';
endfunction

function [qd] = f_arm2GeomInv(xd)

ARM_L1 = 1.0;
ARM_L2 = 1.0;

k1 = (xd(1)*xd(1) + xd(2)*xd(2) - ARM_L1*ARM_L1  - ARM_L2*ARM_L2)/(2*ARM_L1*ARM_L2);
qd(2) = acos(k1);
if (qd(2) > 0) then
 qd(1) = - qd(1);
end;

k1 = atan(xd(2),xd(1));
k2 = (xd(1)*xd(1) + xd(2)*xd(2) + ARM_L1*ARM_L1 - ARM_L2*ARM_L2);
k2 = k2/(2*ARM_L1*ARM_L1*sqrt(xd(1)*xd(1) + xd(2)*xd(2)));
k2 = acos(k2);
	 
if (k2 < 0) then
  qd(1) = k1 + k2;
else
  qd(1) = k1 - k2;
end;

endfunction
