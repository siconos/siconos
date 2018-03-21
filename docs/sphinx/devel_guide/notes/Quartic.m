M=[15,-1,-2;-1,14,1;-2,1,16];
q=[-1;2;3];
q2=[q(2);q(3)];
M2=[M(2,2),M(2,3);M(3,2),M(3,3)];
M1_=[M(2,1);M(3,1)];
mu=0.1;
NormD2=M(1,2)*M(1,2)+M(1,3)*M(1,3);
NormD=sqrt(NormD2);
e=mu*NormD/M(1,1);
d=abs(q(1))/NormD;
p=e*d;
OD=[-(M(1,2)*q(1))/NormD2;-((-M(1,2)*M(1,2)*q(1))/NormD2+q(1))/M(1,3)];
D_Dir=[1;-M(1,2)/M(1,3)];

[Vaux,D] = eig(M2);
d1=D(1,1);
d2=D(2,2);
Vaux2=Vaux';
V=[V(2,1),V(2,2);V(1,1),V(1,2)];
d1=D(2,2);
d2=D(1,1);
OD2=V*OD;
D_Dir2=V*D_Dir;
phi=atan(OD2(2)/OD2(1));
if(OD2(1)<0)
    phi=phi+acos(-1);
end
cosphi=cos(phi);
sinphi=sin(phi);

q2b=V*q2;
M1_b=V*M1_;
a1=-M1_b(1)/mu;
a2=-M1_b(2)/mu;
AA=-e*q2b(2)*cosphi;
BB=e*q2b(1)*sinphi;
CC=(d1-d2)*p+e*cosphi*q2b(1)-e*sinphi*q2b(2);
DD=q2b(1)-p*a1;
EE=-q2b(2)+p*a2;


P4=AA-EE;
P3=-2*CC+2*DD;
P2=4*BB-2*AA;
P1=2*CC+2*DD;
P0=AA+EE;

poly = [P4 P3 P2 P1 P0];
racines=roots(poly);

r1=racines(4);
theta=2*atan(r1);
r=p/(1+e*cos(theta-phi));
RTb=[r*cos(theta);r*sin(theta)];
alpha1=(-q2b(2)+a2*r)/RTb(2)-d2;
alpha2=(-q2b(1)+a1*r)/RTb(1)-d1;
alpha=alpha1;

zero1=(d1-d2)*r*cos(theta)*sin(theta)+q2b(1)*sin(theta)-q2b(2)*cos(theta)-r*(a1*sin(theta)-a2*cos(theta));
zero2=(d1-d2)*RTb(1)*RTb(2)+q2b(1)*RTb(2)-q2b(2)*RTb(1)-(a1*RTb(2)-a2*RTb(1))*sqrt(RTb(1)*RTb(1)+RTb(2)*RTb(2));
RT=V'*RTb
RESIDUb=(r/mu)*M1_b+([alpha+d1,0;0,alpha+d2])*RTb+q2b
RESIDU=(r/mu)*M1_+(M2+[alpha,0;0,alpha])*RT+q2
ONELIPSE2=norm(RTb)/sqrt((D_Dir2(2)*RTb(1)-D_Dir2(1)*RTb(2)-OD2(1)*D_Dir2(2)+D_Dir2(1)*OD2(2))^2/(D_Dir2(1)*D_Dir2(1)+D_Dir2(2)*D_Dir2(2)))-e
ONELIPSE=(M(1,1)*M(1,1)/(mu*mu))*r*r-(q(1)+M(1,2)*RT(1)+M(1,3)*RT(2))^2
s1=M(1,1)*norm(RT)/mu
s2=-q(1)-M(1,2)*RT(1)-M(1,3)*RT(2)
R=[norm(RT)/mu;RT]
M*R+q