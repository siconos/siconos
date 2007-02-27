global GLOBAL_U;
global GLOBAL_DT;
global GLOBAL_HIST;


getf 'lib/modules.sci'
getf 'lib/simulation.sci'
getf 'lib/plot.sci'


getf 'lib/pa10/modele.sci'
getf 'pa10Control.sci'
getf 'pa10Jac.sci'
getf 'pa10Plot.sci'

pa10ModelInit();

// Test du modèle cinématique
function []= kinTest(q,i)

M1=pa10Kin(q);
M2 = pa10Tags(q);

printf("%f , %f \n",M1(i)(1,4),M2(i));
printf("%f , %f \n",M1(i)(2,4),M2(i+8));
printf("%f , %f \n",M1(i)(3,4),M2(i+16));

endfunction


q0=[0,0,0,0,0,0,0];
kinTest(q0,7);

q0=[%pi/2,0,0,0,0,0,0];
kinTest(q0,7);

q0=[0,%pi/2,0,0,0,0,0];
kinTest(q0,7);

q0=[0,0,%pi/2,0,0,0,0];
kinTest(q0,7);

q0=[0,0,0,%pi/2,0,0,0];
kinTest(q0,7);

q0=[0,0,0,0,%pi/2,0,0];
kinTest(q0,7);

q0=[0,0,0,0,0,%pi/2,0];
kinTest(q0,7);

q0=[0,0,0,0,0,0,%pi/2];
kinTest(q0,7);
