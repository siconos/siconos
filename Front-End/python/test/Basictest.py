# SICONOS - PYTHON : bouncing ball example

import pySiconos;
model = pySiconos.Model(0, 10);
model.display();

nsds = model.createNonSmoothDynamicalSystem( False );
nsds.display();
model.setNonSmoothDynamicalSystemPtr( nsds );
