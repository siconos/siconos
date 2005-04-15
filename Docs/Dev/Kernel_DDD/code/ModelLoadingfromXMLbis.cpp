LagrangianLinearTIDS *lltids1;
// Constructor of a LagrangianLinearTIDS with minimal data
lltids1 = new LagrangianLinearTIDS(1, 3, &q0, &v0, &mass,
                                   "BasicPlugin:FExt", &K, &C);
nsds1.addDS(lltids1);
