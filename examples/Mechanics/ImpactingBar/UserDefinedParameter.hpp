// User-defined main parameters
unsigned int nDof = 10;// degrees of freedom for the beam
double t0 = 1e-8;                   // initial computation time
double T = 0.0015;                  // final computation time
double h = 1e-7;                // time step
double position_init = 0.00005;      // initial position
double velocity_init =  -.1;      // initial velocity
double epsilon = 0.0;//1e-1;
double theta = 1/2.0 + epsilon;              // theta for MoreauJeanOSI integrator
//theta = 1.0;
double E = 210e9; // young Modulus
double S = 0.000314; //  Beam Section 1 cm  for the diameter
//S=0.1;
double L = 1.0; // length of the  beam
double rho = 7800.0 ; // specific mass
//rho=1.0;
double g = 9.81; // Gravity
g=0.0;

