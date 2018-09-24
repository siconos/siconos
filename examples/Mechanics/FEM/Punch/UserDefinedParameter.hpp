// User-defined main parameters
double t0 = 1e-8;                   // initial computation time
double T = 0.0015;                  // final computation time
double h = 1e-4;                // time step

T=10000*h;

double position_init = 0.00000;      // initial position
double velocity_init =  -.1;      // initial velocity

double epsilon = 0.0;//1e-1;
double theta = 1/2.0 + epsilon;              // theta for MoreauJeanOSI integrator
//theta = 1.0;

double E = 210e9; // young Modulus
double nu = 0.25; // young Modulus
double rho = 7800.0 ; // specific mass

//rho=1.0;
double g = 9.81; // Gravity
g=0.0;

