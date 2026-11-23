.. _beads_column_example:

Tutorial: a column of beads
==========================================

.. highlight:: c++

This tutorial deals with a column of beads subject to
gravity. It introduces linear Lagrangian systems with linear Lagrangian
relations and the Newton impact law.

Building the nonsmooth dynamical system
---------------------------------------

As shown in the figure below, we consider a ball of mass :math:`m` and radius :math:`R`, described by 3 generalized coordinates :math:`q=(z,x,\theta)`.
The ball is subject to gravity :math:`g`. The system also includes a rigid plane, defined by its position :math:`h` along
the axis Oz. We assume that the plane is fixed and not modeled as a dynamical system.

.. figure:: /figures/mechanics/BouncingBall/BouncingBall.*
   :width: 8cm
   :align: center

   fig 1: Coordinate system

The equation of motion of a ball is

.. math::

   M\ddot q = F_{ext}(t) + p \ \ with \ \ M=\left[\begin{array}{ccc}
   m &0 &0 \\
   0 & m & 0 \\
   0 & 0 & I
   \end{array}\right] \ \ , \ \ I = \frac{3}{5}mR^2 \ \ and \ \ F_{ext} = \left[\begin{array}{c}
   -m g \\
   0  \\
   0
   \end{array}\right]

with

* :math:`M` the inertia term, an :math:`n\times{}n` matrix;
* :math:`p` the force due to the nonsmooth law, i.e. the reaction at impact;
* :math:`F_{ext}(t):  \mathcal R \rightarrow \mathcal R^{n}` the given external force.

This is a linear time-invariant Lagrangian dynamical system, represented by the LagrangianLinearTIDS class (see :ref:`dynamical_systems`).

We now consider a column of "dsNumber" balls like the one above (mass and radius may differ; if necessary, variables are indexed by :math:`i`, the number of the ball). Each ball is governed by a linear system like the one written for a single ball, and may be in contact with the balls above and below it.

Let us now write the input file. As in the first tutorial, we create a new directory, multiBeads, and save the :ref:`template <template_siconos_driver>` as multiBeads.cpp.

We start by setting some parameters, like the number of balls, their initial positions and velocities and so on::

  // User-defined main parameters
  unsigned int dsNumber = 10;      // the number of dynamical systems
  int nDof = 3;           // degrees of freedom for beads
  double increment_position = 1;   // initial position increment from one DS to the following
  double increment_velocity = 0;   // initial velocity increment from one DS to the following
  double t0 = 0;                   // initial computation time
  double T = 10;                   // final computation time
  double h = 0.005;                // time step
  double position_init = 10.5;     // initial position of the lowest bead
  double velocity_init = 0.0;      // initial velocity of the lowest bead
  double R = 0.1;                  // radius of the balls

Then, we define the initial conditions of the balls, and create the
corresponding Dynamical Systems, all linear time-invariant Lagrangian systems.

All the systems are inserted in a container, a DynamicalSystemsSet,
named allDS.

From now on, to simplify, we assume that all
balls have the same mass, :math:`m = 1`, and the same radius, :math:`R=0.1`::

  // -------------------------
  // --- Dynamical systems ---
  // -------------------------

  // mass matrix, set to identity
  auto Mass = std::make_shared<siconos::algebra::SimpleMatrix>(nDof,nDof);
  Mass->setIdentity();
  (*Mass)(2,2) = 3.0/5*R*R;

  // -- Initial positions and velocities --
  // q0[i] and v0[i] correspond to position and velocity of ball i.
  vector<SimpleVector *> q0;
  vector<SimpleVector *> v0;
  q0.resize(dsNumber,NULL);
  v0.resize(dsNumber,NULL);

  for (unsigned i = 0; i < dsNumber; ++i)
  {
      // Memory allocation for q0[i] and v0[i]
      q0[i] = new SimpleVector(nDof);
      v0[i] = new SimpleVector(nDof);
      // set values
      (*(q0[i]))(0) = position_init;
      (*(v0[i]))(0) = velocity_init;
      // Create a new linear Lagrangian dynamical system, with q0[i] and v0[i] as initial conditions,
      // Mass as mass matrix and i as identification number.
      // The system is then inserted in allDS.
      allDS.insert( new LagrangianLinearTIDS(i,nDof,*(q0[i]),*(v0[i]),*Mass));
      // Increment values for next system
      position_init+= increment_position;
      velocity_init+= increment_velocity;
  }

Next, we define the external forces, i.e. gravity, applied to each ball. As explained in :ref:`dsPlugins`, a plugin function is available for these forces (for details on plugin functions, see :ref:`siconos_plugins`). Its signature (the types of its arguments) is given in DefaultPlugin.cpp. We copy it into a new file, say BeadsPlugin.cpp, and define an extern function, gravity::

  const double m = 1; // bead mass
  const double g = 9.81; // gravity
  extern "C" void gravity(unsigned int sizeOfq, double time,	double * fExt, double *param)
  {
      // set fExt components to 0
      for (unsigned int i = 0; i < sizeOfq; i++)
      fExt[i] = 0.0;
      // apply gravity
      fExt[0] = -m*g;
  }

Warning:

* gravity must be an *extern "C"* function, and the code is C, not C++;
* the name of the plugin file, BeadsPlugin.cpp here, must be xxxPlugin.cpp, xxx being anything you want.


Now we have to say "use gravity from BeadsPlugin.cpp to compute the external forces of my systems".
This is done with the "setComputeFExtFunction" function, in multiBeads.cpp::

   //
   CheckInsertDS checkDS;
   for (i=0;i<dsNumber;i++)
      {
        // Memory allocation for q0[i] and v0[i]
	q0[i] = new SimpleVector(nDof);
	v0[i] = new SimpleVector(nDof);
        // set values
	(*(q0[i]))(0) = position_init;
	(*(v0[i]))(0) = velocity_init;
        // Create and insert in allDS a new Lagrangian Linear Dynamical System ...
	checkDS = allDS.insert(new LagrangianLinearTIDS(i,nDof,*(q0[i]),*(v0[i]),*Mass));
        // Note that we now use a CheckInsertDS object: checkDS.first is
	// an iterator that points to the DS inserted above.
        //
	// Set the external forces of the last created system.
	(static_cast<LagrangianDS*>(*(checkDS.first)))->setComputeFExtFunction("BeadsPlugin.so", "gravity");
	// A cast is required, since allDS handles DynamicalSystem*,
	// not LagrangianLinearTIDS*.
	// Increment values for next system
	position_init+= increment_position;
	velocity_init+= increment_velocity;
      }

From now on, any computation of the external forces of a system in allDS calls the function gravity defined in BeadsPlugin.cpp.

*Remark:* :math:`m` and :math:`R` are set inside the BeadsPlugin file,
but it would also be possible, and maybe better, to pass them as
parameters of the gravity function
(see :ref:`siconos_plugins`).

The DynamicalSystems are now defined and all saved in allDS. Let us turn to the Interactions. They are also handled by a container, an InteractionsSet, named allInteractions. The potential interactions are the contacts between beads and the impact on the ground. Thus, for dsNumber systems, there are dsNumber-1 "bead-bead" Interactions plus one between the "bottom bead" and the floor.

We start with the bead-floor Interaction: the bottom ball bounces on the rigid plane, which introduces a constraint on the position of the ball:
:math:`z-R-h\geq 0`.
To define an Interaction, we first need a relation between the local variables at contact and the global coordinates.
As local variables of the Interaction, we introduce :math:`y`, the distance between the ball and the floor, and :math:`\lambda`, the multiplier that corresponds to
the reaction at contact. The relation is then:

.. math::

   y =& Hq + b = [1 \ 0 \ 0] q - R - h \\
   p =& H^t \lambda

(below, we set h=0).

.. compound::

   Finally, we need a nonsmooth law that defines the behaviour of the ball at impact.
   The unilateral constraint reads

   .. math:: 0 \leq y \perp \lambda \geq 0

   completed with a Newton impact law, with a restitution coefficient :math:`e` set to 0.9:

   .. math:: \textrm{if} \ y=0, \ \dot y(t^+) = -e \dot y(t^-)

   :math:`t^+` and :math:`t^-` being the post- and pre-impact times.

The first Interaction can then be constructed::

  // -------------------
  // --- Interactions---
  // -------------------
  InteractionsSet allInteractions;
  // The total number of Interactions
  int interactionNumber = dsNumber;
  // Interaction first bead and floor
  // A set for the systems handled by the "current" Interaction
  DynamicalSystemsSet dsConcerned;
  // Only the "bottom" bead is involved in this first Interaction,
  // therefore DynamicalSystem number 0.
  dsConcerned.insert(allDS.getDynamicalSystemPtr(0));
  // -- Newton impact law --
  double e = 0.9;
  NonSmoothLaw * nslaw0 = new NewtonImpactNSL(e);
  // Lagrangian Relation
  unsigned int interactionSize = 1; // y vector size
  SiconosMatrix *H = new SimpleMatrix(interactionSize,nDof);
  (*H)(0,0) = 1.0;
  SiconosVector *b = new SimpleVector(interactionSize);
  (*b)(0) = -R;
  Relation * relation0 = new LagrangianLinearR(*H,*b);
  // Interaction
  unsigned int num = 0 ; // an id number for the Interaction
  Interaction * inter0 = new Interaction("bead-floor", dsConcerned,num,interactionSize, nslaw0, relation0);
  allInteractions.insert(inter0);

In the same way, the potential contact between two balls introduces new constraints:

:math:`(z_i-R_i)-(z_j-R_j)-h \geq 0`, if ball :math:`i` is on top of ball :math:`j`.

So, for the Interaction between balls :math:`i` and :math:`j`, with :math:`y` the distance between the two balls and :math:`\lambda` the multiplier, we get:

.. math::

   y =& HQ + b = [-1 \ 0 \ 0 \ 1 \ 0 \ 0]Q + R_j-R_i-h \\
   p =& H^t \lambda \\
   Q =& \left[\begin{array}{c}
   q_j \\
   q_i
   \end{array}\right]

With the same nonsmooth law as for the first Interaction::

  // A list of names for the Interactions
  vector<string> id;
  id.resize(interactionNumber-1);
  CheckInsertInteraction checkInter;
  // A vector that will handle all the relations
  vector<Relation*> LLR(interactionNumber-1);
  //
  SiconosMatrix *H1 = new SimpleMatrix(1,2*nDof);
  if (dsNumber>1)
  {
      (*H1)(0,0) = -1.0;
      (*H1)(0,3) = 1.0;
      // Since Ri=Rj and h=0, we do not need to set b.
      Relation * relation = new LagrangianLinearR(*H1);
      for (i=1;(int)i<interactionNumber;i++)
      {
          // The systems handled by the current Interaction ...
	  dsConcerned.clear();
	  dsConcerned.insert(allDS.getDynamicalSystemPtr(i-1));
	  dsConcerned.insert(allDS.getDynamicalSystemPtr(i));
	  // The id: "i"
	  ostringstream ostr;
	  ostr << i;
	  id[i-1]= ostr.str();
	  // The relations
	  LLR[i-1] = new LagrangianLinearR(*relation); // the copy constructor builds all the relations
	  checkInter = allInteractions.insert( new Interaction(id[i-1], dsConcerned,i,interactionSize, nslaw0, LLR[i-1]));
      }
      delete relation;
  }

Note that each Relation belongs to one and only one Interaction (unlike a NonSmoothLaw); that is why we build a new Relation LLR[i-1] for each Interaction.

Everything is now ready to build the NonSmoothDynamicalSystem and the related Model::

    // --------------------------------
    // --- NonSmoothDynamicalSystem ---
    // --------------------------------
    NonSmoothDynamicalSystem * nsds = new NonSmoothDynamicalSystem(allDS, allInteractions);
    // -------------
    // --- Model ---
    // -------------
    Model * multiBeads = new Model(t0,T);
    multiBeads->setNonSmoothDynamicalSystemPtr(nsds); // set NonSmoothDynamicalSystem of this model

The Simulation
--------------

Time-Stepping scheme
""""""""""""""""""""

As a first example, we use Moreau's time-stepping scheme, where the nonsmooth problem is written as an LCP. The process is more or less the same as for the diode bridge, so we do not detail it. The only difference is that the OneStepIntegrator now handles several DynamicalSystems::

  string solverName = "Lemke";      // solver used for the nonsmooth problem
  Simulation* s = new TimeStepping(multiBeads);
  // -- Time discretisation --
  TimeDiscretisation * t = new TimeDiscretisation(h,s);
  // -- OneStepIntegrators --
  double theta = 0.5000001;
  OneStepIntegrator * OSI = new Moreau(allDS , theta ,s);
  // That means that all systems in allDS have the same theta value.
  // -- OneStepNsProblem --
  OneStepNSProblem * osnspb = new LCP(s,"LCP",solverName,10001, 0.001);

Event-Driven algorithm
""""""""""""""""""""""

In this second part, an event-driven algorithm is used to solve the problem. The principle of event-driven simulations is detailed in :ref:`event_driven`.

The dynamics are decomposed into "modes", time intervals where the dynamics are smooth, and discrete events, where the dynamics are nonsmooth.

Here, nonsmooth events are impacts between balls. Each time such an event is detected, a nonsmooth problem is formulated and solved (here as an LCP), while between events the systems are integrated with Lsodar, an ODE solver with root finding.

As for time-stepping, we first build the simulation and then its time discretisation::

  // The simulation belongs to Model multiBeads
  EventDriven* s = new EventDriven(multiBeads);
  TimeDiscretisation * t = new TimeDiscretisation(h,s);

The next step is the declaration of the integrators of the dynamical systems.
The integrator handles all the DynamicalSystems of the Model. During the integration, Lsodar searches for the roots of some equations (the constraints, i.e. the Interactions of the NonSmoothDynamicalSystem). The required OSI type is Lsodar, applied to allDS::

  OneStepIntegrator * OSI = new Lsodar(allDS,s);

Each time a root is found, a new NonSmoothEvent is created and a nonsmooth problem must be formulated and solved. We do not detail this here; just remember that this requires two LCPs, one at the "velocity" level, named impact, and another at the "acceleration" level, named acceleration.
The whole event-driven algorithm for Lagrangian systems is described in :ref:`event_driven_lagrange`::

  OneStepNSProblem * impact = new LCP(s, "impact",solverName,101, 0.0001,"max",0.6);
  OneStepNSProblem * acceleration = new LCP(s, "acceleration",solverName,101, 0.0001,"max",0.6);

The Model is now complete, and we can start the simulation process.

Simulation process
------------------

Time-Stepping
"""""""""""""

Once again, the process is the same as in the first tutorial and is not detailed.
For the output, we save the position and velocity of all the balls::

  s->initialize();
  int k = 0;
  int N = t->getNSteps(); // Number of time steps
  // Prepare the output and save the values at the initial time
  unsigned int outputSize = dsNumber*2+1;
  SimpleMatrix dataPlot(N+1,outputSize ); // Output data matrix
  // time
  dataPlot(k, 0) = multiBeads->getT0();
  // Positions and velocities
  i = 0; // Remember that the DS are sorted in increasing order of their number.
  DSIterator it;
  for(it = allDS.begin();it!=allDS.end();++it)
  {
      dataPlot(k,(int)i*2+1) = static_cast<LagrangianLinearTIDS*>(*it)->getQ()(0);
      dataPlot(k,(int)i*2+2) = static_cast<LagrangianLinearTIDS*>(*it)->getVelocity()(0);
      i++;
  }

Note that we use a "DSIterator", which is simply a pointer into a set of DynamicalSystems: allDS.begin() points to the first object handled by allDS and allDS.end() "just after" the last one. The system currently pointed to is \*it ("content of the pointer"). The loop above thus goes through all the DynamicalSystems and gets the corresponding :math:`q` and :math:`v`.
A static_cast is required, since allDS contains DynamicalSystem objects whereas we need functions specific to LagrangianDS (getQ...).

Next, we write::

  while(k < N)
  {
      k++;
      // solve ...
      s->computeOneStep();
      dataPlot(k, 0) = s->getNextTime();
      //
      i = 0;
      for(it = allDS.begin();it!=allDS.end();++it)
      {
          dataPlot(k,(int)i*2+1) = static_cast<LagrangianLinearTIDS*>(*it)->getQ()(0);
	  dataPlot(k,(int)i*2+2) = static_cast<LagrangianLinearTIDS*>(*it)->getVelocity()(0);
	  i++;
	  s->nextStep();
      }
  }

and to save the output file::

  ioMatrix io("result.dat", "ascii");
  io.write(dataPlot,"noDim");

Event-Driven
""""""""""""

An EventDriven simulation roughly integrates between events, stops at each event and gives it a special treatment. Events are handled by a specific object, the EventsManager, a kind of stack where events are saved in chronological order. It belongs to the Simulation object and can be accessed with::

  EventsManager * eventsManager = s->getEventsManagerPtr();

The manager is built during the initialization, which is still the first step of any simulation process::

  s->initialize();

Among other things, this initialization schedules the time events of the TimeDiscretisation object in the manager. Each time step is saved as a TimeDiscretisationEvent.

The simulation process then consists in:

* checking whether there is a "future" event;
* integrating the system until this future event is reached or until a nonsmooth event is found;
* scheduling the new event, if any;
* processing the system at the event (for example, for a nonsmooth event, formulating and solving one or more LCPs);
* going to the next step.

Once again, this is only a summary; see :ref:`event_driven` for more details on the event-driven strategy.

The resulting code is::

   // While there are events in the manager ...
    while(eventsManager->hasNextEvent())
      {
	eventDriven->computeOneStep();
      }

For the output, we first save the displacements and velocities at each time step::

    while(eventsManager->hasNextEvent())
      {
	k++;
	eventDriven->advanceToEvent();

	eventDriven->processEvents();
        // Positions and velocities for user time steps
	i = 0; // Remember that the DS are sorted in increasing order of their number.
	DSIterator it;
	dataPlot(k, 0) = eventDriven->getStartingTime();
	for(it = allDS.begin();it!=allDS.end();++it)
	  {
	    dataPlot(k,(int)i*2+1) = static_cast<LagrangianLinearTIDS*>(*it)->getQ()(0);
	    dataPlot(k,(int)i*2+2) = static_cast<LagrangianLinearTIDS*>(*it)->getVelocity()(0);
	    i++;
	  }
      }

But when a nonsmooth event occurs, it may be interesting to get the pre- and post-impact values.
In Siconos, the values saved in objects are usually the last computed ones, here the post-impact values.
The next-to-last values are saved in "memory" objects; we get them for a "NonSmoothEvent"::

    while(eventsManager->hasNextEvent())
      {
	k++;
	eventDriven->advanceToEvent();

	eventDriven->processEvents();
	if(eventsManager->getStartingEventPtr()->getType() == "NonSmoothEvent")
	  {
	    i = 0; // Remember that the DS are sorted in increasing order of their number.
	    DSIterator it;
	    dataPlot(k, 0) = eventDriven->getStartingTime();
	    for(it = allDS.begin();it!=allDS.end();++it)
	      {
		dataPlot(k,(int)i*2+1) = (*static_cast<LagrangianLinearTIDS*>(*it)->getQMemoryPtr()->getSiconosVector(1))(0);
		dataPlot(k,(int)i*2+2) = (*static_cast<LagrangianLinearTIDS*>(*it)->getVelocityMemoryPtr()->getSiconosVector(1))(0);
		i++;
	      }
	    k++;
	  }
        // Positions and velocities for user time steps
	i = 0; // Remember that the DS are sorted in increasing order of their number.
	DSIterator it;
	dataPlot(k, 0) = eventDriven->getStartingTime();
	for(it = allDS.begin();it!=allDS.end();++it)
	  {
	    dataPlot(k,(int)i*2+1) = static_cast<LagrangianLinearTIDS*>(*it)->getQ()(0);
	    dataPlot(k,(int)i*2+2) = static_cast<LagrangianLinearTIDS*>(*it)->getVelocity()(0);
	    i++;
	  }
      }

    // Output written in result.dat
    ioMatrix io("result.dat", "ascii");
    io.write(dataPlot,"noDim");

The simulation is now ready. The input file is completed with the required headers and delete instructions at the end.
The complete input files are:

* BeadsColumnTS.cpp for the time-stepping version;
* BeadsColumnED.cpp for the event-driven version;
* BeadsPlugin.cpp, the file that contains the external plugin.


Results
-------

You can now run in a terminal::

  siconos multiBeadsTS.cpp

and then plot the results, for example with gnuplot::

  gnuplot -persist result.gp

result.gp being a command file (see the example in mechanics/MultiBeadsColumn).

The results are shown in fig 2 below:

.. figure:: /figures/mechanics/MultiBeads/MultiBeads.*
   :align: center

   fig 2: Results of the MultiBeads simulation

.. highlight:: python
