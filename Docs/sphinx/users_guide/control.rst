.. _control_toolbox::

Control Toolbox
===============

THIS PAGE IS OUTDATED AND MUST BE REVIEWED

Control Manager
---------------

Rules:

* define a control manager linked to an EXISTING model::

  SP::ControlManager cm(new ControlManager(myModel));

* add Sensors and Actuators to this manager::

    SP::TimeDiscretisation t1(t0,h);	
    SP::Sensor s1 = cm->addSensor(typeS1,t1);
    SP::TimeDiscretisation t2(t0,h);	
    SP::Actuator a1 = cm->addActuator(typeA1,t2);
    // ... 

typeS1 and typeA1 are integers which represent the type of the Sensor/Actuator. See corresponding
sections for details and various types. \n
Important: each actuator/sensor must have its own TimeDiscretisation object. \n
Why: after each process of an event, its TimeDiscretisation is increment (ie tk->tk+1 etc )
and if several events share the same TimeDiscretisation, it will be incremented too many times.

* initialize the manager::

    cm->initialize();

This result in the scheduling of events corresponding to each Sensor/Actuator into the EventsManager of the simulation.
It must be called BEFORE simulation->initialize()

It is also possible to insert a new Sensor or Actuator at any time during the simulation::

   cm->addAndRecordSensor(typeS1,t3);
   cm->addAndRecordActuator(typeA1,t4);

Sensors
-------

* tk
* capture
* map of vectors. Save values for all events? 

Actuators
---------

link to one DS. 

* tk
* a way to compute the value of z
* setZ in DS
