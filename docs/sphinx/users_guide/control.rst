.. SPDX-FileCopyrightText: INRIA
..
.. SPDX-License-Identifier: Apache-2.0

.. _control_toolbox:

Control Toolbox
===============

THIS PAGE IS OUTDATED AND MUST BE REVIEWED

Control Manager
---------------

Rules:

* define a control manager linked to an EXISTING model::

    auto cm = std::make_shared<siconos::control::ControlManager>(myModel);

* add Sensors and Actuators to this manager::

    auto t1 = std::make_shared<siconos::simulation::TimeDiscretisation>(t0,h);
    auto s1 = cm->addSensor(typeS1,t1);
    auto t2 = std::make_shared<siconos::simulation::TimeDiscretisation>(t0,h);
    auto a1 = cm->addActuator(typeA1,t2);
    // ...

typeS1 and typeA1 are integers that give the type of the Sensor/Actuator. See the corresponding
sections for details and the available types.

Important: each actuator/sensor must have its own TimeDiscretisation object,
because its TimeDiscretisation is incremented (tk -> tk+1) each time one of its events is processed:
if several events shared the same TimeDiscretisation, it would be incremented too many times.

* initialize the manager::

    cm->initialize();

This schedules the events of each Sensor/Actuator in the EventsManager of the simulation.
It must be called BEFORE simulation->initialize().

It is also possible to insert a new Sensor or Actuator at any time during the simulation::

   cm->addAndRecordSensor(typeS1,t3);
   cm->addAndRecordActuator(typeA1,t4);

Sensors
-------

* tk
* capture
* map of vectors. Save the values for all events?

Actuators
---------

Linked to one DS.

* tk
* a way to compute the value of z
* setZ in DS
