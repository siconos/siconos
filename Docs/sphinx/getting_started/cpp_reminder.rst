.. _cpp_reminder::


C++ Reminder
============

C++ is the basic language of siconos input files, thus you need at least a basic knowledge of c++.
This topic is covered by many books or tutorials, try for example :cite:`Eckel_cpp`
However, you can find in this first section the main C++ commands you need to know to write a first input file. \n


However, this page presents some basic commands and things you need to know to write a first cpp driver for yout simulation.

Building/destroying and using objects
-------------------------------------

.. highlight:: c++
	       
To describe your problem, you will mainly need to build pointers to objects. 
Siconos use smart pointers so you do not have to bother with deletion.
The namespace SP enclose typedefs of all Siconos classes as smart pointers.

So to build an object of class *CHILD* derived from base class *BASE*, the syntax will be::

  SP::BASE my_object_name(new CHILD(some args ...)

For example, suppose you want to build a LagrangianDS, which belongs
to the base class DynamicalSystem, then do::

  SP::DynamicalSystem nameOfMyDS(new LagrangianDS(someParameters));

*new* call will reserve memory and build the object.

For each object, different types of constructors may be available. For an
exhaustive list see :ref:`siconos_api_reference`.

The access to the methods of the objects is done thanks to "*" or "->"::
  
  nameOfMyDS->getType(); // return the type of the DynamicalSystem
  (*nameOfMyDS).getType(); // Same thing. 

