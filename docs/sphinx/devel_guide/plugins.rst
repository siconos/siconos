.. _siconos_devel_plugins:


About plugins in classes
------------------------

An attempt to write some rules when implementing a plugin for a class attribute ...


Rules
^^^^^
Consider a class with an attribute '_xx', for which you want to propose a plugin-mechanism.


Attributes to add to the class:

- _pluginxx, a :class:`PluggedObject`,
  
Methods to add to the class:

- getPluginXx() returns the _pluginxx object
- setComputeXxFunction(args), args=(path,name) or (fptr) to connect fPtr of the PluggedObject to a user-defined function 'name' in file 'path'.
- computeXx(...) call the plugged function and update _xx content.

Others

- _zeroPlugin() -->  set to 'zero' all PluggedObject of the class
- updatePlugins(time) --> call all plugged functions of the class, for time and current state

Plugged functions can be set either with the construtor or with the setComputeXxFunction method.


Example/template
^^^^^^^^^^^^^^^^

computation of :math:`f(x,t,z)` in :class:`FirstOrderNonLinearDS`.

Attributes and methods

- :class:`FirstOrderNonLinearDS::_f`
- :class:`FirstOrderNonLinearDS::_pluginf`
- :class:`FirstOrderNonLinearDS::computef()`
- :class:`FirstOrderNonLinearDS::setComputeFFunction`
- :class:`FirstOrderNonLinearDS::_zeroPlugin`
- :class:`FirstOrderNonLinearDS::updatePlugins`

  
