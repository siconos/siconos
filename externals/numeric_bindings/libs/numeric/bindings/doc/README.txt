
Creating Documentation
----
To be able to create the documentation, please consult the quickbook documentation for the needed tools, and make sure you have a Boost repository copy, and updater the user-config in the Boost repository to reflect your local settings.

Then, on a posix system, from within the documentation directory, do a 

$ BOOST_ROOT=/path/to/full/boost/install bjam --v2 html

to generate the documentation. The correct resulting documentation will end up at doc/html/index.html.

