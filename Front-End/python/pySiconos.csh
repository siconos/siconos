#! /bin/csh

if ($?LD_LIBRARY_PATH) then
	setenv LD_LIBRARY_PATH $SICONOSPATH/python:$LD_LIBRARY_PATH
else
	setenv LD_LIBRARY_PATH $SICONOSPATH/python
endif
