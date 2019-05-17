#!/usr/bin/env python
# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
  
# Copyright 2018 INRIA.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

# http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
from __future__ import print_function

import siconos.numerics as N
import siconos.kernel as K

from siconos.kernel import \
    FrictionContact,\
    GlobalFrictionContact

from siconos.numerics import \
    FrictionContactProblem,\
    GlobalFrictionContactProblem

try:
    import siconos.fclib as F
    has_fclib=True
except :
    has_fclib=False
    print("No module named siconos.fclib.")

import imp
import sys
import time
import getopt
import random
import h5py
import os 

class FrictionContactTraceParams():
    def __init__(self, dump_itermax=10, dump_probability=None, fileName="filename", title = "title", description = "description", mathInfo= "mathInfo"):
        self._dump_itermax = dump_itermax
        self._dump_probability = dump_probability
        self._fileName = fileName
        self._title = title
        self._description = description
        self._mathInfo = mathInfo

    def display(self):
        print('title',self._title)
    



class FrictionContactTrace(FrictionContact):

    def __init__(self, dim, solver, params=None, nsds=None):
        if params == None:
            self._params = FrictionContactTraceParams()
        else:
            self._params = params
        proba= params._dump_probability
        maxiter =  params._dump_itermax
        if proba is not None:
            self._proba = 1. - proba
            self.condition = self.random_condition
        if maxiter is not None:
            self._maxiter = maxiter
            if proba is None:
                self.condition = self.maxiter_condition
            else:
                self.condition = self.random_and_maxiter_condition
        self._counter = 0
        self._stepcounter = 0
        self._nsds=nsds
        super(FrictionContactTrace, self).__init__(dim, solver)

    def maxiter_condition(self, SO):
        return SO.iparam[N.SICONOS_IPARAM_ITER_DONE] >= self._maxiter

    def random_condition(self, SO):
        return random.random() > self._proba

    def random_and_maxiter_condition(self, SO):
        return self.maxiter_condition(SO) and self.random_condition(SO)

    def compute(self,time):
        info = 0
        self.setKeepLambdaAndYState(False)
        cont = self.preCompute(time)

        if (not cont):
            return 0

        if (self.indexSetLevel() == 999):
            return 0

        self.updateMu()

        if self.getSizeOutput() != 0:

            #            M = BlockCSRMatrix()
            #M.fillM(model.nonSmoothDynamicalSystem().topology().indexSet(1))
            #M.convert()

            #            H = BlockCSRMatrix()


            #t = GlobalFrictionContactProblem()

            #t.M = M.getNumericsMatSparse()

            w_backup = self.w().copy()
            z_backup = self.z().copy()
            SO = self.numericsSolverOptions()
            fclib_written = False
            if self.condition(SO) and has_fclib: 
           
                # problem = self.getNumericsProblemPtr()
                # print(problem, type(problem))
                    
                problem = self.frictionContactProblemPtr()
                #print(problem, type(problem))
                solver_maxiter=SO.iparam[0]
                n_format_string=len(str(solver_maxiter))
                format_string = "{0}-i{1:0"+str(n_format_string)+"d}-{2}-{3}.hdf5"
                filename = format_string.format(self._params._fileName,
                                                              SO.iparam[N.SICONOS_IPARAM_ITER_DONE],
                                                              problem.numberOfContacts,
                                                              self._counter)

                print('filename =', filename)
                if os.path.exists(filename):
                    os.remove(filename)
                    print('WARNING: file '+filename+ ' was existing and has been replaced')
                
                self._counter += 1
                N.frictionContact_fclib_write(problem,
                                              self._params._title,
                                              self._params._description,
                                              self._params._mathInfo,
                                              filename,
                                              -1)
                guess = F.fclib_solution()
                guess.u = w_backup
                guess.r = z_backup
                F.fclib_write_guesses(1, guess, filename)

                with h5py.File(filename, 'r+') as fclib_file:
                    attrs = fclib_file['fclib_local']['info'].attrs
                    attrs.create('numberOfInvolvedDS',
                                 self._nsds.topology().numberOfInvolvedDS(1))

                fclib_written =True

                    
            info = self.solve()

            if fclib_written:
                
                solution = F.fclib_solution()
                solution.u = self.w()
                solution.z = self.z()
                F.fclib_write_solution(solution, filename)



                
            self.postCompute()

        return info
    
class GlobalFrictionContactTraceParams():
    def __init__(self, dump_itermax=10, dump_probability=None, fileName="filename", title = "title", description = "description", mathInfo= "mathInfo"):
        self._dump_itermax = dump_itermax
        self._dump_probability = dump_probability
        self._fileName = fileName
        self._title = title
        self._description = description
        self._mathInfo = mathInfo

    def display(self):
        print('title',self._title)
    



class GlobalFrictionContactTrace(GlobalFrictionContact):

    def __init__(self, dim, solver, params=None, nsds=None):
        if params == None:
            self._params = GlobalFrictionContactTraceParams()
        else:
            self._params = params
        proba= params._dump_probability
        maxiter =  params._dump_itermax
        if proba is not None:
            self._proba = 1. - proba
            self.condition = self.random_condition
        if maxiter is not None:
            self._maxiter = maxiter
            if proba is None:
                self.condition = self.maxiter_condition
            else:
                self.condition = self.random_and_maxiter_condition
        self._counter = 0
        self._stepcounter = 0
        self._nsds=nsds
        super(GlobalFrictionContactTrace, self).__init__(dim, solver)

    def maxiter_condition(self, SO):
        return SO.iparam[N.SICONOS_IPARAM_ITER_DONE] >= self._maxiter

    def random_condition(self, SO):
        return random.random() > self._proba

    def random_and_maxiter_condition(self, SO):
        return self.maxiter_condition(SO) and self.random_condition(SO)

    def compute(self,time):
        info = 0

        cont = self.preCompute(time)

        if (not cont):
            return 0

        if (self.indexSetLevel() == 999):
            return 0

        self.updateMu()

        w_backup = self.w().copy()
        z_backup = self.z().copy()
        SO = self.numericsSolverOptions()

        info = self.solve()
        problem = self.globalFrictionContactProblemPtr()
        if problem.numberOfContacts >0 :
            if self.condition(SO) and has_fclib:
                # problem = self.getNumericsProblemPtr()
                # print(problem, type(problem))

                

                solver_maxiter=SO.iparam[0]
                n_format_string=len(str(solver_maxiter))
                format_string = "{0}-i{1:0"+str(n_format_string)+"d}-{2}-{3}.hdf5"
                filename = format_string.format(self._params._fileName,
                                                SO.iparam[N.SICONOS_IPARAM_ITER_DONE],
                                                problem.numberOfContacts,
                                                self._counter)

                if os.path.exists(filename):
                    os.remove(filename)
                    print('WARNING: file '+filename+ ' was existing and has been replaced')

                self._counter += 1
                N.globalFrictionContact_fclib_write(problem,
                                                    self._params._title,
                                                    self._params._description,
                                                    self._params._mathInfo,
                                                    filename)
                guess = F.fclib_solution()
                guess.u = w_backup
                guess.r = z_backup
                F.fclib_write_guesses(1, guess, filename)




                solution = F.fclib_solution()
                solution.u = self.w()
                solution.z = self.z()
                F.fclib_write_solution(solution, filename)
            
                
        self.postCompute()

        return info
