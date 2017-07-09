// -*- C++ -*-
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.
//
// Copyright 2016 INRIA.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
//

// List here the functions that return a newly allocated object so that swig knows that it has to delete them

// FCLIB stuff
%newobject frictionContact_fclib_read;
%newobject globalFrictionContact_fclib_read;
%newobject from_fclib_local;
%newobject from_fclib_global;

// Problems
%newobject variationalInequality_new;
%newobject frictionContactProblem_new;
%newobject genericMechanical_newFromFile;
%newobject buildEmptyGenericMechanicalProblem;

// NumericsMatrix
%newobject newNumericsSparseLinearSolverParams;
%newobject newNumericsSparseMatrix;
%newobject NM_duplicate;
%newobject NM_new;
%newobject NM_create;
%newobject NM_create_from_data;
%newobject NM_new_SBM;
%newobject SBM_new;
%newobject SBCM_to_SBM;
%newobject  SBCM_new_3x3;

// strange things
%newobject create_NMS_data;
%newobject SN_lumod_dense_allocate;

