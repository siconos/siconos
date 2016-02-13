// -*- C++ -*-
// Siconos, Copyright INRIA 2005-2016
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.
// Siconos is a free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
// Siconos is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Siconos; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
//
// Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
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
%newobject duplicateNumericsMatrix;
%newobject newNumericsMatrix;
%newobject createNumericsMatrix;
%newobject createNumericsMatrixFromData;
%newobject newSparseNumericsMatrix;
%newobject newSBM;
%newobject SBCMToSBM;
%newobject newSparseBlockCoordinateMatrix3x3fortran;

// strange things
%newobject create_NMS_data;
%newobject SN_lumod_dense_allocate;

