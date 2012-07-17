#!/usr/bin/env python

# Siconos-sample, Copyright INRIA 2005-2012.
# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
# Siconos is a free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# Siconos is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Siconos; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#
# Contact: Vincent ACARY, siconos-team@lists.gforge.fr
#

import numpy as np
import numpy.linalg as LA
import time
#import Siconos.Kernel as SK
import math
import linear_system
import numpy.fft as fft
import tools
LinearSystem = linear_system.LinearSystem


sqrt = math.sqrt
pi = math.pi

######################################
## Geometry and Physical parameters ##
######################################
# Young modulus
E = 2.1e11
# Poisson
nu = 0.3
# mass per unit vol
rho = 7800
# Length and cross section dim
Length = 1.0
height = 0.03
width = 0.03
section = height*width
g = 9.81
########################
## FEM discretisation ##
########################

ndof = 3
eltLength = Length/ndof

# initial conditions
initialPos = 0.001
initialVel = -1.0
q0 = initialPos*np.ones((ndof))
velocity0= initialVel*np.ones((ndof))
#beam = SK.LagrangianLinearTIDS(q0,velocity0,Mass)
#beam.setKPtr(Stiffness)

# No weight ...

# Mass matrix
coef = rho*section*eltLength/6.0
diag = coef*4.*np.ones((ndof))
diag[0] = diag[-1] = diag[0]/2.
Mass = np.diag(diag)
Mass.flat[1:ndof*(ndof-1):ndof+1] =  coef*np.ones((ndof-1))
Mass.flat[ndof:ndof*ndof:ndof+1] =  coef*np.ones((ndof-1))
# Stiffness matrix
coef = E*section/eltLength
diag = 2.*coef*np.ones((ndof))
diag[0] = diag[-1] = diag[0]/2.
Stiffness = np.diag(diag)
Stiffness.flat[1:ndof*(ndof-1):ndof+1] =  -coef*np.ones((ndof-1))
Stiffness.flat[ndof:ndof*ndof:ndof+1] =  -coef*np.ones((ndof-1))

# Enforce clamped BC
#Stiffness[0,0] = Stiffness[0,0]*1e3

### Boundary conditions
# number of dof where bc are applied
nc = 1
Bc = np.zeros((nc,ndof))
Fbc = np.zeros((nc))
Bc[0,0] = 1

##################
## Interactions ##
##################

nbContacts = 1

H = np.zeros((nbContacts,ndof))
H[0,-1] = 1
b = np.zeros((nbContacts))
b[:] = 0.0
#relation = SK.LagrangianLinearTIR(H,b)
## Non-Smooth Law
e = 0.0
#nslaw = SK.NewtonImpactNSL(e)

#####

system = LinearSystem(Mass,Stiffness,Bc,Fbc)

omega = 1297.5*2.*pi
# sampling
Nfft = 1

bigH = tools.computeBigH(H,Nfft)

nRowMCP = (ndof+nbContacts+nc)*Nfft
nColMCP = nRowMCP
MCPmat = np.zeros((nRowMCP,nColMCP),dtype='complex128')

step=ndof*Nfft
step2=step+nc*Nfft
for i in range(Nfft):
    G = system.computeG(i,omega)
    MCPmat[i*ndof:(i+1)*ndof,i*ndof:(i+1)*ndof] = G
    MCPmat[step+i*nc:step+(i+1)*nc,i*ndof:(i+1)*ndof] = system.Bc
    MCPmat[i*ndof:(i+1)*ndof,step+i*nc:step+(i+1)*nc] = (system.Bc).T

MCPmat[step2:-1,0:step] = bigH
MCPmat[0:step,step2:-1] = bigH.T

print MCPmat    
print MCPmat.shape

bMCP = np.zeros((nRowMCP),dtype='complex128')
# clamped beam + b = 0 in relation ==> bMCP = 0


Ufft = np.zeros((ndof*Nfft),dtype='complex64')
gamma = 12
energyNormalisation(system.mass,U,gamma,Nfft)
