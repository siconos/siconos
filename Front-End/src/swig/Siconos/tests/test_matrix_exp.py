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

from Siconos.Kernel import *
import numpy as np
from scipy.linalg import expm


def compute_phi_psi(A, B, h, TV=False):
    # variable declaration
    t0 = 0.0  # start time
    T = 1  # end time
    n, m = B.shape

     # Matrix declaration
    x0 = np.random.random(n)
    Csurface = np.random.random((m, n))

    # Declaration of the Dynamical System
    if TV:
        processDS = FirstOrderLinearDS(x0, A)
    else:
        processDS = FirstOrderLinearTIDS(x0, A)
    # Model
    process = Model(t0, T)
    process.nonSmoothDynamicalSystem().insertDynamicalSystem(processDS)
    # time discretisation
    processTD = TimeDiscretisation(t0, h)
    # Creation of the Simulation
    processSimulation = TimeStepping(processTD, 0)
    processSimulation.setName("plant simulation")
    # Declaration of the integrator
    processIntegrator = ZeroOrderHold(processDS)
    processSimulation.insertIntegrator(processIntegrator)

    rel = FirstOrderLinearTIR(Csurface, B)
    idi = "interaction for control"
    nslaw = RelayNSL(m)
    inter = Interaction(idi, processDS, m, m, nslaw, rel)

    process.nonSmoothDynamicalSystem().insertInteraction(inter, True)
    #process.nonSmoothDynamicalSystem().link(inter, processDS)
    # Initialization
    process.initialize(processSimulation)

    # Main loop
    processSimulation.computeOneStep()
    phi = processIntegrator.Phi(processDS)
    psi = processIntegrator.Psi(processDS)

    return (phi, psi)


def test_TI():
    h = 1.0e-4
    k = 0
    err = False
    while k < 100:
        n = np.random.randint(2, 100)
        m = np.random.randint(1, n)
        A = np.random.random((n, n))
        B = np.random.random((n, m))
        errInv = 1
        try:
            Ainv = np.linalg.inv(A)
            Aexp = expm(A*h)
            Psi = Ainv.dot((Aexp-np.eye(n)).dot(B))
            k += 1
            errInv = np.linalg.norm(A.dot(np.linalg.inv(A))-np.eye(n))
        except:
            pass
        if errInv < 1.0e-12:
            try:
                (AexpS, PsiS) = compute_phi_psi(A, B, h)
            except:
                print(n, m)
            errPhi = np.linalg.norm((Aexp - AexpS))
            errPsi = np.linalg.norm(Psi-PsiS)
            if errPhi > 5.0e-12 or errPsi > 5.0e-12:
                print(errPhi, errPsi)
                err = True

    assert err == False


def test_SISO():
    h = 1.0e-4
    k = 0
    err = False
    while k < 100:
        n = np.random.randint(2, 100)
        m = 1
        A = np.random.random((n, n))
        B = np.random.random((n, m))
        errInv = 1
        try:
            Ainv = np.linalg.inv(A)
            Aexp = expm(A*h)
            Psi = Ainv.dot((Aexp-np.eye(n)).dot(B))
            k += 1
            errInv = np.linalg.norm(A.dot(np.linalg.inv(A))-np.eye(n))
        except:
            pass
        if errInv < 1.0e-12:
            try:
                (AexpS, PsiS) = compute_phi_psi(A, B, h)
            except:
                print(n, m)
            errPhi = np.linalg.norm((Aexp - AexpS))
            errPsi = np.linalg.norm(Psi-PsiS)
            if errPhi > 5.0e-12 or errPsi > 5.0e-12:
                print(errPhi, errPsi)
                err = True
    assert err == False



def test_TV():
    h = 1.0e-4
    k = 0
    err = False
    while k < 10:
        n = np.random.randint(2, 100)
        m = np.random.randint(1, n)
        A = np.random.random((n, n))
        B = np.random.random((n, m))
        errInv = 1
        try:
            Ainv = np.linalg.inv(A)
            Aexp = expm(A*h)
            Psi = Ainv.dot((Aexp-np.eye(n)).dot(B))
            k += 1
            errInv = np.linalg.norm(A.dot(np.linalg.inv(A))-np.eye(n))
        except:
            pass
        if errInv < 1.0e-12:
            try:
                (AexpS, PsiS) = compute_phi_psi(A, B, h, TV=True)
            except:
                print(n, m)
            errPhi = np.linalg.norm((Aexp - AexpS))
            errPsi = np.linalg.norm(Psi-PsiS)
            if errPhi > 5.0e-12 or errPsi > 5.0e-12:
                print(errPhi, errPsi)
                err = True
    assert err == False

