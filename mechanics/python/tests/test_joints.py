# SPDX-FileCopyrightText: INRIA
#
# SPDX-License-Identifier: Apache-2.0
#
# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.

import numpy as np
import siconos.mechanics.joints as mj
import siconos.modeling as sm


def test_joints():
    joints = []
    joints.append(mj.CylindricalJointR())
    joints.append(mj.FixedJointR())
    joints.append(mj.KneeJointR())
    joints.append(mj.PivotJointR())
    joints.append(mj.PrismaticJointR())

    point0 = np.random.random(3)
    axis0 = np.random.random(3)
    pos = np.random.random(7)

    # Cylindrical joints
    joints[0].setPoint(0, point0)
    joints[0].setAxis(0, axis0)
    joints[0].setBasePositions(pos, None)

    # Fixed
    joints[1].setBasePositions(pos, None)

    # Knee
    joints[2].setPoint(0, point0)
    joints[2].setBasePositions(pos, None)

    # Pivot
    joints[3].setPoint(0, point0)
    joints[3].setAxis(0, axis0)
    joints[3].setBasePositions(pos, None)

    # Prismatic
    joints[4].setAxis(0, axis0)
    joints[4].setBasePositions(pos, None)

    # Coupler
    joints.append(mj.PrismaticJointR())
    joints.append(mj.CylindricalJointR())
    point1 = np.random.random(3)
    axis1 = np.random.random(3)
    pos1 = np.random.random(7)

    joints[5].setAxis(0, axis1)
    joints[5].setBasePositions(pos1, None)
    joints[6].setPoint(0, point1)
    joints[6].setAxis(0, axis1)
    joints[6].setBasePositions(pos1, None)

    q0 = np.ones(7, dtype=np.float64)
    twist0 = np.ones(6, dtype=np.float64)
    mass = 2.0
    inertia = np.zeros((3, 3), dtype=np.float64, order="F")
    inertia[0, 0] = inertia[1, 1] = inertia[2, 2] = 1.0
    neds1 = sm.NewtonEulerDS(q0, twist0, mass, inertia, sm.alias_t)
    neds2 = sm.NewtonEulerDS(q0, twist0, mass, inertia, sm.alias_t)
    cj = mj.CouplerJointR(joints[0], 0, joints[6], 1, 0.5, neds1, 1, neds2, 1)
    ck = mj.CouplerJointR(joints[4], 0, joints[5], 1, 0.5)
    cl = mj.CouplerJointR(joints[0], 0, joints[6], 1, 0.5, neds1)

    q1 = neds1.q()
    q2 = None
    ck.setBasePositions(q1, q2)
    cj.setBasePositions(neds1.q(), neds2.q())
    cl.setBasePositions(neds1.q(), q2)
