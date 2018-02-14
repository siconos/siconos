"""A few tests for classes and functions from kernel/modelingtools

"""
#!/usr/bin/env python

import numpy as np
import siconos.kernel as K

t0 = 0
T = 10
r = 0.1
g = 9.81
m = 1
e = 0.9
theta = 0.5
h = 0.005


q = np.array([1, 0, 0])
v = np.array([0, 0, 0])
mass = np.eye(3)
mass[2, 2] = 3. / 5 * r * r

weight = np.array([-m * g, 0, 0])
tol = np.finfo(np.double).eps


def test_LagrangianLinearTIDS():
    ball = K.LagrangianLinearTIDS(q, v, mass)
    assert np.allclose(ball.q(), q, rtol=tol, atol=tol)
    assert np.allclose(ball.velocity(), v, rtol=tol, atol=tol)
    assert np.allclose(ball.mass(), mass, rtol=tol, atol=tol)
    ball.setFExtPtr(weight)
    assert np.allclose(ball.fExt(), weight, rtol=tol, atol=tol)


def test_NewtonImpactNSL():
    nslaw = K.NewtonImpactNSL(e)
    assert nslaw.e() == e


def test_LagrangianLinearTIR():
    H = np.array([[1, 0, 0]])
    b = np.zeros(1)
    relation = K.LagrangianLinearTIR(H, b)
    assert np.allclose(relation.jachq(), H, rtol=tol, atol=tol)


def test_Nsds():
    bouncing_ball = K.NonSmoothDynamicalSystem(t0, T)
    assert bouncing_ball.t0() == t0

# First order interaction

A = [[1, 2], [3 , 4]]

C = [[1., 2.], [3, 4.], [5., 6.], [7., 8.]]

D = [[1, 1, -1., 0.],
     [1. , 1. , 0., -1.],
     [1., 0., 0., 0.],
     [0., 1., 0., 0.]]

B = [[0., 0., 3, 4],
     [0., 0., 0., 0.]]


def create_first_order_linear_ds_order0():
    #
    # dynamical system
    #
    init_state = [1, 2]
    ds = K.FirstOrderLinearDS(init_state, A)
    return ds

relation_first_order_classes = {
    K.FirstOrderLinearTIR:(C, B)}

def test_firstorder_interaction_interface_order0():
    """Tests methods that should be available
        for all lagrangian_interactions
    """

    ds = create_first_order_linear_ds_order0()

    for args in relation_first_order_classes:
        print('----- ')
        print('Test class ', args)
        class_name = args
        attr = relation_first_order_classes[args]
        relation  = class_name(*attr)
        nslaw = K.ComplementarityConditionNSL(4)
        inter = K.Interaction(nslaw, relation)
        #
        # NSDS
        #
        t0 = 0.      # start time
        tend = 10.   # end time
        nsds = K.NonSmoothDynamicalSystem(t0, tend)

        # add the dynamical system to the non smooth dynamical system
        nsds.insertDynamicalSystem(ds)

        # link the interaction and the dynamical system
        nsds.link(inter, ds)

        # test interface

        ds.x()
        ds.r()

        inter.computeOutput(t0,0)
        inter.computeInput(t0,0)

        print('y(1)[0] = ',inter.y(0)[0])
        print('lambda(0) = ', inter.lambda_(0))
        assert(inter.y(0)[0] == 5.0)
        assert(inter.lambda_(0)[0] == 0.0)

        try:
            print(ds.x(1))
        except Exception as err:
            print(err)

        try:
            inter.computeOutput(t0,1)
            print('y(1) = ', inter.y(1))
        except Exception as err:
            print(err)

        try:
            inter.computeInput(t0,1)
            print('lambda(1) = ', inter.lambda_(1))
        except Exception as err:
            print(err)

def create_lagrangian_linear_tids_order1():
    #
    # dynamical system
    #
    x = np.zeros(3, dtype=np.float64)
    x[0] = 1.
    x[1] = 2.
    x[2] = 3.
    v = np.zeros_like(x)
    v[0] = 4.
    v[1] = 5.
    v[2] = 6.
    # mass matrix
    mass = np.eye(3, dtype=np.float64)
    mass[2, 2] = 3. / 5 * r * r

    # the dynamical system
    ds = K.LagrangianLinearTIDS(x, v, mass)

    # set external forces
    weight = np.zeros_like(x)
    weight[0] = -m * g
    ds.setFExtPtr(weight)
    return ds

def create_lagrangian_linear_tids_order2():
    #
    # dynamical system
    #
    x = np.zeros(3, dtype=np.float64)
    x[0] = 1.
    x[1] = 2.
    x[2] = 3.
    v = np.zeros_like(x)
    v[0] = 4.
    v[1] = 5.
    v[2] = 6.
    # mass matrix
    mass = np.eye(3, dtype=np.float64)
    mass[2, 2] = 3. / 5 * r * r

    # the dynamical system
    ds = K.LagrangianLinearTIDS(x, v, mass, 2)

    # set external forces
    weight = np.zeros_like(x)
    weight[0] = -m * g
    ds.setFExtPtr(weight)
    return ds


H = np.zeros((1, 3), dtype=np.float64)
H[0, 0] = 1.
H[0, 1] = 2.
H[0, 2] = 3.
D = np.zeros((1, 1), dtype=np.float64)
D[0, 0] = 1.


relation_lag_classes = {
    K.LagrangianLinearTIR: (H,),
    K.LagrangianCompliantLinearTIR: (H,D)}


def test_lagrangian_interaction_interface_order1():
    """Tests methods that should be available
        for all lagrangian_interactions
    """

    ds = create_lagrangian_linear_tids_order1()

    for args in relation_lag_classes:
        print('----- ')
        print('Test class ', args)
        nslaw = K.NewtonImpactNSL(e)

        class_name = args
        attr = relation_lag_classes[args]
        relation  = class_name(*attr)
        #print(class_name, attr[0])


        #relation = K.LagrangianLinearTIR(H)
        inter = K.Interaction(nslaw, relation)

        #
        # NSDS
        #
        t0 = 0.      # start time
        tend = 10.   # end time
        nsds = K.NonSmoothDynamicalSystem(t0, tend)

        # add the dynamical system to the non smooth dynamical system
        nsds.insertDynamicalSystem(ds)

        # link the interaction and the dynamical system
        nsds.link(inter, ds)

        # test interface
        inter.computeOutput(t0,0)

        inter.computeInput(t0,1)

        print('y(1)[0] = ',inter.y(0)[0])
        print('lambda(1) = ', inter.lambda_(1))

        assert(inter.y(0)[0] == 14.0)
        assert(inter.lambda_(1)[0] == 0.0)

        inter.computeOutput(t0,1)

        print('y(1)[0] = ', inter.y(1)[0])
        assert(inter.y(1)[0] == 32.0)

        try:
            print('acceleration', ds.acceleration())
        except Exception as err:
            print(err)

        try:
            inter.computeInput(t0,0)
            print('lambda(0) = ', inter.lambda_(0))
        except Exception as err:
            print(err)

        try:
            inter.computeInput(t0,2)
            print('lambda(2) = ', inter.lambda_(2))
        except Exception as err:
            print(err)

        try:
            inter.computeOutput(t0,2)
              #print('y(2)[0] = ', inter.y(2))
        except Exception as err:
            print(err)

def test_lagrangian_interaction_interface_order2():
    """Tests methods that should be available
        for all lagrangian_interactions
    """

    ds = create_lagrangian_linear_tids_order2()

    for args in relation_lag_classes:
        print('----- ')
        print('Test class ', args)
        nslaw = K.NewtonImpactNSL(e)

        class_name = args
        attr = relation_lag_classes[args]
        relation  = class_name(*attr)
        #print(class_name, attr[0])


        #relation = K.LagrangianLinearTIR(H)
        inter = K.Interaction(nslaw, relation, 0, 2)

        #
        # NSDS
        #
        t0 = 0.      # start time
        tend = 10.   # end time
        nsds = K.NonSmoothDynamicalSystem(t0, tend)

        # add the dynamical system to the non smooth dynamical system
        nsds.insertDynamicalSystem(ds)

        # link the interaction and the dynamical system
        nsds.link(inter, ds)

        # test interface
        inter.computeOutput(t0,0)

        inter.computeInput(t0,1)

        print('y(1)[0] = ',inter.y(0)[0])
        print('lambda(1) = ', inter.lambda_(1))

        assert(inter.y(0)[0] == 14.0)
        assert(inter.lambda_(1)[0] == 0.0)

        inter.computeOutput(t0,1)

        print('y(1)[0] = ', inter.y(1)[0])
        assert(inter.y(1)[0] == 32.0)


        print('acceleration', ds.acceleration())
        inter.computeInput(t0,0)
        print('lambda(0) = ', inter.lambda_(0))
        inter.computeInput(t0,2)
        print('lambda(2) = ', inter.lambda_(2))

        inter.computeOutput(t0,2)
        print('y(2)[0] = ', inter.y(2))

def create_newton_euler_order1():
    #
    # dynamical system
    #
    x = np.zeros(7, dtype=np.float64)
    x[0] = 1.
    x[1] = 2.
    x[2] = 3.
    x[3] = 4.
    x[4] = 5.
    x[5] = 6.
    x[6] = 7.
    v = np.zeros(6)
    v[0] = 4.
    v[1] = 5.
    v[2] = 6.
    # mass matrix
    mass = 1.0
    inertia = np.eye(3, dtype=np.float64)

    # the dynamical system
    ds = K.NewtonEulerDS(x, v, mass, inertia)

    return ds
H = np.zeros((1, 7), dtype=np.float64)
H[0, 0] = 1.
H[0, 1] = 2.
H[0, 2] = 3.
H[0, 3] = 4.
H[0, 4] = 5.
H[0, 5] = 6.
H[0, 6] = 7.

def test_newton_interaction_interface_order1():
    """Tests methods that should be available
        for all lagrangian_interactions
    """

    ds = create_newton_euler_order1()
    print('----- ')

    nslaw = K.NewtonImpactNSL(e)

    relation  = K.NewtonEulerR()
    relation.setJachq(H)
    inter = K.Interaction(nslaw, relation)

    #
    # NSDS
    #
    t0 = 0.      # start time
    tend = 10.   # end time
    nsds = K.NonSmoothDynamicalSystem(t0, tend)

    # add the dynamical system to the non smooth dynamical system
    nsds.insertDynamicalSystem(ds)

    # link the interaction and the dynamical system
    nsds.link(inter, ds)


    ds.q()
    ds.twist()

    # test interface
    inter.computeOutput(t0,0)

    inter.computeInput(t0,1)

    inter.display()
    assert(inter.y(0)[0] == 140.0)
    assert(inter.lambda_(1)[0] == 0.0)

    try:
        print('acceleration', ds.acceleration())
    except Exception as err:
        print(err)

    try:
        inter.computeInput(t0,0)
        print('lambda(0) = ', inter.lambda_(0))
    except Exception as err:
        print(err)


if __name__ == '__main__':
    # print('test_lagrangian_interaction_interface')
    # test_lagrangian_interaction_interface_order1()
    # print('test_lagrangian_interaction_interface')
    # test_lagrangian_interaction_interface_order2()
    print('test_firstorder_interaction_interface_order0')
    test_firstorder_interaction_interface_order0()
    print('test_newton_interaction_interface_order1')
    test_newton_interaction_interface_order1()
