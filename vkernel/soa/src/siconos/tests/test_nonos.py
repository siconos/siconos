from __future__ import annotations
import numpy as np
import nonos as m


def test_version():
    assert m.disks.__version__ == "dev"


def test_make():
    data = m.disks.make_storage()


def test_add():
    data = m.disks.make_storage()
    disk = m.disks.add_body(data)


def test_q():
    import numpy

    data = m.disks.make_storage()
    disk = m.disks.add_body(data)
    disk.set_q(np.array([1, 2, 3], dtype=np.float64))
    q = disk.q()
    assert list(q) == [1, 2, 3]

    q[0] = 4
    assert list(disk.q()) == [4, 2, 3]


def test_nslaw():
    data = m.disks.make_storage()
    nslaw = m.disks.add_nslaw(data)
    nslaw.set_e(0.7)
    assert nslaw.e() == 0.7


def test_disk_shape():
    # the shape is the body's: diskdisk_r binds any pair of bodies
    data = m.disks.make_storage()
    shape = m.disks.add_disk_shape(data)
    shape.set_radius(2.0)
    body = m.disks.add_body(data)
    assert not body.is_disk()
    body.set_disk_shape(shape)
    assert body.is_disk()
    assert shape.radius() == 2.0


def test_time_stepping():

    data = m.disks.make_storage()
    simul = m.disks.add_simulation(data)
    simul.time_discretization().set_t0(0.1)
    simul.time_discretization().set_tmax(10)
    simul.time_discretization().set_h(0.005)

    timedisc = m.disks.add_time_discretization(data)
    timedisc.set_t0(0.3)

    print(f"t0={simul.time_discretization().t0()}")
    assert simul.time_discretization().t0() == 0.1
    assert simul.time_discretization().tmax() == 10
    assert simul.time_discretization().h() == 0.005


def test_interaction_manager():
    data = m.disks.make_storage()
    interman = m.disks.add_interaction_manager(data)
    nslaw = m.disks.add_nslaw(data)
    interman.insert_nonsmooth_law(nslaw, 0, 0).insert_nonsmooth_law(nslaw, 0, 1)


# def test_io():
#    data = m.disks.make_storage()
#    io = m.disks.add_io(data)
#    assert io.positions(0).shape[0] == 0
#    assert io.positions(0).shape[1] == 5
#    disk = m.disks.add_body(data)
#    assert io.positions(0).shape[0] == 1

# test_time_stepping()
