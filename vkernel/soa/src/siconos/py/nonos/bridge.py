import nonos as vkernel
import numpy as np
from math import sqrt
import siconos.numerics as sn

from math import pi

def array(l):
    return np.array(l, dtype=np.float64)

class Stored():
    _data = None

    @classmethod
    def setStorage(cls, data):
        cls._data = data

    @classmethod
    def data(cls):
        return cls._data

    def handle(self):
        return self._handle

class SpaceFilter(Stored):

    _static_shape_counter = 1 # in cf_info 0 -> shape associated to some ds

    def __init__(self, options):
        self._options = options
        self._initialized = False
        self._ngbh = vkernel.disks.add_neighborhood(self.data())
        self._ngbh.create(options.neighborhood_radius)
        self._interman = vkernel.disks.add_interaction_manager(self.data())
        self._handle = vkernel.disks.add_space_filter(self.data())
        self._handle.set_neighborhood(self._ngbh)
        self._handle.set_diskdisk_r(vkernel.disks.add_diskdisk_r(self.data()))
        self._fdisks = {}

    # Fixed segment
    def insertSegment(self, x1, y1, x2, y2):
        segment = vkernel.disks.add_segment_shape(self.data())
        segment.set_points(array([x1, y1, 0, x2, y2, 0]), 0)

        mp = int(max(3, sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)) / self._options.min_radius))

        segment.set_maxpoints(mp) # fix / size of smallest disk
        segment.initialize(0)

        # ident is attached in io.hpp
        # it is negative for static shape
        segment.set_ident(- self._static_shape_counter)
        self._static_shape_counter = self._static_shape_counter + 1

        diskfsegment = vkernel.disks.add_diskfsegment_r(self.data())
        diskfsegment.set_segment(segment)

        self.handle().insert_diskfsegment_r(diskfsegment)

        return segment

    def insertMesh(self, points):
        mesh

    # Fixed line (to be removed)
    def insertLine(self, a, b , c):
        line = vkernel.disks.add_line_shape(self.data())
        line.set_a(a)
        line.set_b(b)
        line.set_c(c)
        line.set_maxpoints(10000)
        line.initialize()

         # negative for static shape
        line.set_ident(- self._static_shape_counter)
        self._static_shape_counter = self._static_shape_counter + 1

        print("new line,  p0:", line.p0())
        print("--------, dir:", line.direction())

        diskline = vkernel.disks.add_diskline_r(self.data())
        diskline.set_line(line)

        self.handle().insert_line(diskline)

        return line

    # Fixed disk
    def insertTranslatedDisk(self, radius, translation):

        disk_shape = None
        if radius in self._fdisks:
            disk_shape = self._fdisks[radius]
        else:
            disk_shape = vkernel.disks.add_disk_shape(self.data())
            disk_shape.set_radius(radius)
            self._fdisks[radius] = disk_shape

        mp = int(2 * pi * radius / self._options.min_radius)

        translated_disk_shape = \
            vkernel.disks.add_translated_disk_shape(self.data())
        translated_disk_shape.set_translated(disk_shape)
        translated_disk_shape.set_translation(array([translation[0],translation[1],0]))
        translated_disk_shape.translated().set_maxpoints(mp)

        # negative for static shape
        translated_disk_shape.set_ident(- self._static_shape_counter)
        self._static_shape_counter = self._static_shape_counter + 1

        diskfdisk = vkernel.disks.add_diskfdisk_r(self.data())
        diskfdisk.set_translated_disk_shape(translated_disk_shape)

        self.handle().insert_diskfdisk_r(diskfdisk)

        return translated_disk_shape

    def insertNonSmoothLaw(self, nslaw, gid1, gid2):
        self._handle.set_nslaw(nslaw.handle()) # one nslaw !
        self._interman.insert_nonsmooth_law(nslaw.handle(), gid1, gid2)

    def updateInteractions(self, step):
        if not self._initialized:
            self._handle.make_points()
            self._ngbh.add_point_sets(step)
            self._ngbh.set_active(0, 0, True)       # disk - disk
            self._ngbh.set_active(0, 1, True)       # disk - segment
            self._ngbh.set_active(0, 2, True)       # disk - fixed disk
            self._ngbh.set_active(1, 1, False)
            self._ngbh.set_active(2, 2, False)
            self._ngbh.set_active(1, 0, False)
            self._ngbh.set_active(2, 0, False)
            self._ngbh.set_active(1, 2, False)
            self._ngbh.set_active(2, 1, False)

            self._initialized = True

        self._ngbh.update(step)

        self._ngbh.search()
        self._handle.update_index_set0(step);
        self._ngbh.sort()

    def removeStaticBody(self, body):
        self._ngbh.search() # needed after sort

        # only segments are removable
        if type(body) == type([]):
            # box2d
            self.handle().remove_static_segment(0, body[0])
            self.handle().remove_static_segment(0, body[1])
            self.handle().remove_static_segment(0, body[2])
            self.handle().remove_static_segment(0, body[3])
        else:
            self.handle().remove_static_segment(0, body)


class NewtonImpactFrictionNSL(Stored):

    def __init__(self, e, not_used, mu, dimension):
        self._handle = vkernel.disks.add_nslaw(self.data())
        self._handle.set_e(e)
        self._handle.set_mu(mu)
        #self._handle.set_dimension(dimension)

class Osi(Stored):

    def __init__(self, theta):
        self._handle = vkernel.disks.add_osi(self.data())
        self._handle.assembled_osi().set_theta(theta)

    def setTheta(self, theta):
        self.handle().assembled_osi().set_theta(theta)

    def setConstraintActivationThreshold(self, cat):
        self.handle().assembled_osi().set_constraint_activation_threshold(cat)

    def setGamma(self, gamma):
        self.handle().assembled_osi().set_gamma(gamma)


class Topology(Stored):

    def __init__(self):
        self._handle = vkernel.disks.add_topology(self.data())

    def indexSetsSize(self):
        return 2

class NonSmoothDynamicalSystem(Stored):

    def __init__(self, t0, T):
        self._t0 = t0
        self._T = T
        self._topology = Topology()
        self._mapid = {}

    def topology(self):
        return self._topology

    def insertDynamicalSystem(self, body):
        self._mapid[int(body.number())] = body

    def dynamicalSystem(self, ds_id):
        return self._mapid[int(ds_id)]



class TimeDiscretisation(Stored):

    def __init__(self, t0, h):
        self._t0 = t0
        self._h = h
        self._handle = \
            vkernel.disks.add_time_discretization(self.data())
        self.handle().set_t0(t0)
        self.handle().set_h(h)

class Simulation(Stored):

    def __init__(self, nsds, timedisc):
        self._need_init = True
        self._nsds = nsds
        self._timedisc = timedisc
        t0 = timedisc.handle().t0()
        h = timedisc.handle().h()
        self._handle = vkernel.disks.add_simulation(self.data())

        # a time discretization is added during add_simulation: replace it
        self._handle.set_time_discretization(self._timedisc.handle())

        self._timedisc.handle().set_tmax(self._nsds._T) # vkernel does not have nsds
        self._timedisc.handle().set_h(h)
        self._timedisc.handle().set_t0(t0)

        self.handle().initialize()
        self.handle().one_step_integrator().assembled_osi().set_theta(0.50001)

    def insertIntegrator(self, osi):
        pass # unimplemented

    def insertNonSmoothProblem(self, osnspb):
        self._osnspb = osnspb

    def insertInteractionManager(self, interman):
        self._interman = interman

    def setNewtonOptions(self, nopts):
        pass # unimplemented

    def setNewtonMaxIteration(self, maxiter):
        pass # unimplemented

    def setNewtonTolerance(self, newtontol):
        pass # unimplemented

    def setNewtonWarningOnNonConvergence(self, opt):
        pass # unimplemented

    def setSkipLastUpdateOutput(self, skipluo):
        pass # unimplemented

    def setSkipLastUpdateInput(self, skiplui):
        pass # unimplemented

    def setSkipResetLambdas(self, skipresetlbds):
        pass # unimplemented

    def setDisplayNewtonConvergence(self, dnc):
        pass # unimplemented

    def setWarningNonsmoothSolver(self, opt):
        pass #unimplemented

    def startingTime(self):
        return self.handle().time_discretization().t0()

    def nextTime(self):
        return self.startingTime() + self.handle().current_step() *\
            self.handle().time_discretization().h()

    def hasNextEvent(self):
        return self.handle().has_next_event()

    def updateInteractions(self):
        self._interman.updateInteractions(self.handle().current_step())

    def computeOneStep(self):
        return self.handle().compute_one_step()

    def clearNSDSChangeLog(self):
        pass # unimplemented

    def nextStep(self):
        pass # unimplemented

    def oneStepNSProblem(self, idx):
        return self._osnspb

# Allow for enumerated bodies with associated disk shape.
class BodyBase(Stored):

    __count = 0
    __shapes = {}

    @classmethod
    def count(cls):
        return cls.__count

    @classmethod
    def set_count(cls, newcount):
        cls.__count = newcount

    @classmethod
    def shapes(cls):
        return cls.__shapes

class Body(BodyBase):

    def __init__(self):
        pass

    def init_fem(self, mesh_data):

        self.set_count(self.count() + 1)
        self._ident = self.count()

        body = vkernel.disks.add_fem(self.data())
        self._handle = body
        body.set_id(self._ident)
        self._mesh_data = mesh_data

    def init_disk(self, radius, mass, position, velocity):

        self.set_count(self.count() + 1)
        self._ident = self.count()

        body = vkernel.disks.add_disk(self.data())
        self._handle = body
        body.set_id(self._ident)
        body.set_q(array(position))
        body.set_velocity(array(velocity))
        body.set_mass_matrix(array([mass, mass, mass*radius*radius/2]))

        disk_shape = None
        if radius in self.shapes():
            disk_shape = self.shapes()[radius]
        else:
            disk_shape = vkernel.disks.add_disk_shape(self.data())
            disk_shape.set_radius(radius)
            self.shapes()[radius] = disk_shape

        body.set_shape(disk_shape)
        body.set_fext(array([0,0,0])) # default

    def getMassValue(self):
        return self.handle().mass_matrix()[0]

    def setConstantFext(self, fext, nargument=None):
        self.handle().set_fext(array(fext))

    def setFExtPtr(self, fext):
        self.handle().set_fext(array(fext))

    def setNumber(self, num):
        # FIX: in mechanics_run insertDynamicalsystem is done before
        # setNumber and as the map is handled in bridge.py by a python
        # map, this is buggy.
        # workaround: setNUmber is disabled
        # self.handle().set_id(num)
        pass

    def number(self):
        return self.handle().id()

    def setQ0Ptr(self, pos, step=0):

        # FIX: next step not necessary 1
        self.handle().set_q_at_step(array(pos), 0)
        self.handle().set_q_at_step(array(pos), 1)

    def setVelocity0Ptr(self, vel, step=0):

        # FIX: next step not necessary 1
        self.handle().set_velocity_at_step(array(vel), step)
        self.handle().set_velocity_at_step(array([0., 0., 0.]), step+1)

    def resetToInitialState(self):
        pass # compatibility

    def swapInMemory(self):
        pass # compatibility

class BodyWrap():

    def __init__(self, bdy):

        self._bdy = bdy

    def getMassValue(self):
        return self._bdy.mass_matrix()[0]

    def setConstantFext(self, fext, nargument=None):
        self._bdy.set_fext(array(fext))

    def number(self):
        return self._bdy.id()

    def setQ0Ptr(self, pos, step=0):

        # FIX: next step not necessary 1
        self._bdy.set_q_at_step(array(pos), 0)
        self._bdy.set_q_at_step(array(pos), 1)

    def setVelocity0Ptr(self, vel, step=0):

        # FIX: next step not necessary 1
        self._bdy.set_velocity_at_step(array(vel), step)
        self._bdy.set_velocity_at_step(array([0., 0., 0.]), step+1)

    def resetToInitialState(self):
        pass # compatibility

    def swapInMemory(self):
        pass # compatibility

# many bodies
class Bodies(BodyBase):

    def __init__(self, radius, mass, positions, velocities):

        start = self.count()
        self.set_count(self.count() + len(positions))
        self._ident = np.arange(start, self.count())

        bodies = vkernel.disks.multiple_add_disk(self.data(), len(positions))
        self._handle = bodies
        bodies.multiple_set_id(self._ident)
        bodies.multiple_set_q(array(positions))
        bodies.multiple_set_velocity(array(velocities))
        bodies.set_mass_matrix(array([mass, mass, mass*radius*radius/2]))

        disk_shape = None
        if radius in self.shapes():
            disk_shape = self.shapes()[radius]
        else:
            disk_shape = vkernel.disks.add_disk_shape(self.data())
            disk_shape.set_radius(radius)
            self.shapes()[radius] = disk_shape

        bodies.set_shape(disk_shape)
        bodies.set_fext(array([0,0,0])) # default

    def get(self):
        return [BodyWrap(bdy) for bdy in self.handle().get()]

    def getMassValue(self):
        return self.handle().mass_matrix()[0]

    def setConstantFext(self, fext, nargument=None):
        self.handle().set_fext(array(fext))

    def setFExtPtr(self, fext):
        self.handle().set_fext(array(fext))

    def setNumber(self, num):
        self.handle().set_id(num)

    def number(self):
        return self.handle().id()

    def setQ0Ptr(self, pos, step=0):

        # FIX: next step not necessary 1
        self.handle().set_q_at_step(array(pos), 0)
        self.handle().set_q_at_step(array(pos), 1)

    def setVelocity0Ptr(self, vel, step=0):

        # FIX: next step not necessary 1
        self.handle().set_velocity_at_step(array(vel), step)
        self.handle().set_velocity_at_step(array([0., 0., 0.]), step+1)

    def resetToInitialState(self):
        pass # compatibility

    def swapInMemory(self):
        pass # compatibility

class TraceParams(Stored):
    def __init__(self, tp_args):
        self._args = tp_args
        self._handle = vkernel.disks.add_trace_params(self.data())
        self.handle().set_filename(tp_args._fileName)
        self.handle().set_title(tp_args._title)
        self.handle().set_description(tp_args._description)
        self.handle().set_math_info(tp_args._mathInfo)


class OSNSPB(Stored):
    def __init__(self, dim, solvopts):

        self._so = vkernel.disks.add_solver_options(self.data())

        if solvopts is not None:
            self._so.create(solvopts.solverId)
            for i,v in enumerate(solvopts.iparam):
                self._so.set_iparam(i, v)

            for i,v in enumerate(solvopts.dparam):
                self._so.set_dparam(i, v)

        else:
            # default
            self._so.create(sn.solver_ids.SICONOS_FRICTION_2D_NSGS)
            self._so.set_iparam(sn.params.SICONOS_IPARAM_MAX_ITER, 100)
            self._so.set_dparam(sn.params.SICONOS_DPARAM_TOL, 1e-3)
            self._so.set_dparam(sn.params.SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT, 10)

        self._fc2d = vkernel.disks.add_fc2d(self.data())
        self._handle = vkernel.disks.add_osnspb(self.data())
        self._handle.set_options(self._so)
        self._fc2d.create(self._so.solver_id())
        self.handle().set_problem(self._fc2d)
        self.handle().set_verbose(False)
        self.handle().set_trace(False)
#        self._fc2d.instance().dimension = 2

    def setMaxSize(self, maxs):
        self._maxSize = maxs

    def setMStorageType(self, stype):
        self._storageType = stype

    def setNumericsVerboseMode(self, vm):
        self._numericsVerboseMode = vm

    def setKeepLambdaAndYState(self, kly):
        self._keepLambdaAndYState = kly

    def setAssemblyType(self, at):
        pass

    def getSizeOutput(self):
        return 0 # unimplemented

    def numericsSolverOptions(self):
        return self.handle().options()


SICONOS_TS_NONLINEAR = "SICONOS_TS_NONLINEAR"


MoreauJeanOSI = Osi

class Unimplemented():

    def __init__(self, *args):
        assert False

MoreauJeanGOSI = Unimplemented
TimeSteppingDirectProjection = Unimplemented

FrictionContact = OSNSPB

class MechanicsIO(Stored):

    def __init__(self):
        self._handle = vkernel.disks.add_io(self.data())
        self._simul = None

    # current step needed for output
    def setSimulation(self, simul):
        self._simulation = simul

    def p0s(self, nsds):
        return self.handle().p0s(
            self._simulation.handle().current_step())

    def radii(self, nsds):
        return self.handle().radii(
            self._simulation.handle().current_step())

    def displacements(self, nsds):
        return self.handle().displacements(
            self._simulation.handle().current_step())

    def positions(self, nsds):
        return self.handle().positions(
            self._simulation.handle().current_step())

    def velocities(self, nsds):
        return self.handle().velocities(
            self._simulation.handle().current_step())

    def contactPoints(self, nsds, output_contact_index_set):
        return self.handle().contact_points(
            self._simulation.handle().current_step())

    def contactInfo(self, nsds, output_contact_index_set):
        return self.handle().contact_info(
            self._simulation.handle().current_step())

    def contactContactWork(self, nsds, output_contact_index_set, omega,
                           tol=1e-8):
        return self.handle().contact_work(
            self._simulation.handle().current_step(), omega, tol)

class SpaceFilterOptions():

    neighborhood_radius = 2.1
    min_radius = 0.5

    def toJson(self):
        return json.dumps(self, default=lambda o: o.__dict__)
