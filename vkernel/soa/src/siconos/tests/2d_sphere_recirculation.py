#!/usr/bin/env python
#
# Example of one object under gravity with one contactor and a ground
# using the Siconos proposed mechanics API
#

# rebuild the work of Thierry Faug in paper 2009 : "Mean steady flow on a wall overflow by free-surface gravity-driven dense flow"
# the grains have form of SPHERES

import sys
#from siconos.mechanics.collision.convexhull import ConvexHull2d
from siconos.mechanics.collision.tools import Contactor
import nonos
from nonos.mechanics_run import MechanicsHdf5Runner, MechanicsHdf5Runner_run_options
from nonos.bridge import SpaceFilterOptions
import siconos.numerics as sn
import siconos.kernel as sk
from siconos.mechanics.collision.bullet import SiconosBulletOptions

import numpy
import random
import math

backend = str(sys.argv[1])
nonos.mechanics_run.set_backend(backend)

bullet_options = SiconosBulletOptions()
bullet_options.worldScale = 100000
bullet_options.contactBreakingThreshold = 0.04
bullet_options.dimension = 1 # siconos.mechanics.collision.bullet.SICONOS_BULLET_2D
bullet_options.perturbationIterations = 0.
bullet_options.minimumPointsPerturbationThreshold = 0.


#inclination = math.pi*4/45. #16°
#inclination = math.pi*0.1 	 #18°
#inclination = math.pi/9     #20°
#inclination = math.pi*11/90 #22°
#inclination = math.pi*2/15  #24 °
#inclination = math.pi*13/90 #26°
#inclination = math.pi*7/45  #28°
inclination = 30      #30°
#inclination = math.pi*8/45  #32°

#default
en = 0.5
mu = 0.5
angle = inclination* math.pi/180.
n_row= 250
T=10.0

test = False
test_performance = False
performance_verbose= True

if (len(sys.argv) > 2):
    en = sys.argv[1]
    mu = float(sys.argv[2])
    inclination = float(sys.argv[3])
    angle = inclination * math.pi/180.
    if len(sys.argv) >= 5:
        n_row = int(sys.argv[4])
    if len(sys.argv) >= 6:
        T = float(sys.argv[5])

#if dist not in ['uniform', 'double', 'exp']:
#    print("dist = [uniform | double | exp]")
#    sys.exit(1)
if float(mu) < 0.1 or float(mu) > 2.0:
    print("mu = [0.1 .. 2.0]")
    sys.exit(1)


distribution = ('uniform', 0.15)


if test:
    n_row = 300
    n_col = 10
    N = n_row*n_col
    hstep = 1e-4
    #T=10*hstep
    T=0.5
    output_period = 0.01
    output_frequency=int (output_period / hstep)
    output_sample_period = 1.0
    output_sample_number = 100
    with_wall=True
    wall_x_position_number=700
    restart=False
else:
    n_col = 3
    n_row = 3
    N = n_row*n_col
    hstep = 1e-2
    output_period = 0.01
    output_frequency=int(output_period/ hstep)
    output_sample_period = 0.5
    output_sample_number = 100

    with_wall=True
    wall_x_position_number=700
    restart=False


if test_performance:
    en=0.0
    n_row = 200
    n_col = 100
    N = n_row*n_col
    T = 1.0
    hstep = 1e-4
    output_period = 0.01
    output_frequency=int(output_period/ hstep)
    output_sample_period = 4.0
    output_sample_number = 100

    with_wall=False
    wall_x_position_number=900
    restart=False

print('output_frequency', output_frequency)

# hdf5 file name
if with_wall:
    fn = '2d_sphere_flow_wall_{0}_N-{1}-e-{2}-mu-{3}-angle-{4:2.2g}.hdf5'.format(wall_x_position_number,N, en, mu, inclination)
else:
    fn = '2d_sphere_flow_N-{0}-e-{1}-mu-{2}-angle-{3:2.2g}.hdf5'.format(N, en, mu, inclination)

# mechanical parameters.

density = 2450.0


#grain and ground
grain_size = 0.001  # diameter of sphere



tank_width_number =n_col+2
tank_width = tank_width_number*grain_size*(1.0+distribution[1])
if with_wall:
    ground_size_number = tank_width_number+wall_x_position_number+200
else:
    ground_size_number = tank_width_number+wall_x_position_number+200



line_grain_size = grain_size

ground_size = ground_size_number*grain_size
ground_thickness = 10 * grain_size

ground_y_shift = 0.09


tank_thickness = 10 * grain_size

#gate : exit of tank
gate_size = 35*grain_size
gate_thickness = tank_thickness/2.


#tank front and behind
tank_behind_size = 1.5*n_row*grain_size + gate_size

tank_behind_x_position = -tank_thickness/2.0
tank_behind_y_position = tank_behind_size *0.5 - ground_y_shift - ground_thickness



tank_front_size = 1.5*n_row*grain_size
# compute the number of part of tank_front

tank_front_size_part = 20*grain_size

n_part_tank_front_size = math.ceil(tank_front_size/tank_front_size_part)
tank_front_size = tank_front_size_part * n_part_tank_front_size

tank_front_x_position = tank_width+tank_thickness+tank_thickness/2.0
tank_front_y_position = gate_size + tank_front_size_part *0.5 - ground_y_shift - ground_thickness

tank_behind_size = tank_front_size #+ gate_size

tank_behind_x_position = tank_thickness/2.0
tank_behind_y_position = tank_behind_size *0.5 - ground_y_shift - ground_thickness

gate_x_position = tank_front_x_position
gate_y_position = gate_size*0.5 - ground_y_shift - ground_thickness

grain_x_start=tank_width  + tank_thickness
grain_y_start=-100.0*grain_size

#wall of obstacle
wall_size_ratio = {}
if with_wall:
    wall_size_ratio[16] = 34 # for 16°
    wall_size_ratio[18] = 27 # for 18°
    wall_size_ratio[20] = 23 # for 20°
    wall_size_ratio[22] = 20 # for 22°
    wall_size_ratio[24] = 18 # for 24°
    wall_size_ratio[26] = 17 # for 26°
    wall_size_ratio[28] = 16 # for 28°
    wall_size_ratio[30] = 13 # for 30°
    wall_size_ratio[32] = 12 # for 32°

    if wall_size_ratio.get(inclination):
        wall_size_ratio_value =  wall_size_ratio.get(inclination)
        print('Wall size ratio is set to', wall_size_ratio_value)
    else:
        wall_size_ratio_value = 15
        print('Warning: wall size ration is set to', wall_size_ratio_value)

    wall_size = wall_size_ratio_value*grain_size
    wall_thickness = 10 * grain_size

    wall_x_position = wall_x_position_number*grain_size + tank_front_x_position +tank_thickness/2.0
    wall_y_position = wall_size/2. - ground_y_shift - ground_thickness



margin_ratio = 0.0  # 1e-05 # pourquoi ?

g  = 9.81
def apply_gravity(body):
	#Fonction define inclination of the gratvity
    g = 9.81
    weight = [body.scalarMass() * g * math.sin(angle), -
              body.scalarMass() * g * math.cos(angle), 0.]
    body.setFExtPtr(weight)  # scalMass() dans quel bibli ?


def create_grain(io, name, cname, grain_size=grain_size, density=1, trans=None,
                 tob=None, ground_sphere=True):
    #Fonction create one grain
    # Definition of a sphere as a primitive shape disk 2d
    io.add_primitive_shape(cname, 'Disk', (grain_size*0.5,), insideMargin=margin_ratio *
                           grain_size, outsideMargin=margin_ratio * grain_size)

    # computation of inertia and volume
    mass = density*(grain_size*1/2.)**2 * math.pi
    inertia = 1/2. * mass * (grain_size*1/2.)**2

    # print('geometric inertia:', inertia)
    # print('volume:', volume)
    # print('mass:', volume*density)
    # print('inertia:', inertia*density)
    io.add_object(name,
                  [Contactor(cname)],
                  translation=trans,
                  # velocity=veloci,
                  mass=mass,
                  time_of_birth=tob,
                  inertia=inertia)

# name : name of object
# cname : name of shape


def grains_locations(n_row=5, n_col=5, shift_ratio=3.0,
                     x_start=0.0, y_start=0.0, grain_size=0.05):
    grains_locations = []
    for j in range(n_row):
        for i in range(n_col):
            x_loc = -(i+1) * shift_ratio * grain_size + x_start
            y_loc = (j+1) * shift_ratio * grain_size + y_start
            # print(x_loc,y_loc)
            trans = [x_loc, y_loc]
            grains_locations.append(trans)

    return grains_locations


def progressbar(it, prefix="", size=60, file=sys.stdout):
    count = len(it)
    def show(j):
        x = int(size*j/count)
        file.write("%s[%s%s] %i/%i\r" % (prefix, "#"*x, "."*(size-x), j, count))
        file.flush()
    show(0)
    for i, item in enumerate(it):
        yield item
        show(i+1)
    file.write("\n")
    file.flush()

def create_grains(io, n_row=5, n_col=5, shift_ratio=3.0,
                  x_start=0.0, y_start=0.0,
                  grain_size=0.05, top=0, rate=0.01, density=1,
                  distribution=('uniform', 0.1)):
	#Fonction creat bulk of grains in tank
    N = n_row * n_col
    sizes = numpy.random.uniform(low=grain_size*(1.0-distribution[1]),
                                 high=grain_size*(1.0+distribution[1]), size=N)

    k = 0
    max_x = -1e24   # le size minimum ?
    max_y = -1e24
    min_x = 1e24   # le size minimum ?
    min_y = 1e24
    grains_locations = []
    for j in progressbar(range(n_row), "Creating grains: ", 60):
        for i in range(n_col):
            # initial translation
            # if (k % 100 == 0):
            #     print('.', end='', flush=True)
            x_loc = -(i+1) * shift_ratio * grain_size + x_start
            y_loc = (j+1) * shift_ratio * grain_size + y_start

            # print(x_loc,y_loc)
            trans = [x_loc, y_loc]
            grains_locations.append(trans)
            name = 'sphere' + '_' + str(i) + '_' + str(j)
            cname = 'SphereCS' + '_' + str(i) + '_' + str(j)
            create_grain(io, name, cname, sizes[k], density, trans,
                         tob=random.random() * rate)
            k += 1
            max_x = max(max_x, x_loc)
            max_y = max(max_y, y_loc)
            min_x = min(min_x, x_loc)
            min_y = min(min_y, y_loc)
    return max_x, max_y, min_x, min_y, grains_locations

# create a ground with line of spheres rigides


def create_grain_ground(io, name, cname, grain_size=grain_size, trans=None):
	#Fonction creat one grain of ground line
    io.add_primitive_shape(cname, 'Disk', (grain_size*0.5,), insideMargin=margin_ratio *
                           grain_size, outsideMargin=margin_ratio * grain_size)

    io.add_object(name, [Contactor(cname)], translation=trans)


def create_ground_line(ground_size_ratio=1500, line_grain_size=0.001):
	#Fonciton creat ground line made out of grain
    for s in range(ground_size_ratio):
        trans_ground = [(s+1)*line_grain_size - line_grain_size*0.5,-line_grain_size*0.5 - ground_y_shift - ground_thickness]
        create_grain_ground(io, 'sphere' +'_' + str(s),
                            'Sphere' + '_' + str(s), line_grain_size, trans_ground)

    return  - ground_y_shift - ground_thickness


if not restart:
    # Creation of the hdf5 file for input/output
    with MechanicsHdf5Runner(mode='w', io_filename=fn) as io:

        max_x, max_y, min_x, min_y, initial_grains_locations  = create_grains(io, n_row=n_row, n_col=n_col,
                                                                              x_start=grain_x_start, y_start=grain_y_start,
                                                                              shift_ratio=1.0+distribution[1], grain_size=grain_size, top=3,
                                                                              rate=0.0, density=density, distribution=distribution)

        # the ground object made with the ground shape. As the mass is
        # not given, it is a static object only involved in contact
        # detection.
        # Definition of the ground shape

        y_max_ground_line= create_ground_line(ground_size_number, line_grain_size)

        #io.add_primitive_shape('Ground', 'Box2d', (ground_size, ground_thickness),
        #insideMargin=0.0, outsideMargin=0.0)

        #io.add_object('ground', [Contactor('Ground')],
        #translation=[ground_size*0.5, - ground_y_shift - ground_thickness])
        if with_wall:
            io.add_primitive_shape('Wall', 'Box2d', (wall_thickness, wall_size),
                                   insideMargin=0.0, outsideMargin=0.0)

            io.add_object('wall', [Contactor('Wall')],
                          translation=[wall_x_position, wall_y_position])

        io.add_primitive_shape('Tank_front_part', 'Box2d', (tank_thickness, tank_front_size_part),
                               insideMargin=0.0, outsideMargin=0.0)

        for p in range(n_part_tank_front_size):
            io.add_object('tank_front_'+str(p), [Contactor('Tank_front_part')],
                          translation=[tank_front_x_position, tank_front_y_position + p*tank_front_size_part ])

        io.add_primitive_shape('Tank_behind_part', 'Box2d', (tank_thickness, tank_front_size_part),
                               insideMargin=0.0, outsideMargin=0.0)

        #io.add_object('tank_behind', [Contactor('Tank_behind')],
        #              translation=[tank_behind_x_position, tank_behind_y_position])

        for p in range(n_part_tank_front_size):
            io.add_object('tank_behind_'+str(p), [Contactor('Tank_behind_part')],
                          translation=[tank_behind_x_position, tank_front_y_position + p*tank_front_size_part ])

        io.add_primitive_shape('Gate', 'Box2d', (gate_thickness, gate_size),
                               insideMargin=0.0, outsideMargin=0.0)

        io.add_object('gate', [Contactor('Gate')],
                     translation=[gate_x_position, gate_y_position], time_of_death=0.1)

        io.add_primitive_shape('Gate_behind', 'Box2d', (tank_thickness, gate_size),
                               insideMargin=0.0, outsideMargin=0.0)

        io.add_object('gate_behon', [Contactor('Gate_behind')],
                     translation=[tank_behind_x_position, gate_y_position])

        # Definition of a non smooth law. As no group ids are specified it
        # is between contactors of group id 0.
        io.add_Newton_impact_friction_nsl('contact', mu=mu, e=en)

else:
    initial_grains_locations  = grains_locations(n_row=n_row, n_col=n_col,
                                                 x_start=grain_x_start, y_start=grain_y_start,
                                                 shift_ratio=1.1, grain_size=grain_size)

    import shutil
    shutil.copy2(fn, 'restart_from_'+ fn )


# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.

options = sk.solver_options_create(sn.SICONOS_FRICTION_2D_NSGS)
options.iparam[sn.SICONOS_IPARAM_MAX_ITER] = 100
options.dparam[sn.SICONOS_DPARAM_TOL] = 1e-3
options.iparam[sn.SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT] = 10

class output_end_run_iteration_hook():

    def __init__(self):
        self._io= None

        self._output_sample_number = output_sample_number
        self._output_sample_period = output_sample_period
        self._output_sample_step_period = int(self._output_sample_period / hstep)

        print('output every',  output_sample_period,' (s) that is every ', self._output_sample_step_period, 'steps with ', self._output_sample_number  , 'samples' )
        pass

    def initialize(self, io):
        self._io= io
        print('end_run_iteration_hook initialize')
        pass


    def call(self, step):
        if self._output_sample_step_period - self._output_sample_number <= (step+1) % self._output_sample_step_period < self._output_sample_step_period :
            print('end_run_iteration_hook: output results at step :', step+1, ' time =', self._io.current_time())
            self._io.output_results(with_timer=True)


class recirculation_start_run_iteration_hook():

    def __init__(self, initial_grains_locations):
        self._io= None
        self._initial_grains_locations= numpy.array(initial_grains_locations)
        self._initial_grains_locations_min = numpy.min(self._initial_grains_locations[:,1])
        self._initial_grains_locations_max = numpy.max(self._initial_grains_locations[:,1])
        self._period = 1e-01 # we recirculate grains every period
        self._frequency = 1/self._period
        self._step_period = int(self._period / hstep)
        self._max_recirculated_rows=20
        self._initial_velocity = self._max_recirculated_rows* grain_size / self._period #/(math.cos(angle))

        print('recirculation every', self._period,' (s) that is every ', self._step_period, 'steps' )
        print('              initial  recirculation velocity {:2.2g}'.format(self._initial_velocity),' (m/s) for travelling',self._max_recirculated_rows, '*grain_size per period'  )


        pass

    def initialize(self, io):
        self._io= io
        self._y_max_actual_recirculation = numpy.inf

        self._current_location_idx=0
        print('start_run_iteration_hook initialize')
        pass


    def recirculation(self,step):

        current_step = self._io._simulation.handle().current_step()
        positions  = self._io._io.positions(self._io._nsds)

        if positions is not None:
#            y = positions[:,2]

            # We search for the ds index that are below a given criteria
#            ds_idx = numpy.nonzero(y < y_max_ground_line)[0]
            x = positions[:,1]
            ds_idx  = numpy.nonzero(x > ground_size / 20 )[0]
            # when we start to recirculate grains, we compute the actual height
            # of the pile in the tank
            if len(ds_idx) > 0 :
                self._y_max_actual_recirculation = numpy.max(positions[:,2])
                print("positions:")
                print(positions)
                #print('self._y_max_actual_recirculation', self._y_max_actual_recirculation)

            y_shift = - self._initial_grains_locations_min + self._y_max_actual_recirculation + grain_size
            self._current_location_idx = 0 # we recirculate at the initial point

            nsds = self._io._nsds

            print('ds_idx = {}'.format(ds_idx))
            for i in ds_idx :

                print('i={}, type(i)={}'.format(i, type(i)))
                n_ds = int(positions[i,0])

                print('n_ds={}'.format(n_ds))

                ds = nsds.dynamicalSystem(n_ds)

                print('ds.number()={}, n_ds={}'.format(ds.number(), n_ds))
                assert (ds.number() == n_ds)
                # we reset the grain location on the shifted initial grid
                x_loc =  self._initial_grains_locations[self._current_location_idx,0]
                y_loc =  self._initial_grains_locations[self._current_location_idx,1] + y_shift

                print('reset initial position ds', n_ds, 'with position', positions[n_ds-1,1], positions[n_ds-1,2], 'at position', x_loc,y_loc)
                ds.setQ0Ptr([x_loc,y_loc,0])

                ds.setVelocity0Ptr([0.0, 0.01, 0.], current_step)
#                ds.setVelocity0Ptr([- self._initial_velocity,- self._initial_velocity,0], current_step+1)
                ds.resetToInitialState()
                ds.swapInMemory()

                # we increment the index on the initial grid
                self._current_location_idx = self._current_location_idx +1

            print('recirculation number of bodies recirculated :', len(ds_idx), ' flow rate',len(ds_idx)*self._frequency ,'grain/s')

            if (len(ds_idx) /n_col >= self._max_recirculated_rows):
                print('Warning: too many bodies are recirculated')
        pass


    def call(self, step):
        #print('step', step)
        #print('step_period', self._step_period)

        if step % self._step_period == 0:
            print('start_run_iteration_hook: recirculation at step :', step)
            self.recirculation(step)



srih = recirculation_start_run_iteration_hook(initial_grains_locations)
erih = output_end_run_iteration_hook()



run_options=MechanicsHdf5Runner_run_options()
run_options['t0']=0
run_options['T']= T
run_options['h']=hstep

run_options['bullet_options']=bullet_options
run_options['solver_options']=options
run_options['constraint_activation_threshold']=1e-05
run_options['start_run_iteration_hook']=srih
run_options['end_run_iteration_hook'] = erih

run_options['Newton_options']=sk.SICONOS_TS_LINEAR

run_options['skip_last_update_output']=True
run_options['skip_reset_lambdas']=True
run_options['osns_assembly_type']= sk.REDUCED_DIRECT

#run_options['osns_assembly_type']= sk.GLOBAL_REDUCED
#run_options['osi']= sk.MoreauJeanGOSI
#run_options['skip_last_update_input']=True
if performance_verbose:
    run_options['verbose']=True
    run_options['with_timer']=True
    run_options['explode_Newton_solve']=False
    run_options['explode_computeOneStep']=False

run_options['violation_verbose'] = True
run_options['output_frequency']=output_frequency

run_options['time_stepping']=None
run_options['set_external_forces']=apply_gravity

vnative_options = SpaceFilterOptions()
vnative_options.neighborhood_radius = grain_size * 1.01
vnative_options.min_radius = grain_size / 2

run_options['vnative_options'] = vnative_options

with MechanicsHdf5Runner(mode='r+',  io_filename=fn) as io:

    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    io.run(run_options)
