#!/usr/bin/env @Python_EXECUTABLE@
# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
#
# Copyright 2025 INRIA.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Run a pre-configured Siconos "mechanics" HDF5 file."""

from math import cos, sin, acos, sqrt, atan, pi
import scipy.constants
import numpy as np

import bisect
import numbers
import shutil
import json

# Siconos imports
import siconos.io.mechanics_hdf5
import siconos.nonsmooth_formulations
import siconos.numerics as sn
import siconos.modeling as sm
import siconos.nonsmooth_formulations as nsf
import siconos.simulation as simu
import siconos.integrators as integrators

# Siconos Mechanics imports
# import siconos.mechanics.collision.tools as smc_tools
import siconos.mechanics
import siconos.mechanics.collision
import siconos.mechanics.joints
import siconos.mechanics.quaternions
import siconos.io as sio
from siconos.io.FrictionContactTrace import GlobalFrictionContactTrace as GFCTrace
from siconos.io.FrictionContactTrace import FrictionContactTrace as FCTrace
from siconos.io.FrictionContactTrace import (
    GlobalRollingFrictionContactTrace as GRFCTrace,
)
import siconos.io.shape_collection


# rotation around the origin
def rotate_point(x, y, alpha):
    x_new = x * cos(alpha) - y * sin(alpha)
    y_new = x * sin(alpha) + y * cos(alpha)
    return x_new, y_new


class RunnerConfig:
    """A class to manage the backend and the default classes to be used

    This class allows you to configure the backend for the simulation.
    It supports three backends:
    - 'bullet' for Bullet Physics Engine
    - 'occ' for OpenCascade Geometry (occ)
    - 'native' for a native simulation
    - 'vnative' for a vectorized native simulation

    """

    def __init__(self, backend="bullet"):
        self.default_manager_class = None
        self.default_simulation_class = None
        self.default_body_class = None
        self.use_bullet = False
        self.backend = backend
        self.occ = None
        self.bullet = None
        self.set_backend(backend)

    def setup_default_classes(self):
        import siconos.mechanics

        have_bullet = siconos.mechanics.have_bullet
        have_occ = siconos.mechanics.have_occ
        if self.backend == "bullet" and have_bullet:
            import siconos.mechanics.collision.bullet

            self.default_simulation_class = siconos.simulation.TimeStepping
            self.default_body_class = siconos.mechanics.collision.RigidBodyDS
            self.bullet = siconos.mechanics.collision.bullet

            def bullet_manager(bullet_options):
                if bullet_options is None:
                    bullet_options = self.bullet.SiconosBulletOptions()
                return self.bullet.SiconosBulletCollisionManager(bullet_options)

            self.default_manager_class = bullet_manager
            self.use_bullet = have_bullet
        elif self.backend == "occ" and have_occ:
            import siconos.mechanics.occ as occ

            self.occ = occ
            self.default_manager_class = lambda options: occ.OccSpaceFilter()
            self.default_simulation_class = occ.OccTimeStepping
            self.default_body_class = occ.OccBody
        elif self.backend == "native":

            def native_manager(options):
                sp = siconos.mechanics.collision.SpaceFilter()
                sp.setBBoxfactor(3)
                sp.setCellsize(6)
                return sp

            self.default_manager_class = native_manager
            self.default_simulation_class = siconos.simulation.TimeStepping
            self.default_body_class = siconos.mechanics.collision.RigidBodyDS

        elif self.backend == "vnative":
            global sm
            global sio
            global nsf
            # global simu
            import nonos as vkernel
            import nonos.bridge as sm
            import nonos.bridge as sio
            import nonos.bridge as nsf

            self.use_bullet = False
            sm.Stored.setStorage(vkernel.disks.make_storage())
            self.default_manager_class = sio.SpaceFilter
            self.default_simulation_class = sio.Simulation
            self.default_body_class = sio.Body
            self.default_bodies_class = sio.Bodies

        else:
            raise RuntimeError(f"Unavailable chosen backend {self.backend}")

    def set_backend(self, b):
        self.backend = b
        self.setup_default_classes()


class MechanicsHdf5Runner_run_options(dict):
    def __init__(self):
        d = {}
        self._d_comment = {}

        self._valid_options_keys = []

        def create_option(d, key, type_t, default, info):
            self._valid_options_keys.append(key)

            d[key] = default
            self._d_comment[key] = {}
            self._d_comment[key]["type"] = type_t
            self._d_comment[key]["default"] = default
            self._d_comment[key]["info"] = info

        create_option(
            d, "t0", "real, optional", 0.0, """initial time of the simulation"""
        )
        create_option(
            d, "T", "real, optional", 10.0, """final time of the simulation"""
        )

        create_option(d, "h", "real, optional", 5e-4, """time-step size""")
        create_option(
            d,
            "theta",
            "real, optional",
            0.50001,
            """theta parameter for Moreau-Jean one-step integrator in [0,1] """,
        )
        create_option(
            d,
            "gamma",
            "real, optional",
            None,
            """gamma parameter for Moreau-Jean one-step integrator OSI  """,
        )
        create_option(
            d,
            "set_external_forces",
            "python function, optional",
            None,
            """ if None, set to self.apply_gravity
                      function used to apply forces onto the body with signature :
                      def funcname(body):
                      ...
                      with body  a siconos Body (siconos.SecondOrderDS) """,
        )
        create_option(
            d,
            "gravity_scale",
            "real, optional",
            1.0,
            """ scaling factor for the gravity.
                      1.     for meters (default).
                      1./100 for centimeters.
                      This parameter may be needed for small
                      objects because of Bullet collision margin (0.04). """,
        )

        create_option(
            d,
            "bullet_options",
            "?, optional",
            None,
            """set of options for the interaction manager (e.g. SiconosBulletOptions)""",
        )

        create_option(
            d,
            "vnative_options",
            "?, optional",
            None,
            """set of options for the vnative backend""",
        )

        create_option(
            d,
            "multipoints_iterations",
            "boolean, optional",
            None,
            """if true use bullet "multipoint iterations (Obsolete)""",
        )

        create_option(
            d,
            "time_stepping",
            "siconos.kernel.Simulation, optional",
            None,
            "a siconos simulation type instance",
        )
        create_option(
            d,
            "interaction_manager",
            "SiconosCollisionManager, optional",
            None,
            """ user-defined interaction handler (e.g. from Bullet)
                      (depends on the backend, e.g. Bullet or OCC).
                      Warning: overwrite the value
                      provided during MechanicsHdf5Runner init.""",
        )

        # default hook options
        create_option(d, "start_run_iteration_hook", "", None, """ """)
        create_option(d, "before_next_step_iteration_hook", "", None, """ """)
        create_option(d, "end_run_iteration_hook", "", None, """ """)
        create_option(d, "controller", "", None, """ """)

        # create_optionaln(d,
        #                      "body_class",
        #                      "siconos.mechanics.RigidBodyDS, optional",
        #                      None,
        #                      """ class used for body definition (e.g. RigidBodyDS)""")
        # create_option(d,
        #                      "shape_class",
        #                      "?, optional",
        #                      None,
        #                      """class used for shape definition (e.g. occ.OccContactShape)""")
        # create_option(d,
        #                      "face_class",
        #                      "?, optional",
        #                      None,
        #                      """ class used for face definition (e.g. occ.OccContactFace) occ only?""")
        # create_option(d,
        #                      "edge_class",
        #                      "?, optional",
        #                      None,
        #                      """ class used for edge definition (e.g. occ.OccContactEdge) occ only?""")

        # default osi options
        create_option(
            d,
            "osi",
            "integrators.OneStepIntegrator, optional",
            integrators.MoreauJeanOSI,
            "class type used to describe one-step integration",
        )
        create_option(
            d,
            "constraint_activation_threshold",
            "real, optional",
            None,
            """threshold under which constraint is assume to be active.
                      if None, default value is taken as default value of osi class""",
        )
        create_option(
            d,
            "constraint_activation_threshold_velocity",
            "real, optional",
            None,
            """threshold under which constraint is assume to be active at the veloicity level.
                      if None, default value is taken as default value of osi class """,
        )
        create_option(
            d,
            "activate_with_negative_relative_velocity",
            "real, optional",
            None,
            """ activate constraints at the velocity level only if the relative velocity is negative
                      if None, default value is taken as default value of osi class """,
        )
        create_option(
            d,
            "projection_itermax",
            "int, optional",
            20,
            """ max number of iteration for projection
                      (only for TimeSteppingDirectProjection) """,
        )
        create_option(
            d,
            "projection_tolerance",
            "real, optional",
            1e-8,
            """ tolerance for the violation of the equality constraints at the  position level
                      (only for TimeSteppingDirectProjection)""",
        )
        create_option(
            d,
            "projection_tolerance_unilateral",
            "real, optional",
            1e-8,
            """ tolerance for the violation of the unilateral constraints at the  position level
                      (only for TimeSteppingDirectProjection) """,
        )

        # default Newton solve options
        create_option(
            d,
            "Newton_options",
            "siconos.simulation.TYPE, optional",
            siconos.simulation.NONLINEAR,
            """simu.TimeStepping options to control the Newton loop
                      possible values : LINEAR, LINEAR_IMPLICIT, NONLINEAR, NONLINEAR_FULL""",
        )
        create_option(
            d,
            "Newton_max_iter",
            "int, optional",
            20,
            """maximum number of iterations allowed for the Newton method""",
        )
        create_option(
            d,
            "Newton_tolerance",
            "real, optional",
            1e-10,
            """ required tolerance for the Newton method""",
        )
        create_option(
            d,
            "Newton_warning_on_nonconvergence",
            "boolean, optional",
            True,
            """ display a warning if the Newton method does not converge""",
        )
        create_option(
            d,
            "Warning_nonsmooth_solver",
            "boolean, optional",
            True,
            """ display a warning if the nonsmooth does not converge""",
        )
        create_option(
            d,
            "display_Newton_convergence",
            "boolean, optional",
            False,
            """ display the information about the convergence of the Newton method""",
        )
        create_option(
            d,
            "skip_last_update_output",
            "boolean, optional",
            False,
            """ Skip the computation of the last update of the output (kinematic contact variable)
                      in order to save time""",
        )
        create_option(
            d,
            "skip_last_update_input",
            "boolean, optional",
            False,
            """ Skip the computation of the last update of the input (contact reaction variable)
                      in order to save time""",
        )
        create_option(
            d,
            "skip_reset_lambdas",
            "boolean, optional",
            False,
            """ Skip the reset to 0.0 of the input (contact reaction variable) """,
        )

        # default osnpb options
        create_option(
            d,
            "solver_options",
            "numerics.SolverOptions, optional",
            None,
            """ SolverOptions for siconos.numerics solvers
                      we advice to create it with  siconos.numerics.solver_options_create
                      if solver_option is None, we leave siconos/numerics choosing the default option
                      (see numerics solvers documentation for details)""",
        )
        create_option(
            d,
            "solver_options_pos",
            "numerics.SolverOptions, optional",
            None,
            """ SolverOptions for siconos.numerics solvers at the poistion level
                      we advice to create it with  siconos.numerics.solver_options_create
                      if solver_option is None, we leave siconos/numerics choosing the default option
                      (see numerics solvers documentation for details)""",
        )

        create_option(
            d,
            "osnspb_max_size",
            "int, optional",
            0,
            """  estimation of the maximum number of dynamical systems taken into account.
                      Useful for memory pre-allocations and optimisation.
                      if equal to 0, it will be set to
                      simulation().nonSmoothDynamicalSystem().topology().numberOfConstraints() """,
        )
        create_option(
            d,
            "osns_assembly_type",
            "siconos.nonsmooth_formulations.TYPE, optional",
            None,
            """  Assembly method for the Delassus operator. Possible value are:
                      siconos.nonsmooth_formulations.REDUCED_BLOCK, GLOBAL, REDUCED_DIRECT, GLOBAL_REDUCED.
                      if None, the osnspb uses thd default assembly type defined in LinearONSS.hpp""",
        )

        d["friction_contact_trace_params"] = None
        create_option(
            d,
            "friction_contact_trace_params",
            "FrictionContactTraceParams instance, optional",
            None,
            """  """,
        )

        # default output options
        create_option(
            d,
            "output_frequency",
            "int, optional",
            None,
            """ log and screen outputs frequency
                      The initial step (k=0) and the first step (k=1) are always written in hdf5 file.
                      if None set equal to 1. set equal to 0 to cancel output""",
        )
        create_option(
            d,
            "output_backup",
            "boolean, optional",
            False,
            """ if True, make a backup of the output (hdf5 file)""",
        )
        create_option(
            d, "output_backup_frequency", "int, optional", None, """ backup frequency"""
        )
        create_option(
            d,
            "output_contact_index_set",
            "int, optional",
            1,
            """ index set level for outputting contact info
                      1: only the contact active at the velocity level are written
                      0: all contact are written""",
        )
        create_option(
            d,
            "output_contact_forces",
            "boolean, optional",
            True,
            """if True, the contact forces are written in the hdf5 file  """,
        )
        create_option(
            d,
            "output_contact_info",
            "boolean, optional",
            True,
            """if True, the contact information are written in the hdf5 file.
                      it mainly contains for each contact, the bodies involved in the contact.""",
        )

        create_option(
            d,
            "output_contact_work",
            "boolean, optional",
            True,
            """if True, the contact work is computed and are written in the hdf5 file""",
        )
        create_option(
            d,
            "output_energy_work",
            "boolean, optional",
            False,
            """if True, the kinetic and work for each bodies are computed and written in the hdf5 file""",
        )

        # default verbose and debug  options
        create_option(
            d,
            "verbose",
            "boolean, optional",
            True,
            """if true, print current step information""",
        )
        create_option(
            d,
            "verbose_progress",
            "boolean, optional",
            True,
            """if true, print the number of step""",
        )
        create_option(
            d,
            "numerics_verbose",
            "boolean, optional",
            False,
            """if true, activate numerics verbosity """,
        )
        create_option(
            d,
            "numerics_verbose_level",
            "int, optional",
            0,
            """numerics verbosity level""",
        )
        create_option(
            d,
            "violation_verbose",
            "boolean, optional",
            False,
            """if true, print info about contact violation """,
        )

        create_option(
            d,
            "with_timer",
            "boolean, optional",
            False,
            """if true, use a timer for log output in std::output and hdf5""",
        )
        create_option(
            d,
            "with_timer_output_at_the_end",
            "boolean, optional",
            True,
            """ if True, we store the timers into a dict and output it in hdf5 at the end
                      fastest method but need to complete the simulation""",
        )

        create_option(
            d,
            "explode_computeOneStep_in_python",
            "boolean, optional",
            False,
            """ if True, the ComputeOneStep function of siconos.TimeStepping is exploded in a python version
                      in order to add more log/trace during Newton loop. """,
        )

        create_option(
            d,
            "explode_computeOneStepNSProblem_in_python",
            "boolean, optional",
            False,
            """ if True, the ComputeOneStepNSProblem function of siconos.OneStepProblem is exploded in a python version
                      in order to add more log/trace during Newton loop. """,
        )
        create_option(
            d,
            "exit_tolerance",
            "boolean, optional",
            False,
            """ if True, the simulation run exits if the tolerance is not reached """,
        )

        super(self.__class__, self).__init__(d)

    def display(self):
        def print_comment(d_comment_item):
            if d_comment_item is None:
                print("  | no info on this option ")
            else:
                for k in d_comment_item.keys():
                    print("  |   {0}: {1}".format(k, d_comment_item[k]))

        print("display run options")
        print("{0} = {1}".format("option", "value"))
        for k in self.keys():
            print("{0} = {1}".format(k, self[k]))
            print_comment(self._d_comment.get(k))

    def check_valid_run_options(self):

        if self.get("explode_Newton_solve") is not None:
            msg = "run_options.check_valid_run_options():"
            msg += "explode_Newton_solve option is obsolete."
            msg += " Use instead explode_computeOneStep_in_python"
            raise RuntimeError(msg)

        if self.get("explode_computeOneStep") is not None:
            msg = "run_options.check_valid_run_options():"
            msg += "explode_computeOneStep option is obsolete."
            msg += "Use instead explode_computeOneStepNSProblem_in_python"
            raise RuntimeError(msg)

        for k in self.keys():
            if k not in self._valid_options_keys:
                msg = (
                    "run_options.check_valid_run_options() : the key "
                    + str(k)
                    + " in run_options dictionnary is not a valid key"
                )
                raise RuntimeError(msg)


class MechanicsHdf5Runner(siconos.io.mechanics_hdf5.MechanicsHdf5):
    """a Hdf5 context manager that reads the translations and
    orientations of collision objects from hdf5 file

    It provides functions to output translations and orientations in
    the same file during simulation (output is done by default in
    pos.dat)

    Parameters
    ----------
    io_filename: string, optional
         hdf5 file name, default = <caller>.hdf5, caller being the name
         without ext of the file that instanciates the Runner.
     mode: string, optional
         h5 mode (w, r, append), default = 'w'
     interaction_manager: SiconosCollisionManager, optional
         user-defined interaction handler (e.g. from Bullet), default=None
     nsds: siconos::modeling::NonSmoothDynamicalSystem, optional
         default = None
     simulation: Simulation, optional
         default = None
     osi: OneStepIntegrator, optional
         default = None
     shape_filename: string
         vtk file describing a mesh, default = None
     set_external_forces: python function, optional
         function used to apply forces onto the body.
         Must be :

         def funcname(body)

         body being a siconos Body (DS)
         Default : apply gravity forces.
     gravity_scale: real, optional
         multiplication factor for the gravity.
         1.     for meters (default).
         1./100 for centimeters.
         This parameter may be needed for small
         objects because of Bullet collision margin.
     collision_margin: real number, optional
         tolerance for collision, default = None (0.04 in Shape builder)
     use_compression: boolean, optional
         true to use compression for h5 file, default=False
     output_domains: boolean, optional
         if trueoutputs info regarding contact point domains
         default=False
     verbose: boolean, optional
        default=True

    """

    def __init__(
        self,
        config=None,
        io_filename=None,
        io_filename_backup=None,
        mode="w",
        interaction_manager=None,
        nsds=None,
        simulation=None,
        osi=None,
        shape_filename=None,
        set_external_forces=None,
        gravity_scale=None,
        collision_margin=None,
        use_compression=False,
        output_domains=False,
        verbose=True,
    ):
        super(MechanicsHdf5Runner, self).__init__(
            io_filename,
            mode,
            io_filename_backup,
            use_compression,
            output_domains,
            verbose,
        )
        if config is None:
            config = RunnerConfig()  # default backend is Bullet
        self.config = config
        self._interman = interaction_manager
        self._nsds = nsds
        self._simulation = simulation
        self._osi = osi
        self._osnspb = None
        self._static = {}
        self._shape = None
        self._occ_contactors = dict()
        self._io = sio.MechanicsIO()
        self._set_external_forces = set_external_forces
        self._shape_filename = shape_filename
        self._number_of_shapes = 0
        self._number_of_permanent_interactions = 0
        self._number_of_dynamic_objects = 0
        self._number_of_static_objects = 0
        self._gravity_scale = gravity_scale
        self._collision_margin = collision_margin
        self._output_frequency = 1
        self._output_backup_frequency = 1
        self._keep = []
        self._scheduled_births = []
        self._scheduled_deaths = []
        self._births = dict()
        self._deaths = dict()
        self._initializing = True
        self._output_contact_index_set = 1
        self._start_run_iteration_hook = None
        self._end_run_iteration_hook = None
        self._before_next_step_iteration_hook = None
        self._ds_positions = None
        self._ds_boundary_conditions = {}
        self._time_stepping_class = None
        self._k = None
        self._k0 = None
        self._q0 = []
        self._v0 = []
        self.fext = []
        self.weight = []
        self._timing = {}

        if self.config.backend == "vnative":
            self.get_io_array = lambda array: array
        else:
            self.get_io_array = lambda array: array.transpose()

    def __enter__(self):
        super(MechanicsHdf5Runner, self).__enter__()

        if self._gravity_scale is None:
            self._gravity_scale = 1  # 1 => m, 1/100. => cm

        if self._set_external_forces is None:
            self._set_external_forces = self.apply_gravity

        # -- Build shape collection and save it in self._shape
        if self._shape_filename is None:
            shape_collection_ref = self
        else:
            shape_collection_ref = self._shape_filename

        # Note FP: I can't find an example using vtk file for shape_filename and I can not
        # understant how it can work. To be reviewed.
        self._shape = siconos.io.shape_collection.ShapeCollection(
            io=shape_collection_ref,
            collision_margin=self._collision_margin,
            backend=self.config.backend,
        )
        self.print_verbose(
            f"Start Mechanics Runner - Backend: {self.config.backend} - Mode: {self._mode}\n"
        )
        return self

    def log(self, fun, with_timer=False, after=True):
        if with_timer:
            t = siconos.io.tools.Timer()

            def logged(*args):
                t.update()
                if not after:
                    print(
                        "[io.mechanics] |-->start {0:42s} ...".format(fun.__name__),
                        flush=True,
                    )

                output = fun(*args)
                endt = t.elapsed()

                if not after:
                    print(
                        "[io.mechanics] |-->end {0:44s} .... {1:6.2e} s".format(
                            fun.__name__, endt
                        ),
                        flush=True,
                    )
                else:
                    print(
                        "[io.mechanics] | {0:50s} .... {1:6.2e} s".format(
                            fun.__name__, endt
                        )
                    )

                # timing in hdf5
                if self._run_options["with_timer_output_at_the_end"]:
                    # we store timing in a dictionnary
                    if self._timing.get(fun.__name__) is None:
                        self._timing[fun.__name__] = [endt]
                    else:
                        self._timing[fun.__name__].append(endt)
                else:
                    # we store timing in hdf5
                    siconos.io.mechanics_hdf5.group(self.log_data(), fun.__name__)
                    siconos.io.mechanics_hdf5.add_line(
                        siconos.io.mechanics_hdf5.data(
                            self.log_data()[fun.__name__], "timing", 1
                        ),
                        endt,
                    )
                    if isinstance(output, numbers.Number):
                        siconos.io.mechanics_hdf5.add_line(
                            siconos.io.mechanics_hdf5.data(
                                self.log_data()[fun.__name__], "value", 1
                            ),
                            float(output),
                        )
                return output

            return logged
        else:

            def silent(*args):
                output = fun(*args)
                return output

            return silent

    def apply_gravity(self, body):
        g = scipy.constants.g / self._gravity_scale
        if hasattr(body, "scalarMass"):
            scalar_mass = body.scalarMass
        elif hasattr(body, "getMassValue"):
            scalar_mass = body.getMassValue()
        else:
            raise RuntimeError(f"Can't get mass value for this kind of DS {type(body)}")
        if self._dimension == 3:
            weight = [0, 0, -scalar_mass * g]
        elif self._dimension == 2:
            weight = [0, -scalar_mass * g, 0.0]
        self.weight.append(np.array(weight, dtype=np.float64))
        body.setConstantFext(self.weight[-1], siconos.modeling.alias_t)

    def import_nonsmooth_law(self, name):
        if self._interman is not None:
            nslawClass = getattr(sm, self._nslaws_data[name].attrs["type"])
            if nslawClass == sm.NewtonImpactFrictionNSL:
                nslaw = nslawClass(
                    float(self._nslaws_data[name].attrs["e"]),
                    0.0,
                    float(self._nslaws_data[name].attrs["mu"]),
                    self._dimension,
                )
            elif nslawClass == sm.FremondImpactFrictionNSL:
                nslaw = nslawClass(
                    float(self._nslaws_data[name].attrs["e"]),
                    0.0,
                    float(self._nslaws_data[name].attrs["mu"]),
                    self._dimension,
                )
            elif nslawClass == sm.NewtonImpactRollingFrictionNSL:
                if self._dimension == 3:
                    nslaw = nslawClass(
                        float(self._nslaws_data[name].attrs["e"]),
                        0.0,
                        float(self._nslaws_data[name].attrs["mu"]),
                        float(self._nslaws_data[name].attrs["mu_r"]),
                        5,
                    )
                elif self._dimension == 2:
                    nslaw = nslawClass(
                        float(self._nslaws_data[name].attrs["e"]),
                        0.0,
                        float(self._nslaws_data[name].attrs["mu"]),
                        float(self._nslaws_data[name].attrs["mu_r"]),
                        3,
                    )
            elif nslawClass == sm.NewtonImpactNSL:
                nslaw = nslawClass(float(self._nslaws_data[name].attrs["e"]))
            elif nslawClass == sm.RelayNSL:
                nslaw = nslawClass(
                    int(self._nslaws_data[name].attrs["size"]),
                    float(self._nslaws_data[name].attrs["lb"]),
                    float(self._nslaws_data[name].attrs["ub"]),
                )
            if not nslaw:
                raise AssertionError("no nslaw")
            # assert(nslaw)
            self._nslaws[name] = nslaw
            gid1 = int(self._nslaws_data[name].attrs["gid1"])
            gid2 = int(self._nslaws_data[name].attrs["gid2"])
            if gid1 >= 0 and gid2 >= 0:
                self._interman.insertNonSmoothLaw(nslaw, gid1, gid2)

    def import_native_object(
        self,
        name,
        translation,
        orientation,
        velocity,
        contactors,
        mass,
        given_inertia,
        body_class,
        shape_class,
        birth=False,
        number=None,
    ):
        if mass is None:
            # a static object
            # for native only plans
            self._static[name] = {
                "number": number,
                "origin": translation,
                "orientation": orientation,
            }

            ctor = contactors[0]

            bdy = None
            # only one contactor
            if self._shape.attributes(ctor.shape_name)["primitive"] == "Line":
                a = self._shape._io.shapes()[ctor.shape_name][:][0][0]
                b = self._shape._io.shapes()[ctor.shape_name][:][0][1]
                c = self._shape._io.shapes()[ctor.shape_name][:][0][2]

                bdy = self._interman.insertLine(a, b, c)
            elif self._shape.attributes(ctor.shape_name)["primitive"] == "Segment":
                x10 = self._shape._io.shapes()[ctor.shape_name][:][0][0]
                y10 = self._shape._io.shapes()[ctor.shape_name][:][0][1]
                x20 = self._shape._io.shapes()[ctor.shape_name][:][0][2]
                y20 = self._shape._io.shapes()[ctor.shape_name][:][0][3]

                alpha = orientation[0]  # 2D

                x1a, y1a = rotate_point(x10, y10, alpha)
                x2a, y2a = rotate_point(x20, y20, alpha)

                tx = translation[0]
                ty = translation[1]

                x1 = x1a + tx
                y1 = y1a + ty
                x2 = x2a + tx
                y2 = y2a + ty

                bdy = self._interman.insertSegment(x1, y1, x2, y2)
            elif self._shape.attributes(ctor.shape_name)["primitive"] == "Box2d":
                thickness = self._shape._io.shapes()[ctor.shape_name][:][0][0]
                size = self._shape._io.shapes()[ctor.shape_name][:][0][1]

                # insert 4 segments
                xa0 = -thickness / 2
                xb0 = thickness / 2
                ya0 = -size / 2
                yb0 = size / 2

                alpha = orientation[0]  # 2D

                x1a = xa0
                y1a = ya0
                x2a = xa0
                y2a = yb0
                x3a = xb0
                y3a = yb0
                x4a = xb0
                y4a = yb0

                x1r, y1r = rotate_point(x1a, y1a, alpha)
                x2r, y2r = rotate_point(x2a, y2a, alpha)
                x3r, y3r = rotate_point(x3a, y3a, alpha)
                x4r, y4r = rotate_point(x4a, y4a, alpha)

                tx = translation[0]
                ty = translation[1]

                x1 = x1r + tx
                y1 = y1r + ty
                x2 = x2r + tx
                y2 = y2r + ty
                x3 = x3r + tx
                y3 = y3r + ty
                x4 = x4r + tx
                y4 = y4r + ty

                s1 = self._interman.insertSegment(x1, y1, x2, y2)
                s2 = self._interman.insertSegment(x2, y2, x3, y3)
                s3 = self._interman.insertSegment(x3, y3, x4, y4)
                s4 = self._interman.insertSegment(x4, y4, x1, y1)

                bdy = [s1, s2, s3, s4]

            elif self._shape.attributes(ctor.shape_name)["primitive"] == "Disk":
                r = self._shape._io.shapes()[ctor.shape_name][:][0][0]
                bdy = self._interman.insertTranslatedDisk(r, translation)
            else:
                self.print_verbose(
                    "unknown primitive:{}".format(
                        self._shape.attributes(ctor.shape_name)["primitive"]
                    )
                )
                raise RuntimeError("unknown primitive")
            body = bdy
            flag = "static"
        else:
            # a dynamic object
            ctor = contactors[0]
            # shape = self._shape.get(ctor.shape_name)
            # attrs = self._shape.attributes(ctor.shape_name)
            initial_pos = np.concatenate([translation, orientation], axis=0)
            self._q0.append(initial_pos.copy())
            self._v0.append(velocity)
            if self.config.backend == "vnative":
                body_class = self.config.default_body_class
            else:
                if self._shape.attributes(ctor.shape_name)["primitive"] == "Disk":
                    body_class = siconos.mechanics.collision.Disk
                elif self._shape.attributes(ctor.shape_name)["primitive"] == "Circle":
                    body_class = siconos.mechanics.collision.Circle
            radius = self._shape._io.shapes()[ctor.shape_name][:][0][0]
            body = body_class(radius, mass, self._q0[-1], self._v0[-1])
            self._set_external_forces(body)
            self._nsds.insertDynamicalSystem(body)
            if birth and self._verbose:
                self.print_verbose(
                    "birth of body named {0}, translation {1}, orientation {2}".format(
                        name, translation, orientation
                    )
                )
            flag = "dynamic"

            if number is not None:
                body.setNumber(int(number))

        return body, flag

    def import_occ_object(
        self,
        name,
        translation,
        orientation,
        velocity,
        contactors,
        mass,
        given_inertia,
        body_class,
        shape_class,
        face_class,
        edge_class,
        birth=False,
        number=None,
    ):

        assert self.config.backend == "occ"

        if mass is None:
            # a static object
            body = None

            self._static[name] = {
                "number": number,
                "origin": translation,
                "orientation": orientation,
            }
            flag = "static"
        else:
            if not np.isscalar(mass) or mass <= 0:
                raise RuntimeError("Warning mass must be a positive scalar")

            if body_class is None:
                body_class = self.config.default_body_class

            # assert (given_inertia is not None)
            if given_inertia is None:
                raise AssertionError("given_inertia is  None")
            else:
                if np.shape(given_inertia) == (3,):
                    given_inertia = np.diag(given_inertia)
                elif np.shape(given_inertia) != (3, 3):
                    self.print_verbose("Wrong shape of inertia")
                inertia = np.asfortranarray(given_inertia, dtype=np.float64)

            body = body_class(
                np.concatenate([translation, orientation], axis=0),
                velocity,
                mass,
                inertia,
            )

            if number is not None:
                body.setNumber(int(number))
            flag = "dynamic"

        ref_shape = {
            ctor.instance_name: self.config.occ.OccContactShape(
                self._shape.get(
                    ctor.shape_name,
                    shape_class,
                    face_class,
                    edge_class,
                    new_instance=True,
                )
            )
            for ctor in contactors
        }

        ref_added = dict()
        for contactor in contactors:
            contact_shape = None
            reference_shape = ref_shape[contactor.instance_name]
            self._keep.append(reference_shape)

            if hasattr(contactor, "contact_type"):
                if contactor.contact_type == "Face":
                    contact_shape = self.config.occ.OccContactFace(
                        reference_shape, contactor.contact_index
                    )

                elif contactor.contact_type == "Edge":
                    contact_shape = self.config.occ.OccContactEdge(
                        reference_shape, contactor.contact_index
                    )

            if contact_shape is not None:
                if name not in self._occ_contactors:
                    self._occ_contactors[name] = dict()

                self._occ_contactors[name][contactor.instance_name] = contact_shape

                if body is not None:
                    body.addContactShape(
                        contact_shape,
                        contactor.translation,
                        contactor.orientation,
                        contactor.collision_group,
                    )

            if reference_shape not in ref_added:
                if body is not None:
                    body.addShape(
                        reference_shape,
                        contactor.translation,
                        contactor.orientation,
                    )
                else:

                    self.config.occ.occ_move(
                        reference_shape,
                        list(contactor.translation) + list(contactor.orientation),
                    )
                ref_added[reference_shape] = True

        if body is not None:
            self._set_external_forces(body)

            # add the dynamical system to the non smooth
            # dynamical system
            self._nsds.insertDynamicalSystem(body)
            self._nsds.setName(body, str(name))

            if birth and self._verbose:
                self.print_verbose(
                    "birth of body named {0}, translation {1}, orientation {2}".format(
                        name, translation, orientation
                    )
                )

        return body, flag

    def import_bullet_object(
        self,
        name,
        translation,
        orientation,
        velocity,
        contactors,
        mass,
        inertia,
        body_class,
        shape_class,
        birth=False,
        number=None,
    ):
        if body_class is None:
            body_class = self.config.default_body_class

        if self._dimension == 2:
            body_class = siconos.mechanics.collision.RigidBody2dDS

        if self._interman is not None and "input" in self._data:
            body = None
            # ---------------
            # a static object
            # ---------------
            if mass is None:
                cset = siconos.mechanics.collision.SiconosContactorSet()
                csetpos = np.concatenate([translation, orientation], axis=0)
                for c in contactors:
                    shp = self._shape.get(c.shape_name)
                    if type(shp) is siconos.io.shape_collection.NativeSegmentShape:
                        # segment -> box2d
                        # to ==> in mechanics_hdf5 at the declaration level
                        x1 = shp.params[0]
                        y1 = shp.params[1]
                        x2 = shp.params[2]
                        y2 = shp.params[3]
                        d = sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)
                        shp = siconos.mechanics.collision.SiconosBox2d(d, 0.0)
                        shp.setInsideMargin(0.0)
                        shp.setOutsideMargin(0.01)
                        trans = [(x1 + x2) / 2, (y1 + y2) / 2, 0.0]
                        c.translation = trans
                        if x2 != x1:
                            orien = atan((y2 - y1) / (x2 - x1))
                        else:
                            orien = pi / 2
                        c.orientation = [cos(orien / 2), 0, 0, sin(orien / 2)]
                        pos = np.concatenate([c.translation, c.orientation], axis=0)
                        cset.append(
                            siconos.mechanics.collision.SiconosContactor(
                                shp, pos, c.collision_group
                            )
                        )
                    else:
                        pos = np.concatenate([c.translation, c.orientation], axis=0)
                        cset.append(
                            siconos.mechanics.collision.SiconosContactor(
                                shp, pos, c.collision_group
                            )
                        )
                        self.print_verbose(
                            "              Adding shape %s to static contactor"
                            % c.shape_name,
                            "at relative position",
                            pos,
                        )

                staticBody = self._interman.addStaticBody(cset, csetpos, number)

                self._static[name] = {
                    "number": number,
                    "origin": translation,
                    "orientation": orientation,
                    "shape": shp,
                }
                # In the case of a static object, we return the staticBody and a flag
                # this will be used to remove the contactor if we have a time of death
                return staticBody, "static"

            # ---------------
            # a dynamic object
            # ---------------
            else:
                if not np.isscalar(mass) or mass <= 0:
                    self.print_verbose("Warning mass must be a positive scalar")
                inertia_ok = False
                initial_pos = np.concatenate([translation, orientation], axis=0)
                if self._dimension == 3:
                    if inertia is not None:
                        if np.shape(inertia) == (3,):
                            inertia = np.diag(inertia)
                        elif np.shape(inertia) != (3, 3):
                            self.print_verbose("Wrong shape of inertia")
                        inertia = np.asfortranarray(inertia, dtype=np.float64)
                        inertia_ok = True

                elif self._dimension == 2:
                    if inertia is not None:
                        if (
                            np.shape(inertia) == (1, 1)
                            or np.shape(inertia) == (1,)
                            or np.isscalar(inertia)
                        ):
                            inertia_ok = True

                if inertia_ok:
                    if not inertia.flags["F_CONTIGUOUS"]:
                        inertia = inertia.copy(order="F")
                    assert inertia.flags["F_CONTIGUOUS"]
                    # create the dynamics object with mass and inertia
                    self._q0.append(initial_pos.copy())
                    self._v0.append(velocity)
                    body = body_class(self._q0[-1], self._v0[-1], mass, inertia)
                    body.setUseContactorInertia(False)
                else:
                    if inertia is not None:
                        if self._dimension == 3:
                            self.print_verbose(
                                "**** Warning inertia for object named {0} does not have the"
                                " correct shape: {1} instead of (3, 3) or (3,)".format(
                                    name, np.shape(inertia)
                                )
                            )
                            self.print_verbose(
                                "**** Inertia will be computed with the shape"
                                " of the first contactor"
                            )
                        elif self._dimension == 2:
                            self.print_verbose(
                                "**** Warning inertia for object named {0} does not have the"
                                " correct shape: {1} instead of (1, 1) or (1,)"
                                "or scalar".format(name, np.shape(inertia))
                            )
                            self.print_verbose(
                                "**** Inertia will be computed with the shape of"
                                " the first contactor"
                            )
                    if self._dimension == 3:
                        inertia = np.zeros((3, 3), order="F", dtype=np.float64)
                    elif self._dimension == 2:
                        inertia = np.zeros((1, 1), order="F", dtype=np.float64)
                    np.fill_diagonal(inertia, 1)
                    self._q0.append(initial_pos.copy())
                    self._v0.append(velocity)
                    body = body_class(self._q0[-1], self._v0[-1], mass, inertia)
                    body.setUseContactorInertia(True)

                self.fext.append(self._input[name].get("allow_self_collide", None))
                if self.fext[-1] is not None:
                    body.setConstantFext(self.fext[-1], siconos.modeling.alias_t)

                self_collide = self._input[name].get("allow_self_collide", None)
                if self_collide is not None:
                    body.setAllowSelfCollide(not not self_collide)

                cset = siconos.mechanics.collision.SiconosContactorSet()
                for c in contactors:
                    shp = self._shape.get(c.shape_name)
                    pos = np.concatenate([c.translation, c.orientation], axis=0)
                    cset.append(
                        siconos.mechanics.collision.SiconosContactor(shp, pos, c.group)
                    )
                    self.print_verbose(
                        "              Adding shape %s to dynamic contactor"
                        % c.shape_name,
                        "at relative position",
                        pos,
                    )

                body.setContactors(cset)

            if body:
                # set id number
                if number is not None:
                    body.setNumber(int(number))

                # set external forces
                self._set_external_forces(body)

                # add the dynamical system to the non smooth
                # dynamical system
                self._nsds.insertDynamicalSystem(body)
                self._nsds.setName(body, str(name))
            return body, "dynamic"

    def make_coupler_jointr(self, ds1_name, ds2_name, coupled, references):
        topo = self._nsds.topology()
        dof1, dof2, ratio = coupled[0, :]
        refds_name = None
        if len(references) > 0:
            joint1_name = references[0]
            joint1 = topo.getInteraction(str(joint1_name)).relation()
            joint2 = joint1
            joint1_ds1 = self.joints()[joint1_name].attrs["object1"]
            joint1_ds2 = self.joints()[joint1_name].attrs.get("object2", None)
        if len(references) > 1:
            joint2_name = references[1]
            joint2 = topo.getInteraction(str(joint2_name)).relation()
            joint2_ds1 = self.joints()[joint2_name].attrs["object1"]
            joint2_ds2 = self.joints()[joint2_name].attrs.get("object2", None)

            if len(references) > 2:
                refds_name = references[2]
            else:
                # if second joint provided but no reference, then
                # infer refds_name from the reference joints
                dss = set([joint1_ds1, joint1_ds2, joint2_ds1, joint2_ds2])
                diff = list(dss.difference(set([ds1_name, ds2_name])))
                # there must be exactly one reference in common that
                # is not either of the DSs
                refds_name = diff[0]

        if refds_name:
            refds = topo.getDynamicalSystem(str(refds_name))

            # Determine reference indexes:
            # Assert if neither ds in reference joints is the
            # ref ds and that the other ds is ds1/ds2 as
            # appropriate.
            refidx1, refidx2 = 0, 0
            if joint1_ds1 == ds1_name:
                refidx1 = 1
                if joint1_ds2 != refds_name:
                    raise AssertionError("joint1_ds2 != refds_name")
                # assert(joint1_ds2 == refds_name)
            else:
                if joint1_ds1 != refds_name:
                    raise AssertionError("joint1_ds1 != refds_name")
                # assert(joint1_ds1 == refds_name)
                if joint1_ds2 != ds1_name:
                    raise AssertionError("joint1_ds2 != ds1_name")
                # assert(joint1_ds2 == ds1_name)
            if joint2_ds1 == ds2_name:
                refidx2 = 1
                if joint2_ds2 != refds_name:
                    raise AssertionError("joint2_ds2 != refds_name")
                # assert(joint2_ds2 == refds_name)
            else:
                # assert(joint2_ds1 == refds_name)
                if joint2_ds1 != refds_name:
                    raise AssertionError("joint2_ds1 != refds_name")
                # assert(joint2_ds2 == ds2_name)
                if joint2_ds2 != ds2_name:
                    raise AssertionError("joint2_ds2 != ds1_name")

            joint = siconos.mechanics.joints.CouplerJointR(
                joint1,
                int(dof1),
                joint2,
                int(dof2),
                ratio,
                refds,
                refidx1,
                refds,
                refidx2,
            )
        else:
            # otherwise we assume no reference, coupler is directly
            # between ds1 and ds2
            joint = siconos.mechanics.joints.CouplerJointR(
                joint1, int(dof1), joint2, int(dof2), ratio
            )
        return joint

    def import_joint(self, name):
        if self._interman is None:
            self.print_verbose(
                "Try to import joints but no interaction manager has been defined."
                " Nothing will be done."
            )
            return

        nsds = self._nsds
        topo = nsds.topology()

        # Read hdf file: ["data"]["joints"][name]
        joint_type = self.joints()[name].attrs["type"]
        joint_class = getattr(siconos.mechanics.joints, joint_type)
        absolute = bool(self.joints()[name].attrs.get("absolute", True))
        allow_self_collide = self.joints()[name].attrs.get("allow_self_collide", None)
        stops = self.joints()[name].attrs.get("stops", None)
        nslaws = self.joints()[name].attrs.get("nslaws", None)
        friction = self.joints()[name].attrs.get("friction", None)
        coupled = self.joints()[name].attrs.get("coupled", None)
        references = self.joints()[name].attrs.get("references", None)

        # work-around h5py unicode bug
        # https://github.com/h5py/h5py/issues/379
        if references is not None:
            references = [
                r.decode("utf-8") if isinstance(r, bytes) else r for r in references
            ]

        points = self.joints()[name].attrs.get("points", [])
        axes = self.joints()[name].attrs.get("axes", [])
        siconos.io.mechanics_hdf5.check_points_axes(name, joint_class, points, axes)

        ds1_name = self.joints()[name].attrs["object1"]
        ds1 = topo.getDynamicalSystem(ds1_name)
        ds2 = None

        if "object2" in self.joints()[name].attrs:
            ds2_name = self.joints()[name].attrs["object2"]
            ds2 = topo.getDynamicalSystem(ds2_name)

        if joint_class == siconos.mechanics.joints.CouplerJointR:
            # This case is a little different, handle it specially
            # assert(references is not None)
            # assert(np.shape(coupled) == (1, 3))
            if references is None:
                raise AssertionError("references is None")
            if np.shape(coupled) != (1, 3):
                raise AssertionError("np.shape(coupled) != (1, 3)")

            joint = self.make_coupler_jointr(ds1_name, ds2_name, coupled, references)
            coupled = None  # Skip the use for "coupled" below, to
            # install joint-local couplings

        else:
            # Generic NewtonEulerJointR interface
            joint = joint_class()
            for n, p in enumerate(points):
                p = np.asarray(p, dtype=np.float64)
                joint.setPoint(n, p)
            for n, a in enumerate(axes):
                a = np.asarray(a, dtype=np.float64)
                joint.setAxis(n, a)
            joint.setAbsolute(absolute)

        q1 = ds1.q()
        q2 = None if ds2 is None else ds2.q()

        joint.setBasePositions(q1, q2)

        if allow_self_collide is not None:
            joint.setAllowSelfCollide(not not allow_self_collide)
        joint_nslaw = sm.EqualityConditionNSL(joint.numberOfConstraints())
        joint_inter = sm.Interaction(joint_nslaw, joint)
        self._nsds.link(joint_inter, ds1, ds2)
        nsds.setName(joint_inter, str(name))

        # Add a e=0 joint by default, otherwise user can specify
        # the impact law by name or a list of names for each axis.
        if stops is not None:
            if np.shape(stops)[1] != 3:
                raise AssertionError("np.shape(stops)[1] != 3")
            # assert np.shape(stops)[1] == 3, 'Joint stops shape must be (?, 3)'
            if nslaws is None:
                nslaws = [sm.NewtonImpactNSL(0.0)] * np.shape(stops)[0]
            elif isinstance(nslaws, bytes):
                nslaws = [self._nslaws[nslaws.decode("utf-8")]] * np.shape(stops)[0]
            elif isinstance(nslaws, str):
                nslaws = [self._nslaws[nslaws]] * np.shape(stops)[0]
            else:
                if np.shape(nslaws)[0] != np.shape(stops)[0]:
                    raise AssertionError("np.shape(nslaws)[0] != np.shape(stops)[0]")
                # assert(np.shape(nslaws)[0] == np.shape(stops)[0])
                nslaws = [self._nslaws[nsl] for nsl in nslaws]
            for n, (nsl, (axis, pos, _dir)) in enumerate(zip(nslaws, stops)):
                # "bool()" is needed because type of _dir is
                # numpy.bool_, which SWIG doesn't handle well.
                stop = siconos.mechanics.joints.JointStopR(
                    joint, pos, bool(_dir < 0), int(axis)
                )
                stop_inter = sm.Interaction(nsl, stop)
                self._nsds.link(stop_inter, ds1, ds2)
                nsds.setName(stop_inter, "%s_stop%d" % (str(name), n))

        # The per-axis friction NSL, can be ''
        if friction is not None:
            if isinstance(friction, str):
                friction = [friction]
            elif isinstance(friction, bytes):
                friction = [friction.decode("utf-8")]
            else:
                friction = [
                    (f.decode("utf-8") if isinstance(f, bytes) else f) for f in friction
                ]
            for ax, fr_nslaw in enumerate(friction):
                if fr_nslaw == "":  # no None in hdf5, use empty string
                    continue  # instead for no NSL on an axis
                nslaw = self._nslaws[fr_nslaw]
                fr = siconos.mechanics.joints.JointFrictionR(joint, ax)
                fr_inter = sm.Interaction(nslaw, fr)
                self._nsds.link(fr_inter, ds1, ds2)
                nsds.setName(fr_inter, "%s_friction%d" % (str(name), ax))

        # An array of tuples (dof1, dof2, ratio) specifies
        # coupling between a joint's DoFs (e.g., to turn a
        # cylindrical joint into a screw joint)
        if coupled is not None:
            if len(coupled.shape) == 1:
                coupled = np.array([coupled])
            if coupled.shape[1] != 3:
                raise AssertionError("coupled.shape[1] != 3")
            # assert(coupled.shape[1] == 3)
            for n, (dof1, dof2, ratio) in enumerate(coupled):
                cpl = siconos.mechanics.joints.CouplerJointR(
                    joint, int(dof1), joint, int(dof2), ratio
                )
                cpl.setBasePositions(q1, q2)
                cpl_inter = sm.Interaction(sm.EqualityConditionNSL(1), cpl)
                self._nsds.link(cpl_inter, ds1, ds2)
                nsds.setName(cpl_inter, "%s_coupler%d" % (str(name), n))

    def import_boundary_conditions(self, name):
        if self._interman is not None:
            topo = self._nsds.topology()

            bc_type = self.boundary_conditions()[name].attrs["type"]
            bc_class = getattr(sm, bc_type)

            ds1_name = self.boundary_conditions()[name].attrs["object1"]
            ds1 = topo.getDynamicalSystem(ds1_name)

            if bc_type == "HarmonicBC":
                bc = bc_class(
                    self.boundary_conditions()[name].attrs["indices"],
                    self.boundary_conditions()[name].attrs["a"],
                    self.boundary_conditions()[name].attrs["b"],
                    self.boundary_conditions()[name].attrs["omega"],
                    self.boundary_conditions()[name].attrs["phi"],
                )

            elif bc_type == "BoundaryCondition":
                velocity_values = self.boundary_conditions()[name].attrs.get("v", None)
                if velocity_values is None:
                    # fixed bc with zero prescribed velocities
                    bc = bc_class(self.boundary_conditions()[name].attrs["indices"])
                else:
                    # fixed bc with zero prescribed velocities
                    bc = bc_class(
                        self.boundary_conditions()[name].attrs["indices"],
                        velocity_values,
                    )

            # set bc to the ds1

            ds1.setBoundaryConditions(bc)

            # joint_inter = sm.Interaction(joint_nslaw, joint)
            #    self._nsds.\
            #        link(joint_inter, ds1)

    def import_permanent_interactions(self, name):
        """ """
        if (
            self._interman is not None
            and "input" in self._data
            and self.permanent_interactions() is not None
        ):
            topo = self._nsds.topology()
            pinter = self.permanent_interactions()[name]
            body1_name = pinter.attrs["body1_name"]
            body2_name = pinter.attrs["body2_name"]

            try:
                ds1 = topo.getDynamicalSystem(body1_name)
            except Exception:
                ds1 = None

            try:
                ds2 = topo.getDynamicalSystem(body2_name)
            except Exception:
                ds2 = None

            # static object in second
            if ds1 is None:
                ds2, ds1 = ds1, ds2

            contactor1_name = pinter.attrs.get("contactor1_name", None)
            contactor2_name = pinter.attrs.get("contactor2_name", None)

            if contactor1_name is None:
                contactor1_names = self._occ_contactors[body1_name].keys()
            else:
                contactor1_names = [contactor1_name]

            if contactor2_name is None:
                contactor2_names = self._occ_contactors[body2_name].keys()
            else:
                contactor2_names = [contactor2_name]

            distance_calculator = pinter.attrs["distance_calculator"]
            offset1 = pinter.attrs["offset1"]
            offset2 = pinter.attrs["offset2"]

            body1 = self._input[body1_name]
            body2 = self._input[body2_name]

            real_dist_calc = {
                "cadmbtb": self.config.occ.CadmbtbDistanceType,
                "occ": self.config.occ.OccDistanceType,
            }

            for contactor1_name in contactor1_names:
                for contactor2_name in contactor2_names:
                    ctr1 = body1[contactor1_name]
                    ctr2 = body2[contactor2_name]

                    cg1 = int(ctr1.attrs["group"])
                    cg2 = int(ctr2.attrs["group"])
                    nslaw = self._interman.nonSmoothLaw(cg1, cg2)

                    cocs1 = self._occ_contactors[body1_name][contactor1_name]
                    cocs2 = self._occ_contactors[body2_name][contactor2_name]

                    if ds2 is None:
                        self.print_verbose(
                            "moving contactor {0} of static object {1} to {2}".format(
                                contactor2_name,
                                body2_name,
                                list(
                                    np.array(body2.attrs["translation"])
                                    + np.array(ctr2.attrs["translation"])
                                )
                                + list(
                                    siconos.mechanics.quaternions.quaternion_multiply(
                                        ctr2.attrs["orientation"],
                                        body2.attrs["orientation"],
                                    )
                                ),
                            )
                        )
                        self.config.occ.occ_move(
                            cocs2,
                            list(
                                np.array(body2.attrs["translation"])
                                + np.array(ctr2.attrs["translation"])
                            )
                            + list(
                                siconos.mechanics.quaternions.quaternion_multiply(
                                    ctr2.attrs["orientation"],
                                    body2.attrs["orientation"],
                                )
                            ),
                        )

                    cp1 = self.config.occ.ContactPoint(cocs1)
                    cp2 = self.config.occ.ContactPoint(cocs2)

                    relation = self.config.occ.OccR(
                        cp1, cp2, real_dist_calc[distance_calculator]()
                    )

                    relation.setOffset1(offset1)
                    relation.setOffset2(offset2)

                    inter = sm.Interaction(nslaw, relation)

                    if ds2 is not None:
                        self._nsds.link(inter, ds1, ds2)
                    else:
                        self._nsds.link(inter, ds1)

                    # keep pointers
                    self._keep.append([cocs1, cocs2, cp1, cp2, relation])

    def import_object(
        self,
        name,
        body_class=None,
        shape_class=None,
        face_class=None,
        edge_class=None,
        birth=False,
        translation=None,
        orientation=None,
        velocity=None,
    ):
        """Import an object by name,
        possibly overriding initial position and velocity.
        """
        obj = self._input[name]
        self.print_verbose("Import object name:", name)
        self.print_verbose("              number (id): {0} ".format(obj.attrs["id"]))
        if translation is None:
            translation = obj.attrs["translation"]
        if orientation is None:
            orientation = obj.attrs["orientation"]
        if velocity is None:
            velocity = obj.attrs["velocity"]

        velocity = np.asarray(velocity, dtype=np.float64)
        translation = np.asarray(translation, dtype=np.float64)
        orientation = np.asarray(orientation, dtype=np.float64)

        # bodyframe center of mass
        center_of_mass = np.asarray(
            obj.attrs.get("center_of_mass", [0, 0, 0]), dtype=np.float64
        )

        mass = obj.attrs.get("mass", None)
        inertia = obj.attrs.get("inertia", None)

        if mass is None:
            self.print_verbose("              static object")
            self.print_verbose(
                "              position", list(translation) + list(orientation)
            )
        else:
            self.print_verbose("              dynamic object")
            self.print_verbose(
                "              position", list(translation) + list(orientation)
            )
            self.print_verbose("              velocity", velocity)

        # self.print_verbose('mass = ', mass)
        # self.print_verbose('inertia = ', inertia)

        input_ctrs = [ctr for _n_, ctr in obj.items()]

        contactors = []
        occ_type = False

        for ctr in input_ctrs:
            if "type" in ctr.attrs:
                # occ contact
                occ_type = True
                contactors.append(
                    siconos.mechanics.collision.tools.Contactor(
                        instance_name=ctr.attrs["instance_name"],
                        shape_name=ctr.attrs["shape_name"],
                        collision_group=ctr.attrs["group"].astype(int),
                        contact_type=ctr.attrs["type"],
                        contact_index=ctr.attrs["contact_index"].astype(int),
                        relative_translation=np.subtract(
                            ctr.attrs["translation"].astype(float), center_of_mass
                        ),
                        relative_orientation=ctr.attrs["orientation"].astype(float),
                    )
                )
            elif "group" in ctr.attrs:
                # bullet contact
                if occ_type:
                    raise AssertionError("occ type is found")
                # assert not occ_type
                contactors.append(
                    siconos.mechanics.collision.tools.Contactor(
                        instance_name=ctr.attrs["instance_name"],
                        shape_name=ctr.attrs["shape_name"],
                        collision_group=ctr.attrs["group"].astype(int),
                        relative_translation=np.subtract(
                            ctr.attrs["translation"].astype(float), center_of_mass
                        ),
                        relative_orientation=ctr.attrs["orientation"].astype(float),
                    )
                )
            else:
                # occ shape
                occ_type = True
                # fix: not only contactors here
                contactors.append(
                    siconos.mechanics.collision.tools.Shape(
                        instance_name=ctr.attrs["instance_name"],
                        shape_name=ctr.attrs["shape_name"],
                        relative_translation=np.subtract(
                            ctr.attrs["translation"].astype(float), center_of_mass
                        ),
                        relative_orientation=ctr.attrs["orientation"].astype(float),
                    )
                )

        if occ_type:
            # Occ object
            body, flag = self.import_occ_object(
                name,
                translation,
                orientation,
                velocity,
                contactors,
                mass,
                inertia,
                body_class,
                shape_class,
                face_class,
                edge_class,
                birth=birth,
                number=self.instances()[name].attrs["id"],
            )
        elif self.config.backend == "native" or self.config.backend == "vnative":
            body, flag = self.import_native_object(
                name,
                translation,
                orientation,
                velocity,
                contactors,
                mass,
                inertia,
                body_class,
                shape_class,
                birth=birth,
                number=self.instances()[name].attrs["id"],
            )
        else:
            # Bullet object
            body, flag = self.import_bullet_object(
                name,
                translation,
                orientation,
                velocity,
                contactors,
                mass,
                inertia,
                body_class,
                shape_class,
                birth=birth,
                number=self.instances()[name].attrs["id"],
            )

        # import boundary conditions
        bc_name = self._ds_boundary_conditions.get(name, None)
        if bc_name is not None:
            self.import_boundary_conditions(bc_name)

        # schedule its death immediately
        time_of_death = obj.attrs.get("time_of_death", None)

        if time_of_death is not None:
            bisect.insort_left(self._scheduled_deaths, time_of_death)
            if time_of_death in self._deaths:
                self._deaths[time_of_death].append((name, obj, body, flag))
            else:
                self._deaths[time_of_death] = [(name, obj, body, flag)]

    def import_objects(self, name, body_class=None, shape_class=None,
                       face_class=None, edge_class=None, birth=False):
        """
        Import several objects at once.
        """

        # This is implemented only for the vnative backend
        assert self.config.backend == "vnative"

        obj = self._input[name]

        translations = np.asarray(obj["translations"][:], dtype=np.float64)
        orientations = np.asarray(obj["orientations"][:], dtype=np.float64)
        positions = np.concatenate([translations, orientations], axis=1)
        velocities = np.asarray(obj["velocities"][:], dtype=np.float64)
        mass = obj.attrs.get("mass", None)
        inertia = obj.attrs.get("inertia", None)
        base_id = obj.attrs["id"]
        count = obj.attrs["count"]

        # Create contactor (same for all grains)
        import siconos.mechanics.collision.tools as smct_tools
        contactor = smct_tools.Contactor(
            shape_name=obj.attrs["shape_name"],
            collision_group=obj.attrs.get("group", 0)
        )

        self.print_verbose(f"Importing aggregate {name} with {count} grains...")

        # For vnative backend, create bodies efficiently
        ctor = contactor
        radius = self._shape._io.shapes()[ctor.shape_name][:][0][0]

        # Only dynamic objects
        assert mass is not None

        bodies = self.config.default_bodies_class(radius, mass, positions, velocities)
        np.vectorize(self._nsds.insertDynamicalSystem)(bodies.get())
        np.vectorize(self._set_external_forces)(bodies.get())

        self.print_verbose(f"Finished importing aggregate {name}")
        return bodies

    def import_scene(self, time, body_class, shape_class, face_class, edge_class):
        """
        From the specification given in the hdf5 file with the help of
        add* functions, import into the NSDS:
          - the static objects
          - the dynamic objects
          - the joints
        and into the interaction_manager:
          - the nonsmooth laws
        that have a specified time of birth <= current time.
        """

        # Ensure we count up from zero for implicit DS numbering
        if self.config.backend != "vnative":
            sm.DynamicalSystem.resetCount(0)

        for _ in self._ref:
            self._number_of_shapes += 1

        # import dynamical systems
        if self._interman is not None and "input" in self._data:
            # We map the boundary conditions with the object into
            #  self._ds_boundary_conditions.
            # In that way, boundary conditions will be imported once the
            # object is created in import_object and not necessarily at the
            # beginning of the simulation.
            for name in self.boundary_conditions():
                ds1_name = self.boundary_conditions()[name].attrs["object1"]
                self._ds_boundary_conditions[ds1_name] = name
                # only one boundary condition per object. list is not needed
                # ds_bc = self._ds_boundary_conditions.get(ds1_name, None)
                # if ds_bc is None:
                #     self._ds_boundary_conditions[ds1_name] = [name]
                # else:
                #     self._ds_boundary_conditions[ds1_name].append(name)

            # get pointers
            dpos_data = self.dynamic_data()
            velocities = self.velocities_data()

            # some dict to prefetch values in order to
            # speedup cold start in the case of many objects
            xdpos_data = dict()
            xvelocities = dict()

            if dpos_data is not None and len(dpos_data) > 0:
                max_time = max(dpos_data[:, 0])
                id_last = np.where(abs(dpos_data[:, 0] - max_time) < 1e-9)[0]
                id_vlast = np.where(abs(velocities[:, 0] - max_time) < 1e-9)[0]

                # prefetch positions and velocities in dictionaries
                # this avoids calls to np.where for each object
                for ldpos_data in dpos_data[id_last, :]:
                    xdpos_data[ldpos_data[1]] = ldpos_data
                for lvelocities in velocities[id_vlast, :]:
                    xvelocities[lvelocities[1]] = lvelocities

            else:
                # should not be used
                max_time = None
                id_last = None
            self.print_verbose("import dynamical systems ...")
            for name, obj in sorted(self._input.items(), key=lambda x: x[0]):
                 # Check if this is an aggregate object first
                obj_type = obj.attrs.get("type", "")

                if obj_type in ("dynamic_aggregate", "static_aggregate"):
                    # Import aggregate as individual bodies
                    self.import_objects(name, body_class, shape_class, face_class, edge_class)
                    continue

                mass = obj.attrs.get("mass", None)
                time_of_birth = obj.attrs.get("time_of_birth", -1)
                time_of_death = obj.attrs.get("time_of_death", float("inf"))

                if time_of_birth >= time:
                    #
                    # in the future
                    #
                    bisect.insort_left(self._scheduled_births, time_of_birth)
                    if time_of_birth in self._births:
                        self._births[time_of_birth].append((name, obj))
                    else:
                        self._births[time_of_birth] = [(name, obj)]
                elif time_of_death <= time:
                    # object already dead do not import
                    self.print_verbose("object", name, "already dead do not import")
                    # input()
                    # pass
                else:
                    #
                    # this is for now
                    #
                    # cold restart if output previously done

                    if (
                        mass is not None
                        and dpos_data is not None
                        and len(dpos_data) > 0
                    ):
                        xpos = xdpos_data[obj.attrs["id"]]
                        translation = (xpos[2], xpos[3], xpos[4])
                        orientation = (xpos[5], xpos[6], xpos[7], xpos[8])
                        xvel = xvelocities[obj.attrs["id"]]
                        velocity = (
                            xvel[2],
                            xvel[3],
                            xvel[4],
                            xvel[5],
                            xvel[6],
                            xvel[7],
                        )
                        if self._dimension == 2:
                            angle = 2.0 * acos(translation[2])
                            orientation = (angle,)
                            translation = (translation[0], translation[1])
                            velocity = (velocity[0], velocity[1], velocity[5])

                    else:
                        # start from initial conditions
                        translation = obj.attrs["translation"]
                        orientation = obj.attrs["orientation"]
                        velocity = obj.attrs["velocity"]

                    self.import_object(
                        name=name,
                        body_class=body_class,
                        shape_class=shape_class,
                        face_class=face_class,
                        edge_class=edge_class,
                        translation=translation,
                        orientation=orientation,
                        velocity=velocity,
                        birth=False,
                    )

            # import nslaws
            # note: no time of birth for nslaws and joints
            self.print_verbose("import nslaws ...")
            for name in self._nslaws_data:
                self.import_nonsmooth_law(name)

            # As for the  boundary conditions, the joints should be imported
            # after the importation of the object, or the two objects. If one
            # of the object has a time of birth this will fail.

            self.print_verbose("import joints ...")
            for name in self.joints():
                self.import_joint(name)

            self.print_verbose("import permanent interactions ...")
            for name in self.permanent_interactions():
                self.import_permanent_interactions(name)

    def current_time(self):
        if self._initializing:
            return self._simulation.startingTime()
        else:
            return self._simulation.nextTime()

    def import_births(
        self,
        body_class=None,
        shape_class=None,
        face_class=None,
        edge_class=None,
    ):
        """
        Import new objects into the NSDS.
        """
        time = self.current_time()

        ind_time = bisect.bisect_right(self._scheduled_births, time)

        current_times_of_births = set(self._scheduled_births[:ind_time])
        self._scheduled_births = self._scheduled_births[ind_time:]

        for time_of_birth in current_times_of_births:
            for name, _ in self._births[time_of_birth]:
                self.import_object(
                    name, body_class, shape_class, face_class, edge_class, birth=True
                )

    def execute_deaths(self):
        """
        Remove objects from the NSDS
        """
        time = self.current_time()

        ind_time = bisect.bisect_right(self._scheduled_deaths, time)
        # print('ind_time', ind_time)
        # print('self._scheduled_deaths', self._scheduled_deaths )
        # print('self._deaths', self._deaths )

        current_times_of_deaths = set(self._scheduled_deaths[:ind_time])
        # print('current_times_of_deaths', current_times_of_deaths )
        self._scheduled_deaths = self._scheduled_deaths[ind_time:]
        # print(self._scheduled_deaths)
        for time_of_death in current_times_of_deaths:
            # print(self._deaths[time_of_death])
            for _, _, body, flag in self._deaths[time_of_death]:
                if flag == "static":
                    self._interman.removeStaticBody(body)
                elif flag == "dynamic":
                    self._interman.removeBody(body)
                    self._nsds.removeDynamicalSystem(body)
                else:
                    msg = "execute_deaths : unknown object type"
                    msg += "It should static or dynamic"
                    raise RuntimeError(msg)
        # input()

    def output_static_objects(self):
        """
        Outputs translations and orientations of static objects
        """
        time = self.current_time()
        p = 0

        # append new static position
        current_line = self._static_data.shape[0]
        self._static_data.resize(current_line + len(self._static), 0)

        p = current_line
        for static in self._static.values():
            translation = static["origin"]
            rotation = static["orientation"]
            if self._dimension == 3:
                self._static_data[p, :] = [
                    time,
                    static["number"],
                    translation[0],
                    translation[1],
                    translation[2],
                    rotation[0],
                    rotation[1],
                    rotation[2],
                    rotation[3],
                ]
            elif self._dimension == 2:
                # VA. change the position such that is corresponds to a 3D object
                self._static_data[p, :] = [
                    time,
                    static["number"],
                    translation[0],
                    translation[1],
                    0.0,
                    cos(rotation[0] / 2.0),
                    0.0,
                    0.0,
                    sin(rotation[0] / 2.0),
                ]

            p += 1

        # print('current_line , self._static_data', current_line, self._static_data)
        # input()

    def output_radii(self):
        """
        Outputs radii.
        """

        current_line = self._radii_data.shape[0]

        radii = self.get_io_array(self._io.radii(self._nsds))

        if radii is not None:
            self._radii_data.resize(current_line + radii.shape[0], 0)
            self._radii_data[current_line:, :] = radii

    def output_p0s(self):
        """
        Outputs p0 vectors
        """

        current_line = self._p0s_data.shape[0]

        p0s = self.get_io_array(self._io.p0s(self._nsds))

        if p0s is not None:
            self._p0s_data.resize(current_line + p0s.shape[0], 0)
            self._p0s_data[current_line:, :] = p0s

    def output_dynamic_objects(self, initial=False):
        """
        Outputs translations and orientations of dynamic objects.
        """

        current_line = self._dynamic_data.shape[0]

        time = self.current_time()

        # io.positions returns ds state vectors in columns.
        # Each column corresponds to one DS. First value in the column
        # is the ds number.
        positions = self.get_io_array(self._io.positions(self._nsds))

        if positions.shape[0] > 0:
            number_of_ds = positions.shape[0]
            self._dynamic_data.resize(current_line + number_of_ds, 0)

            times = np.empty((number_of_ds, 1))
            times.fill(time)
            if self._dimension == 3:
                self._dynamic_data[current_line:, :] = np.concatenate(
                    (times, positions), axis=1
                )
            elif self._dimension == 2:
                # VA. change the position such that is corresponds to a 3D object
                new_positions = np.zeros((number_of_ds, 8))

                new_positions[:, 0] = positions[:, 0]  # ds number
                new_positions[:, 1] = positions[:, 1]  # x position
                new_positions[:, 2] = positions[:, 2]  # y position

                new_positions[:, 4] = np.cos(positions[:, 3] / 2.0)
                new_positions[:, 7] = np.sin(positions[:, 3] / 2.0)

                self._dynamic_data[current_line:, :] = np.concatenate(
                    (times, new_positions), axis=1
                )

    def output_velocities(self):
        """
        Output velocities of dynamic objects
        """
        current_line = self._velocities_data.shape[0]

        time = self.current_time()

        velocities = self.get_io_array(self._io.velocities(self._nsds))

        if velocities.shape[0] > 0:
            number_of_ds = velocities.shape[0]
            self._velocities_data.resize(current_line + number_of_ds, 0)

            times = np.empty((number_of_ds, 1))
            times.fill(time)
            if self._dimension == 3:
                self._velocities_data[current_line:, :] = np.concatenate(
                    (times, velocities[:, :-1]), axis=1
                )
            elif self._dimension == 2:
                # VA. change the position such that is corresponds to a 3D object
                new_velocities = np.zeros((number_of_ds, 7))

                new_velocities[:, 0] = velocities[:, 0]  # ds number
                new_velocities[:, 1] = velocities[:, 1]  # x velocity
                new_velocities[:, 2] = velocities[:, 2]  # y velocity

                new_velocities[:, 6] = velocities[:, 3]  # theta velocity
                self._velocities_data[current_line:, :] = np.concatenate(
                    (times, new_velocities), axis=1
                )

    def output_contact_forces(self):
        """
        Outputs contact forces
        _output_contact_index_set default value is 1.
        """
        if self._nsds.topology().indexSetsSize() > 1:
            time = self.current_time()

            contact_points = self._io.contactPoints(
                self._nsds, self._output_contact_index_set
            )

            if contact_points.shape[0] > 0:
                current_line = self._cf_data.shape[0]
                # Increase the number of lines in cf_data
                # (h5 dataset with chunks)
                self._cf_data.resize(current_line + contact_points.shape[0], 0)
                times = np.empty((contact_points.shape[0], 1))
                times.fill(time)

                if self._dimension == 3 or self.config.backend == "vnative":
                    self._cf_data[current_line:, :] = np.concatenate(
                        (times, contact_points), axis=1
                    )

                elif self._dimension == 2:

                    # VA. change the contact info such that
                    # this corresponds to a 3D object
                    new_contact_points = np.zeros((contact_points.shape[0], 25))
                    new_contact_points[:, 0] = contact_points[:, 0]  # mu
                    new_contact_points[:, 1] = contact_points[:, 1]  # posa
                    new_contact_points[:, 2] = contact_points[:, 2]
                    # new_contact_points[:, 3]
                    new_contact_points[:, 4] = contact_points[:, 3]  # posb
                    new_contact_points[:, 5] = contact_points[:, 4]
                    # new_contact_points[:, 6]
                    new_contact_points[:, 7] = contact_points[:, 5]  # nc
                    new_contact_points[:, 8] = contact_points[:, 6]
                    # new_contact_points[:, 9]
                    # cf
                    new_contact_points[:, 10] = contact_points[:, 7]
                    new_contact_points[:, 11] = contact_points[:, 8]
                    # new_contact_points[:, 12]
                    new_contact_points[:, 13] = contact_points[:, 9]  # gap
                    new_contact_points[:, 14] = contact_points[:, 10]
                    # new_contact_points[:, 15]
                    # relative velocity
                    new_contact_points[:, 16] = contact_points[:, 11]
                    new_contact_points[:, 17] = contact_points[:, 12]
                    # new_contact_points[:, 18]
                    # reaction impulse
                    new_contact_points[:, 19] = contact_points[:, 13]
                    new_contact_points[:, 20] = contact_points[:, 14]
                    # new_contact_points[:, 21]
                    # inter id
                    new_contact_points[:, 22] = contact_points[:, 15]
                    new_contact_points[:, 23] = contact_points[:, 16]  # ds 1
                    new_contact_points[:, 24] = contact_points[:, 17]  # ds 2
                    self._cf_data[current_line:, :] = np.concatenate(
                        (times, new_contact_points), axis=1
                    )

                # return the number of contacts
                return len(contact_points)
            return 0
        return 0

    def output_contact_info(self):
        """
        Outputs contact infos
        _output_contact_index_set default value is 1.
        """
        if self._nsds.topology().indexSetsSize() > 1:
            time = self.current_time()
            contact_info = self._io.contactInfo(
                self._nsds, self._output_contact_index_set
            )

            if contact_info is not None:
                current_line = self._cf_info.shape[0]
                # Increase the number of lines in cf_data
                # (h5 dataset with chunks)
                self._cf_info.resize(current_line + contact_info.shape[0], 0)
                times = np.empty((contact_info.shape[0], 1))
                times.fill(time)

                self._cf_info[current_line:, :] = np.concatenate(
                    (times, contact_info), axis=1
                )
                # return the number of contacts
                return len(contact_info)
            return 0
        return 0

    def output_contact_work(self):
        """
        Outputs contact contact_work
        _output_contact_index_set default value is 1.
        """
        if self._nsds.topology().indexSetsSize() > 1:
            time = self.current_time()
            contact_work = self._io.contactContactWork(
                self._nsds,
                self._output_contact_index_set,
                omega=self._run_options.get("theta"),
            )

            # print(contact_work)
            if contact_work is not None:
                current_line = self._cf_work.shape[0]
                # Increase the number of lines in cf_data
                # (h5 dataset with chunks)
                self._cf_work.resize(current_line + contact_work.shape[0], 0)
                times = np.empty((contact_work.shape[0], 1))
                times.fill(time)

                self._cf_work[current_line:, :] = np.concatenate(
                    (times, contact_work), axis=1
                )
                # return the number of contacts
                # print(self._cf_contact_work[:,:])
                return len(contact_work)
            return 0
        return 0

    def output_energy_and_work(self):
        """
        Outputs energy_and_work
        _output_contact_index_set default value is 1.
        """
        if self._nsds.topology().indexSetsSize() > 1:
            time = self.current_time()

            kinetic_sum = 0.0
            # external_forces_work_sum_old = 0.0
            external_forces_work_sum = 0.0
            energy_contact_work = np.zeros(8)

            positions = self.get_io_array(self._io.positions(self._nsds))
            nsds = self._nsds
            if positions is not None:
                ds_idx = positions[:, 0]
                for i in ds_idx:
                    n_ds = int(i)
                    neds = nsds.dynamicalSystem(n_ds)
                    kinetic = neds.computeKineticEnergy()
                    kinetic_sum = kinetic + kinetic_sum

            work_forces = self._osi.computeWorkForces()
            if work_forces is not None:
                # print(work_forces)
                external_forces_work_sum = np.sum(work_forces[:, 1])

            cf_work = self._io.contactContactWork(
                self._nsds,
                self._output_contact_index_set,
                omega=self._run_options.get("theta"),
            )
            # print('cf_work', cf_work)

            if cf_work.shape[0] > 0:
                # print('cf_work', cf_work)

                normal_work = cf_work[:, 1]
                normal_work_negative = np.where(normal_work < 0, normal_work, 0)
                # print('normal_work_negative', normal_work_negative)

                tangent_work = cf_work[:, 2]
                tangent_work_negative = np.where(tangent_work < 0, tangent_work, 0)
                # print('tangent_work_negative', tangent_work_negative)

                normal_work_negative_sum = np.sum(normal_work_negative)
                tangent_work_negative_sum = np.sum(tangent_work_negative)

                normal_work_sum = np.sum(normal_work)
                tangent_work_sum = np.sum(tangent_work)

                normal_work_sum_theta = np.sum(cf_work[:, 3])
                tangent_work_sum_theta = np.sum(cf_work[:, 4])

            else:
                normal_work_negative_sum = 0.0
                tangent_work_negative_sum = 0.0

                normal_work_sum = 0.0
                tangent_work_sum = 0.0
                normal_work_sum_theta = 0.0
                tangent_work_sum_theta = 0.0

            energy_contact_work[0] = kinetic_sum
            energy_contact_work[1] = external_forces_work_sum

            energy_contact_work[2] = normal_work_sum
            energy_contact_work[3] = tangent_work_sum

            energy_contact_work[4] = normal_work_sum_theta
            energy_contact_work[5] = tangent_work_sum_theta

            energy_contact_work[6] = normal_work_negative_sum
            energy_contact_work[7] = tangent_work_negative_sum

            current_line = self._energy_work.shape[0]
            self._energy_work.resize(current_line + 1, 0)
            times = np.empty((1, 1))
            times.fill(time)
            energy_contact_work = np.array(energy_contact_work).reshape(
                (1, len(energy_contact_work))
            )

            self._energy_work[current_line:, :] = np.concatenate(
                (times, energy_contact_work), axis=1
            )

        return 0

    def output_domains(self):
        """
        Outputs domains of contact points
        """
        if self._nsds.topology().indexSetsSize() > 1:
            time = self.current_time()
            domains = self.get_io_array(self._io.domains(self._nsds))

            if domains is not None:
                current_line = self._domain_data.shape[0]
                self._domain_data.resize(current_line + domains.shape[0], 0)
                times = np.empty((domains.shape[0], 1))
                times.fill(time)

                self._domain_data[current_line:, :] = np.concatenate(
                    (times, domains), axis=1
                )

    def output_solver_infos(self):
        """
        Outputs solver #iterations & precision reached
        """

        time = self.current_time()
        so = self._simulation.oneStepNSProblem(0).numericsSolverOptions()

        current_line = self._solv_data.shape[0]
        self._solv_data.resize(current_line + 1, 0)
        if self.config.backend == "vnative":
            # need a fix in bridge.py
            iterations = so.iparam(sn.params.SICONOS_IPARAM_ITER_DONE)
            precision = so.dparam(sn.params.SICONOS_DPARAM_RESIDU)
        else:
            iterations = so.iparam[sn.params.SICONOS_IPARAM_ITER_DONE]
            precision = so.dparam[sn.params.SICONOS_DPARAM_RESIDU]
        if so.solverId == sn.solver_ids.SICONOS_GENERIC_MECHANICAL_NSGS:
            local_precision = so.dparam[3]  # Check this !
        elif so.solverId == sn.solver_ids.SICONOS_FRICTION_3D_NSGS:
            local_precision = 0.0
        else:
            local_precision = precision

        self._solv_data[current_line, :] = [
            time,
            iterations,
            precision,
            local_precision,
        ]

    def output_results(self, with_timer=False):
        self.log(self.output_static_objects, with_timer)()

        self.log(self.output_dynamic_objects, with_timer)()

        self.log(self.output_velocities, with_timer)()

        if self.config.backend == "vnative":
            self.log(self.output_p0s, with_timer)()

        if self._output_contact_forces:
            self.log(self.output_contact_forces, with_timer)()

        if self._output_contact_info and (
            self.config.backend == "bullet" or self.config.backend == "vnative"
        ):
            self.log(self.output_contact_info, with_timer)()
        else:
            self.print_verbose(
                "[warning] output_contact_info is only "
                "available with bullet backend for the moment"
            )
            self.print_verbose(
                "          to remove this message set output_contact_info options to False"
            )

        if self._output_contact_work:
            self.log(self.output_contact_work, with_timer)()
            if (
                self._run_options["skip_last_update_output"]
                or self._run_options["skip_last_update_input"]
            ):
                self.print_verbose(
                    "[warning] output_contact_work with "
                    "the options skip_last_update_output=True\n"
                    " or skip_last_update_input=True\n"
                    " could result in wrong output"
                )
        else:
            self.print_verbose(
                "[warning] output_contact_work is only "
                "available with bullet backend for the moment"
            )
            self.print_verbose(
                "          to remove this message set output_contact_info options to False"
            )

        if self._output_energy_work:
            self.log(self.output_energy_and_work, with_timer)()

        if self._should_output_domains:
            self.log(self.output_domains, with_timer)()

        if self.config.backend != "vnative":
            self.log(self.output_solver_infos, with_timer)()

        if self.config.backend == "vnative":
            self.log(self.output_radii, with_timer)()

        self.log(self._out.flush)()

    def output_run_options(self):
        """
        Outputs run_options
        """
        d = self._run_options.copy()

        so = d["solver_options"]
        if so:
            d["solver_options"] = {}
            d["solver_options"]["solverId"] = so.solverId
            d["solver_options"]["solver name"] = so.solverId
            d["solver_options"]["iparam size"] = so.iSize
            for i in range(so.iSize):
                d["solver_options"]["iparam[" + str(i) + "]"] = int(so.iparam[i])
            for i in range(so.dSize):
                if so.dparam[i] <= 1e24:
                    d["solver_options"]["dparam[" + str(i) + "]"] = float(so.dparam[i])
            # d['solver_options']['numberOfInternalSolvers']=so.numberOfInternalSolvers
            # fix it

        sop = d["solver_options_pos"]
        if sop:
            d["solver_options_pos"] = {}
            d["solver_options_pos"]["solverId"] = so.solverId
            d["solver_options_pos"]["solver name"] = so.solverId
            d["solver_options_pos"]["iparam size"] = so.iSize
            for i in range(so.iSize):
                d["solver_options_pos"]["iparam[" + str(i) + "]"] = int(so.iparam[i])
            for i in range(so.dSize):
                if so.dparam[i] <= 1e24:
                    d["solver_options_pos"]["dparam[" + str(i) + "]"] = float(
                        so.dparam[i]
                    )
            # d['solver_options_pos']['numberOfInternalSolvers']=so.numberOfInternalSolvers
            # fix it

        bo = d["bullet_options"]
        if bo:
            opt = [
                "clearOverlappingPairCache",
                "contactBreakingThreshold",
                "contactProcessingThreshold",
                "dimension",
                "enablePolyhedralContactClipping",
                "enableSatConvex",
                "minimumPointsPerturbationThreshold",
                "perturbationIterations",
                "useAxisSweep3",
                "worldScale",
            ]  # fix it
            opt = [
                "contactProcessingThreshold",
                "dimension",
                "enablePolyhedralContactClipping",
                "enableSatConvex",
                "minimumPointsPerturbationThreshold",
                "perturbationIterations",
                "useAxisSweep3",
                "worldScale",
            ]
            d["bullet_options"] = {}
            for e in opt:
                #                print('getattr(bo, e)', getattr(bo, e))
                d["bullet_options"][e] = getattr(bo, e)

        # to fix the serialization of run_options, we should use pickle
        # which is able to serialize python object.

        d["friction_contact_trace_params"] = "not serialized"  # fix it
        d["osi"] = "not serialized"  # fix it
        d["time_stepping"] = "not serialized"  # fix it

        d["start_run_iteration_hook"] = "not serialized"  # fix it
        d["before_next_step_iteration_hook"] = "not serialized"  # fix it
        d["end_run_iteration_hook"] = "not serialized"  # fix it

        if d["set_external_forces"] is not None:
            try:
                d["set_external_forces"] = (
                    type(d["set_external_forces"]).__name__ + "(name serialized)"
                )
            except KeyError:
                d["set_external_forces"] = "not serialized"

        if d["controller"] is not None:
            try:
                d["controller"] = (
                    getattr(d["controller"], "__name__", type(d["controller"]).__name__)
                    + " (serialized)"
                )
            except Exception:
                d["controller"] = "not serialized"

        # Special care for enum, to make them json-complient
        def serialize_enum(obj):
            enum_types = (
                siconos.nonsmooth_formulations.LinearOSNSAssemblyType,
                siconos.simulation.TimeSteppingType,
            )
            if self.config.use_bullet:
                enum_types += (
                    siconos.mechanics.collision.bullet.SiconosBulletDimension,
                )
            if isinstance(obj, enum_types):
                return obj.value
            elif isinstance(obj, (str, int, float, bool, list, dict, type(None))):
                return obj
            elif hasattr(obj, "__class__") and hasattr(obj, "__dir__"):
                return {
                    key: serialize_enum(getattr(obj, key))
                    for key in dir(obj)
                    if not key.startswith("_") and hasattr(obj, key)
                }

            raise TypeError(
                f"Object of type {obj.__class__.__name__} is not JSON serializable"
            )

        dict_json = json.dumps(d)
        self._run_options_data.attrs["options"] = dict_json

    def print_solver_infos(self):
        """
        Outputs solver #iterations & precision reached
        """
        time = self.current_time()
        so = self._simulation.oneStepNSProblem(0).numericsSolverOptions()
        iterations = so.iparam[sn.params.SICONOS_IPARAM_ITER_DONE]
        precision = so.dparam[sn.params.SICONOS_DPARAM_RESIDU]
        msg = "Numerics solver info at time : {0:10.6f}".format(time)
        msg += " iterations = {0:8d}".format(iterations)
        msg += " precision = {0:5.3e}".format(precision)
        self.print_verbose(msg)

    def import_external_functions(self):
        topo = self._nsds.topology()

        for name in self._external_functions:
            ext_fun = self._external_functions[name]
            plugin_name = ext_fun.attrs["plugin_name"]
            plugin_function_name = ext_fun.attrs["plugin_function_name"]
            body_name = ext_fun.attrs["body_name"]

            ds = topo.getDynamicalSystem(body_name)

            if "function_name" in ext_fun.attrs:
                function_name = ext_fun.attrs["function_name"]

                getattr(ds, function_name)(plugin_name, plugin_function_name)
            else:
                bc_indices = ext_fun.attrs["bc_indices"]
                # a boundary condition
                bc = sm.BoundaryCondition(bc_indices)
                bc.setComputePrescribedVelocityFunction(
                    plugin_name, plugin_function_name
                )
                ds.setBoundaryConditions(bc)

    def computeOneStep_python(self, with_timer):
        s = self._simulation

        # 1 - s.initialize:
        # self.log(s.initialize, with_timer)()
        self.log(s.initializeOSIAssociations, with_timer)()
        self.log(s.initializeIndexSets, with_timer)()
        self.log(s.applyNSDSChangelogForDS, with_timer)()
        self.log(s.updateWorldFromDS, with_timer)()
        self.log(s.updateInteractions, with_timer)()
        self.log(s.initializeNSDSChangelog, with_timer)()

        self.log(s.firstInitialize, with_timer)()

        # 2  advanceToEvent:
        if not s.skipResetLambdas:
            self.log(s.resetLambdas, with_timer)()

        # 2.1   newtonSolve:
        # Again the access to s._newtonTolerance generates a segfault due to director
        newtonTolerance = s.newtonTolerance
        newtonMaxIteration = s.newtonMaxIteration

        # return _kernel.TimeStepping_newtonSolve(self, criterion, maxStep)
        # RuntimeError: accessing protected member newtonSolve
        # s.newtonSolve(newtonTolerance, newtonMaxIteration);

        newtonNbIterations = 0
        isNewtonConverge = False
        explode_computeOneStepNSProblem_in_python = self._run_options.get(
            "explode_computeOneStepNSProblem_in_python"
        )

        # self.log(s.initializeNewtonSolve, with_timer)()
        # explode version

        self.log(s.updateAndSwapAllOutput, with_timer)()
        self.log(s.updateIndexSets, with_timer)()
        self.log(s.initializeOneStepNSProblem, with_timer)()
        self.log(s.computeInitialStateOfTheStep, with_timer)()
        self.log(s.updateDSPlugins, with_timer)(s.nextTime())
        self.log(s.computeResidu, with_timer)()

        # missing computeResiduY
        # self.log(s.computeResiduY, with_timer)()

        if s.newtonOptions == simu.LINEAR or s.newtonOptions == simu.LINEAR_IMPLICIT:
            if s.newtonOptions == simu.LINEAR_IMPLICIT:
                self.log(s.prepareNewtonIteration, with_timer)()
            self.log(s.computeFreeState, with_timer)()
            info = 0
            if s.numberOfOSNSProblems > 0:
                if explode_computeOneStepNSProblem_in_python:
                    fc = self._osnspb
                    # self.log(fc.updateInteractionBlocks, with_timer)()
                    self.log(fc.preCompute, with_timer, after=False)(s.nextTime())
                    self.log(fc.updateMu, with_timer)()
                    if self._run_options.get("osi") == integrators.MoreauJeanOSI:
                        if fc.getSizeOutput != 0:
                            info = self.log(fc.solve, with_timer)()
                            self.log(fc.postCompute, with_timer)()
                    else:
                        info = self.log(fc.solve, with_timer)()
                        self.log(fc.postCompute, with_timer)()
                else:
                    info = self.log(s.computeOneStepNSProblem, with_timer)(
                        simu.constants.SICONOS_OSNSP_TS_VELOCITY
                    )
                self.log(s.DefaultCheckSolverOutput, with_timer)(info)
                if not s.skipLastUpdateInput():
                    self.log(s.updateAllInput, with_timer)()
                self.log(s.computeIteration, with_timer)()

        else:
            while (not isNewtonConverge) and newtonNbIterations < newtonMaxIteration:
                # self.print_verbose('newtonNbIterations',newtonNbIterations)
                info = 0
                newtonNbIterations = newtonNbIterations + 1
                self.log(s.prepareNewtonIteration, with_timer)()
                self.log(s.computeFreeState, with_timer)()
                if s.numberOfOSNSProblems > 0:
                    if explode_computeOneStepNSProblem_in_python:
                        fc = self._osnspb
                        self.log(fc.preCompute, with_timer)(s.nextTime())
                        self.log(fc.updateMu, with_timer)()
                        if fc.getSizeOutput != 0:
                            info = self.log(fc.solve, with_timer)()
                            self.log(fc.postCompute, with_timer)()
                    else:
                        info = self.log(s.computeOneStepNSProblem, with_timer)(
                            simu.constants.SICONOS_OSNSP_TS_VELOCITY
                        )
                self.log(s.DefaultCheckSolverOutput, with_timer)(info)
                self.log(s.updateAllInput, with_timer)()
                self.log(s.computeIteration, with_timer)()
                if (not isNewtonConverge) and (newtonNbIterations < newtonMaxIteration):
                    self.log(s.updateOutput, with_timer)()
                isNewtonConverge = self.log(s.newtonCheckConvergence, with_timer)(
                    newtonTolerance
                )
                if s.displayNewtonConvergence():
                    s.displayNewtonConvergenceInTheLoop()

            if s.displayNewtonConvergence():
                s.displayNewtonConvergenceAtTheEnd(info, newtonMaxIteration)

        self.log(s.updateState, with_timer)()
        if not s.skipLastUpdateOutput():
            self.log(s.updateOutput, with_timer)()

    def build_run_options_from_old_arguments_in_kwargs(
        self,
        run_options=None,
        with_timer=False,
        time_stepping=None,
        interaction_manager=None,
        bullet_options=None,
        vnative_options=None,
        controller=None,
        gravity_scale=1.0,
        t0=0,
        T=10,
        h=0.0005,
        multipoints_iterations=None,
        theta=0.50001,
        gamma=0.0,
        Newton_options=siconos.simulation.NONLINEAR,
        Newton_max_iter=20,
        set_external_forces=None,
        solver_options=None,
        solver_options_pos=None,
        osnspb_max_size=0,
        exit_tolerance=None,
        projection_itermax=20,
        projection_tolerance=1e-8,
        projection_tolerance_unilateral=1e-8,
        numerics_verbose=False,
        numerics_verbose_level=0,
        violation_verbose=False,
        verbose=True,
        verbose_progress=True,
        output_frequency=None,
        output_backup=False,
        output_backup_frequency=None,
        output_contact_forces=True,
        output_contact_info=True,
        output_contact_work=True,
        output_energy_work=False,
        friction_contact_trace_params=None,
        output_contact_index_set=1,
        osi=integrators.MoreauJeanOSI,
        constraint_activation_threshold=0.0,
        explode_Newton_solve=False,
        explode_computeOneStep=False,
        explode_computeOneStep_in_python=False,
        explode_computeOneStepNSProblem_in_python=False,
        display_Newton_convergence=False,
        start_run_iteration_hook=None,
        end_run_iteration_hook=None,
        before_next_step_iteration_hook=None,
        skip_last_update_output=False,
    ):
        """Run a simulation from a set of parameters described in a hdf5 file."""
        build_from_kwargs = False
        if run_options is None:
            self.print_verbose("\nWARNING: no run_options given.")
            self.print_verbose("Please, consider to use a run options dictionnary.")
            self.print_verbose("Otherwise, some new options may not be available,")
            self.print_verbose("or becomes obsolete.")
            self.print_verbose("for instance:\n")

            run_options = MechanicsHdf5Runner_run_options()

            run_options["with_timer"] = with_timer
            run_options["time_stepping"] = time_stepping
            run_options["interaction_manager"] = interaction_manager
            run_options["bullet_options"] = bullet_options
            run_options["controller"] = controller
            run_options["gravity_scale"] = gravity_scale
            run_options["t0"] = t0
            run_options["T"] = T
            run_options["h"] = h
            run_options["multipoints_iterations"] = multipoints_iterations
            run_options["theta"] = theta
            run_options["gamma"] = gamma
            run_options["Newton_options"] = Newton_options
            run_options["Newton_max_iter"] = Newton_max_iter
            run_options["interaction_manager"] = None
            run_options["set_external_forces"] = set_external_forces
            run_options["solver_options"] = solver_options
            run_options["solver_options_pos"] = solver_options_pos
            run_options["osnspb_max_size"] = osnspb_max_size
            run_options["exit_tolerance"] = exit_tolerance
            run_options["projection_itermax"] = projection_itermax
            run_options["projection_tolerance"] = projection_tolerance
            run_options["projection_tolerance_unilateral"] = (
                projection_tolerance_unilateral
            )
            run_options["numerics_verbose"] = numerics_verbose
            run_options["numerics_verbose_level"] = numerics_verbose_level
            run_options["violation_verbose"] = violation_verbose
            run_options["verbose"] = verbose
            run_options["verbose_progress"] = verbose_progress
            run_options["output_frequency"] = output_frequency
            run_options["output_backup"] = output_backup
            run_options["output_backup_frequency"] = output_backup_frequency
            run_options["friction_contact_trace_params"] = friction_contact_trace_params
            run_options["output_contact_index_set"] = output_contact_index_set
            run_options["osi"] = osi
            run_options["constraint_activation_threshold"] = (
                constraint_activation_threshold
            )
            run_options["explode_computeOneStep_in_python"] = explode_Newton_solve
            run_options["explode_computeOneStepNSProblem_in_python"] = (
                explode_computeOneStepNSProblem_in_python
            )
            run_options["display_Newton_convergence"] = display_Newton_convergence
            run_options["start_run_iteration_hook"] = start_run_iteration_hook
            run_options["before_next_step_iteration_hook"] = (
                before_next_step_iteration_hook
            )
            run_options["end_run_iteration_hook"] = end_run_iteration_hook
            run_options["skip_last_update_output"] = skip_last_update_output
            run_options["output_contact_forces"] = output_contact_forces
            run_options["output_contact_info"] = output_contact_info
            run_options["output_contact_work"] = output_contact_work
            run_options["output_energy_work"] = output_energy_work

            build_from_kwargs = True

        self._run_options = run_options

        return build_from_kwargs

    def run_initialize(self):

        run_options = self._run_options

        run_options.check_valid_run_options()

        if run_options["verbose"]:
            run_options.display()

        self.print_verbose("setup model simulation ...")
        if run_options["set_external_forces"] is not None:
            self._set_external_forces = run_options["set_external_forces"]

        interaction_manager = run_options["interaction_manager"]
        if interaction_manager is None:
            interaction_manager = self.config.default_manager_class

        if run_options["time_stepping"] is None:
            self._time_stepping_class = self.config.default_simulation_class
        else:
            self._time_stepping_class = run_options["time_stepping"]

        if run_options["output_frequency"] is not None:
            self._output_frequency = run_options["output_frequency"]

        if run_options["output_backup_frequency"] is not None:
            self._output_backup_frequency = run_options["output_backup_frequency"]

        if run_options["output_backup"] is not None:
            self._output_backup = run_options["output_backup"]

        if run_options["output_contact_forces"] is not None:
            self._output_contact_forces = run_options["output_contact_forces"]

        if run_options["output_contact_info"] is not None:
            self._output_contact_info = run_options["output_contact_info"]

        if run_options["output_contact_work"] is not None:
            self._output_contact_work = run_options["output_contact_work"]

        if run_options["output_energy_work"] is not None:
            self._output_energy_work = run_options["output_energy_work"]

        if run_options["gravity_scale"] is not None:
            self._gravity_scale = run_options["gravity_scale"]

        t0 = run_options["t0"]
        T = run_options["T"]
        h = run_options["h"]

        # cold restart
        times = set()
        if self.dynamic_data() is not None and len(self.dynamic_data()) > 0:
            dpos_data = self.dynamic_data()
            times = set(dpos_data[:, 0])
            t0 = float(max(times))

        # Time-related parameters for this simulation run
        self._k0 = 1 + int(t0 / h)
        self._k = self._k0
        kT = self._k0 + int((T - t0) / h)
        if T > t0:
            self.print_verbose("")
            msg = "Simulation will run from time {0:.4f} ".format(t0)
            msg += "to {0:.4f}s, ".format(T)
            msg += "step {} to step {} (h={}, ".format(self._k0, kT, h)
            msg += "times=[{},{}])".format(
                min(times) if len(times) > 0 else "?",
                max(times) if len(times) > 0 else "?",
            )
            self.print_verbose(msg)
            self.print_verbose("")
        else:
            msg = "Simulation time {0} >= T={1}, exiting.".format(t0, T)
            self.print_verbose(msg)
            exit(0)

        # Respect run() parameter for multipoints_iterations for
        # backwards compatibility, but this is overridden by
        # SiconosBulletOptions if one is provided.
        multipoints_iterations = run_options.get("multipoints_iterations")
        bullet_options = run_options.get("bullet_options")
        vnative_options = run_options.get("vnative_options")
        if (multipoints_iterations is not None) and (bullet_options is not None):
            msg = "[io.mechanics] run(): one cannot give multipoints_iterations"
            msg += " and bullet_options simultaneously. \n"
            msg += "                             multipoints_iterations will be marked"
            msg += " as obsolete. use preferably bullet_options\n"
            msg += "                             with"
            msg += " bullet_options.perturbationIterations and"
            msg += "bullet_options.minimumPointsPerturbationThreshold."
            raise RuntimeError(msg)

        if (
            bullet_options is None and self.config.use_bullet
        ):  # siconos.mechanics.have_bullet:
            bullet_options = self.config.bullet.SiconosBulletOptions()
            if multipoints_iterations:
                bullet_options.perturbationIterations = 3 * multipoints_iterations
                bullet_options.minimumPointsPerturbationThreshold = (
                    3 * multipoints_iterations
                )

               
        # MB: this may be in conflict with 'dimension' attribute
        if bullet_options is not None and self.config.bullet is not None:
            # we are using bullet
            if self._dimension == 2:
                if bullet_options.dimension != self.config.bullet.TwoD:
                    self.print_verbose(
                        """Warning. The infered dimention in attrs["dimension"] is 2D
                        but the bullet_options are not consistent
                        we impose bullet_options.dimension == self.config.bullet.TwoD""")
                    bullet_options.dimension == self.config.bullet.TwoD
        else:
            if self._out.attrs.get("dimension", None) is None:
                # this is a second place to set the default
                self._dimension = 3

                
        if self.config.backend == "vnative":
            if vnative_options is None:
                vnative_options = sio.SpaceFilterOptions()
            self._interman = interaction_manager(vnative_options)
        else:
            self._interman = interaction_manager(bullet_options)

        joints = list(self.joints())
        if hasattr(self._interman, "useEqualityConstraints"):
            if len(joints) == 0:
                self._interman.useEqualityConstraints(False)
            else:
                self._interman.useEqualityConstraints(True)

        # (0) NonSmooth Dynamical Systems definition
        self._nsds = sm.NonSmoothDynamicalSystem(t0, T)
        nsds = self._nsds

        self.print_verbose("import scene ...")
        body_class = run_options.get("body_class")
        shape_class = run_options.get("shape_class")
        face_class = run_options.get("face_class")
        edge_class = run_options.get("egde_class")

        self.import_scene(t0, body_class, shape_class, face_class, edge_class)

        self._output_contact_index_set = run_options.get("output_contact_index_set")

        # (1) OneStepIntegrators
        osi = run_options.get("osi")
        self._osi = osi(run_options.get("theta"))
        if run_options.get("gamma"):
            self._osi.setGamma(run_options.get("gamma"))
        if run_options.get("constraint_activation_threshold"):
            self._osi.setConstraintActivationThreshold(
                run_options["constraint_activation_threshold"]
            )

        if run_options.get("activate_with_negative_relative_velocity"):
            self._osi.setActivateWithNegativeRelativeVelocity(
                run_options["activate_with_negative_relative_velocity"]
            )

        if run_options.get("constraint_activation_threshold_velocity"):
            self._osi.setConstraintActivationThresholdVelocity(
                run_options["constraint_activation_threshold_velocity"]
            )

        # (2) Time discretisation --
        if self.config.backend != "vnative":
            timedisc = simu.TimeDiscretisation(t0, h)
        else:
            # cannot override simu, why ?
            timedisc = sio.TimeDiscretisation(t0, h)

        # (3) choice of default OneStepNonSmoothProblem
        # w.r.t the type of nslaws
        nslaw_type_list = []
        for name in self._nslaws_data:
            nslaw_type_list.append(self._nslaws_data[name].attrs["type"])

        # print(set(nslaw_type_list))

        # This trick is used to add the EqualityConditionNSL
        # to the list of nslaw type
        # this must be improved by adding the EqualityConditionNSL
        # in self._nslaws_data
        # when a joint is imported.
        # For the moment, the nslaw is implicitely added
        # when we import_joint but is not stored
        # self._nslaws_data

        if len(joints) > 0:
            nslaw_type_list.append("EqualityConditionNSL")

        nb_of_nslaw_type = len(set(nslaw_type_list))

        # Check nslaw/OSI compatibility
        if osi == integrators.MoreauJeanGOSI:
            checknsl = "NewtonImpactFrictionNSL" in set(
                nslaw_type_list
            ) or "NewtonImpactRollingFrictionNSL" in set(nslaw_type_list)
            if not checknsl:
                msg = "MoreauJeanGOSI can only deal "
                msg += "with NewtonImpactFrictionNSL or NewtonImpactRollingFrictionNSL."
                raise RuntimeError(msg)
            if nb_of_nslaw_type > 1:
                msg = "MoreauJeanGOSI cannot only deal with multiple"
                msg += "impact laws at the same time "
                raise RuntimeError(msg)

        # Creates and initialises the one-step nonsmooth problem.
        # The OSI and/or the nonsmooth law drives the choice.

        # rationale for choosing numerics solver options
        # if solver_option is None --> we leave Siconos/kernel choosing the default option
        # else we use the user solver_options
        solver_options = run_options.get("solver_options")
        osnspb_max_size = run_options.get("osnspb_max_size")
        osnspb_assembly_type = run_options.get("osns_assembly_type")
        friction_contact_trace_params = run_options.get("friction_contact_trace_params")
        if friction_contact_trace_params is None:
            # Global friction contact.
            if osi == integrators.MoreauJeanGOSI:
                if "NewtonImpactFrictionNSL" in set(nslaw_type_list):
                    if solver_options is None:
                        osnspb = nsf.GlobalFrictionContact(self._dimension)
                    else:
                        osnspb = nsf.GlobalFrictionContact(
                            self._dimension, solver_options
                        )
                elif "NewtonImpactRollingFrictionNSL" in set(nslaw_type_list):
                    if self._dimension == 3:
                        dimension_contact = 5
                    elif self._dimension == 2:
                        dimension_contact = 3
                    if solver_options is None:
                        osnspb = nsf.GlobalRollingFrictionContact(dimension_contact)
                    else:
                        osnspb = nsf.GlobalRollingFrictionContact(
                            dimension_contact, solver_options
                        )
                osnspb.setMStorageType(sn.params.NM_SPARSE)
                # if sid == sn.solver_ids.SICONOS_GLOBAL_FRICTION_3D_ADMM:
                #     osnspb.setMStorageType(sn.params.NM_SPARSE)
                #     # which is the default for gfc in kernel, is it?
                # else:
                #     osnspb.setMStorageType(sn.params.NM_SPARSE_BLOCK)
                osnspb.setMaxSize(osnspb_max_size)
            else:
                if "EqualityConditionNSL" in set(nslaw_type_list):
                    if solver_options is None:
                        osnspb = nsf.GenericMechanical()
                    else:
                        osnspb = nsf.GenericMechanical(solver_options)
                else:
                    if (
                        ("NewtonImpactFrictionNSL" in set(nslaw_type_list))
                        or ("FremondImpactFrictionNSL" in set(nslaw_type_list))
                        or (len(set(nslaw_type_list)) == 0)
                    ):
                        if solver_options is None:
                            osnspb = nsf.FrictionContact(self._dimension)
                        else:
                            osnspb = nsf.FrictionContact(
                                self._dimension, solver_options
                            )
                        osnspb.setMaxSize(osnspb_max_size)
                        osnspb.setMStorageType(sn.params.NM_SPARSE_BLOCK)

                    elif "NewtonImpactRollingFrictionNSL" in set(nslaw_type_list):
                        if self._dimension == 3:
                            dimension_contact = 5
                        elif self._dimension == 2:
                            dimension_contact = 3
                        if solver_options is None:
                            osnspb = nsf.RollingFrictionContact(dimension_contact)
                        else:
                            osnspb = nsf.RollingFrictionContact(
                                dimension_contact, solver_options
                            )
                    else:
                        msg = "Unknown nslaw type"
                        msg += str(set(nslaw_type_list))
                        raise RuntimeError(msg)

            if osnspb_assembly_type:
                osnspb.setMStorageType(sn.params.NM_SPARSE)
                osnspb.setAssemblyType(osnspb_assembly_type)

        else:  # With trace
            pass
            if self.config.backend == "vnative":
                osnspb = nsf.FrictionContact(self._dimension, solver_options)

                vnative_tp = nsf.TraceParams(friction_contact_trace_params)
                osnspb.handle().set_trace_params(vnative_tp.handle())
                osnspb.handle().set_trace(True)
            else:
                if solver_options is None:
                    solver_options = sn.solver_options_create(
                        sn.solver_ids.SICONOS_FRICTION_3D_NSGS
                    )
                # sid = solver_options.solverId
                if osi == integrators.MoreauJeanGOSI:
                    if "NewtonImpactFrictionNSL" in set(nslaw_type_list):
                        osnspb = GFCTrace(
                            3, solver_options, friction_contact_trace_params, nsds
                        )
                        osnspb.setMStorageType(sn.params.NM_SPARSE)
                        osnspb.setMaxSize(osnspb_max_size)
                    elif "NewtonImpactRollingFrictionNSL" in set(nslaw_type_list):
                        osnspb = GRFCTrace(
                            5, solver_options, friction_contact_trace_params, nsds
                        )
                        osnspb.setMStorageType(sn.params.NM_SPARSE)
                        osnspb.setMaxSize(osnspb_max_size)
                    else:
                        msg = "Unknown nslaw type"
                        msg += str(set(nslaw_type_list))
                        raise RuntimeError(msg)
                else:
                    osnspb = FCTrace(
                        3, solver_options, friction_contact_trace_params, nsds
                    )
                    osnspb.setMaxSize(osnspb_max_size)
                    osnspb.setMStorageType(sn.params.NM_SPARSE_BLOCK)

        numerics_verbose = run_options.get("numerics_verbose")
        osnspb.setNumericsVerboseMode(numerics_verbose)
        if numerics_verbose:
            sn.numerics_set_verbose(run_options.get("numerics_verbose_level"))

        # keep previous solution
        osnspb.setKeepLambdaAndYState(True)

        self._osnspb = osnspb

        # (6) Simulation setup with (1) (2) (3) (4) (5)
        if self._time_stepping_class == simu.TimeSteppingDirectProjection:
            if run_options["solver_options_pos"] is None:
                osnspb_pos = nsf.MLCPProjectOnConstraints(
                    sn.solver_ids.SICONOS_MLCP_ENUM, 1.0
                )
            else:
                osnspb_pos = nsf.MLCPProjectOnConstraints(
                    run_options["solver_options_pos"], 1.0
                )

            osnspb_pos.setMaxSize(osnspb_max_size)
            osnspb_pos.setMStorageType(sn.params.NM_DENSE)
            # "not yet implemented for sparse storage"
            osnspb_pos.setNumericsVerboseMode(numerics_verbose)
            osnspb_pos.setKeepLambdaAndYState(True)
            simulation = self._time_stepping_class(
                nsds, timedisc, self._osi, osnspb, osnspb_pos
            )
            simulation.setProjectionMaxIteration(run_options["projection_itermax"])
            simulation.setConstraintTolUnilateral(
                run_options["projection_tolerance_unilateral"]
            )
            simulation.setConstraintTol(run_options["projection_tolerance"])
        else:
            simulation = self._time_stepping_class(nsds, timedisc)
            simulation.insertIntegrator(self._osi)
            simulation.insertNonSmoothProblem(osnspb)

        simulation.insertInteractionManager(self._interman)

        simulation.setNewtonOptions(run_options["Newton_options"])
        simulation.setNewtonMaxIteration(run_options["Newton_max_iter"])
        simulation.setNewtonTolerance(run_options["Newton_tolerance"])
        simulation.setNewtonWarningOnNonConvergence(
            run_options["Newton_warning_on_nonconvergence"]
        )
        simulation.setWarningNonsmoothSolver(run_options["Warning_nonsmooth_solver"])

        simulation.setSkipLastUpdateOutput(run_options.get("skip_last_update_output"))
        simulation.setSkipLastUpdateInput(run_options.get("skip_last_update_input"))
        simulation.setSkipResetLambdas(run_options.get("skip_reset_lambdas"))

        verbose = run_options.get("verbose")
        if verbose:
            simulation.setDisplayNewtonConvergence(
                run_options.get("display_Newton_convergence")
            )

        self._simulation = simulation

        # time step necessary for output
        if self.config.backend == "vnative":
            self._io.setSimulation(self._simulation)

        if len(self._plugins) > 0:
            self.print_verbose("import plugins ...")
            self.import_plugins()

        if len(self._external_functions) > 0:
            self.print_verbose("import external functions ...")
            self.import_external_functions()
        controller = run_options.get("controller")
        if controller is not None:
            controller.initialize(self)

        self._start_run_iteration_hook = run_options.get("start_run_iteration_hook")
        if self._start_run_iteration_hook is not None:
            self._start_run_iteration_hook.initialize(self)
        self._end_run_iteration_hook = run_options.get("end_run_iteration_hook")
        if self._end_run_iteration_hook is not None:
            self._end_run_iteration_hook.initialize(self)
        self._before_next_step_iteration_hook = run_options.get(
            "before_next_step_iteration_hook"
        )
        if self._before_next_step_iteration_hook is not None:
            self._before_next_step_iteration_hook.initialize(self)

        self.print_verbose("first output static and dynamic objects ...")
        self.output_static_objects()
        self.output_dynamic_objects()
        self.output_velocities()

        if self._should_output_domains:
            self.log(self.output_domains, run_options["with_timer"])()

        # self.output_run_options()

        self.print_verbose("start simulation ...")
        self._initializing = False

    def solver_verbose(self, number_of_contacts):

        so = self._simulation.oneStepNSProblem(0).numericsSolverOptions()
        if self.config.backend == "vnative":
            # need a fix in bridge.py
            iterations = so.iparam(sn.params.SICONOS_IPARAM_ITER_DONE)
            precision = so.dparam(sn.params.SICONOS_DPARAM_RESIDU)
        else:
            iterations = so.iparam[sn.params.SICONOS_IPARAM_ITER_DONE]
            precision = so.dparam[sn.params.SICONOS_DPARAM_RESIDU]

        solver_output = {}
        mask = "|      "
        solver_output["solver iter"] = [iterations, mask + "{:<10d}"]
        solver_output["solver error"] = [precision, mask + "{:8.4e}"]

        print_violation = None

        if self._run_options.get("violation_verbose") and number_of_contacts > 0:
            print_violation = {}
            if len(self._simulation.y_output(0, 0)) > 0:
                y = self._simulation.y_output(0, 0)
                yplus = np.zeros((2, len(y)))
                yplus[0, :] = y
                y = np.min(yplus, axis=1)
                violation_max = np.max(-y)
                if self._collision_margin is not None:
                    if violation_max >= self._collision_margin:
                        self.print_verbose(
                            "  violation max is larger than the collision_margin"
                        )
                lam = self._simulation.lambda_input(1, 0)
                print_violation["violation max"] = [violation_max, mask + "{:8.4e}"]
                print_violation["reaction max"] = [np.max(lam), mask + "{:8.4e}"]

            if len(self._simulation.y_output(1, 0)) > 0:
                v = self._simulation.y_output(1, 0)
                vplus = np.zeros((2, len(v)))
                vplus[0, :] = v
                v = np.max(vplus, axis=1)
                print_violation["velocity max"] = [np.max(v), mask + "{:8.4e}"]
                print_violation["velocity min"] = [np.min(v), mask + "{:8.4e}"]

        if print_violation is not None:
            print_solver_verbose = {**solver_output, **print_violation}
        else:
            print_solver_verbose = solver_output

        # print banner
        ll = ["| {:<14} ".format(k) for k in print_solver_verbose.keys()]
        ll.append("|")
        self.print_verbose(" ".join(ll))

        # print results
        ll = []
        for k in print_solver_verbose.keys():
            fmt = print_solver_verbose[k][1]
            value = print_solver_verbose[k][0]
            ll.append(fmt.format(value))
        ll.append("|")
        self.print_verbose(" ".join(ll))

    def contact_statistics_verbose(self):
        # Note these are not the same and neither is correct.
        # "_interman.statistics" gives the number of contacts
        # collected by the collision engine, but it's possible some
        # are not in indexset1.  Meanwhile checking the size of
        # the non-smooth problem is wrong when there are joints.
        if self.config.use_bullet:
            number_of_contacts = self._interman.statistics().new_interactions_created
            number_of_contacts += (
                self._interman.statistics().existing_interactions_processed
            )
            if self._verbose and number_of_contacts > 0:
                bullet_statistics = self._interman.statistics()
                self.print_verbose(
                    "bullet_statistics:",
                    "new_interactions_created :",
                    bullet_statistics.new_interactions_created,
                    "existing_interactions_processed :",
                    bullet_statistics.existing_interactions_processed,
                    "interaction_warnings :",
                    bullet_statistics.interaction_warnings,
                )
                self.print_verbose(
                    "number of contacts",
                    number_of_contacts,
                    "(detected)",
                    self._osnspb.getSizeOutput // self._dimension,
                    "(active at velocity level. approx)",
                )
                # self.print_solver_infos()

        else:
            if self.config.backend != "vnative":
                number_of_contacts = self._osnspb.getSizeOutput // self._dimension
            else:
                number_of_contacts = self._osnspb.getSizeOutput() // self._dimension

            if self._run_options["verbose"] and number_of_contacts > 0:
                msg = "number of active contacts at the velocity level (approx)"
                self.print_verbose(msg, number_of_contacts)
                self.print_solver_infos()

        return number_of_contacts

    def run_loop(self):
        verbose = self._run_options.get("verbose")
        self._verbose = verbose
        with_timer = self._run_options.get("with_timer")
        body_class = self._run_options.get("body_class")
        shape_class = self._run_options.get("shape_class")
        face_class = self._run_options.get("face_class")
        edge_class = self._run_options.get("edge_class")
        controller = self._run_options.get("controller")
        friction_contact_trace_params = self._run_options.get(
            "friction_contact_trace_params"
        )
        exit_tolerance = self._run_options.get("exit_tolerance")
        t0 = self._run_options.get("t0")
        T = self._run_options.get("T")
        h = self._run_options.get("h")

        while self._simulation.hasNextEvent():
            if self._run_options.get("verbose_progress"):
                self.print_verbose(
                    "step",
                    self._k,
                    "of",
                    self._k0 + int((T - t0) / h) - 1,
                    " time : {0:12.8f}".format(self.current_time()),
                )

            if self._start_run_iteration_hook is not None:
                if (
                    self.log(
                        self._start_run_iteration_hook.call, with_timer, after=False
                    )(self._k)
                    is False
                ):
                    break

            self.log(self.import_births, with_timer)(
                body_class, shape_class, face_class, edge_class
            )

            self.log(self.execute_deaths, with_timer)()

            if controller is not None:
                controller.step()

            if friction_contact_trace_params is not None:
                self._osnspb._stepcounter = self._k

            if self._run_options.get("explode_computeOneStep_in_python"):
                if self._time_stepping_class == simu.TimeStepping:
                    self.log(self.computeOneStep_python, with_timer, after=False)(
                        with_timer
                    )
                else:
                    self.print_verbose(
                        "| [warning]. simulation of type",
                        self._time_stepping_clasimu.__name__,
                        " has no exploded version",
                    )
                    self.log(self._simulation.computeOneStep, with_timer)()
            else:
                if self.config.backend == "vnative":
                    self.log(self._simulation.updateInteractions, with_timer)()
                self.log(self._simulation.computeOneStep, with_timer)()

            number_of_contacts = self.log(self.contact_statistics_verbose, with_timer)()

            self.log(self.solver_verbose, with_timer)(number_of_contacts)

            cond = self._output_frequency and (self._k % self._output_frequency == 0)
            if cond or self._k == 1:
                if verbose:
                    self.print_verbose(
                        "output results in hdf5 file at step ",
                        self._k,
                    )

                self.log(self.output_results, with_timer)()

            if self._output_backup:
                if (self._k % self._output_backup_frequency == 0) or (self._k == 1):
                    # close io file, hdf5 memory is cleaned
                    self._out.close()
                    try:
                        shutil.copyfile(self._io_filename, self._io_filename_backup)
                    except shutil.Error as e:
                        siconos.io.tools.warn(str(e))
                    # open the file again
                    finally:
                        self.__enter__()

            self.log(self._simulation.clearNSDSChangeLog, with_timer)()

            if exit_tolerance:
                solver_options = self._osnspb.numericsSolverOptions()
                precision = solver_options.dparam[sn.params.SICONOS_DPARAM_RESIDU]
                if precision > exit_tolerance:
                    print("precision is larger exit_tolerance")
                    return False

            if self._before_next_step_iteration_hook is not None:
                if (
                    self.log(self._before_next_step_iteration_hook.call, with_timer)(
                        self._k
                    )
                    is False
                ):
                    break

            self.log(self._simulation.nextStep, with_timer)()

            if self._end_run_iteration_hook is not None:
                if (
                    self.log(self._end_run_iteration_hook.call, with_timer)(self._k)
                    is False
                ):
                    break

            self.print_verbose("")
            self._k += 1

        return True

    def output_timer_at_the_end(self):
        if len(self._timing) > 0:
            for k in self._timing.keys():
                siconos.io.mechanics_hdf5.group(self.log_data(), k)
                timing_data = np.array(self._timing[k])
                data_set = siconos.io.mechanics_hdf5.data(
                    self.log_data()[k], "timing", 1
                )
                current_line = data_set.shape[0]
                data_set.resize(current_line + timing_data.shape[0], 0)
                data_set[current_line : current_line + timing_data.shape[0]] = (
                    timing_data[:].reshape(timing_data.shape[0], 1)
                )

    def run(self, *args, **kwargs):

        build_from_kwargs = self.build_run_options_from_old_arguments_in_kwargs(
            *args, **kwargs
        )

        if build_from_kwargs:
            run_options_default = MechanicsHdf5Runner_run_options()
            print("run_options = MechanicsHdf5Runner_run_options()")
            for k in self._run_options.keys():
                if k in kwargs.keys():

                    # print('arg', kwargs[k],run_options_default[k] )
                    if kwargs[k] is not run_options_default[k]:
                        # print('diff', kwargs[k],run_options_default[k] )
                        print('run_options["{0}"]={1}'.format(k, kwargs[k]))

            # input('Enter a key to continue')

        with_timer = self._run_options.get("with_timer")

        info = self.log(self.run_initialize, with_timer)()
        info = self.log(self.run_loop, with_timer)()

        if with_timer and self._run_options["with_timer_output_at_the_end"]:
            self.output_timer_at_the_end()

        return info
