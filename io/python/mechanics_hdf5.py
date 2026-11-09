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
import os
import sys
import numpy as np
import pickle
import siconos.io.tools
import siconos.io.format_hdf5 as format_hdf5
from siconos.mechanics.quaternions import quaternion_get
import siconos.mechanics.collision.tools as smct
from siconos.mechanics import have_occ  # , have_bullet

# Mapping from joint type to expected number of points and axes.
# Format: joint_type -> (number_of_points, number_of_axes)
joint_points_axes = {
    "KneeJointR": (1, 0),
    "PivotJointR": (1, 1),
    "PrismaticJointR": (0, 1),
    "CylindricalJointR": (1, 1),
    "FixedJointR": (0, 0),
}


def check_points_axes(name, joint_class, points, axes):
    """
    Check that the number and shape of provided `points` and `axes`
    match what is expected for the given joint type.

    Parameters
    ----------
    name : str
        The name of the joint
    joint_class : str
        The type of joint. Must be one of the keys in `joint_points_axes`.
    points : np.ndarray, optional
        array of points
    axes : np.ndarray, optional
        array of axes
    """

    def check(x, idx):
        def er():
            n = joint_points_axes[joint_class][idx]
            item_type = "point" if idx == 0 else "axis"
            raise ValueError(f"{joint_class} ({name}) expects {n} {item_type}(got {x})")

        if np.shape(x) == (0,) or np.shape(x) == ():
            num = 0
        else:
            if len(np.shape(x)) != 2 or np.shape(x)[1] != 3:
                er()
            num = np.shape(x)[0]
        if (
            joint_class in joint_points_axes
            and joint_points_axes[joint_class][idx] != num
        ):
            er()

    check(points, 0)
    check(axes, 1)


def group(h, name, must_exist=True):
    try:
        return h[name]
    except KeyError:
        if must_exist:
            return h.create_group(name)
        else:
            try:
                return h.create_group(name)
            except ValueError:
                # could not create group, return None
                # (file is probably in read-only mode)
                return None


def data(h, name, nbcolumns, use_compression=False):
    try:
        return h[name]
    except KeyError:
        comp = use_compression and nbcolumns > 0
        return h.create_dataset(
            name,
            (0, nbcolumns),
            maxshape=(None, nbcolumns),
            chunks=[None, (4000, nbcolumns)][comp],
            compression=[None, "gzip"][comp],
            compression_opts=[None, 9][comp],
            dtype='f4'
        )


def text_data(text):
    """A variable-length string stored as a (1,) dataset (shapes, plugins)."""
    return np.array([text], dtype=format_hdf5.vlen_dtype)


def add_line(dataset, line):
    dataset.resize(dataset.shape[0] + 1, 0)
    dataset[dataset.shape[0] - 1, :] = line


#
# misc fixes
#
# fix ctr.'name' in old hdf5 files
#
def upgrade_io_format(filename):

    with MechanicsHdf5(filename, mode="a") as io:

        for instance_name in io.instances():
            for contactor_instance_name in io.instances()[instance_name]:
                contactor = io.instances()[instance_name][contactor_instance_name]
                if "name" in contactor.attrs:
                    siconos.io.tools.warn(
                        """
contactor {0} attribute 'name': renamed in 'shape_name'
                    """
                    )
                    contactor.attrs["shape_name"] = contactor["name"]
                    del contactor["name"]


class MechanicsHdf5(object):
    """a MechanicsHdf5 context manager, used to prepare a simulation description
    to be executed by MechanicsRunner.

    Parameters
    ----------
    io_filename: string, optional
        hdf5 file name, default = <caller>.hdf5, caller being the name
        without ext of the file that instanciates the Runner.
    mode: string, optional
        h5 mode (w, r, append), default = 'w'
    io_filename_backup: string, optional
        name of a backup (copy) file for hdf5 outputs.
        Backup every <output_frequency> step. Default = <caller>_last.hdf5
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
        io_filename=None,
        mode="w",
        io_filename_backup=None,
        use_compression=False,
        output_domains=False,
        verbose=True,
    ):
        if io_filename is None:
            self._io_filename = "{0}.hdf5".format(
                os.path.splitext(os.path.basename(sys.argv[0]))[0]
            )
        else:
            self._io_filename = io_filename

        if io_filename_backup is None:
            self._io_filename_backup = "{0}_last.hdf5".format(
                os.path.splitext(self._io_filename)[0]
            )
        else:
            self._io_filename_backup = io_filename_backup

        self._output_backup = False
        self._mode = mode
        self._run_options = None
        self._nslaws = dict()
        self._adaptor = None
        self._out = None
        self._data = None
        self._number_of_shapes = 0
        self._number_of_permanent_interactions = 0
        self._number_of_dynamic_objects = 0
        self._number_of_static_objects = 0
        self._use_compression = use_compression
        self._should_output_domains = output_domains
        self._verbose = verbose
        self._fem_dof_mappings = {}  # ds_id -> {'dof_indices': array, 'spatial_order': array, 'coords': array}

    def __enter__(self):
        """Reminder: this function will be called when a 'with'
        statement
        will be executed with the present class.

        Warning : it means that this class must be called inside
        a with statement to be properly initialized!
        """
        # -- Opens the hdf5 file and its content, as described by the format --
        # created if missing (write modes); in read mode, a missing entry is None
        # if the format says it is optional, else the reading fails.
        # domain is only created when the domains are output.
        self._adaptor = format_hdf5.FormatAdaptor(
            self._io_filename, self._mode, compression=self._use_compression
        )
        self._adaptor.open(exclude=() if self._should_output_domains else ("domain",))
        self._out, self._data, self._dimension = self._adaptor.hdf5, self._adaptor.data, self._adaptor.dimension
        return self

    def __exit__(self, type_, value, traceback):
        self._adaptor.close()

    def print_verbose(self, *args, **kwargs):
        if self._verbose:
            print("[io.mechanics]", *args, **kwargs)

    def print_io_mechanics(self, *args, **kwargs):
        print("[io.mechanics]", *args, **kwargs)

    # hdf5 structure
    def dimension(self):
        """
        dimension : get the dimension (2 or 3) of the scene
        """
        return self._dimension

    def shapes(self):
        """
        Shapes : parameterized primitives or user defined
                 (convex set or meshes)
        """
        return self._adaptor["ref"]

    def permanent_interactions(self):
        """
        Permanent interactions.
        """
        return self._adaptor["permanent_interactions"]

    def radii_data(self):
        """
        Radii of objects
        """
        return self._adaptor["radii"]

    def p0s_data(self):
        """
        Vector p0 for each object
        """
        return self._adaptor["p0s"]

    def static_data(self):
        """
        Coordinates and orientations of static objects.
        """
        return self._adaptor["static"]

    def dynamic_data(self):
        """
        Coordinates and orientations of dynamic objects.
        """
        return self._adaptor["dynamic"]

    def displacement_data(self, fem_id):
        """
        Displacements
        """
        nbcolumns = self._data[f"fem_displacements_{fem_id}"].attrs["dof_count"]
        return data(self._data, f"fem_displacements_{fem_id}", nbcolumns,
                    use_compression=False)

    def velocities_data(self):
        """
        Velocities of dynamic objects
        """
        return self._adaptor["velocities"]

    def contact_forces_data(self):
        """
        Contact points information.
        """
        return self._adaptor["cf"]

    def contact_internal_variable_data(self):
        """
        Contact internal variables information.
        """
        return self._adaptor["cf_internal_variables"]

    def contact_info_data(self):
        """
        Contact points information.
        """
        return self._adaptor["cf_info"]

    def contact_work_data(self):
        """
        Contact points information.
        """
        return self._adaptor["cf_work"]

    def domains_data(self):
        """
        Contact point domain information.
        """
        return self._adaptor["domain"]

    def solver_data(self):
        """
        Solver output
        """
        return self._adaptor["solv"]

    def log_data(self):
        """
        log output
        """
        return self._adaptor["log"]

    def instances(self):
        """
        Scene objects.
        """
        return self._adaptor["input"]

    def nonsmooth_laws(self):
        """
        Non smooth laws between group of contactors.
        """
        return self._adaptor["nslaws"]

    def joints(self):
        """
        Joints between dynamic objects or between an object and the scenery.
        """
        return self._adaptor["joints"]

    def boundary_conditions(self):
        """
        Boundary conditions applied to  dynamic objects
        """
        return self._adaptor["boundary_conditions"]

    def p0s(self, step):
        """Get p0 data for a specific time step (fixed-dof systems only)."""
        if self._adaptor["p0s"] is None:
            return None
        # Get dynamic objects at this step to know count
        dyn_data = self.dynamic_data()
        if dyn_data is None or len(dyn_data) == 0:
            return None
        # Find rows for this step
        step_times = dyn_data[:, 0]
        # Find indices for this step (assuming step corresponds to index in unique times)
        unique_times, indices = np.unique(step_times, return_index=True)
        if step >= len(indices):
            return None
        start_idx = indices[step]
        if step + 1 < len(indices):
            end_idx = indices[step + 1]
        else:
            end_idx = len(dyn_data)
        # Count dynamic objects at this step
        n_dyn = end_idx - start_idx
        # p0s_data is appended per step, so calculate offset
        # This assumes constant number of fixed-dof objects per step
        # Better: store count per step, but for now assume fixed
        p0s_all = self._adaptor["p0s"][:]
        if len(p0s_all) == 0:
            return None
        # Estimate entries per step
        entries_per_step = len(p0s_all) // len(unique_times)
        p0_start = step * entries_per_step
        p0_end = min(p0_start + entries_per_step, len(p0s_all))
        return p0s_all[p0_start:p0_end]

    def radii(self, step):
        """Get radii for a specific time step (fixed-dof systems only)."""
        if self._adaptor["radii"] is None:
            return None
        radii_all = self._adaptor["radii"][:]
        if len(radii_all) == 0:
            return None
        dyn_data = self.dynamic_data()
        if dyn_data is None:
            return None
        unique_times = np.unique(dyn_data[:, 0])
        entries_per_step = len(radii_all) // len(unique_times)
        r_start = step * entries_per_step
        r_end = min(r_start + entries_per_step, len(radii_all))
        return radii_all[r_start:r_end]


    def _add_shape(self, name, data, **attrs):
        """Add a shape to ref with the next shape id (nothing if it exists)."""
        if name not in self._adaptor["ref"]:
            self._adaptor.add_into_group("ref", name, data=data, id=self._number_of_shapes, **attrs)
            self._number_of_shapes += 1

    def add_plugin_source(self, name, filename):
        """
        Add C source plugin
        """
        if name not in self._adaptor["plugins"]:
            self._adaptor.add_into_group("plugins", name, filename=filename,
                                     data=text_data(siconos.io.tools.str_of_file(filename)))

    def add_external_function(
        self, name, body_name, function_name, plugin_name, plugin_function_name
    ):
        self._adaptor.add_into_group("external_functions", name, body_name=body_name,
                                 function_name=function_name, plugin_name=plugin_name,
                                 plugin_function_name=plugin_function_name)

    def add_external_bc_function(
        self, name, body_name, bc_indices, plugin_name, plugin_function_name
    ):
        self._adaptor.add_into_group("external_functions", name, body_name=body_name,
                                 bc_indices=bc_indices, plugin_name=plugin_name,
                                 plugin_function_name=plugin_function_name)

    def add_mesh_from_string(
        self, name, shape_data, scale=None, insideMargin=None, outsideMargin=None
    ):
        """
        Add a mesh shape from a string.
        Accepted format : mesh encoded in VTK .vtp format
        """
        self._add_shape(name, text_data(shape_data), type="vtp", scale=scale,
                        insideMargin=insideMargin, outsideMargin=outsideMargin)

    def add_mesh_from_file(
        self, name, filename, scale=None, insideMargin=None, outsideMargin=None
    ):
        """Add a mesh shape from a file.
        Accepted format : .stl or mesh encoded in VTK .vtp format
        """

        import vtk

        if filename[0] != os.path.sep:
            filename = os.path.join(
                os.path.split(os.path.abspath(sys.argv[0]))[0], filename
            )
        if name not in self._adaptor["ref"]:

            if os.path.splitext(filename)[-1][1:] == "stl":
                reader = vtk.vtkSTLReader()
                reader.SetFileName(filename)
                reader.Update()

                if reader.GetErrorCode() != 0:
                    print("vtkSTLReader error", reader.GetErrorCode())
                    sys.exit(1)

                with siconos.io.tools.tmpfile() as tmpf:
                    writer = vtk.vtkXMLPolyDataWriter()
                    writer.SetInputData(reader.GetOutput())
                    writer.SetFileName(tmpf[1])
                    writer.Write()

                    shape_data = siconos.io.tools.str_of_file(tmpf[1])

            else:
                assert os.path.splitext(filename)[-1][1:] == "vtp"
                shape_data = siconos.io.tools.str_of_file(filename)

            self.add_mesh_from_string(
                name,
                shape_data,
                scale=scale,
                insideMargin=insideMargin,
                outsideMargin=outsideMargin,
            )

    def add_height_map(
        self, name, heightmap, rectangle, insideMargin=None, outsideMargin=None
    ):
        """
        Add a heightmap represented as a SiconosMatrix
        """
        assert heightmap.shape[0] >= 2
        assert heightmap.shape[1] >= 2
        # rect: lengths (x, y) of the rectangle where the heightmap is placed
        self._add_shape(name, heightmap, type="heightmap", rect=rectangle,
                        insideMargin=insideMargin, outsideMargin=outsideMargin)

    def add_brep_from_string(self, name, shape_data):
        """
        Add a brep contained in a string.
        """
        # raw str, or __getstate__ as with pythonocc
        text, occ_indx = (shape_data, None) if type(shape_data) is str else shape_data[:2]
        self._add_shape(name, text_data(text), type="brep", occ_indx=occ_indx)

    def add_occ_shape(self, name, occ_shape):
        """
        Adds an OpenCascade TopoDS_Shape to the hdf5 group 'ref'

        It takes an OpenCascade shape (e.g., a sphere, box, etc.),
        converts it into the STEP format,
        and stores it as a dataset

        Parameters:
        -----------
        name : str
            The name associated with the shape. This name will be used as the
            identifier for the shape in the internal storage.
            If the name already exists in the storage, the shape will not be added again.

        occ_shape : TopoDS_Shape
            The OpenCascade shape to be added, expected to be an instance
            of the `TopoDS_Shape` class.
        """
        if name in self._adaptor["ref"]:
            return

        from OCC.Core.STEPControl import STEPControl_Writer, STEPControl_AsIs

        # step format is used for the storage.
        step_writer = STEPControl_Writer()
        step_writer.Transfer(occ_shape, STEPControl_AsIs)

        shape_data = None

        with siconos.io.tools.tmpfile() as tmpf:

            step_writer.Write(tmpf[1])

            tmpf[0].flush()
            shape_data = siconos.io.tools.str_of_file(tmpf[1])
        self._add_shape(name, text_data(shape_data), type="step")

    def add_shape_data_from_file(self, name, filename):
        """
        Add shape data from a file.
        """
        if name not in self._adaptor["ref"]:
            self._add_shape(name, text_data(siconos.io.tools.str_of_file(filename)),
                            type=os.path.splitext(filename)[1][1:].lower() or "unknown")

    def add_interaction(
        self,
        name,
        body1_name,
        contactor1_name=None,
        body2_name=None,
        contactor2_name=None,
        distance_calculator="cadmbtb",
        offset1=0.0,
        offset2=0.0,
    ):
        """
        Add permanent interactions between two objects contactors.
        """
        if name not in self._adaptor["permanent_interactions"]:
            self._adaptor.add_into_group(
                "permanent_interactions", name, id=self._number_of_permanent_interactions,
                type="permanent_interaction", body1_name=body1_name, body2_name=body2_name,
                contactor1_name=contactor1_name, contactor2_name=contactor2_name,
                distance_calculator=distance_calculator, offset1=offset1, offset2=offset2)
            self._number_of_permanent_interactions += 1

    def add_convex_shape(
        self,
        name,
        points,
        insideMargin=None,
        outsideMargin=None,
        avoid_internal_edge_contact=False,
    ):
        """
        Add a convex shape defined by a list of points.

        outsideMargin is the value of margin that substract
        from the actual contact distance

        If insideMargin is positive, the convex hull is shrunken by that amount
        (each face is moved by "shrink" length units towards
        the center along its normal).
        This value is then added to outsideMargin to compensate the shrink.
        The convex hull appears for the user as in its original size.



        """
        # infer the dimension of the problem
        if np.shape(points)[1] == 2:
            self._dimension = 2
        else:
            if self._dimension == 2:
                raise ValueError(
                    "It is not yet possible to mix 2D and 3D primitives shapes"
                )
            self._dimension == 3
        self._out.attrs["dimension"] = self._dimension

        self._add_shape(name, np.asarray(points, dtype="f4"), type="convex",
                        insideMargin=insideMargin, outsideMargin=outsideMargin,
                        avoid_internal_edge_contact=True if avoid_internal_edge_contact else None)

    def add_primitive_shape(
        self, name, primitive, params, insideMargin=None, outsideMargin=None
    ):
        """
        Add a primitive shape.

        Todo: Comments on insideMargin and outsideMargin that depends on the
              primitive shape type.
        """
        # infer the dimension of the problem
        if (
            primitive == "Disk"
            or primitive == "Box2d"
            or primitive == "Line"
            or primitive == "Circle"
            or primitive == "Segment"
        ):
            self._dimension = 2
        else:
            if self._dimension == 2:
                raise ValueError(
                    "It is not yet possible to mix 2D and 3D primitives shapes"
                )
            self._dimension == 3
        self._out.attrs["dimension"] = self._dimension

        self._add_shape(name, np.asarray([params], dtype="f4"), type="primitive", primitive=primitive,
                        insideMargin=insideMargin, outsideMargin=outsideMargin)

    def add_object(
        self,
        name,
        shapes,
        translation=None,
        orientation=None,
        velocity=None,
        use_volume_centroid_as_initial_translation=False,
        material=None,
        boundary_conditions=None,
        nodal_forces=None,
        mass=None,
        center_of_mass=None,
        inertia=None,
        time_of_birth=-1,
        time_of_death=-1,
        fext=None,
        allow_self_collide=False,
    ):
        """Add an object with associated shapes as a list of Volume or
        Contactor objects. Contact detection and processing is
        defined by the Contactor objects. The Volume objects are used for
        the computation of inertia and center of mass if not provided.

        The  body-fixed frame is assumed to be the global inertial
        frame. This means that
        1. By default, the center of mass is located at the origin.
        The initial translation is applied from this point, so that x_g(0) = translation
        2. the  orientation is identical to the inertial frame.
        The initial orientation is applied to the inertial frame to obtain
        the body-fixed frame.


        Each Contactor and Volume object may have a relative
        translation and a relative orientation expressed in the bodyframe
        coordinates.

        Parameters
        ----------
        name: string
            The name of the object.

        shapes: iterable
            The list of associated Contactor or Volume objects.

        translation: array_like of length 3 or 2 (dimension =2)
            Initial translation of the object (mandatory)

        orientation: array_like of length 3 (Euler Angles) or 4
            (unit quaternion), or 1 (dimension =2)
            Initial orientiation of the object. By default, identity

        velocity: array_like of length 6,  or 3 (dimension =2)
            Initial velocity of the object. The default velocity is zero.
            dimension =3 :
               The components are those of the translation velocity along
               x, y and z axis and the rotation velocity around x, y and
               z axis.
            dimension =2 :
               The components are those of the translation velocity along
               x, y and the rotation velocity z axis.

        material: array_like of length 3.
            (density, elastic Young modulus, poisson ratio)
            The material of the finite elements object. If material is
            given, mass is directly computed and the mass parameter is not
            used.

        mass: float
            The mass of the object, if it is None the object is defined as
            a static object involved only in contact detection.
            The default value is None.

        center_of_mass: array_like of length 3
            The position of the center of mass expressed in the body frame
            coordinates.

        inertia: array_like of length 3 or 3x3 matrix.
            The principal moments of inertia (array of length 3) or
            a full 3x3 inertia matrix

        use_volume_centroid_as_initial_translation: boolean.
            if True and if a Volume is given is the list of shape,
            the position of
            the volume centroid is used as initial translation.

        """

        # Does nothing if an object of the same name is already registered.
        if name in self._adaptor["input"]:
            return

        # Ensure all input parameters have proper values
        # (set default values, convert tuple or list to numpy arrays ...)

        rigid = material is None  # material only for the elastic case

        translation = np.asarray(translation, dtype=np.float64)

        if center_of_mass is None:
            center_of_mass = np.zeros(3, dtype=np.float64)

        if rigid and self._dimension == 3:
            if orientation is None:
                orientation = np.array([1, 0, 0, 0], dtype=np.float64)
            if velocity is None:
                velocity = np.zeros(6, dtype=np.float64)

            ori = quaternion_get(orientation)
            assert len(translation) == 3
            assert len(ori) == 4

        elif rigid and self._dimension == 2:
            if orientation is None:
                orientation = np.zeros(1, dtype=np.float64)
            if velocity is None:
                velocity = np.zeros(3, dtype=np.float64)
            assert len(translation) == 2
            ori = np.asarray(orientation, dtype=np.float64)

        if not rigid:
            if orientation is None:
                orientation = 0.0
            if velocity is None:
                velocity = 0.0
            if translation is None:
                translation = 0.0
            ori = np.atleast_1d(np.asarray(orientation, dtype=np.float64))
            translation = np.atleast_1d(np.asarray(translation, dtype=np.float64))
            velocity = np.atleast_1d(np.asarray(velocity, dtype=np.float64))

        is_center_of_mass_computed = False

        if (inertia is None) or (mass is None):
            # if :
            # - we need to compute mass or inertia
            # - and a Volume is present in the shapes list
            # - and occ is available
            if any(
                map(
                    lambda s: isinstance(s, smct.Volume),
                    shapes,
                )
            ):

                # --- We compute inertia and center of mass ---
                # This part requires occ.
                assert have_occ
                import siconos.io.occ_tools

                # get the list of volumes (as an iterator)
                volumes = filter(
                    lambda s: isinstance(s, smct.Volume),
                    shapes,
                )

                # compute mass, inertia for these volumes
                computed_mass, com, computed_inertia, computed_inertia_matrix = (
                    siconos.io.occ_tools.compute_inertia_and_center_of_mass(
                        volumes, self
                    )
                )

                self.print_verbose("{0}: computed mass from Volume".format(name))
                self.print_verbose(
                    "{0}: computed center of mass:".format(name),
                    com[0],
                    com[1],
                    com[2],
                )
                self.print_verbose("{0}: computed mass:".format(name), computed_mass)
                self.print_verbose(
                    "{0}: computed inertia:".format(name),
                    computed_inertia[0],
                    computed_inertia[1],
                    computed_inertia[2],
                )
                self.print_verbose(
                    "{0}: computed inertia matrix:".format(name),
                    computed_inertia_matrix,
                )
                is_center_of_mass_computed = True
                if mass is None:
                    mass = computed_mass

                if inertia is None:
                    inertia = computed_inertia_matrix

        if use_volume_centroid_as_initial_translation and is_center_of_mass_computed:
            translation = com
            for s in shapes:
                s.translation = s.translation - com

        if mass is not None and np.isscalar(mass) and mass <= 0.0:
            self.print_verbose(
                "The use of a mass equal to zero"
                " to define a static object is deprecated."
            )
            self.print_verbose("Set mass=None to define a static object")

        static = (mass is None or mass == 0) and material is None
        obj = self._adaptor.add_into_group(
            "input", name,
            id=-(self._number_of_static_objects + 1) if static else self._number_of_dynamic_objects + 1,
            type="dynamic" if mass is not None else "static",
            time_of_birth=time_of_birth if time_of_birth >= 0 else None,
            time_of_death=time_of_death if time_of_death >= 0 else None,
            material=None if material is None else np.asarray(material, dtype=np.float64),
            boundary_conditions=None if boundary_conditions is None
            else np.asarray(boundary_conditions, dtype=np.int32),
            nodal_forces=None if nodal_forces is None else np.asarray(nodal_forces, dtype=np.float64),
            mass=mass,
            translation=translation,
            orientation=ori,
            velocity=np.asarray(velocity, dtype=np.float64),
            center_of_mass=np.asarray(center_of_mass, dtype=np.float64),
            inertia=inertia,
            fext=None if fext is None else np.asarray(fext, dtype=np.float64),
            allow_self_collide=allow_self_collide,
        )

        for num, ctor in enumerate(shapes):
            # default name for a contactor
            instance_name = ctor.instance_name or "{0}-{1}".format(ctor.shape_name, num)
            parameters = getattr(ctor, "parameters", None)
            self._adaptor.add_into_group(
                "input/" + name, instance_name,
                instance_name=instance_name,
                shape_name=ctor.shape_name,
                group=ctor.collision_group if hasattr(ctor, "group") else None,
                # np.void to write the pickled parameters, see http://docs.h5py.org/en/latest/strings.html
                parameters=None if parameters is None else np.void(pickle.dumps(parameters)),
                type=getattr(ctor, "contact_type", None),
                contact_index=getattr(ctor, "contact_index", None),
                translation=ctor.translation,
                orientation=quaternion_get(ctor.orientation),
            )

        if static:
            self._number_of_static_objects += 1
        else:
            self._number_of_dynamic_objects += 1

        return obj

    def add_objects(
        self,
        name,
        shapes,
        translations,
        orientations=None,
        velocities=None,
        mass=None,
        inertia=None,
        time_of_birth=-1,
        time_of_death=-1,
    ):
        """
        Add an aggregate of objects (e.g., granular material) with shared properties.

        All objects share the same shape, mass and inertia. Only translations,
        orientations and velocities are specific to each grain.

        Parameters
        ----------
        name: string
            The name of the aggregate.
        shapes: Contactor or list of Contactors
            The shape(s) for all grains (only first contactor used for now).
        translations: array_like (N, 2) or (N, 3)
            Initial translations for all N grains.
        orientations: array_like (N,) for 2D or (N, 4) for 3D, optional
            Initial orientations. Default: zero rotation.
        velocities: array_like (N, 3) for 2D or (N, 6) for 3D, optional
            Initial velocities. Default: zero.
        mass: float, optional
            Mass for each grain (all grains have same mass).
        inertia: array_like, optional
            Inertia for each grain (all grains have same inertia).
        """
        if name in self._adaptor["input"]:
            return

        # shape info, the same for all grains
        shape_info = shapes[0] if isinstance(shapes, (list, tuple)) else shapes

        translations = np.asarray(translations, dtype=np.float64)
        if self._dimension == 2 and translations.shape[-1] != 2:
            raise ValueError(
                f"Expected Nx2 translations for 2D, got {translations.shape}"
            )
        if self._dimension == 3 and translations.shape[-1] != 3:
            raise ValueError(
                f"Expected Nx3 translations for 3D, got {translations.shape}"
            )

        if orientations is None:
            if self._dimension == 2:
                orientations = np.zeros(len(translations), dtype=np.float64)
            else:
                orientations = np.tile(
                    np.array([1.0, 0.0, 0.0, 0.0], dtype=np.float64),
                    (len(translations), 1),
                )
        else:
            orientations = np.asarray(orientations, dtype=np.float64)
            if self._dimension == 2 and orientations.ndim == 1:
                orientations = orientations.reshape(-1, 1)

        if velocities is None:
            velocities = np.zeros((len(translations), 3 if self._dimension == 2 else 6),
                                  dtype=np.float64)
        else:
            velocities = np.asarray(velocities, dtype=np.float64)

        obj = self._adaptor.add_into_group(
            "input", name,
            type="dynamic_aggregate" if mass is not None else "static_aggregate",
            id=self._number_of_dynamic_objects + 1 if mass is not None
            else -(self._number_of_static_objects + 1),
            shape_name=shape_info.shape_name,
            group=getattr(shape_info, "collision_group", getattr(shape_info, "group", None)),
            mass=mass,
            inertia=None if inertia is None else np.asarray(inertia, dtype=np.float64),
            count=len(translations),
            time_of_birth=time_of_birth if time_of_birth >= 0 else None,
            time_of_death=time_of_death if time_of_death >= 0 else None,
        )
        for key, values in (("translations", translations), ("orientations", orientations),
                            ("velocities", velocities)):
            self._adaptor.add_into_group("input/" + name, key, data=values,
                                     compression=self._use_compression)

        if mass is not None:
            self._number_of_dynamic_objects += len(translations)
        else:
            self._number_of_static_objects += len(translations)

        return obj

    def _add_nslaw(self, name, nslaw_type, **attrs):
        """Add a nonsmooth law, or update the parameters of the existing one."""
        if name in self._adaptor["nslaws"] and self._adaptor["nslaws"][name].attrs["type"] != nslaw_type:
            self.print_verbose(
                "[warning] a nslaw is already existing with the same name ",
                name,
                " but not the same type",
            )
        if name in self._adaptor["nslaws"]:  # update the parameters (the type is kept), checked against the format
            spec = self._adaptor._spec("nslaws/" + name)["attrs"]
            for k, v in attrs.items():
                format_hdf5.check_datatype(v, spec[k], self._adaptor.symbols)
            self._adaptor["nslaws"][name].attrs.update(attrs)
        else:
            self._adaptor.add_into_group("nslaws", name, type=nslaw_type, **attrs)

    def add_Newton_impact_rolling_friction_nsl(
        self, name, mu, mu_r, e=0, collision_group1=0, collision_group2=0
    ):
        """
        Add a nonsmooth law for contact between 2 groups.
        Only NewtonImpactRollingFrictionNSL are supported.
        name is an user identifiant and must be unique,
        mu is the coefficient of friction,
        e is the coefficient of restitution on the contact normal,
        gid1 and gid2 define the group identifiants.

        """
        self._add_nslaw(name, "NewtonImpactRollingFrictionNSL", mu=mu, mu_r=mu_r, e=e, gid1=collision_group1, gid2=collision_group2)

    def add_binary_cohesive_nsl(self, name, mu,  e=0, sigma_c=0, delta_c=0,
                                gamma=0.0,
                                collision_group1=0,
                                collision_group2=0):
        """
        Add a nonsmooth law for contact between 2 groups.
        Only BinaryCohesionNSL are supported.
        name is an user identifiant and must be unique,
        mu is the coefficient of friction,
        e is the coefficient of restitution on the contact normal,
        gamma is the ratio of tangent to normal cohesive strength,
        gid1 and gid2 define the group identifiants.

        """
        self._adaptor.add_into_group("nslaws", name, type="BinaryCohesiveNSL", mu=mu, e=e,
                                 sigma_c=sigma_c, delta_c=delta_c, gamma=gamma,
                                 gid1=collision_group1, gid2=collision_group2)

    def add_Fremond_impact_friction_nsl(
        self, name, mu, e=0, collision_group1=0, collision_group2=0
    ):
        """
        Add a nonsmooth law for contact between 2 groups.
        Only NewtonImpactFrictionNSL are supported.
        name is an user identifiant and must be unique,
        mu is the coefficient of friction,
        e is the coefficient of restitution on the contact normal,
        gid1 and gid2 define the group identifiants.

        """
        self._add_nslaw(name, "FremondImpactFrictionNSL", mu=mu, e=e, gid1=collision_group1, gid2=collision_group2)

    def add_Newton_impact_friction_nsl(
        self, name, mu, e=0, collision_group1=0, collision_group2=0
    ):
        """
        Add a nonsmooth law for contact between 2 groups.
        Only NewtonImpactFrictionNSL are supported.
        name is an user identifiant and must be unique,
        mu is the coefficient of friction,
        e is the coefficient of restitution on the contact normal,
        gid1 and gid2 define the group identifiants.

        """
        self._add_nslaw(name, "NewtonImpactFrictionNSL", mu=mu, e=e, gid1=collision_group1, gid2=collision_group2)

    # Note, default groups are -1 here, indicating not to add them to
    # the nslaw lookup table for contacts, since 1D impacts are
    # useless in this case.  They are however useful for joint stops.
    def add_Newton_impact_nsl(
        self, name, e=0, collision_group1=-1, collision_group2=-1
    ):
        """
        Add a nonsmooth law for contact between 2 groups.
        Only NewtonImpactNSL are supported.
        name is a user identifier and must be unique,
        e is the coefficient of restitution on the contact normal,
        gid1 and gid2 define the group identifiers.

        As opposed to add_Newton_impact_friction_nsl, the default groups are
        -1, making the NSL unassociated with point contacts.  It can
        by used for joint stops however.
        """
        self._add_nslaw(name, "NewtonImpactNSL", e=e, gid1=collision_group1, gid2=collision_group2)

    # Note, default groups are -1 here, indicating not to add them to
    # the nslaw lookup table for contacts, since 1D impacts are
    # useless in this case.  They are however useful for joint friction.
    def add_relay_nsl(
        self, name, lb, ub, size=1, collision_group1=-1, collision_group2=-1
    ):
        """
        Add a nonsmooth law for contact between 2 groups.
        Only NewtonImpactNSL are supported.
        name is a user identifier and must be unique,
        e is the coefficient of restitution on the contact normal,
        gid1 and gid2 define the group identifiers.

        As opposed to add_Newton_impact_friction_nsl, the default groups are
        -1, making the NSL unassociated with point contacts.  It can
        by used for joint stops however.
        """
        self._add_nslaw(name, "RelayNSL", size=size, lb=lb, ub=ub, gid1=collision_group1, gid2=collision_group2)


    def add_joint(
        self,
        name,
        object1,
        object2=None,
        points=[[0, 0, 0]],
        axes=[[0, 1, 0]],
        joint_class="PivotJointR",
        absolute=None,
        allow_self_collide=None,
        nslaws=None,
        stops=None,
        friction=None,
        coupled=None,
        references=None,
    ):
        """
        add a joint between two objects
        """
        check_points_axes(name, joint_class, points, axes)
        if coupled is not None:
            # must be a list of tuples of two integers (DoF indexes) and a float (ratio)
            for c in coupled:
                assert len(c) == 3
        if references is not None:
            # must be a list of two joint names and one DS name
            assert len(references) == 2 or len(references) == 3
        attrs = dict(
            object1=object1, object2=object2, type=joint_class,
            points=points, axes=axes,
            absolute=absolute if absolute in [True, False] else None,
            allow_self_collide=allow_self_collide if allow_self_collide in [True, False] else None,
            # nslaws: name of one nslaw, or a list of names same length as stops
            nslaws=None if nslaws is None else np.array(nslaws, dtype="S"),
            stops=stops,  # table of [[axis, pos, dir]..]
            # friction: an NSL name (e.g. RelayNSL), or list of same
            friction=None if friction is None else np.array(friction, dtype="S"),
            coupled=None if coupled is None else np.array(coupled),
            references=None if references is None else np.array(references, dtype="S"),
        )
        if name in self._adaptor["joints"]:  # update the existing joint, checked against the format
            spec = self._adaptor._spec("joints/" + name)["attrs"]
            attrs = {k: v for k, v in attrs.items() if v is not None}
            for k, v in attrs.items():
                format_hdf5.check_datatype(v, spec[k], self._adaptor.symbols)
            self._adaptor["joints"][name].attrs.update(attrs)
        else:
            self._adaptor.add_into_group("joints", name, **attrs)

    def add_boundary_condition(
        self,
        name,
        object1,
        indices=None,
        bc_class="HarmonicBC",
        v=None,
        a=None,
        b=None,
        omega=None,
        phi=None,
    ):
        """
        add boundarycondition to the object object1

        implementation only works for HarmonicBC for the moment
        """
        if name not in self.boundary_conditions():
            if bc_class == "HarmonicBC":
                params = dict(a=a, b=b, omega=omega, phi=phi)
            elif bc_class == "BoundaryCondition":
                params = dict(v=v)  # v might be None for "fixed" (val=0) BC
            else:
                raise NotImplementedError
            self._adaptor.add_into_group("boundary_conditions", name, object1=object1,
                                     indices=indices, type=bc_class, **params)
