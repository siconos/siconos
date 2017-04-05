from math import cos, sin
from numpy.linalg import norm

class Material(object):
    """Some material properties that may be associated to shapes.

    Parameter
    ---------

    density: float
        The material density.
    """

    def __init__(self,
                 density):
        self.density = density


class MovedShape(object):
    """A shape with a reference to some shape data and that is moved by a
    given relative translation and relative orientation. The
    orientation may be given as a quaternion or in the form (axis,
    angle) where len(axis) == 3 and angle is a scalar.

    The default value for the translation is
    [0, 0, 0] and [1, 0, 0, 0] for the orientation.

    Parameters
    ----------

    shape_data
      some instance of an implementation of the shape.

    relative_translation: array_like of length 3
      translation in the bodyframe coordinates.

    relative_orientation: array_like of length 4 (quaternion) or (axis, angle)
      The orientation of the shape in the bodyframe coordinates.
      It may be expressed with a quaternion [w, x, y, z] or a couple
      ([x, y, z], angle)
    """
    def __init__(self,
                 shape_data,
                 relative_translation=[0, 0, 0],
                 relative_orientation=[1, 0, 0, 0]):

        if len(relative_orientation) == 2:
            # axis + angle
            axis = relative_orientation[0]
            assert len(axis) == 3
            angle = relative_orientation[1]
            assert type(angle) is float
            n = sin(angle / 2) / norm(axis)

            ori = [cos(angle / 2.), axis[0] * n, axis[1] * n, axis[2] * n]
        else:
            # a given quaternion
            assert len(relative_orientation) == 4
            ori = relative_orientation

        self.data = shape_data
        self.translation = relative_translation
        self.orientation = ori


class Shape(MovedShape):
    """
    A MovedShape with optional instance name.

    Parameters
    ----------

    shape_data
        Some instance of an implementation of the shape.

    instance_name: string, optional
        The name of the instance for further reference.

    relative_translation: array_like of length 3
        Translation in the bodyframe coordinates.

    relative_orientation: array_like of length 4 (quaternion) or (axis, angle)
        The orientation of the shape in the bodyframe coordinates.
        It may be expressed with a quaternion [w, x, y, z] or a couple
        ([x, y, z], angle)
    """

    def __init__(self,
                 shape_data,
                 instance_name=None,
                 relative_translation=[0, 0, 0],
                 relative_orientation=[1, 0, 0, 0]):

        self.instance_name = instance_name
        super(Shape, self).__init__(shape_data,
                                    relative_translation,
                                    relative_orientation)


class Volume(Shape):
    """A Shape with associated parameters.

    Parameters
    ----------

    shape_data
        some instance of an implementation of the shape.

    instance_name: string, optional
        The name of the instance for further reference.

    parameters
        The parameters associated to the volume.
        This may be informations about the material and its density.

    relative_translation: array_like of length 3
        translation in the bodyframe coordinates.

    relative_orientation: array_like of length 4 (quaternion) or (axis, angle)
        The orientation of the shape in the bodyframe coordinates.
        It may be expressed with a quaternion [w, x, y, z] or a couple
        ([x, y, z], angle)
    """

    def __init__(self,
                 shape_data,
                 instance_name=None,
                 parameters=Material(density=1),
                 relative_translation=[0, 0, 0],
                 relative_orientation=[1, 0, 0, 0]):
        self.instance_name = instance_name
        self.parameters = parameters
        super(Shape, self).__init__(shape_data,
                                    relative_translation,
                                    relative_orientation)


class Contactor(Shape):
    """A Contactor is a Shape that belongs to a collision group and may
    have associated parameters. Depending on the geometrical engine
    used, some informations may be added to the contactor, such as the
    kind of contact and the concerned part of the shape.  Note that
    contact laws must then be defined between collision groups.

    Parameters
    ----------

    shape_data
        Some instance of an implementation of the shape.

    instance_name: string, optional
        The name of the instance for further reference.

    collision_group: int
        The collision group, an integer >= 0.

    parameters: optional
        Parameters associated to the contactor.

    relative_translation: array_like of length 3
        Translation of in the bodyframe coordinates.

    relative_orientation: array_like of length 4 (quaternion) or (axis, angle)
        The orientation of the shape in the bodyframe coordinates.
        It may be expressed with a quaternion [w, x, y, z] or a couple
        ([x, y, z], angle)

    """

    def __init__(self,
                 shape_data,
                 instance_name=None,
                 collision_group=0,
                 parameters=None,
                 contact_type=None,
                 contact_index=None,
                 relative_translation=[0, 0, 0],
                 relative_orientation=[1, 0, 0, 0]):

        self.group = collision_group
        self.parameters = parameters
        self.contact_type = contact_type
        self.contact_index = contact_index

        super(Contactor, self).__init__(shape_data,
                                        instance_name,
                                        relative_translation,
                                        relative_orientation)


