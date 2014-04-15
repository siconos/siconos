

class Contactor:
    def __init__(self,
                 shape_name,
                 collision_group=0,
                 relative_position=[0, 0, 0],
                 relative_orientation=[1, 0, 0, 0]):
        self.name = shape_name
        self.group = collision_group
        self.position = relative_position
        self.orientation = relative_orientation

