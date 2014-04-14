

class Contactor:
    def __init__(self,
                 shape_name,
                 collision_group,
                 relative_position,
                 relative_orientation):
        self.name = shape_name
        self.group = collision_group
        self.position = relative_position
        self.orientation = relative_orientation

