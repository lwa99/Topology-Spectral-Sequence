class Scalar:
    """
    In out setting, scalars are elements in a fintie prime Field of order (and characteristic) c
    """
    def __init__(self, c: int, val):
        self.c = c
        if val is not None:
            self.val = val % self.c
        else:
            self.val = None

    # Basic Functionalities

    def __eq__(self, other):
        return self.c == other.c and (self.val == other.val)

    def __str__(self):
        return str(self.val) + " (mod " + str(self.c) + ")"

    # Modifiers
    def increase_by(self, s: 'Scalar'):
        self.val = (self.val + s.val) % self.c

    def multiply_by(self, s: 'Scalar'):
        self.val = (self.val * s.val) % self.c

    def take_inverse(self):
        self.val = pow(self.val, -1, self.c)

    def update(self, val):
        self.val = val % self.c

    # Operations that yields new objects
    def __add__(self, other: 'Scalar'):
        return Scalar(self.c, self.val + other.val)

    def __sub__(self, other):
        return Scalar(self.c, self.val - other.val)

    def __mul__(self, other):
        return Scalar(self.c, self.val * other.val)

    def get_inverse(self):
        return Scalar(self.c, pow(self.val, -1, self.c))

    def __truediv__(self, other: 'Scalar'):
        return Scalar(self.c, self.val * pow(other.val, -1, self.c))

    # Static method to produce empty scalar wrapper
    @staticmethod
    def get_empty_scalar(c: int):
        return Scalar(c, None)


if __name__ == "__main__":
    pass
