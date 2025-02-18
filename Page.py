import Element
import Module


class Page:
    def __init__(self, generators: list[str], generator_bigrades: list[tuple[int, int]]):
        self.generators = generators
        self.generator_bigrades = generator_bigrades
        self.h_spaces: list[list[Module]] = []