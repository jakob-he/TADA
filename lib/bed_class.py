import enum

from . import tad, gene, enhancer, cnv, bed


class BedClass(enum.Enum):
    TAD = enum.auto()
    GENE = enum.auto()
    ENHANCER = enum.auto()
    CNV = enum.auto()
    BED = enum.auto()

    @classmethod
    def from_str(cls, name):
        """Get enum type from a string as a case-insensitive operation."""
        if isinstance(name, cls):
            return name
        name = name.upper()
        return cls[name]

    def get_class(self):
        """Get the associated class for the given value."""
        if self == self.TAD:
            return tad.Tad
        elif self == self.GENE:
            return gene.Gene
        elif self == self.ENHANCER:
            return enhancer.Enhancer
        elif self == self.CNV:
            return cnv.CNV
        elif self == self.BED:
            return bed.Bed
        else:
            raise TypeError(f"Type {self} has no associated class.")
