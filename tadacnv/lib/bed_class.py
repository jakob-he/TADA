import enum

from . import tad, cnv, bed, gene


class BedClass(enum.Enum):
    TAD = enum.auto()
    CNV = enum.auto()
    BED = enum.auto()
    GENE = enum.auto()

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
        elif self == self.CNV:
            return cnv.CNV
        elif self == self.BED:
            return bed.Bed
        elif self == self.GENE:
            return gene.Gene
        else:
            raise TypeError(f"Type {self} has no associated class.")
