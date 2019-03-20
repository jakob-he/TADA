"""
The CNV class is an extension of the Bed class.
"""
from .bed import Bed


class CNV(Bed):
    def __init__(self, line, column_names):
        """A Class containing the information about a copy number variant.
        Args:
            line: The line of bed file (This is just a string).
            column_names: Additional columns except chr,start and end in the bed file sorted by appearence. (e.g. sequence)
        Returns:
            A new Gene object.
        """
        super().__init__(line,column_names)
        self.tads = []
        self.boundary_spanning = False

    def __str__(self):
        """CNV Objects are representated by all the columns entered in the bed file"""
        tads = "\n".join([str(tad) for tad in self.tads])
        return f'{self.chr}\t{self.start}\t{self.end}\nTADS\n{tads}'
