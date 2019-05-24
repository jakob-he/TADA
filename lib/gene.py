"""
The Gene class is an extension of the Bed class.
"""
from .bed import Bed


class Gene(Bed):
    def __init__(self, line, column_names):
        """A Class containing the information about one gene from a Bed file.
        Args:
            line: The line of bed file (This is just a string).
            column_names: Additional columns except chr,start and end in the bed file sorted by appearence. (e.g. gene name)
        Returns:
            A new Gene object.
        """
        super().__init__(line,column_names)
