"""
The Enhancer class is an extension of the Bed class.
"""
from .bed import Bed


class Enhancer(Bed):
    def __init__(self, line, column_names):
        """A Class containing the information about one enhancer from a Bed file.
        Args:
            line: The line of bed file (This is just a string).
            column_names: Additional columns except chr,start and end in the bed file sorted by appearence. (e.g. enhancer name)
        Returns:
            A new Enhancer object.
        """
        super().__init__(line,column_names)
