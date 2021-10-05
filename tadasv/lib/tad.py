"""
The Tad class is an extension of the bed class.
"""
from .bed import Bed


class Tad(Bed):
    def __init__(self, line, column_names, **kwargs):
        """Basic class containing the information from a BED file element.
        Args:
            line: The line of bed file (This is just a string).
            column_names: Additional columns except chr,start and end in the bed file sorted by appearence. (e.g. enhancer name)
        Returns:
            A new Tad object.
        """
        super().__init__(line, column_names, **kwargs)
        self.annotations = {}

    def __str__(self):
        """An object of the Tad class is representeted by the chromsome, start and end position
        as well as the recursive representation of its annotations"""
        representation = f'{self.chr}\t{self.start}\t{self.end}\n'
        for name, annotations in self.annotations.items():
            annotation_repr = '\n'.join([str(annotation) for annotation in annotations])
            representation += f'\n{name}:\n{annotation_repr}'
        return representation
