"""
The Tad class is an extension of the bed class.
"""
from .bed import Bed


class Tad(Bed):
    def __init__(self, line, column_names):
        """Basic class containing the information from a BED file element.
        Args:
            line: The line of bed file (This is just a string).
            column_names: Additional columns except chr,start and end in the bed file sorted by appearence. (e.g. enhancer name)
        Returns:
            A new Tad object.
        """
        super().__init__(line, column_names)
        self.genes = []
        self.enhancer = []

    def __str__(self):
        """An object of the Tad class is representeted by the chromsome, start and end position
        as well as the recursive representation of its genes and enhancer"""
        genes = "\n".join([str(gene) for gene in self.genes])
        enhancer = "\n".join([str(enh) for enh in self.enhancer])
        return f'{self.chr}\t{self.start}\t{self.end}\nGenes:\n{genes}\nEnhancer:\n{enhancer}'

    def count_genes(self):
        return len(self.genes)

    def count_enhancer(self):
        return len(self.enhancer)
