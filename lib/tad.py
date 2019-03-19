"""
The Tad class is an extension of the bed class.
"""
from .bed import Bed

class Tad(Bed):
    def __init__(self, line):
        """Basic class containing the information from a BED file element.
        Args:
            line: The line of bed file (This is just a string).
        Returns:
            A new Bed object.
        """
        self.data = line.strip().split('\t')
        self.chr = self.data[0]
        self.start = int(self.data[1])
        self.end = int(self.data[2])
        self.data = {}
        self.genes = []
        self.enhancer = []

    def count_genes(self):
        return len(self.genes)

    def count_enhancer(self):
        return len(self.enhancer)
