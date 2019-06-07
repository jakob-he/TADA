"""
The CNV class is an extension of the Bed class.
"""
from .bed import Bed
from . import utils
import numpy as np

class CNV(Bed):
    def __init__(self, line, column_names, vcf=False):
        """A Class containing the information about a copy number variant.
        Args:
            line: The line of bed file (This is just a string).
            column_names: Additional columns except chr,start and end in the bed file sorted by appearence. (e.g. sequence)
        Returns:
            A new Gene object.
        """
        if vcf:
            self.data = line.strip().split('\t')
            location = self.data[1].split(':')
            self.chr = location[0]
            position = location[1].split('-')
            self.start = int(position[0])
            self.end = int(position[1])
        else:
            super().__init__(line,column_names)
        self.tads = []
        self.boundary_spanning = False
        self.annotation_distances = {}

    def __str__(self):
        """CNV Objects are representated by all the columns entered in the bed file"""
        tads = "\n".join([str(tad) for tad in self.tads])
        return f'{self.chr}\t{self.start}\t{self.end}\nTADS\n{tads}'

    def calculate_overlap_and_distances(self):
        """Calculates the distance and overlap for each annotation in the same TAD as the CNV. Currenlty modifed for binary overlap."""
        if self.tads:
            for tad in self.tads:
                for annotation_name, annotations in tad.annotations.items():
                    self.annotation_distances[annotation_name] = []
                    for annotation in annotations:
                        overlap = utils.getOverlap([self.start,self.end],[annotation.start,annotation.end])
                        if overlap == 0:
                            distance = 0
                            #distance = utils.getDistance([self.start,self.end],[annotation.start,annotation.end])
                        else:
                            distance = 1
                        self.annotation_distances[annotation_name].append(distance)

    def get_binary_features(self):
        """Returns binary features which are either directly derived from the TADs or based on the CNV itself.
        The output is a numpy boolean feature vector."""
        features = [self.boundary_spanning,any(overlap for overlap in self.gene_distances),any(overlap for overlap in self.enhancer_distances),any(tad.high_pLI for tad in self.tads),any(tad.high_Phastcon for tad in self.tads)]
        return np.array(features)
