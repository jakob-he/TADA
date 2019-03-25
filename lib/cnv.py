"""
The CNV class is an extension of the Bed class.
"""
from .bed import Bed
from . import utils


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

    def get_genes(self):
        genes = []
        for tad in self.tads:
            genes.extend(tad.genes)
        return genes

    def get_enhancer(self):
        enhancer = []
        for tad in self.tads:
            enhancer.extend(tad.enhancer)
        return enhancer

    def calculate_overlap_and_distances(self):
        """Calculates the distance and overlap for each gene and enhancer in the same TAD as the CNV"""
        genes = self.get_genes()
        enhancers = self.get_enhancer()
        self.gene_distances = []
        self.gene_overlaps = []
        if self.tads:
            for gene in genes:
                gene_overlap = utils.getOverlap([self.start,self.end],[gene.start,gene.end])
                if gene_overlap == 0:
                    gene_distance = utils.getDistance([self.start,self.end],[gene.start,gene.end])
                else:
                    gene_distance = 0
                self.gene_distances.append(gene_distance)
                self.gene_overlaps.append(gene_overlap)


            self.enhancer_distances = []
            self.enhancer_overlaps = []
            for enhancer in enhancers:
                enhancer_overlap = utils.getOverlap([self.start,self.end],[enhancer.start,enhancer.end])
                if enhancer_overlap == 0:
                    enhancer_distance = utils.getDistance([self.start,self.end],[enhancer.start,enhancer.end])
                else:
                    enhancer_distance = 0
                self.enhancer_distances.append(enhancer_distance)
                self.enhancer_overlaps.append(enhancer_overlap)
